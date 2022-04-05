function predict_label = CNN_SVM_prediction(imdb,OP)
%
% predict_label:
%   if OP.class = 1 then: predict_label = 1:Wake, 2:REM sleep, 3:NREM light sleep, 4:NREM deep sleep
%   if OP.class = 2 then: predict_label = 1:NREM sleep (light + deep), 2:REM sleep, 3: wake
%   if OP.class = 3 then: predict_label = 1:NREM deep sleep, 2:NREM light sleep, 3:Wake + REM sleep
%   if OP.class = 4 then: predict_label = 1:NREM sleep (light + deep), 2: Wake + REM sleep

imdb.meta.features(isinf(imdb.meta.features))=NaN;
sqi=nanmean(imdb.meta.sqi,1);
sqi(find(isnan(sqi)))=nanmean(sqi);

% load model
load(OP.model);

net=stat.net;
net.layers{1,end}.type='softmax';
mean_features=stat.mean_features;
std_features=stat.std_features;

mean_sqi=stat.mean_sqi;
std_sqi=stat.std_sqi;

model=stat.model;

% standardization
for i=8:13
    for j=1:size(imdb.meta.features,1)
        imdb.meta.features(j,i)=log(imdb.meta.features(j,i));
    end
end
imdb.meta.features(find(isinf(imdb.meta.features)))=NaN;

% subtract data_mean
for k=1:size(imdb.images.data,4)
    imdb.images.data(:,:,1,k)=imdb.images.data(:,:,1,k)-stat.data_mean;
end

for i=1:size(imdb.meta.features,2)
    imdb.meta.features(find(isnan(imdb.meta.features(:,i))),i)=mean_features(i);
    imdb.meta.features(:,i)=bsxfun(@minus, imdb.meta.features(:,i),mean_features(i));
    imdb.meta.features(:,i)=bsxfun(@rdivide, imdb.meta.features(:,i),std_features(i));
end
sqi=bsxfun(@minus,sqi,mean_sqi);
sqi=bsxfun(@rdivide,sqi,std_sqi);

% CNN prediction
if size(imdb.images.data,4)>1000
    split=round(size(imdb.images.data,4)/1000);
else
    split=1;
end

retall=[];
for kk=1:split
    each=ceil(size(imdb.images.data,4)/split);
    res=vl_simplenn(net,imdb.images.data(:,:,:,(kk-1)*each+1:min(kk*each,size(imdb.images.data,4))));
    
    ret=[];
    for i=1:OP.nclass
        for j=1:size(res(end).x,4)
            ret(i,j)=res(end).x(1,1,i,j);
        end
    end
    retall=[retall,ret];
end
% probability of CNN output
ret=retall;

xtest=[];
for i=1:size(ret,1)
    xtest(:,end+1)=ret(i,:);
end
for i=1:13
    if i~=6 % remove DFA feature
        xtest(:,end+1)=imdb.meta.features(:,i);
    end
end
xtest(:,end+1)=sqi;

% SVM prediction
[predict_label, accuracy_test, ytest_out] = svmpredict(ones(size(xtest,1),1), xtest, model,'-b 1');
end


    
    
    


    



