function [ v1, tv1] = RR_Preprocessing_for_MSE( NN, tNN )
%RR_PREPROCESSING_FOR_MSE Summary of this function goes here
%   

%[NN, tNN, flagged_beats] = RRIntervalPreprocess(RR,tRR,annotations,HRVparams);

tv1 = [];
v1 = [];


th = 10;


% NN = flipud(NN(:));
% tNN = flipud(tNN(:));
v = abs(diff(tNN));
ii = 0;
for idx = 1:length(v)
    try
        medianValue = median(NN(idx-th: idx+th));
    catch 
        try
            medianValue = median(NN(idx-th: idx));
        catch
            try
                medianValue = median(NN(idx:idx+th));
            catch
                medianValue = median(NN(1:idx));
            end
        end
    end
    
    thresh = 2 * medianValue;
    
       if(min(v1)<0)
           5
       end        
            
    if v(idx) > thresh
       
       samps = floor(v(idx)/medianValue); 
       
%        tNN(idx+1)
       
       v1 = [v1 ; NN(idx)];
       tv1 = [tv1; tNN(idx)];
       ii=ii+1;
       % adding 'samps' points
       tv1 = [tv1; ones(samps,1)];
       v1 = [v1; ones(samps,1)];
       ii=ii+samps;
       tv1(end) = tNN(idx+1) - NN(idx+1);
       v1(end) = medianValue;
       
       for k = 1: samps-1
       
           tv1(end-k) = tv1(end-k+1) - v1(end-k+1);
           v1(end-k) = medianValue;
       
       end
       
       % chack if fisrt value in the gap is correct
       
       if tv1(ii-samps+1) - tNN(idx) < 0.5*medianValue
           tv1(ii-samps+1) = [];
           v1(ii-samps+1) = [];
           ii=ii-1;
           samps = samps -1;
       end   
       v1(ii-samps+1) = tv1(ii-samps+1) - tNN(idx);
           


       
%         v1 = [v1; NN(idx) ];
%         tv1 = [tv1; tNN(idx); tNN(idx)-medianValue];
%         
%         for k = 1:samps-1          
%             v1 = [v1; medianValue];
%             tv1 = [tv1; tv1(end)-medianValue];
%         end
%         
%         if abs(tNN(idx+1)-tv1(end)) <= 0.5*medianValue
%             tv1(end)=[];
%             v1(end) = abs(tNN(idx+1)-tv1(end));
%         else
%             v1 =[v1; abs(tNN(idx+1)-tv1(end))];
%         end
%    
        

    else
        v1 = [v1 ; NN(idx)];
        tv1 = [tv1; tNN(idx)];
        ii=ii+1;
    end

end

v1 = [v1; NN(end)];
tv1 = [tv1; tNN(end)];


