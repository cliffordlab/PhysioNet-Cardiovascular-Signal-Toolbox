function mf = meanfilt1(data,n)

for i=1:length(data)
    mf(i)=mean(data(max(1,i-round((n-1)/2)):min(i+round((n-1)/2),length(data))));
end
    