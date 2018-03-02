function COSEn = COSEn(data, m, r,fs)
% COSEn	
%           Ref.    Douglas E. Lake and J. Randall Moorman
%                   Accurate estimation of entropy in very short physiological time series: the problem of atrial fibrillation detection in implanted ventricular devices
%                   Am J Physiol Heart Circ Physiol 300: H319¨CH325, 2011
% 
% $ Modified:2015.10.10 by Chengyu Liu
%           @ Emory
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information


% data=[179,185,184,183,184,184,185,185,180,178,177,176]*4;
% data=[159   157   160   160   159   160   160   157   160   158   157   157]*4;
% data=[240   240   239   240   239   240   239   239   239   239   239
% 240]*4;
% data=[109,108,108,109,109,109,108,109,108,108,109,108,108,109,108,108,108,109,108,109,109,109,108,108,109,108,108,109,109,109]*4;
% m=1;
% r=30;

N=length(data);
for i=1:N-m
    x1(i,:)=data(i:i+m-1);
    x2(i,:)=data(i:i+m);
end

ratio=0.4;
if N>20
%     Thr=ratio*N;
    Thr=5;
else
    Thr=5;
end

Min_numerator=0;
kk=0;
while Min_numerator<Thr
    [Min_numerator,r]=ReGet_min_numerator(x1,x2,N,m,r,fs);
    kk=kk+1;    
end
if kk==1 || kk>20
    if N<20
        r=r-1;
    else
        r=r-1000/fs;
    end
end

if Min_numerator==N-m
    while Min_numerator>=Thr
        [Min_numerator,r]=ReGet_min_numerator2(x1,x2,N,m,r,fs);
    end
    if  r<0
        r=r+2;
    else 
        r=r+1;
    end
end

for i=1:N-m
    t1=0;
    t2=0;
    for j=1:N-m
        d1max(i,j)=max(abs(x1(i,:)-x1(j,:)));
        d2max(i,j)=max(abs(x2(i,:)-x2(j,:)));
        if  d1max(i,j)<r
            t1=t1+1;
        end
        if  d2max(i,j)<r
                t2=t2+1;
        end        
        c1(i)=(t1-1)/(N-m-1);
        c2(i)=(t2-1)/(N-m-1);
    end
end
b1=0;
b2=0;
for i=1:N-m
    b1=b1+c1(i);
    b2=b2+c2(i);
end
B1=b1/(N-m);
B2=b2/(N-m);
COSEn=log(B1)/log(B2)+log(2*r/1000)-log(mean(data)/1000);



