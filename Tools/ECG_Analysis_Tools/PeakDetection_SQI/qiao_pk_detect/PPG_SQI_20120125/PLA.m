% Piecewise linear approximation (PLA) 
% See: HJLM. Vullings, Automated ECG segmentation with Dynamic Time Warping 

function [y1 PLA] = PLA(input,s,th)

% input:    source data
% s:        step
% th:       threshold

if nargin<3
    th=10;
end
if nargin<2
    s=10;
end
    
n=length(input);
% Low pass filter 
% [b,a] = cheby1(3,0.5,40/62.5);
% y1=filter(b,a,input);
% b = fir1(31,40/62.5);
% y=conv(b,input);
% y1=y(int16((length(y)-length(input))/2):(int16((length(y)-length(input))/2+length(input)-1)));

s1=s;
pp=1;
PLA(pp)=1;
pp=pp+1;
i=1;
while i<n
    if (i+s1>=n)
        i_plus_s=n;
    else
        i_plus_s=i+s1;
    end
    interrupt=0;
    while (interrupt==0)
    j=i+1;
    while j<=i_plus_s
        distance=input(i_plus_s)-input(i);
        dcur=input(j)-input(i)-(distance*(j-i)/(i_plus_s-i));
        % distance is larger than threshold, truncate old line
        if abs(dcur)>th
            s1=j-i;
            i_plus_s=i+s1;
            j=i+1;
            interrupt=1;
            continue;
        end
        j=j+1;
    end
    if interrupt==1
        PLA(pp)=j-1;
        pp=pp+1;
        i=j-2;
        s1=s;
    else
        % distance is smaller than threshold, expand old line
        if (i_plus_s>=n)
            i_plus_s=n;
            break;
        else
            if ((i_plus_s+s1)>=n)
                i_plus_s=n;
            else
                i_plus_s=i_plus_s+s1;
            end
        end
    end
    end
    i=i+1;
end
PLA(pp)=n;
y1=input;    
% PLA=y1(int16((length(y1)-length(input))/2):(int16((length(y1)-length(input))/2+length(input)-1)));
