function [ymodify y2modify r] = draw_dtw(y1,pla1,p,y2,pla2,q)

l=length(p);
i=1;
j=1;
la=pla1(p(1));
lb=pla2(q(1));
while p(i+1)==p(1)
    i=i+1;
end
while q(j+1)==q(1)
    j=j+1;
end
point=max(i,j)+1;
outputx=[];
outputy=[];
while (point<=l)
    tp=1;
    tq=1;
    while (point+1<=l && ((p(point+1)==p(point)) || (q(point+1)==q(point))))
        point=point+1;
    end
    if point<=l
        laold=la;
        lbold=lb;
        la=pla1(p(point));
        lb=pla2(q(point));
        %intvb=(la-laold)/(lb-lbold);
        intvb=(lb-lbold)/(la-laold)/10;
        xx=pla1(p(i)):intvb:pla1(p(point));
        yy=y2(pla2(q(j)):pla2(q(point)));
        x1=xx(1):(xx(length(xx))-xx(1))/(length(yy)-1):xx(length(xx));
        if length(xx)<=1
            y=yy(end);
        else
            y = griddedInterpolant(x1,yy,'spline');
            y = y(xx);
            % y = interp1(x1,yy,xx,'spline'); % 12-19-2017 Modified by Giulia Da Poian
            % replaced interp1 with griddedInterpolant for speed
        end
%        plot(xx,y,'k');
        outputx=[outputx,xx];
        outputy=[outputy,y];
        i=point;
        j=point;
    end
    point=point+1;
end
% outputxx=unique(outputx);
j=1;
i=1;    
while i<=length(outputx)-1 
    while i<length(outputx)-1 && outputx(i)==outputx(i+1)
        i=i+1;
    end
    outputxx(j)=outputx(i);
    outputyy(j)=outputy(i);
    i=i+1;    
    j=j+1;
end
outputxx(j)=outputx(i);
outputyy(j)=outputy(i);
% 12-01-2017 added by Giulia Da Poian to solve Error: The grid vectors must contain unique points
[~, uidx] = unique(outputxx);
y=interp1(outputxx(uidx),outputyy(uidx),1:(outputxx(length(outputxx(uidx)))-1)/(length(y1)-1):outputxx(length(outputxx(uidx))),'spline');

y(isnan(y))=0;
%plot(y,'m');
y2modify=y;
% ylength=min(y,y1);
[S,R] = size(y);
[R2,Q] = size(y1);
if R ~= R2
    difference=dist(y,y1');
else
    difference=dist(y,y1);
end
    
meany1=sqrt(sum(y1.^2));
r=difference/meany1;
if r<0.0010
    ymodify=y1.*0.9+y'.*0.1;
else
    ymodify=y1;
end