% function [w ta tb] = simmx_dtw(y1,pla1,y2,pla2)
%
% calculate a sim matrix between y1 and y2.
%

function [w ta tb] = simmx_dtw(y1,pla1,y2,pla2)

% slope1(1)=y1(pla1(1));
% slope2(1)=y2(pla2(1));

slope1(1)=0;
slope2(1)=0;

ta(1)=1;
tb(1)=1;
for i=1:length(pla1)-1
    slope1(i+1)=(y1(pla1(i+1))-y1(pla1(i)))/(pla1(i+1)-pla1(i));
    ta(i+1)=pla1(i+1)-pla1(i);
end
for i=1:length(pla2)-1
    slope2(i+1)=(y2(pla2(i+1))-y2(pla2(i)))/(pla2(i+1)-pla2(i));
    tb(i+1)=pla2(i+1)-pla2(i);
end
al=length(slope1);
bl=length(slope2);
for i=1:bl
    A1(:,i)=slope1;
end
for j=1:al
    B1(j,:)=slope2';
end

w = abs(A1-B1);