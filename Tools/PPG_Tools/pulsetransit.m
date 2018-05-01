function ptt = pulsetransit(ecgann, pulseann)
%   ptt = pulsetransit(ecgann, pulseann)
%
%   OVERVIEW:   Calculate PTT
%               match pulse annotations with corresponding ecg annotations
%
%   INPUT:      
%
%   OUTPUT:     
%
%   DEPENDENCIES & LIBRARIES:
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Main script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%

ptt = NaN((length(ecgann)-1),3);
for i = 1:(length(ecgann)-1)
    x = find(pulseann < ecgann(i+1) & pulseann > ecgann(i));
    ptt(i,1) = ecgann(i);
    if x > 0
        if length(x)>1
            ptt(i,2) = pulseann(x(1));
        else
            ptt(i,2) = pulseann(x);
        end
    end
    ptt(i,3) = ptt(i,2) - ptt(i,1); 
end




end