function idxRRtoBeRemoved = FindSpikesInRR(RR, th)
% 
% idxRRtoBeRemoved = FindSpikesInRR(RR, th)
%
% OVERVIEW : Code used to clean RR intervals that change more than a given 
% threshold (eg., th = 0.2 = 20%) with respect to the median value of the 
% previous 5 and next 5 RR intervals (using a forward-backward approach).
%
% INPUTS:
%       RR : a single row of rr interval data in seconds
%       th : threshold percent limit of change from one interval to the next
% OUTPUTS:
%       idxRRtoBeRemoved : a single vector of indexes related to RR
%                          intervals corresponding to a change > th
%
%   DEPENDENCIES & LIBRARIES:
%       PhysioNet Cardiovascular Signal Toolbox
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   REFERENCE: 
%       Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%       Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   Written by Giulia Da Poian (giulia.dap@gmail.com), 09-13-2017
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%


if size(RR,1)>size(RR,2)
    RR = RR';
end

% Forward search 
FiveRR_MedianVal = medfilt1(RR,5); % compute as median RR(-i-2: i+2)
% shift of three position to aligne with to corresponding RR 
FiveRR_MedianVal = [RR(1:5) FiveRR_MedianVal(3:end-3)];
rr_above_th = find((abs(RR-FiveRR_MedianVal)./FiveRR_MedianVal)>=th);

RR_forward = RR;
RR_forward(rr_above_th) = NaN;


% Backward search 
RRfilpped = fliplr(RR);

FiveRR_MedianVal = medfilt1(RRfilpped,5); % compute as median RR(-i-2: i+2)
% shift of three position to aligne with to corresponding RR 
FiveRR_MedianVal = [RRfilpped(1:5) FiveRR_MedianVal(3:end-3)];
rr_above_th = find(abs(RRfilpped-FiveRR_MedianVal)./FiveRR_MedianVal>=th);
rr_above_th = sort(length(RR)-rr_above_th+1);


RR_backward = RRfilpped;
RR_backward(rr_above_th) = NaN;


% Combine 

idxRRtoBeRemoved = (isnan(RR_forward) & isnan(RR_backward));


end


