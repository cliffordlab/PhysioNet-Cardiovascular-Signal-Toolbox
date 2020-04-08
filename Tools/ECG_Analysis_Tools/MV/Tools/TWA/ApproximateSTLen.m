function stlen = ApproximateSTLen(s, stend)
% OVERVIEW, This function gives an estimate of the stlen length for the
% analysis window.
% 
% INPUT         MANDATORY       DESCRIPTION
%               s               Array of s point locations in analysis window.
%
%               stend           Array of s-toff intervals in current analysis window.
% 
% OUTPUT        stlen           median length of st segment for
%                               the current analysis window.
%
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Shamim Nemati, 
%       editted by Ismail Sadiq on 10/26/2019. 
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  

if (max(stend) == 0)
    % if no st segments detected consider it to be 0.5 times the rr
    % interval
    rr = mean(s(2:length(s)) - s(1:length(s) - 1));
    stlen = ceil(0.5 * rr);
else
    % otherwise consider it as median st segment.
    stlen = floor(nanmedian(stend(find(stend))));
end;
    
return;