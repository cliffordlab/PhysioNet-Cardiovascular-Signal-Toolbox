function [cleanRRAnalysisWindows,idx_af] = RemoveAFsegments(RRAnalysisWindows,AfAnalysisWindows, AFtest, HRVparams)

%   cleanRRAnalysisWindows = RemoveAFsegments(RRAnalysisWindows,AfAnalysisWindows, AFtest, HRVparams )
%
%	OVERVIEW:
%       Remove windows containing AF segments 
%
%   INPUT:
%       RRAnalysisWindows : vector containing windows indexes corresponding 
%                           to RR interval segmentation for HRV analysis    
%       AfAnalysisWindows : vector containing windows indexes corresponding 
%                           to RR interval segmentation for AF analysis
%       AFtest            : vector containing flags (1) for AF windows and 
%                           (0) for non-AF segemnts 
%       HRVparams         : 
%   OUTPUT:
%       cleanRRAnalysisWindows : vector containing windows indexes corresponding 
%                                to RR interval segmentation for HRV analysis
%                                where windows containing AF are identified
%                                by a NaN index and therfore exluded during
%                                analysis
%       idx_af                 : indexes of the AF windows
%
%   DEPENDENCIES & LIBRARIES:
%       PhysioNet Cardiovascular Signal Toolbox
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   REFERENCE: 
%   Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%   Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Giulia Da Poian (giulia.dap@gmail.com) on Sep 6, 2017.
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


cleanRRAnalysisWindows = RRAnalysisWindows;

% Exclude AF Segments
afsegs = AfAnalysisWindows(AFtest == 1); % afsegs is in seconds
idx_af = [];

if ~isempty(afsegs)
    for k = 1:length(afsegs)
        try
            idx_af = [idx_af find(RRAnalysisWindows <= afsegs(k) & HRVparams.windowlength + RRAnalysisWindows > afsegs(k))];
        catch
        end
    end
    idx_af = unique (idx_af); % Do not consider duplicate positions
    cleanRRAnalysisWindows(idx_af) = NaN;  % Remove all the windows with AF   
end

