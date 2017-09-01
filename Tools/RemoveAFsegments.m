function cleanRRAnalysisWindows = RomoveAFsegments(RRAnalysisWindows,AfAnalysisWindows, AFtest, HRVparams)

%   cleanRRAnalysisWindows = RomoveAFsegments(RRAnalysisWindows,AfAnalysisWindows, AFtest, HRVparams )
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
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       This script written by Giulia Da Poian
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
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
    cleanRRAnalysisWindows(idx_af) = NaN;   
end

