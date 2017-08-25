function cleanRRAnalysisWindows = RomoveAFsegments(RRAnalysisWindows, AFtest)

%   cleanRRAnalysisWindows = RomoveAFsegments(RRAnalysisWindows, AFtest)
%
%	OVERVIEW:
%       Perform Atrial Fibrillation (AF) detection 
%
%   INPUT:
%       RRAnalysisWindows : string containing the identifier of the subject to be analyze      
%       AFtest            :
%   OUTPUT:
%       cleanRRAnalysisWindows :  
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
idx_afsegs = (find(AFtest == 1));
if ~isempty(idx_afsegs)
    afsegs = AfAnalysisWindows(idx_afsegs);	% afsegs is in seconds
    for k = 1:length(afsegs)
        try
            idx_af(k) = find(RRAnalysisWindows <= afsegs(k) & HRVparams.increment + RRAnalysisWindows > afsegs(k));
        catch
        end
    end
    cleanRRAnalysisWindows(idx_af) = NaN;
end