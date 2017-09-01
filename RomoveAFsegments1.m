function cleanRRAnalysisWindows = RomoveAFsegments1(tNN,NN,AfAnalysisWindows, AFtest, HRVparams)

%   cleanRRAnalysisWindows = RomoveAFsegments(RRAnalysisWindows,AfAnalysisWindows, AFtest, HRVparams )
%
%	OVERVIEW:
%       Remove windows containing AF segments 
%
%   INPUT:
%       tNN               : 
%       NN                :
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
%       Script written by Giulia Da Poian
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 





% Exclude AF Segments
idx_afsegs = AfAnalysisWindows(AFtest == 1);
cleanRRAnalysisWindows = [];

if ~isempty(idx_afsegs) 
    
    for kk = 1:length(idx_afsegs)+1
        if kk==1
            tstart = 0;
            tstop = idx_afsegs(kk);
        elseif kk==length(idx_afsegs)+1
            tstart = idx_afsegs(kk-1)+HRVparams.af.increment;
            tstop = tNN(end);
        else
            tstart = idx_afsegs(kk-1)+HRVparams.af.increment;
            tstop = idx_afsegs(kk);
        end   
        
        tmp_NN = NN(tNN>=tstart & tNN <= tstop);
        tmp_tNN = tNN(tNN>=tstart & tNN <= tstop)-tstart;
        if ~isempty(tmp_NN) 
           cleanRRAnalysisWindows = [cleanRRAnalysisWindows CreateWindowRRintervals(tmp_tNN, tmp_NN, HRVparams)+tstart];
        end
        
    end
end