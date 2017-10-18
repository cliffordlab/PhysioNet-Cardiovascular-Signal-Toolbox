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
%
%   Written by Giulia Da Poian (giulia.dap@gmail.com) 


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
    cleanRRAnalysisWindows(idx_af) = NaN;  % Remove all the windows with AF   
end

