function [cleanNN, cleantNN, flagged_beats] = RRIntervalPreprocess(rr,t_rr,annotations,HRVparams)
%
%   [cleanNN, cleantNN, flagged_beats] = RRIntervalPreprocess(rr, time, annotations, HRVparams)
%
%   OVERVIEW:   This function preprocesses RR interval data in preparation
%               for running HRV analyses. Noise and non-normal beats are
%               removed from the RR interval data and replaced with
%               interpolated data using a interpolation method of choice.
%
%   INPUT:      rr          - a single row of rr interval data in seconds
%               t_rr        - time of the rr interval data 
%                             (seconds)
%               annotations - (optional) annotations of the RR data at each
%                              point indicating the quality of the beat
%               HRVparams   - struct of settings for hrv_toolbox analysis
%
%   OUTPUT:     cleanNN       - normal normal interval data
%               cleantNN      - time of NN interval
%               flagged_beats - percent of original data removed
%   REFERENCE: 
%       Clifford, G. (2002). "Characterizing Artefact in the Normal 
%       Human 24-Hour RR Time Series to Aid Identification and Artificial 
%       Replication of Circadian Variations in Human Beat to Beat Heart
%       Rate using a Simple Threshold."
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)    
%
%   09-013-2017
%   Edited by Giulia Da Poian
%   -Updated the method for measuring changes in the current RR interval 
%   from the last N (N=5) intervals and excluding intervals that 
%   change by more than a certain percentage define by 
%   HRVparams.preprocess.per_limit
%   - Removed SQI, windows containing low quality SQI segments are removed  
%     in EvalTimeDomainHRVstats and EvalFrequencyDomainHRVstats
%
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%


% CINC 2002 - Cite for rationale for % good data for cutoff

% 2004 paper - var goes up with numbers of beats dropped out
% don't put phantom beats in for the lomb

% check input
if isempty(t_rr)
    t_rr = cumsum(rr);
end

if isempty(annotations)
    annotations = repmat('N',[length(rr) 1]);
end

if nargin < 4
	error('Not enough input arguments!')
end

figures = HRVparams.gen_figs;

% 1. Identify data that is too close together and Remove
% These are not counted towards the total signal removed
idx_remove = find(diff(t_rr) < 1/HRVparams.Fs); 
% could make this the refractory period - and have the variable in the settings
% document
rr(idx_remove+1) = [];
t_rr(idx_remove) = [];
annotations(idx_remove) = [];
clear idx_remove;

% 2. Find Artifact Annotations and Remove Those Points
% These are not counted towards the total signal removed
idx_remove = strcmp('|', annotations); % find artifacts
rr(idx_remove) = [];
annotations(idx_remove) = [];
t_rr(idx_remove) = [];

idx_remove2 = strcmp('~', annotations); % find artifacts
rr(idx_remove2) = [];
annotations(idx_remove2) = [];
t_rr(idx_remove2) = [];
clear idx_remove idx_remove2

% 3. Remove Large Gaps in Data At Beginning of Dataset
%if t(1) > settings.preprocess.forward_gap
%    t = t - t(1);
%end

% 4. Remove Large RR intervals Caused by Gaps
% These are not counted towards the total signal removed
idx_remove = find(rr >= HRVparams.preprocess.gaplimit);
rr(idx_remove) = [];
t_rr(idx_remove) = [];
annotations(idx_remove) = [];
clear idx_remove;

% Detect PVCs and Implement Phantom Beat
% Look for premature ectopic beats
% replace value of ectopic beat with a linear interpolated point
% replace sample number
% now next beat becomes normal
% with more than two Vent ectopic beats, cut off data.
% pvcthresh = .3;
% rrdiff = diff(rr);
% idx = find(rrdiff > pvcthresh);
% for i = 1:length(idx)
%     ppvc = rr(idx(i)+1);
%     if rr(idx(i)+1+1) > rr(idx(i)+1)
%         if (rr(idx(i)+1+2) - rr(idx(i))) < pvcthresh
%             % If all is true, we can assume we have a PVC
%             gap = rr(idx(i)+1+2)-rr(idx(i)+1-1);
%             even = gap/3;
%             rr(i) = NaN;
%             rr(i+1) = NaN;
%             avgrr = mean(rr(i-10:i-1));
%             rr(i) = even + rr(i-1);
%             rr(i+2) = even + rr(i);
%         end
%     end
% end

%hold on; plot(t,rr);

% 5. Find Non-'N' Annotations
[~,~,ann_outliers] = AnnotationConversion(annotations);
if ~(length(annotations) == length(rr))
    ann_outliers = zeros(length(rr),1);
end


% 6. Find RR Over Given Percentage Change 
perLimit = HRVparams.preprocess.per_limit;
idxRRtoBeRemoved = FindSpikesInRR(rr, perLimit); 

% 7. Find Long Running Outliers

    

% 8. Combine Annotations and Percentage Outliers
% Combination of methods
outliers_combined = idxRRtoBeRemoved(:) + ann_outliers(:);
outliers = logical(outliers_combined); 

% remove extra points that may have been introduced at the end of the file
if (length(outliers) > length(ann_outliers)) 
    fprintf('Too many Outliers - Removing last Point\n');
    z = length(outliers);
    outliers(z) = [];
end

% timestamp = [find(diff([-1; outliers; -1]) ~= 0)]; % where does outliers change?
% runlength = diff(timestamp); % how long is each segment
% % to split them for 0's and 1's
% runlength0 = runlength(1+(outliers(1)==1):2:end);
% runlength1 = runlength(1+(outliers(1)==0):2:end);
% 
% r0_rmv_small = (find(runlength0 > 10));
% r1_rmv_small = (find(runlength1 > 10));
% idx_groups_outliers =  find(runlength1 > 10);

% 9. Remove or Interpolate Outliers
idx_outliers = find(outliers == 1);

% Keep count of outliers 
numOutliers = length(idx_outliers);

rr_original = rr;
rr(idx_outliers) = NaN;

switch HRVparams.preprocess.method_outliers
    case 'cub'
        NN_Outliers = interp1(t_rr,rr,t_rr,'spline','extrap');
        t_Outliers = t_rr;
    case 'pchip'
        NN_Outliers = interp1(t_rr,rr,t_rr,'pchip');
        t_Outliers = t_rr;
    case 'lin'
        NN_Outliers = interp1(t_rr,rr,t_rr,'linear','extrap'); 
        t_Outliers = t_rr;
    case 'rem'
        NN_Outliers = rr;
        NN_Outliers(idx_outliers) = [];
        t_Outliers = t_rr;
        t_Outliers(idx_outliers) = []; 
    otherwise % By default remove outliers
        NN_Outliers = rr;
        NN_Outliers(idx_outliers) = [];
        t_Outliers = t_rr;
        t_Outliers(idx_outliers) = [];
end

if figures
    figure;
    plot(t_rr,rr_original,t_Outliers,NN_Outliers);
    legend('raw','interp1(after outliers removed)')
end

% 10. Identify Non-physiologic Beats
toohigh = NN_Outliers > HRVparams.preprocess.upperphysiolim;    % equivalent to RR = 2
toolow = NN_Outliers < HRVparams.preprocess.lowerphysiolim;     % equivalent to RR = .375

idx_toolow = find(toolow == 1);
NN_NonPhysBeats = NN_Outliers;
NN_NonPhysBeats(idx_toolow) = NaN;
numOutliers = numOutliers + length(idx_toolow);


switch HRVparams.preprocess.method_unphysio
    case 'cub'
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'spline','extrap');
        t_NonPhysBeats = t_Outliers;
        flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
    case 'pchip'
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'pchip');
        t_NonPhysBeats = t_Outliers;
        flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
    case 'lin'
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'linear','extrap'); 
        t_NonPhysBeats = t_Outliers;
        flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
    case 'rem'
        NN_NonPhysBeats(idx_toolow) = [];
        t_NonPhysBeats = t_Outliers;
        t_NonPhysBeats(idx_toolow) = []; % Review this line of code for improvement
    otherwise % use cubic spline interpoletion as default
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'pchip');
        t_NonPhysBeats = t_Outliers;
end

if figures
    hold on;
    plot(t_NonPhysBeats,NN_NonPhysBeats+.01);
    hold on; plot(t_NonPhysBeats,toolow,'o')
    legend('raw','interp1(after outliers removed)',...
        'interp2(after too low)','toolow')
end



% 11. Interpolate Through Beats that are Too Fast
toohigh = NN_NonPhysBeats > HRVparams.preprocess.upperphysiolim;    % equivalent to RR = 2

idx_outliers_2ndPass = find(logical(toohigh(:)) ~= 0);
NN_TooFastBeats = NN_NonPhysBeats;
NN_TooFastBeats(idx_outliers_2ndPass) = NaN;
numOutliers = numOutliers + length(idx_outliers_2ndPass);
if strcmp(HRVparams.preprocess.method_unphysio,'rem')
    flagged_beats = numOutliers;
end

switch HRVparams.preprocess.method_outliers
    case 'cub'            
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'spline','extrap');
        t_TooFasyBeats = t_NonPhysBeats;
    case 'pchip'
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'pchip');
        t_TooFasyBeats = t_NonPhysBeats;
    case 'lin'
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'linear','extrap');
        t_TooFasyBeats = t_NonPhysBeats;
    case 'rem'
        NN_TooFastBeats(idx_outliers_2ndPass) = [];
        t_TooFasyBeats = t_NonPhysBeats;
        t_TooFasyBeats(idx_outliers_2ndPass) = []; % Review this line of code for improvement
    otherwise % USe cubic spline interpoletion as default
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'spline','extrap');
        t_TooFasyBeats = t_NonPhysBeats;
end

if figures
    hold on;
    plot(t_TooFasyBeats,NN_TooFastBeats);
    legend('raw','interp1(after outliers removed)',...
        'interp2(after too low)','toolow','interp3 (after too fast removed)')
end

% 12. Remove erroneous data at the end of a record 
%       (i.e. a un-physiologic point caused by removing data at the end of
%       a record)

while NN_TooFastBeats(end) > HRVparams.preprocess.upperphysiolim	% equivalent to RR = 2
    NN_TooFastBeats(end) = [];
    t_TooFasyBeats(end) = [];
end


cleanNN = NN_TooFastBeats;
cleantNN = t_TooFasyBeats;
end

