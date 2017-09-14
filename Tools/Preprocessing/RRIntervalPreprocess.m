function [cleanNN, cleantNN, flagged_beats] = RRIntervalPreprocess(rr,time,annotations,sqi,HRVparams)
%
%   [cleanNN, cleantNN, flagged_beats] = RRIntervalPreprocess(rr, time, annotations, sqi, HRVparams)
%
%   OVERVIEW:   This function preprocesses RR interval data in preparation
%               for running HRV analyses. Noise and non-normal beats are
%               removed from the RR interval data and replaced with
%               interpolated data using a interpolation method of choice.
%
%   INPUT:      rr          - a single row of rr interval data in seconds
%               time        - the time indices of the rr interval data 
%                             (seconds)
%               annotations - (optional) annotations of the RR data at each
%                              point indicating the quality of the beat
%               sqi         - Signal Quality Index; Requires a matrix with
%                             at least two columns. Column 1 should be
%                             timestamps of each sqi measure, and Column 2
%                             should be SQI on a scale from 1 to 100.
%                             Additional columns can be included with
%                             additional sqi at the same timestamps
%               HRVparams   - struct of settings for hrv_toolbox analysis
%
%   OUTPUT:     cleanNN       - normal normal interval data
%               cleantNN      - time indices of NN data
%               flagged_beats - percent of original data removed
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%       Clifford, G. (2002). "Characterizing Artefact in the Normal 
%       Human 24-Hour RR Time Series to Aid Identification and Artificial 
%       Replication of Circadian Variations in Human Beat to Beat Heart
%       Rate using a Simple Threshold."
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)    
%
%   09-013-2017
%   Edited by Giulia Da Poian
%   Updated the method for measuring changes in the current RR interval 
%   from the last N (N=5) intervals and excluding intervals that 
%   change by more than a certain percentage define by 
%   HRVparams.preprocess.per_limit
%
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%


% CINC 2002 - Cite for rationale for % good data for cutoff

% 2004 paper - var goes up with numbers of beats dropped out
% don't put phantom beats in for the lomb

% check input
if isempty(time)
    time = cumsum(rr);
end
if isempty(sqi)
    % Assume SQI is good and no signal needs to be removed
    sqi(:,1) = 1:1:time(end);
    sqi(:,2) = 100*ones(length(sqi(:,1)),1);
end
if isempty(annotations)
    annotations = repmat('N',[length(rr) 1]);
end

if nargin < 5
	error('Not enough input arguments!')
end

figures = HRVparams.gen_figs;

%% 1. Identify data that is too close together and Remove
% These are not counted towards the total signal removed
idx_remove = find(diff(time) < 1/HRVparams.Fs); 
% could make this the refractory period - and have the variable in the settings
% document
rr(idx_remove+1) = [];
time(idx_remove) = [];
annotations(idx_remove) = [];
clear idx_remove;

%% 2. Find Artifact Annotations and Remove Those Points
% These are not counted towards the total signal removed
idx_remove = strcmp('|', annotations); % find artifacts
rr(idx_remove) = [];
annotations(idx_remove) = [];
time(idx_remove) = [];

idx_remove2 = strcmp('~', annotations); % find artifacts
rr(idx_remove2) = [];
annotations(idx_remove2) = [];
time(idx_remove2) = [];
clear idx_remove idx_remove2

%% 3. Remove Large Gaps in Data At Beginning of Dataset
%if t(1) > settings.preprocess.forward_gap
%    t = t - t(1);
%end

%% 4. Remove Large RR intervals Caused by Gaps
% These are not counted towards the total signal removed
idx_remove = find(rr >= HRVparams.preprocess.gaplimit);
rr(idx_remove) = [];
time(idx_remove) = [];
annotations(idx_remove) = [];
clear idx_remove;

%% Detect PVCs and Implement Phantom Beat
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

%% 5. Find Non-'N' Annotations
[~,~,ann_outliers] = AnnotationConversion(annotations);
if ~(length(annotations) == length(rr))
    ann_outliers = zeros(length(rr),1);
end


%% 6. Find RR Over Given Percentage Change 
perLimit = HRVparams.preprocess.per_limit;
idxRRtoBeRemoved = FindSpikesInRR(rr, perLimit); 

%% 7. Find Long Running Outliers

    

%% 8. Combine Annotations and Percentage Outliers
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

%% 9. Remove or Interpolate Outliers
idx_outliers = find(outliers == 1);

% Keep count of outliers 
numOutliers = length(idx_outliers);

rr_before_interp = rr;
rr(idx_outliers) = NaN;
if strcmp(HRVparams.preprocess.method_outliers,'cub')
    NNinterp = interp1(time,rr,time,'spline','extrap');
    tinterp = time;
elseif strcmp(HRVparams.preprocess.method_outliers,'pchip')
    NNinterp = interp1(time,rr,time,'pchip');
    tinterp = time;
elseif strcmp(HRVparams.preprocess.method_outliers,'lin')
    NNinterp = interp1(time,rr,time,'linear','extrap'); 
    tinterp = time;
elseif strcmp(HRVparams.preprocess.method_outliers,'rem')
    NNinterp = rr;
    NNinterp(idx_outliers) = [];
    tinterp = time;
    tinterp(idx_outliers) = []; 
else
    % By default remove outliers
    NNinterp = rr;
    NNinterp(idx_outliers) = [];
    tinterp = time;
    tinterp(idx_outliers) = [];
end

if figures
    figure;
    plot(time,rr_before_interp,tinterp,NNinterp);
    legend('raw','interp1(after outliers removed)')
end

%% 10. Identify Non-physiologic Beats
toohigh = NNinterp > HRVparams.preprocess.upperphysiolim;    % equivalent to RR = 2
toolow = NNinterp < HRVparams.preprocess.lowerphysiolim;     % equivalent to RR = .375

idx_toolow = find(toolow == 1);
NNinterp2 = NNinterp;
NNinterp2(idx_toolow) = NaN;
numOutliers = numOutliers + length(idx_toolow);


if strcmp(HRVparams.preprocess.method_unphysio,'cub')
    NNinterp2 = interp1(tinterp,NNinterp2,tinterp,'spline','extrap');
    tinterp2 = tinterp;
    flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
elseif strcmp(HRVparams.preprocess.method_unphysio,'pchip')
    NNinterp2 = interp1(tinterp,NNinterp2,tinterp,'pchip');
    tinterp2 = tinterp;
    flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
elseif strcmp(HRVparams.preprocess.method_unphysio,'lin')
    NNinterp2 = interp1(tinterp,NNinterp2,tinterp,'linear','extrap'); 
    tinterp2 = tinterp;
    flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
elseif strcmp(HRVparams.preprocess.method_unphysio,'rem')
    NNinterp2(idx_toolow) = [];
    tinterp2 = tinterp;
    tinterp2(idx_toolow) = []; % Review this line of code for improvement
    
else
    % Default is cubic spline
    NNinterp2 = interp1(tinterp,NNinterp2,tinterp,'pchip');
    tinterp2 = tinterp;
    flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
end

if figures
    hold on;
    plot(tinterp2,NNinterp2+.01);
    hold on; plot(time,toolow,'o')
    legend('raw','interp1(after outliers removed)',...
        'interp2(after too low)','toolow')
end



%% 11. Interpolate Through Beats that are Too Fast
toohigh = NNinterp2 > HRVparams.preprocess.upperphysiolim;    % equivalent to RR = 2

idx_outliers_2ndPass = find(logical(toohigh(:)) ~= 0);
NNinterp3 = NNinterp2;
NNinterp3(idx_outliers_2ndPass) = NaN;
numOutliers = numOutliers + length(idx_outliers_2ndPass);
if strcmp(HRVparams.preprocess.method_unphysio,'rem')
    flagged_beats = numOutliers;
end

if strcmp(HRVparams.preprocess.method_outliers,'cub')            
    NNinterp3 = interp1(tinterp2,NNinterp3,tinterp2,'spline','extrap');
    tinterp3 = tinterp2;
elseif strcmp(HRVparams.preprocess.method_outliers,'pchip')
    NNinterp3 = interp1(tinterp2,NNinterp3,tinterp2,'pchip');
    tinterp3 = tinterp2;
elseif strcmp(HRVparams.preprocess.method_outliers,'lin')
    NNinterp3 = interp1(tinterp2,NNinterp3,tinterp2,'linear','extrap');
    tinterp3 = tinterp2;
elseif strcmp(HRVparams.preprocess.method_outliers,'rem')
    NNinterp3(idx_outliers_2ndPass) = [];
    tinterp3 = tinterp2;
    tinterp3(idx_outliers_2ndPass) = []; % Review this line of code for improvement
else
    % Default is cubic spline
    NNinterp3 = interp1(tinterp2,NNinterp3,tinterp2,'spline','extrap');
    tinterp3 = tinterp2;
end

if figures
    hold on;
    plot(tinterp3,NNinterp3);
    legend('raw','interp1(after outliers removed)',...
        'interp2(after too low)','toolow','interp3 (after too fast removed)')
end

%% 12. Remove erroneous data at the end of a record 
%       (i.e. a un-physiologic point caused by removing data at the end of
%       a record)

while NNinterp3(end) > HRVparams.preprocess.upperphysiolim	% equivalent to RR = 2
    NNinterp3(end) = [];
    tinterp3(end) = [];
end


cleanNN = NNinterp3;
cleantNN = tinterp3;
end

