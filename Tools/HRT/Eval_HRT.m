function [TO, TS,nPVCs] = Eval_HRT(RRInts, tRRInts, Labels, sqi, HRVparams, WindIdxs)

%   [TO, TSn, PVCs] = Eval_HRT(RRInts,Labels, HRVparams, WindIdxs)
%   OVERVIEW:
%       This function return TO and TS, i.e., the basic parameters of 
%       heart rate turbulence (HRT) used to quantify the return to 
%       equilibrium of heart rate  after a premature ventricular 
%       contraction (PVC) for each window
%
%   INPUTS:
%       RRInts          : Vector containing RR intervals data (in seconds)
%       tRRInts         : 
%       Labels          : Vector containing annotations of the RR data at 
%                         each point indicating the type of the beat (see 
%                         https://www.physionet.org/physiobank/annotations.shtml)
%       sqi             :
%       HRVparams       :
%       WindIdxs        :
%
%   OUTPUTS:
%       TO         : average turbulence onset (TO) 
%       TS         : turbulence slop (TS) of the average tachogram
%       nPVCs      : numer of PVCs used to build the average tachogram and
%                    compute the mean TO
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
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   ORIGINAL SOURCE AND AUTHORS:     
%       Giulia Da Poian    
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
  

% Verify input arguments
if isempty(WindIdxs)
    WindIdxs = 0;   
end
if isempty(sqi) 
     sqi(:,1) = tRRInts;
     sqi(:,2) = ones(length(tRRInts),1);
end

% Preallocate arrays (all NaN) before entering the loop
TO = nan(length(WindIdxs),1);
TS = nan(length(WindIdxs),1);
nPVCs = nan(length(WindIdxs),1);

BeatsBefore = HRVparams.HRT.BeatsBefore;
BeatsAfter = HRVparams.HRT.BeatsAfter;
GraphOn = HRVparams.HRT.GraphOn;
windowlength = HRVparams.HRT.windowlength*3600; % Convert hours in seconds
SQI_LowQualityThresh = HRVparams.sqi.LowQualityThreshold;
RejectionThreshold = HRVparams.RejectionThreshold;
filterMethod = HRVparams.HRT.filterMethod;

%Analyze by Window

% Loop through each window of RR data
for iWin = 1:length(WindIdxs)
    % Check window for sufficient data
    if ~isnan(WindIdxs(iWin))
        % Isolate data in this window
        idxRRinWin = find(tRRInts >= WindIdxs(iWin) & tRRInts < WindIdxs(iWin) + windowlength);
        
        SQIinWin = sqi(sqi(:,1) >= WindIdxs(iWin) & sqi(:,1) < WindIdxs(iWin) + windowlength,:);
        RRinWin = RRInts(idxRRinWin);
        LabelsinWin = Labels(idxRRinWin);

        % Analysis of SQI for the window
        LowQualityIdxs = find(SQIinWin(:,2) < SQI_LowQualityThresh);

        % If enough data has an adequate SQI, perform the calculations
        if numel(LowQualityIdxs)/length(SQIinWin(:,2)) < RejectionThreshold
            [TO(iWin), TS(iWin), nPVCs(iWin)] = HRT_Analysis(RRinWin,LabelsinWin,BeatsBefore, BeatsAfter, filterMethod, GraphOn);
        end     

    end % end check for sufficient data
    
end % end of loop through window



end % end of Eval_HRT function



