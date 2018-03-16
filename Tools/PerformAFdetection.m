function [afresults, AfAnalysisWindows, AFfile] = PerformAFdetection(subjectID,tNN,NN,sqi,HRVparams)    
%   PerformAFdetection(subjectID,tNN,NN,HRVparams)  
%
%	OVERVIEW:
%       Perform Atrial Fibrillation (AF) detection 
%
%   INPUT:
%       subjectID : string containing the identifier of the subject to be analyze      
%       tNN       : a single row of time indices of the rr interval 
%                   data (seconds)
%       NN        : a single row of NN (normal normal) interval
%                   data in seconds
%       sqi       : Signal Quality Index; Requires a matrix with
%                   at least two columns. Column 1 should be timestamps 
%                   should be timestamps of each sqi measure, and Column 2 
%                   should be SQI on a scale from 0 to 1.
%       HRVparam  : struct of settings for hrv_toolbox analysis
%
%   OUTPUT:
%       afresults : a single row containing a flag (1) when AF is
%                   dettected in a window and 0 if no AF 
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
%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Giulia Da Poian (giulia.dap@gmail.com) on Sep 6, 2017.
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Modified on 02.14.2018 to include SQI check before evaluating AF; if a
%   windows contains low SQI signal do not compute AF and mark the window
%   as NaN


if isempty(sqi) 
     sqi(:,1) = tNN;
     sqi(:,2) = ones(length(tNN),1);
end

% 1. Calculate AF Features
AfAnalysisWindows = CreateWindowRRintervals(tNN,NN,HRVparams,'af');
NNsamps = NN .* HRVparams.Fs;

AFtest = nan(length(AfAnalysisWindows),1);

for idx = 1:length(AfAnalysisWindows)
    tstart = AfAnalysisWindows(idx);
    
    if ~isnan(tstart)  
      
        idxInWin = find(tNN >= tstart & tNN< tstart + HRVparams.af.windowlength);
       
        % Added by Giulia: exclude from the Analysis low quality
        % segments (SQI< SQI_threshold)
        sqiWin = sqi(sqi(:,1) >= tstart & sqi(:,1) < tstart + HRVparams.af.windowlength,2);
        LowQualityIdxs = find(sqiWin < HRVparams.sqi.LowQualityThreshold); 
        
        % If enough data has an adequate SQI, perform the calculations
        cond1 = (numel(LowQualityIdxs)/length(sqiWin)) <  HRVparams.RejectionThreshold;        
        % RR interval time series must be > 12 ans <60 to extract feautures
        cond2 = (length(NNsamps(idxInWin)) > 12 && length(NNsamps(idxInWin)) < 60);
  
        if (cond1 && cond2)
            features_af = AF_features(NNsamps(idxInWin),HRVparams.Fs);
            AFtest(idx) = SVM_AFdetection_withoutTrainingModel(features_af,1);            
        end
    end    
end


% 2. Export AF Data as CSV File

afresults = AFtest(:);
afcol_titles = {'AFtest'};

outputType = 'AF';
AFfile = SaveHRVoutput(subjectID,AfAnalysisWindows,afresults, ...
    afcol_titles, outputType, HRVparams, tNN, NN);
