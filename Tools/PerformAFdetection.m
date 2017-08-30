function [afresults, AfAnalysisWindows] = PerformAFdetection(subjectID,tNN,NN,HRVparams)    
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
%       HRVparam  : struct of settings for hrv_toolbox analysis
%
%   OUTPUT:
%       afresults : a single row containing a flag (1) when AF is
%                   dettected in a window and 0 if no AF 
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

% 1. Calculate AF Features
AfAnalysisWindows = CreateWindowRRintervals(tNN,NN,HRVparams,'af');
NN_afcalc = NN .* HRVparams.Fs;

for i = 1:length(AfAnalysisWindows)
    tstart = features(i);
    if isnan(tstart) 
        features_af = NaN;
        AFtest(i) = NaN;
    else
        idx_af = find(tNN >= tstart & tNN < tstart + HRVparams.af.windowlength);
        features_af = AF_features(NN_afcalc(idx_af),HRVparams.Fs);
        AFtest(i) = SVM_AFdetection_withoutTrainingModel(features_af,1);
    end    
end


% 2. Export AF Data as CSV File

afresults = AFtest(:);
afcol_titles = {'AFtest'};

type = 'AF';
GenerateHRVresultsOutput(subjectID,AfAnalysisWindows,afresults,afcol_titles, type, HRVparams, tNN, NN);
