function [afresults, AfAnalysisWindows, AFfile] = PerformAFdetection(subjectID,tNN,NN,HRVparams)    
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
%
%   Written by Giulia Da Poian (giulia.dap@gmail.com)



% 1. Calculate AF Features
AfAnalysisWindows = CreateWindowRRintervals(tNN,NN,HRVparams,'af');
NN_afcalc = NN .* HRVparams.Fs;

for i = 1:length(AfAnalysisWindows)
    tstart = AfAnalysisWindows(i);
    if isnan(tstart) 
        features_af = NaN;
        AFtest(i) = NaN;
    else
        idx_af = find(tNN >= tstart & tNN < tstart + HRVparams.af.windowlength);
        % When the RR interval time series is < 12 or >60 cannot extract feautures
        if (length(NN_afcalc(idx_af)) < 12 || length(NN_afcalc(idx_af)) > 60 )
            features_af = NaN;
            AFtest(i) = 0;
        else
            features_af = AF_features(NN_afcalc(idx_af),HRVparams.Fs);
            AFtest(i) = SVM_AFdetection_withoutTrainingModel(features_af,1);
        end
    end    
end


% 2. Export AF Data as CSV File

afresults = AFtest(:);
afcol_titles = {'AFtest'};

outputType = 'AF';
AFfile = SaveHRVoutput(subjectID,AfAnalysisWindows,afresults, ...
    afcol_titles, outputType, HRVparams, tNN, NN);
