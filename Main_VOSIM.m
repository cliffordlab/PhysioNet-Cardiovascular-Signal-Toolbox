function [results, ResultsFileName ] = Main_VOSIM(InputSig,t,InputFormat,HRVparams,subjectID,annotations,varargin)
%  ====================== VOSIM Toolbox Main Script ======================
%
%   Main_VOSIM(InputSig,t,annotations,InputFormat,ProjectName,subjectID)
%	OVERVIEW:
%       Main "Validated Open-Source Integrated Matlab" VOSIM Toolbox script
%       Configured to accept RR intervals as well as raw data as input file
%
%   INPUT:
%       InputSig    - Vector containing RR intervals data or ECG waveform  
%       t           - Time indices of the rr interval data (seconds) or
%                     ECG time
%       InputFormat - String that specifiy if the input vector is: 
%                     'RRinetrvals' for RR interval data 
%                     'ECGWaveform' for ECG waveform
%       HRVparams   - struct of settings for hrv_toolbox analysis that can
%                     be obtained using InitializeHRVparams.m function 
%                     HRVparams = InitializeHRVparams();
%       subjectID   - (optional) string to identify current subject
%       annotations - (optional) annotations of the RR data at each point
%                     indicating the quality of the beat 
%   OPTIONAL INPUTS:
%       Use InputSig, Type pairs for additional signals such as ABP 
%       or PPG signal. The input signal must be a vector containing
%       signal waveform and the Type: 'ABP' and\or 'PPG'.
%       
%
%   OUTPUT 
%       results         - HRV time and frequency domain metrics as well
%                         as AC and DC, SDANN and SDNNi
%       ResultsFileName - Name of the file containing the results
%
%       NOTE: before running this script review and modifiy the parameters
%             in "initialize_HRVparams.m" file accordingly with the specific
%             of the new project (see the readme.txt file for further details)   
%   EXAMPLES
%       - rr interval input
%       Main_VOSIM(RR,t,'RRintervals',HRVparams)
%       - ECG wavefrom input
%       Main_VOSIM(ECGsig,t,'ECGWavefrom',HRVparams,'101')
%       - ECG waveform and also ABP and PPG waveforms
%       Main_VOSIM(ECGsig,t,'ECGWaveform',HRVparams,[],[], abpSig, 'ABP', ppgSig, 'PPG')
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

if nargin < 4
    error('Wrong number of input arguments')
end
if nargin < 5
    subjectID = '0000';
end
if nargin < 6
    annotations = [];
end

if length(varargin) == 1 || length(varargin) == 3
    error('Incomplete Signal-Type pair')
elseif length(varargin)  == 2
    extraSigType = varargin(2);
    extraSig = varargin{1};
elseif length(varargin)  == 4
    extraSigType = [varargin(2) varargin(4)];
    extraSig = [varargin{1} varargin{3}];
end

results = [];
col_titles = {};

if isa(subjectID,'cell'); subjectID = string(subjectID); end


% Start HRV analysis
try   
    if strcmp(InputFormat, 'ECGWaveform')
        % Convert ECG waveform in rr intervals
        [t, rr, jqrs_ann, SQIvalue , SQIidx] = ConvertRawDataToRRIntervals(InputSig, HRVparams, subjectID);
        sqi = [SQIidx', SQIvalue'];
        GenerateHRVresultsOutput(subjectID,[], sqi ,{'WinSQI','SQI'},'SQI',HRVparams,[],[]);  
    else
        rr = InputSig; 
        sqi = [];
    end

    % Exlude undesiderable data from RR series (i.e., arrhytmia, low SQI, ectopy, artefact, noise)

    [NN, tNN] = RRIntervalPreprocess(rr,t,annotations, HRVparams);
    RRwindowStartIndices = CreateWindowRRintervals(tNN, NN, HRVparams);
    
    % 1. Atrial Fibrillation Detection
    if HRVparams.af.on == 1
        [AFtest, AfAnalysisWindows] = PerformAFdetection(subjectID,t,rr,HRVparams);
        % Create RRAnalysisWindows contating AF segments
        [RRwindowStartIndices, AFWindows]= RemoveAFsegments(RRwindowStartIndices,AfAnalysisWindows, AFtest,HRVparams);
        fprintf('AF analysis completed for patient %s \n', subjectID);
    end
   
    % 2. Calculate time domain HRV metrics - Using VOSIM Toolbox Functions        
    if HRVparams.timedomain.on == 1
        [NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt, SDNN, NNiqr, ...
        RMSSD,pnn50,btsdet,fdflagTime] = EvalTimeDomainHRVstats(NN,tNN,sqi,HRVparams,RRwindowStartIndices);
        % Export results
        results = [ results, RRwindowStartIndices(:), NNmean(:),NNmedian(:),NNmode(:),NNvariance(:),...
                    NNskew(:),NNkurt(:),SDNN(:), NNiqr(:),RMSSD(:), pnn50(:),...
                    btsdet(:),fdflagTime(:)];
        col_titles = [col_titles {'t_win','NNmean','NNmedian','NNmode','NNvar',...
                      'NNskew','NNkurt','SDNN','NNiqr','RMSSD','pnn50',...
                      'beatsdetected','WinFlagTime'}];
    end
    
    % 3. Frequency domain  metrics (LF HF TotPow) - Using VOSIM Toolbox Functions
    if HRVparams.freq.on == 1
        [ulf, vlf, lf, hf, lfhf, ttlpwr, fdflagFreq] = ...
         EvalFrequencyDomainHRVstats(NN,tNN,sqi,HRVparams,RRwindowStartIndices);
         % Export results
         results = [results, ulf(:),vlf(:),lf(:),hf(:), lfhf(:),ttlpwr(:),fdflagFreq(:)];
          col_titles = [col_titles {'ulf','vlf','lf','hf', 'lfhf','ttlpwr','WinFlagFreq'}];
    end
    
    % 4. PRSA, AC and DC values
    if HRVparams.prsa.on == 1
        try
            [ac,dc,~] = prsa(NN, tNN, HRVparams, sqi, RRwindowStartIndices );
        catch
            ac = NaN; 
            dc = NaN;
        end
        % Export results
        results = [results, ac(:), dc(:)];
        col_titles = [col_titles {'ac' 'dc'}];
    end
    
    
    % Generates Output - Never comment out
    ResultsFileName = GenerateHRVresultsOutput(subjectID,RRwindowStartIndices,results,col_titles, [],HRVparams, tNN, NN);
    fprintf('HRV metrics for file ID %s saved in the output folder \n File name: %s \n', subjectID, ResultsFileName);

    % 5. Multiscale Entropy
    if HRVparams.MSE.on == 1
        try
            mse = ComputeMultiscaleEntropy(NN,HRVparams.MSE.MSEpatternLength, HRVparams.MSE.RadiusOfSimilarity, HRVparams.MSE.maxCoarseGrainings);  
        catch
            mse = NaN;
            fprintf('MSE failed for file ID %s \n', subjectID);
        end
         % Save Results for MSE
        results = mse;
        col_titles = {'MSE'};
        % Generates Output - Never comment out
        GenerateHRVresultsOutput(subjectID,[],results,col_titles, 'MSE', HRVparams, tNN, NN);
    end   
    
    % 6. Analyze additional signals (ABP, PPG or both)
    if ~isempty(varargin)
        fprintf('Analyizing %s \n', extraSigType{:});
        Analyze_ABP_PPG_Waveforms(extraSig,extraSigType,HRVparams,jqrs_ann,subjectID);
    end
        
    % 7. Some statistics on %ages windows removed (poor quality and AF)
    %    save on file  
    RemovedWindowsStats(RRwindowStartIndices,AFWindows,HRVparams,subjectID);
    
    fprintf('HRV Analysis completed for file ID %s \n',subjectID )
  
catch
    % Write subjectID on log file
    fid = fopen(strcat(HRVparams.writedata,filesep,'Error.txt'),'a');
    results = NaN;
    ResultsFileName = '';
    fprintf(fid, '%s \n', subjectID);
    fclose(fid); 
    fprintf('Analysis not performed for file ID %s \n', subjectID);
end % end of HRV analysis







end %== function ================================================================
%




