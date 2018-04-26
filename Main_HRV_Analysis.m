function [HRVout, ResultsFileName ] = Main_HRV_Analysis(InputSig,t,InputFormat,HRVparams,subID,ann,sqi,varargin)

%  ====== HRV Toolbox for PhysioNet Cardiovascular Signal Toolbox =========
%
%   Main_HRV_Analysis(InputSig,t,InputFormat,HRVparams,subID,ann,sqi,varargin)
%	OVERVIEW:
%       Main "Validated Open-Source Integrated Matlab" VOSIM Toolbox script
%       Configured to accept RR intervals as well as raw data as input file
%
%   INPUT:
%       InputSig    - Vector containing RR intervals data (in seconds) 
%                     or ECG/PPG waveform  
%       t           - Time indices of the rr interval data (seconds) or
%                     leave empty for ECG/PPG input
%       InputFormat - String that specifiy if the input vector is: 
%                     'RRIntervals' for RR interval data 
%                     'ECGWaveform' for ECG waveform
%                     'PPGWaveform' for PPG signal
%       HRVparams   - struct of settings for hrv_toolbox analysis that can
%                     be obtained using InitializeHRVparams.m function 
%                     HRVparams = InitializeHRVparams();
%
%      
%   OPTIONAL INPUTS:
%       subID       - (optional) string to identify current subject
%       ann         - (optional) annotations of the RR data at each point
%                     indicating the type of the beat 
%       sqi         - (optional) Signal Quality Index; Requires a 
%                     matrix with at least two columns. Column 1 
%                     should be timestamps of each sqi measure, and 
%                     Column 2 should be SQI on a scale from 0 to 1.
%       Use InputSig, Type pairs for additional signals such as ABP 
%       or PPG signal. The input signal must be a vector containing
%       signal waveform and the Type: 'ABP' and\or 'PPG'.
%       
%
%   OUTPUS:
%       results         - HRV time and frequency domain metrics as well
%                         as AC and DC, SDANN and SDNNi
%       ResultsFileName - Name of the file containing the results
%
%       NOTE: before running this script review and modifiy the parameters
%             in "initialize_HRVparams.m" file accordingly with the specific
%             of the new project (see the readme.txt file for further details)   
%
%   EXAMPLES
%       - rr interval input
%       Main_HRV_Analysis(RR,t,'RRIntervals',HRVparams)
%       - ECG wavefrom input
%       Main_HRV_Analysis(ECGsig,t,'ECGWavefrom',HRVparams,'101')
%       - ECG waveform and also ABP and PPG waveforms
%       Main_HRV_Analysis(ECGsig,t,'ECGWaveform',HRVparams,[],[],[], abpSig, 
%                         'ABP', ppgSig, 'PPG')
%
%   DEPENDENCIES & LIBRARIES:
%       HRV Toolbox for PhysioNet Cardiovascular Signal Toolbox
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   REFERENCE: 
%   Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%   Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       This script written by Giulia Da Poian
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if nargin < 4
    error('Wrong number of input arguments')
end
if nargin < 5
    subID = '0000';
end
if nargin < 6
    ann = [];
end
if nargin < 7
    sqi = [];
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

if isa(subID,'cell'); subID = string(subID); end


% Control on signal length

if (strcmp(InputFormat, 'ECGWaveform') && length(InputSig)/HRVparams.Fs< HRVparams.windowlength) ...
        || (strcmp(InputFormat, 'PPGWaveform') && length(InputSig)/HRVparams.Fs<HRVparams.windowlength) ...
        ||  (strcmp(InputFormat, 'RRIntervals') && t(end)/HRVparams.Fs<HRVparams.windowlength)
   error('Error: imput signal or RR interval time series is shorter than the desired analysis windows (%i s).', HRVparams.windowlength)
end



% Start HRV analysis
try   
    switch InputFormat
        case 'ECGWaveform'
            % Convert ECG waveform in rr intervals
            [t, rr, jqrs_ann, SQIvalue , SQIidx] = ConvertRawDataToRRIntervals(InputSig, HRVparams, subID);
            sqi = [SQIidx', SQIvalue'];            
        case 'PPGWaveform'
            [rr,t] = Analyze_ABP_PPG_Waveforms(InputSig,{'PPG'},HRVparams,[],subID);
        case 'RRIntervals'
            rr = InputSig; 
        otherwise
            error('Wrong Input Type! This function accepts: ECGWaveform, PPGWaveform or RRIntervals')           
    end

    % 1. Preprocess Data, AF detection, create Windows Indexes  
    [NN, tNN, WinIdxs, AFWindows,out] = PreparDataForHRVAnlysis(rr,t,ann,sqi,HRVparams,subID);
        
    HRVout = WinIdxs';
    HRVtitle = {'t_win'};
   
    % 3. Calculate time domain HRV metrics - Using VOSIM Toolbox Functions        
    if HRVparams.timedomain.on 
        TimeMetrics = EvalTimeDomainHRVstats(NN,tNN,sqi,HRVparams,WinIdxs);
        % Export results
        HRVout = [HRVout cell2mat(struct2cell(TimeMetrics))'];
        HRVtitle = [HRVtitle fieldnames(TimeMetrics)'];
    end
    
    % 4. Frequency domain  metrics (LF HF TotPow) - Using VOSIM Toolbox Functions
    if HRVparams.freq.on 
        FreqMetrics = EvalFrequencyDomainHRVstats(NN,tNN,sqi,HRVparams,WinIdxs);
        % Export results
        HRVout = [HRVout cell2mat(struct2cell(FreqMetrics))'];
        HRVtitle = [HRVtitle fieldnames(FreqMetrics)'];
    end
    
    % 5. PRSA, AC and DC values
    if HRVparams.prsa.on 
        [ac,dc,~] = prsa(NN, tNN, HRVparams, sqi, WinIdxs );
        % Export results
        HRVout = [HRVout, ac(:), dc(:)];
        HRVtitle = [HRVtitle {'ac' 'dc'}];
    end
    
    % 6.Poincare Features
    if HRVparams.poincare.on
         [SD1, SD2, SD12Ratio] = EvalPoincareOnWindows(NN, tNN, HRVparams, WinIdxs, sqi);
         % Export results
         HRVout = [HRVout, SD1(:),SD2(:),SD12Ratio(:)];
         HRVtitle = [HRVtitle {'SD1', 'SD2', 'SD1SD2'}];
    end
    
    % 7.Entropy Features
    if HRVparams.Entropy.on
        m = HRVparams.Entropy.patternLength;
        r = HRVparams.Entropy.RadiusOfSimilarity;
        [SampEn, ApEn] = EvalEntropyMetrics(NN, tNN, m ,r, HRVparams, WinIdxs, sqi);
        % Export results
        HRVout = [HRVout, SampEn(:),ApEn(:)];
        HRVtitle = [HRVtitle {'SampEn', 'ApEn'}];
    end
    
    % Generates Output - Never comment out
    ResultsFileName.HRV = SaveHRVoutput(subID,WinIdxs,HRVout,HRVtitle, [],HRVparams, tNN, NN);
    
    % 8. Multiscale Entropy (MSE)
    if HRVparams.MSE.on 
        try
            mse = EvalMSE(out.NN_gapFilled,out.tNN_gapFilled,sqi,HRVparams,out.WinIdxsMSE);
        catch
            mse = NaN;
            fid = fopen([HRVparams.writedata filesep 'AnalysisError.txt','a']);
            fprintf(fid, 'MSE analysis error for subject %s \n',subID );
            fclose(fid);
        end
         % Save Results for MSE
        Scales = 1:HRVparams.MSE.maxCoarseGrainings;
        HRVout = [Scales' mse];
        for i=1:length(out.WinIdxsMSE)
            Windows{i} = strcat('t_', num2str(WindIdxs(i)));
        end
        HRVtitle = {'Scales' Windows{:}'};
        ResultsFileName.MSE = SaveHRVoutput(subID,[],HRVout,HRVtitle, 'MSE', HRVparams, tNN, NN);
    end   

       
    % 9. DetrendedFluctuation Analysis (DFA)
    if HRVparams.DFA.on
        try
            [alpha1, alpha2] = EvalDFA(out.NN_gapFilled,out.tNN_gapFilled,sqi,HRVparams,out.WinIdxsDFA);   
            % Save Results for DFA
            HRVout = [out.WinIdxsDFA' alpha1 alpha2];
            HRVtitle = {'t_win' 'alpha1' 'alpha2'};
            ResultsFileName.DFA = SaveHRVoutput(subID,[],HRVout,HRVtitle, 'DFA', HRVparams, tNN, NN);
        catch
            fid = fopen([HRVparams.writedata filesep 'AnalysisError.txt'],'a');
            fprintf(fid, 'DFA analysis error for subject %s \n',subID );
            fclose(fid);
        end
    end
    
    % 10. Heart Rate Turbulence Analysis (HRT)
    if HRVparams.HRT.on
        try
            % Create analysis windows from original rr intervals
            WinIdxsHRT = CreateWindowRRintervals(t, rr, HRVparams,'HRT');
            [TO, TS, nPVCs] = Eval_HRT(rr,t,ann,sqi, HRVparams, WinIdxsHRT);
            % Save Results for DFA
            HRVout = [WinIdxsHRT' TO TS nPVCs];
            HRVtitle = {'t_win' 'TO' 'TS' 'nPVCs'};
            ResultsFileName.HRT = SaveHRVoutput(subID,[],HRVout,HRVtitle, 'HRT', HRVparams, t, rr);
        catch
            fid = fopen([HRVparams.writedata filesep 'AnalysisError.txt'],'a');
            fprintf(fid, 'HRT analysis error for subject %s \n',subID );
            fclose(fid);
        end
    end
    
    % 11. Analyze additional signals (ABP, PPG or both)
    if ~isempty(varargin)
        try
            fprintf('Analyizing %s \n', extraSigType{:});
            Analyze_ABP_PPG_Waveforms(extraSig,extraSigType,HRVparams,jqrs_ann,subID);
        catch
            fid = fopen([HRVparams.writedata filesep 'AnalysisError.txt'],'a');
            fprintf(fid, 'ABP/PPG analysis error for subject %s \n',subID );
            fclose(fid);
        end
    end
        
    % 12. Some statistics on %ages windows removed (poor quality and AF)
    %    save on file  
    RemovedWindowsStats(WinIdxs,AFWindows,HRVparams,subID);
    
    fprintf('HRV Analysis completed for subject ID %s \n',subID);
    
    fid = fopen([HRVparams.writedata filesep 'FileSuccessfullyAnalyzed.txt'],'a');
    fprintf(fid, '%s \n',subID );
    fclose(fid);
  
catch
    % Write subjectID on log file
    fid = fopen(strcat(HRVparams.writedata,filesep,'AnalysisError.txt'),'a');
    HRVout = NaN;
    ResultsFileName = '';
    fprintf(fid, 'Basic HRV Analysis faild for subject: %s \n', subID);
    fclose(fid); 
    %fprintf('Analysis not performed for file ID %s \n', subjectID);
end % end of HRV analysis


end %== function ================================================================
%

function [NN, tNN, WinIdxs,AFWindows,out] = PreparDataForHRVAnlysis(rr,t,annotations,sqi,HRVparams,subjectID)

    out = []; % Struct used to save DFA and MSE preprocessed data
 
    % Exclude undesiderable data from RR series (i.e., arrhytmia, low SQI, ectopy, artefact, noise)
    [NN, tNN] = RRIntervalPreprocess(rr,t,annotations, HRVparams);  
    WinIdxs = CreateWindowRRintervals(tNN, NN, HRVparams);    % Create Windows for Time and Frequency domain 
    
    % Create Windows for MSE and DFA and preprocess
    if HRVparams.MSE.on || HRVparams.DFA.on
       % Additional pre-processing to deal with missing data for MSE and DFA analysis     
       [out.NN_gapFilled, out.tNN_gapFilled] = RR_Preprocessing_for_MSE_DFA( NN, tNN );
    end
    if HRVparams.MSE.on
       out.WinIdxsMSE = CreateWindowRRintervals(out.tNN_gapFilled, out.NN_gapFilled, HRVparams,'mse');
    end
    if HRVparams.DFA.on
        out.WinIdxsDFA = CreateWindowRRintervals(out.tNN_gapFilled, out.NN_gapFilled, HRVparams,'dfa');
    end    
    
    % 2. Atrial Fibrillation Detection
    if HRVparams.af.on 
        [AFtest, AfAnalysisWindows] = PerformAFdetection(subjectID,t,rr,sqi,HRVparams);
        fprintf('AF analysis completed for subject %s \n', subjectID);
        % Remove RRAnalysisWindows contating AF segments
        [WinIdxs, AFWindows]= RemoveAFsegments(WinIdxs,AfAnalysisWindows, AFtest,HRVparams);
        if HRVparams.MSE.on
            out.WinIdxsMSE = RemoveAFsegments(out.WinIdxsMSE,AfAnalysisWindows, AFtest,HRVparams);
        end
        if HRVparams.DFA.on 
            out.WinIdxsDFA = RemoveAFsegments(out.WinIdxsDFA,AfAnalysisWindows, AFtest,HRVparams);
        end
    else
        AFWindows = [];
    end
    
end

