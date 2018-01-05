function [results, ResultsFileName ] = Main_VOSIM(InputSig,t,InputFormat,HRVparams,subjectID,annotations,sqi,varargin)
%  ====================== VOSIM Toolbox Main Script ======================
%
%   Main_VOSIM(InputSig,t,annotations,InputFormat,ProjectName,subjectID)
%	OVERVIEW:
%       Main "Validated Open-Source Integrated Matlab" VOSIM Toolbox script
%       Configured to accept RR intervals as well as raw data as input file
%
%   INPUT:
%       InputSig    - Vector containing RR intervals data (in seconds) 
%                     or ECG/PPG waveform  
%       t           - Time indices of the rr interval data (seconds) or
%                     live empty for ECG/PPG input
%       InputFormat - String that specifiy if the input vector is: 
%                     'RRIntervals' for RR interval data 
%                     'ECGWaveform' for ECG waveform
%                     'PPGWaveform' for PPG signal
%       HRVparams   - struct of settings for hrv_toolbox analysis that can
%                     be obtained using InitializeHRVparams.m function 
%                     HRVparams = InitializeHRVparams();
%       subjectID   - (optional) string to identify current subject
%       sqi         - (optional) Signal Quality Index; Requires a 
%                     matrix with at least two columns. Column 1 
%                     should be timestamps of each sqi measure, and 
%                     Column 2 should be SQI on a scale from 0 to 1.
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
%       Main_VOSIM(RR,t,'RRIntervals',HRVparams)
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

if isa(subjectID,'cell'); subjectID = string(subjectID); end


% Start HRV analysis
try   
    switch InputFormat
        case 'ECGWaveform'
            % Convert ECG waveform in rr intervals
            [t, rr, jqrs_ann, SQIvalue , SQIidx] = ConvertRawDataToRRIntervals(InputSig, HRVparams, subjectID);
            sqi = [SQIidx', SQIvalue'];            
        case 'PPGWaveform'
            [rr,t] = Analyze_ABP_PPG_Waveforms(InputSig,{'PPG'},HRVparams,[],subjectID);
        case 'RRIntervals'
            rr = InputSig; 
        otherwise
            error('Wrong Input Type! This function accepts: ECGWaveform, PPGWaveform or RRIntervals')           
    end

    
    % 1. Exclude undesiderable data from RR series (i.e., arrhytmia, low SQI, ectopy, artefact, noise)
    [NN, tNN] = RRIntervalPreprocess(rr,t,annotations, HRVparams);
    % Create Windows for Time and Frequency domain 
    RRwindowStartIndices = CreateWindowRRintervals(tNN, NN, HRVparams);
   
    % Additional pre-processing to deal with missing data for MSE and DFA analysis     
    [ NN_gapFilled, tNN_gapFilled] = RR_Preprocessing_for_MSE( NN, tNN );
    % Create Windows for MSE and DFA 
    if HRVparams.MSE.on
        MSEwindowStartIndices = CreateWindowRRintervals(tNN_gapFilled, NN_gapFilled, HRVparams,'mse');
    end
    if HRVparams.DFA.on 
        DFAwindowStartIndices = CreateWindowRRintervals(tNN_gapFilled, NN_gapFilled, HRVparams,'dfa');
    end
    
    % 2. Atrial Fibrillation Detection
    if HRVparams.af.on == 1
        [AFtest, AfAnalysisWindows] = PerformAFdetection(subjectID,t,rr,HRVparams);
        fprintf('AF analysis completed for patient %s \n', subjectID);
        % Remove RRAnalysisWindows contating AF segments
        [RRwindowStartIndices, AFWindows]= RemoveAFsegments(RRwindowStartIndices,AfAnalysisWindows, AFtest,HRVparams);
            if HRVparams.MSE.on
                MSEwindowStartIndices = RemoveAFsegments(MSEwindowStartIndices,AfAnalysisWindows, AFtest,HRVparams);
            end
            if HRVparams.DFA.on 
                DFAwindowStartIndices = RemoveAFsegments(DFAwindowStartIndices,AfAnalysisWindows, AFtest,HRVparams);
            end
    else
        AFWindows = [];
    end
    
    results = RRwindowStartIndices';
    col_titles = {'t_win'};
   
    % 3. Calculate time domain HRV metrics - Using VOSIM Toolbox Functions        
    if HRVparams.timedomain.on == 1
        TimeMetrics = EvalTimeDomainHRVstats(NN,tNN,sqi,HRVparams,RRwindowStartIndices);
        % Export results
        results = [ results cell2mat(struct2cell(TimeMetrics))'];
        col_titles = [col_titles fieldnames(TimeMetrics)'];
    end
    
    % 4. Frequency domain  metrics (LF HF TotPow) - Using VOSIM Toolbox Functions
    if HRVparams.freq.on == 1
        FreqMetrics = EvalFrequencyDomainHRVstats(NN,tNN,sqi,HRVparams,RRwindowStartIndices);
        % Export results
        results = [results cell2mat(struct2cell(FreqMetrics))'];
        col_titles = [col_titles fieldnames(FreqMetrics)'];
    end
    
    % 5. PRSA, AC and DC values
    if HRVparams.prsa.on == 1
        [ac,dc,~] = prsa(NN, tNN, HRVparams, sqi, RRwindowStartIndices );
        % Export results
        results = [results, ac(:), dc(:)];
        col_titles = [col_titles {'ac' 'dc'}];
    end
    
    % 6.Poincare Features
    if HRVparams.poincare.on==1
         [SD1, SD2, SD1_SD2_ratio] = EvalPoincareOnWindows(NN, tNN, HRVparams, RRwindowStartIndices, sqi);
         % Export results
         results = [results, SD1(:),SD2(:),SD1_SD2_ratio(:)];
         col_titles = [col_titles {'SD1', 'SD2', 'SD1SD2'}];
    end
    
    % 7.Entropy Features
    if HRVparams.Entropy.on==1
        m = HRVparams.Entropy.patternLength;
        r = HRVparams.Entropy.RadiusOfSimilarity;
        [SampEn, ApEn] = EvalEntropyMetrics(NN, tNN, m,r, HRVparams, RRwindowStartIndices, sqi);
        % Export results
        results = [results, SampEn(:),ApEn(:)];
        col_titles = [col_titles {'SampEn', 'ApEn'}];
    end
    
    % Generates Output - Never comment out
    ResultsFileName = GenerateHRVresultsOutput(subjectID,RRwindowStartIndices,results,col_titles, [],HRVparams, tNN, NN);
    fprintf('HRV metrics for file ID %s saved in the output folder \n File name: %s \n', subjectID, ResultsFileName);


    
    % 8. Multiscale Entropy (MSE)
    if HRVparams.MSE.on == 1
        try
            mse = EvalMSE(NN_gapFilled,tNN_gapFilled,sqi,HRVparams,MSEwindowStartIndices);
        catch
            mse = NaN;
            fprintf('MSE failed for file ID %s \n', subjectID);
        end
         % Save Results for MSE
        results = [MSEwindowStartIndices' mse'];
        for i=1:HRVparams.MSE.maxCoarseGrainings
            Scales{i}=strcat('Scale', num2str(i));
        end
        col_titles = {'t_win' Scales{:}};
        GenerateHRVresultsOutput(subjectID,[],results,col_titles, 'MSE', HRVparams, tNN, NN);
    end   
    
    
    % 9. DetrendedFluctuation Analysis (DFA)
    if HRVparams.DFA.on == 1
        alpha = EvalDFA(NN_gapFilled,tNN_gapFilled,sqi,HRVparams,DFAwindowStartIndices);   
        % Save Results for DFA
        results = [DFAwindowStartIndices' alpha];
        col_titles = {'t_win' 'alpha'};
        GenerateHRVresultsOutput(subjectID,[],results,col_titles, 'DFA', HRVparams, tNN, NN);
    end
    
    % 10. Analyze additional signals (ABP, PPG or both)
    if ~isempty(varargin)
        fprintf('Analyizing %s \n', extraSigType{:});
        Analyze_ABP_PPG_Waveforms(extraSig,extraSigType,HRVparams,jqrs_ann,subjectID);
    end
        
    % 11. Some statistics on %ages windows removed (poor quality and AF)
    %    save on file  
    RemovedWindowsStats(RRwindowStartIndices,AFWindows,HRVparams,subjectID);
    
    fprintf('HRV Analysis completed for file ID %s \n',subjectID )
  
catch
    % Write subjectID on log file
    fid = fopen(strcat(HRVparams.writedata,filesep,'Error.txt'),'a');
    results = NaN;
    ResultsFileName = '';
    fprintf(fid, 'Analysis faild for subject: %s \n', subjectID);
    fclose(fid); 
    fprintf('Analysis not performed for file ID %s \n', subjectID);
end % end of HRV analysis







end %== function ================================================================
%




