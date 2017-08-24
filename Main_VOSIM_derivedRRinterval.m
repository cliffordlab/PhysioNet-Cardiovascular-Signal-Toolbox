%  ====================== VOSIM Toolbox Main Script ======================
%
%	OVERVIEW:
%       Main "Validated Open-Source Integrated Matlab" VOSIM Toolbox script
%       It can be used  for the anlysis of MIT-B
%
%   INPUT:
%       Configured to accept RR intervals as well as raw data as input file
%       
%       NOTE: before running this script review and modifiy the parameters
%             in "initialize_HRVparams.m" file accordingly with the specific
%             of the new project (see the readme.txt file for further details)       
%   OUTPUT:
%       HRV Metrics 
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Main script written by Giulia Da Poian
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%
clear;
clc;
%% 1. Initizialization file sets up glopal variables and parameters 
%    (i.e., thresholds, window setitngs, noise limits and specral analysis)
HRVparams = InitializeHRVparams('mitarr');

%% 2. Look for all the current project files to be processed 
% Crate a list of all file to be analyze (% filesTBA: files to be analyzed) 
% in the study depending on the extension and on the input data format
% All files in the specified HRVparams.readdata folder will be used in the 
% analysis

[subjectIDs,filesTBA] = GenerateListOfFilesTBA(HRVparams.ext, HRVparams.readdata,[]);
% NOTE: all files with extension HRVparams.ext within the readdata folder 
%       will be analyzed


%% 3. Analyze each file and extract HRV parameters

for file_idx = 1:length(filesTBA)  
    
    try 
        %% 3.1 Load data
        switch HRVparams.input_data_format
            % Raw data input --------------------------------------------------
            case "input_waveform"
                %%% Import raw ECG signal (modify this import based on your format)
                [ECG_time,ECG_RawData] = rdsamp(filesTBA{file_idx},1);
                if size(ECG_time,1) > size(ECG_time,2)  
                    ECG_time =ECG_time';
                    ECG_RawData = ECG_RawData';
                end
               
                [t, rr] = ConvertRawDataToRRIntervals(ECG_RawData, HRVparams,subjectIDs{file_idx});

                % RR intervals input ----------------------------------------------
            case "input_RR_intervals"
                %%% If a *.qrs file is provides otherwise modify the
                %%% following line accordingly with tour RR input
                %%% Must provide both RR seriers and time of each RR
                %%% interval
                [samples, t, rr, ~, ~, annotations, HRVparams.Fs, ~, ~] = read_qrs(filesTBA{file_idx},HRVparams.datatype);
        end

        %% 3.2 Exlude undesiderable data from RR series  
        %      (i.e., arrhytmia, low SQI, ectopy, artefact, noise)

        [NN, tNN, fbeats] = RRIntervalPreprocess(rr,t,[], [], HRVparams);
        RRAnalysisWindows = CreateWindowRRintervals(tNN, NN, HRVparams);

%         %% 3.2 AF and ectopic beats (PVC) detection
%         try
%         [AFtest, AfAnalysisWindows] = PerformAFdetection(subjectIDs{file_idx},tNN,NN,HRVparams);
%        
%         % Exclude AF Segments
% 	    idx_afsegs = (find(AFtest == 1));
%         if ~isempty(idx_afsegs)
%             afsegs = AfAnalysisWindows(idx_afsegs);	% afsegs is in seconds
%             for k = 1:length(afsegs)
%                 try
%                     idx_af(k) = find(RRAnalysisWindows <= afsegs(k) & HRVparams.increment + RRAnalysisWindows > afsegs(k));
%                 catch
%                 end
%             end
%             afsegs_winall = RRAnalysisWindows(idx_af);
%             RRAnalysisWindows(idx_af) = NaN;
%         end
%         catch
%         end

        %% 3.4 Calculate time domain HRV metrics - Using VOSIM Toolbox Functions        
        
        [NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt, SDNN, NNiqr, ...
            RMSSD,pnn50,btsdet,avgsqi,fbeatw, windows_all] = ...
            EvalTimeDomainHRVstats(NN,tNN,[],HRVparams,RRAnalysisWindows,fbeats);

        %% 3.5 Frequency domain  metrics (LF HF TotPow) - Using VOSIM Toolbox Functions

        [ulf, vlf, lf, hf, lfhf, ttlpwr, methods, fdflag, window] = ...
           EvalFrequencyDomainHRVstats(NN,tNN, [],HRVparams,RRAnalysisWindows);

        %% 3.6 PRSA
        try
            [ac,dc,~] = prsa(NN, tNN, [], RRAnalysisWindows, HRVparams);
        catch
            ac = NaN; 
            dc = NaN;
        end

        %% 3.7 SDANN and SDNNi
        [SDANN, SDNNI] = ClalcSDANN(RRAnalysisWindows, tNN, NN(:),HRVparams); 

        %% 3.8 Export HRV Metrics as CSV File
        %Uncomment the following lines for All Results
        results = [RRAnalysisWindows(:), ac(:),dc(:),ulf(:),vlf(:),lf(:),hf(:), ...
                   lfhf(:),ttlpwr(:),fdflag(:), NNmean(:),NNmedian(:), ...
                   NNmode(:),NNvariance(:),NNskew(:),NNkurt(:),SDNN(:),...
                   NNiqr(:),RMSSD(:),pnn50(:),btsdet(:),fbeatw(:)];

        col_titles = {'t_win','ac','dc','ulf','vlf','lf','hf',...
                      'lfhf','ttlpwr','fdflag','NNmean','NNmedian',...
                      'NNmode','NNvar','NNskew','NNkurt','SDNN',...
                      'NNiqr','RMSSD','pnn50','beatsdetected','corrected_beats'};
        
        % Save results
        ResultsFileName = GenerateHRVresultsOutput(subjectIDs(file_idx),RRAnalysisWindows,results,col_titles, [],HRVparams, tNN, NN);

  
    catch
        
        
        
    end
end % end of for loop over TBAfiles

fprintf('Alanlysis Completed. Results saved in: ')






