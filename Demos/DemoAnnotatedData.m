%	OVERVIEW:
%       Demo using annotated data 
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
%       Main script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
HRVparams = InitializeHRVparams('demo');   % Initialize settings for demo



% Check existence of Input\Output data folders and add to search path

if  isempty(HRVparams.readdata) || ~exist([pwd filesep HRVparams.readdata], 'dir')    
    error('Invalid data INPUT folder');    % If folder name is empty
end
addpath(HRVparams.readdata)


addpath(HRVparams.readdata)
HRVparams.writedata = [HRVparams.writedata filesep 'Annotated'];
if ~exist(HRVparams.writedata, 'dir')
   mkdir(HRVparams.writedata)
end
addpath(HRVparams.writedata)


% Check for a list of files to be analyzed in current directory
% in .mat format

[subjectIDs,filesTBA] = GenerateListOfFilesTBA(HRVparams.ext, HRVparams.readdata,[]);

% Prepare for parallel loop by eliminating variables
clear nummatchingfiles x i filename flag match
numsub = length(subjectIDs);

% NOTE: This loop can be run in parallel by changing the loop to a parfor
% loop.

for i_patient = 1:numsub   
    try
    
        %% 1. Import Patient Data
        % HRVparams = initialize_HRVparams('demo'); % enable this when using parfor loops
        RRwindowStartIndices = [];
        tNN = [];
        NN = [];
        [samples, t, rr, ~, ~, annotations, HRVparams.Fs, ~, ~] = read_qrs(filesTBA{i_patient},HRVparams.datatype);
        %%% QUESTION: Does read_qrs accept filename on all platforms, or just
        %%%     LINUX?
        %%% TO DO - Improve functionality of specifying data to read function

        %% 2. Preprocess RR Data - Using HRV Toolbox
        % Remove noise, Remove ectopy, Don't detrend (yet)
        [NN, tNN, fbeats] = RRIntervalPreprocess(rr,t,annotations, HRVparams);

        %% 3. Calculate Windows
        RRwindowStartIndices = CreateWindowRRintervals(tNN, NN, HRVparams);   

        %% 4. Calculate AF Features
        
        [AFtest, AFwindowsStartIndices] = PerformAFdetection(subjectIDs{i_patient},tNN,NN,HRVparams);
        RRwindowStartIndices = RemoveAFsegments(RRwindowStartIndices,AFwindowsStartIndices, AFtest,HRVparams);
        
        %% 5. Calculate time domain HRV metrics - Using HRV Toolbox
        [NNmean,NNmedian,NNmode,NNvariance,NNskew,NNkurt, SDNN, NNiqr, ...
            RMSSD,pnn50,btsdet,avgsqi,fbeatw, RRwindowStartIndices] = EvalTimeDomainHRVstats(NN,tNN,[],HRVparams,RRwindowStartIndices);

        %% 6. Frequency domain HRV metrics (LF HF TotPow)
        %       All Inputs in Seconds
        %%% TO DO: Remove necessity of creating phantom beats with lomb 

        [ulf, vlf, lf, hf, lfhf, ttlpwr, methods, fdflag, window] = ...
            EvalFrequencyDomainHRVstats(NN,tNN, [],HRVparams,RRwindowStartIndices);

        %% 7. PRSA
        try
            [ac,dc,~] = prsa(NN, tNN, [], RRwindowStartIndices, HRVparams);
        catch
            ac = NaN; 
            dc = NaN;
            error_flag(i_patient) = subjectids(i_patient);
        end

        %% 8. SDANN and SDNNi
        [SDANN, SDNNI] = ClalcSDANN(RRwindowStartIndices, tNN, NN(:),HRVparams);

        %% 9. Compare to YALE & Joe Mietus's HRV Toolbox
        % J M Toolbox Results
        % (Requires a lot of processing time because of WFDB functions)
        % jm = compareJMHRV(NN,tNN,samples,annotations,num2str(subjectids(i_patient)),s);
        % 
        % yale = compareYale(subjectids(i_patient),t_win,s);

        % % Plots of Comparisons
        % % (Yale data was scaled on account of input being in ms rather than s)
        % figure; plot(t_win(1:200),vlf(1:200),t_win(1:200),jm.JMvlf,t_win(1:200),yale_mod_vlf(1:200)./1000000); legend('CL HRV TBx','JM HRV TBx','Yale'); title('VLF')
        % figure; plot(t_win(1:200),NNmean(1:200),t_win(1:200),jm.JMNNmean+.001,t_win(1:200),yale_mod_rr(1:200)./1000); legend('CL HRV TBx','JM HRV TBx + .001 offset'); title('Mean RR')
        % %figure; plot(t_win(1:200),pnn50(1:200),t_win(1:200),jm.JMpnn50+.001,t_win(1:200),yale_mod_pnn50(1:200)); legend('CL HRV TBx','JM HRV TBx + .001 offset'); title('pnn50')
        % figure; plot(t_win(1:200),hf(1:200),t_win(1:200),jm.JMhf,t_win(1:200),yale_mod_hf(1:200)./1000000); legend('CL HRV TBx','JM HRV TBx','Yale'); title('HF')
        % figure; plot(t_win(1:200),lf(1:200),t_win(1:200),jm.JMlf,t_win(1:200),yale_mod_lf(1:200)./1000000); legend('CL HRV TBx','JM HRV TBx','Yale'); title('LF')
        % figure; plot(t_win(1:200),RMSSD(1:200),t_win(1:200),jm.JMrmssd,t_win(1:200),yale_mod_rmssd(1:200)./1000); legend('CL HRV TBx','JM HRV TBx','Yale'); title('RMSSD')
        % figure; plot(t_win(1:200),SDNN(1:200),t_win(1:200),jm.JMSDNN+.0001,t_win(1:200),yale_mod_sdnn(1:200)./1000); legend('CL HRV TBx','JM HRV TBx + .0001 offset','Yale'); title('SDNN')
        % figure; plot(t_win(1:200),ttlpwr(1:200),t_win(1:200),jm.JMttlpwr,t_win(1:200), yale_mod_tp(1:200)./1000000); legend('CL HRV TBx','JM HRV TBx','Yale'); title('Total Power')
        % figure; plot(t_win(1:200),lfhf(1:200),t_win(1:200),jm.JMlfhf+.0001,t_win(1:200),yale_mod_lfhf(1:200)); legend('CL HRV TBx','JM HRV TBx + .0001 offset','Yale'); title('LF / HF Ratio')
  
        %% 10. Save Results
        % Uncomment the following lines for All Results
        results = [RRwindowStartIndices(:),ac(:),dc(:),...
            ulf(:),vlf(:),lf(:),hf(:),lfhf(:),ttlpwr(:),fdflag(:),...
            NNmean(:),NNmedian(:),NNmode(:),...
            NNvariance(:),NNskew(:),NNkurt(:),SDNN(:),NNiqr(:),RMSSD(:),pnn50(:),btsdet(:),fbeatw(:)];
        
        col_titles = {'t_win','ac','dc','ulf','vlf','lf','hf','lfhf',...
            'ttlpwr','fdflag','NNmean','NNmedian','NNmode','NNvar','NNskew',...
            'NNkurt','SDNN','NNiqr','RMSSD','pnn50','beatsdetected','corrected_beats'};

        % Uncomment the following lines for just Mean
        % results = [NNmean(:)];
        % col_titles = {'NN Mean'};

        % Generates Output - Never comment out
        resFilenameHRV = GenerateHRVresultsOutput(subjectIDs(i_patient), ...
            RRwindowStartIndices,results,col_titles, [], HRVparams, tNN, NN);
              
        % Compare generated output file with the reference one
        
        currentFile = [HRVparams.writedata filesep resFilenameHRV '.csv'];
        referenceFile = ['ReferenceOutput' filesep 'Annotated_HRV_allwindows.csv'];
        testHRV = CompareOutput(currentFile,referenceFile);
        
        if testHRV
            fprintf('** DemoAnnotatedData: TEST SUCCEEDED ** \n ')
             fprintf('A file named %s.csv \n has been saved in %s \n', ...
            resFilenameHRV, HRVparams.writedata);
        else
            fprintf('** DemoAnnotatedData: TEST FAILED ** \n')
        end
        
    catch
        if isnumeric(subjectIDs(i_patient))
            current_filename = ['error_' num2str(subjectIDs(i_patient))];
        else
            current_filename = ['error_' subjectIDs(i_patient)];
        end
            
        results = NaN;
        col_titles = {'NaN'};
        fprintf('Error on subject %s \n', char(subjectIDs(i_patient)));

        fprintf('** DemoAnnotatedData: TEST FAILED ** \n')
    end
    
end



