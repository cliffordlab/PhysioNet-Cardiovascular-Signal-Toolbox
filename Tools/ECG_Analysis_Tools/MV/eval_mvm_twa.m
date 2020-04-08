function [MVMResult,TWAResult] = eval_mvm_twa(ecg,ann,fs)
%OVERVIEW This function evaluates MVM and TWA for ecg data sampled at 1000
% Hz and the following list of annotations, qrs onset, q, r, s, qrs offsett,
% t offset.
%
%   INPUTS          MANDATORY       DESCRIPTION
%                   ecg             N by M array of ecg data. Contains M
%                                   channels of ecg with N data points.
%
%                   fs              Sampling frequency for ecg, needs to 1000 Hz.
%   
%                   OPTIONAL
%                   ann             Structure containing the following
%                                   annotations, qrs onset, q, r, s, qrs
%                                   offset, t offset. If annotations are
%                                   not available pass an empty variable
%                                   ([]) as input for ann. Annotations will
%                                   be generated.
%
%   OUTPUTS
%                   MVMResult                       Structure containing the results of MVM
%                                                   analysis. The fields in the structure
%                                                   are listed below.
%
%                   Structure.Field                 
%                   MVMResult.energyinband_array    Array containing the energy in the every 2-7
%                                                   beat interval computed in the beatquency domain, for each analysis window.
%
%                   MVMResult.sqi_array             signal quality index
%                                                   for each analysis window.
%
%                   MVMResult.heart_rate_est_arr    average heart rate
%                                                   estimate for each analysis window.
%
%                   TWAResult               Structure containing the
%                                           results of TWA analysis. The fields in the structure
%                                           are listed below.
%                   
%                   TWAResult.HR            Estimate of average heart rate for each analysis window.        
%                   TWAResult.VAlt          Estimate of TWA amplitude for
%                                           each analysis window.
%                   TWAResult.VAlt_Sig      Estimate of TWA amplitude for
%                                           each analysis window which are statistically
%                                           significant compared to the noise threshold.
%                   TWAResult.Noise_Median  The median of the gamma
%                                           distribution used to model the noise.
%                   TWAResult.Noise_95      The 95th percentile of the
%                                           gamma distribution used to model the noise.
%                   TWAResult.VAltPt        The location on the S-T offset
%                                           segment corresponding to the maximum difference between
%                                           the average even and average odd beat.
%                   TWAResult.successful    Flag if set to 1 means TWA
%                                           analysis is successful.
%
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Ismail Sadiq
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  
%

% compute minimum amplitude resolution
amplitude_res = diff(ecg);
amplitude_res = min(amplitude_res(amplitude_res > eps));
if (amplitude_res > 0.006) % if amplitude resolution > 6 uV give warning on accuracy of twa results.
    disp(['Warning: Amplitude resolution for input ecg is ' num2str(amplitude_res) '. This greater than the suggested minimum resolution of 6 uV. Results for MV analysis may not be reliable.']);
end

if (fs == 1000)
    
    % Add dependencies
    addpath(genpath('../../../../PhysioNet-Cardiovascular-Signal-Toolbox-master')); % Cardio vascular toolbox
    addpath(genpath('./Tools/')); % Helper functions for mvtoolbox
    
    % check if ann provided as input
    annfields = {'QRSon', 'Q', 'R', 'S', 'QRSoff', 'Toff'};
    if (~isempty(ann))
        ann_check = zeros(1,length(annfields));
        for ii = 1:length(ann_check)
            if (isfield(ann,annfields{1,ii}) && ~isempty(getfield(ann,annfields{1,ii})))
                ann_check(ii) = 1;
            end
        end
        
        % check if any field missing
        missingi = find(ann_check == 0);
        if (~isempty(missingi))
            anntxt = [];
            for ii = 1:length(missingi)
                anntxt = [anntxt annfields{1,missingi(ii)} ', '];
            end
            disp(['Warning: Missing annotations, ' anntxt(1:end-2)])
            % generate annotations if missing
            disp('Incomplete annotations provided as input. Generating missing annotations.')
            HRVparams = InitializeHRVparams('MVanalysis'); %ecg = double(ECG_signal(:,7))/max(double(ECG_signal(:,7)));%std(double(ECG_signal(:,7)));
            [qrs_pos,sign,en_thres] = jqrs(ecg,HRVparams);
            ECG_header.nsig = 1; ECG_header.freq = fs; ECG_header.nsamp = length(ecg);
            wavedet_config.setup.wavedet.QRS_detection_only = 0;
            % plot(ecg); hold on; scatter(qrs_pos, ecg(qrs_pos)); hold off;
            [ann,~,~] = wavedet_3D_detector(ecg, qrs_pos', ECG_header, wavedet_config );
            
            [detection_flg,ann] = improvfiducials(ann, fs, ecg);
            
            if(~detection_flg)
                disp('Warning: Unable to improve fiducial point detection.')
            end
            
        end
        
    else
        disp('Annotations not provided as input. Generating missing annotations.')
        HRVparams = InitializeHRVparams('MVanalysis'); %ecg = double(ECG_signal(:,7))/max(double(ECG_signal(:,7)));%std(double(ECG_signal(:,7)));
        [qrs_pos,sign,en_thres] = jqrs(ecg,HRVparams);
        ECG_header.nsig = 1; ECG_header.freq = fs; ECG_header.nsamp = length(ecg);
        wavedet_config.setup.wavedet.QRS_detection_only = 0;
        % plot(ecg); hold on; scatter(qrs_pos, ecg(qrs_pos)); hold off;
        [ann,~,~] = wavedet_3D_detector(ecg', qrs_pos', ECG_header, wavedet_config );
        
        [detection_flg,ann] = improvfiducials(ann, fs, ecg);
        
        if(~detection_flg)
            disp('Warning: Unable to improve fiducial point detection.')
        end
        
    end
    
    % figure(1); plot(ecg); hold on; % Check annotation if needed
    % stem(ann.QRSon, ecg(ann.QRSon)); stem(ann.Q, ecg(ann.Q));
    % stem(ann.R, ecg(ann.R)); stem(ann.S, ecg(ann.S));
    % stem(ann.QRSoff, ecg(ann.QRSoff)); stem(ann.Toff, ecg(ann.Toff)); hold off
    % legend('ECG','QRSon','Q','R','S','QRSoff','Toff')
    
    % Variable for storing mvm for qrs
    numofleads = size(ecg,2);
    energyinband_array = cell(1,numofleads);
    sqi_array = cell(1,numofleads);
    heart_rate_est_arr = cell(1,numofleads);
    % Compute Morph Var QRS
    disp(['Evaluating QRS MVM ...']);
    normalize = 1;
    segment_size = 5;
    for lead = 1:size(ecg,2)
        [energyinband,sqi,heart_rate_est] = ComputeMVM(ecg(:,lead),ann,fs,segment_size,normalize);
        energyinband_array(1,lead) = mat2cell(energyinband, size(energyinband,1), size(energyinband,2));
        sqi_array(1,lead) = mat2cell(sqi, size(sqi,1), size(sqi,2));
        heart_rate_est_arr(1,lead) = mat2cell(heart_rate_est, size(heart_rate_est,1), size(heart_rate_est,2));
    end
    MVMResult.energyinband_array = energyinband_array;
    MVMResult.sqi_array = sqi_array;
    MVMResult.heart_rate_est_arr = heart_rate_est_arr;
    
    % Twa analysis
    pause(2)
    disp('Evaluating TWA ...');
    TWAResult = ComputeTWA(ecg,ann,fs);
    
else
    disp(['Input ecg is sampled at ' num2str(fs) ' Hz, needs to be sampled at 1000 Hz. Unable to perform MV analysis.'])
    MVMResult = NaN; TWAResult = NaN;
end

end

