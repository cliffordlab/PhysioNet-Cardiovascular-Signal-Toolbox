function [ gentwaECG_n,gentwaECG ] = gen_twa_ecg(twa_amp,fs,duration,hr,rr, SNR)
%OVERVIEW, This function generates synthetic ecg with TWAs.
%
% INPUTS        MANDATORY               DESCRIPTION
%               twa_amp (uV)            The T wave alternan ampltiude in uV.
%
%               fs (Hz)                 Sampling frequency of the synthetic ecg
%                                       which will be generated. The script
%                                       has been tested for fs = 1000
%                                       Hz. For other fs values the synthetic 
%                                       ecg maybe resampled. 
%
%               duration (min.)         The duration of the ecg to be
%                                       sythesized in minutes.
%
%               hr (bpm)                The heart rate (HR) in beats per minute
%                                       (bpm).
%
%               rr (rpm)                The respiration rate (RR) in
%                                       respirations per minute (rpm).
%
%               SNR (dB)                Signal to noise ratio for sythetic
%                                       ecg signal.
%
% OUTPUTS
%               gentwaECG_n (mV)        Synthetic ecg with additive noise.
%                                       The TWA amplitude is specified by
%                                       the twa_amp variable.
%
%               gentwaECG (mV)          Noise free synthetic ecg. The TWA amplitude 
%                                       is specified by the twa_amp variable. 
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

lead = 1; % ecg lead for TWA fixed to 1 
% if no twa generate normal ECG
[signal] = generate_resp_modulated_ecg(fs, duration, rr, hr, 1, SNR);

% if twa_amp > 0, generate second set of ECG with elevated T wave
% amplitude.
if (twa_amp ~= 0)
    if (twa_amp < 5 || twa_amp > 100)
        disp('Only ECG with TWA amplitude between 5 and 100 uV can be generated.')
    else
        if (twa_amp >= 5 && twa_amp < 5.4)
            twa_amp = 5.4;
        end
        
        if (twa_amp > 99.7 && twa_amp <= 100)
            twa_amp = 99.7;
        end
        
        lambdascalar = [1.007 1.014 1.021 1.027 1.034 1.040 1.046 1.052 1.0585 1.0655 1.0725 1.079 1.086 1.093 1.099 1.106 1.112 1.118 1.124 1.130];
        twa_amp_array = [5.4 10.6 16.1 20.6 26 30.5 35.1 39.7 44.6 50.1 55.4 60.5 65.26 71.2 75.9 81.3 85.9 90.5 94.9 99.7];
        % try to find exact lambda
        lambda = lambdascalar(find(twa_amp_array == twa_amp));
        
        % interpolate if needed
        if (isempty(lambda))
            % upper and lower bound of twa amplitude
            upper = find(twa_amp_array > twa_amp, 1);
            lower = find(twa_amp_array <= twa_amp,1, 'last');
            % interpolate
            lambda = lambdascalar(lower) + (lambdascalar(upper) - lambdascalar(lower))/(twa_amp_array(upper) - twa_amp_array(lower)) * (twa_amp - twa_amp_array(lower));
        end
        % generate high amplitude t wave signal
        [signal0] = generate_resp_modulated_ecg(fs, duration, rr, hr, lambda,SNR);
    end
    
    % merge low and high t wave amplitude series to get TWA ecg
    % get r peaks for each ECG, use lead 1 for TWA generation
    [hrv, R_t, R_amp, R_index, S_t, S_amp]  = rpeakdetect(signal(:,lead),1000);
    [hrv0, R_t0, R_amp0, R_index0, S_t0, S_amp0]  = rpeakdetect(signal0(:,lead),1000);
    
    R_original = round(R_t*1000);
    R_twa = round(R_t0*1000);
    % beat original
    startidx = round((R_original(1)+R_original(2))/2);
    endidx = round((R_original(2)+R_original(3))/2);
    gentwaECG = signal(startidx:endidx,lead)';
    for beatidx = 3:(length(R_original)-1)
        if (mod(beatidx,2)) % even beat, twa
            startidx = round((R_twa(beatidx-1)+R_twa(beatidx))/2);
            endidx = round((R_twa(beatidx)+R_twa(beatidx+1))/2);
            temporarybeat = signal0(startidx:endidx,lead);
            gentwaECG = [gentwaECG temporarybeat'];
        else % odd beat, normal
            startidx = round((R_original(beatidx-1)+R_original(beatidx))/2);
            endidx = round((R_original(beatidx)+R_original(beatidx+1))/2);
            temporarybeat = signal(startidx:endidx,lead);
            gentwaECG = [gentwaECG temporarybeat'];
        end
    end
end

% add noise
[ gentwaECG_n ] = AddNoise(gentwaECG',SNR,fs);% Add noise
    
%figure(1); plot(gentwaECG); hold on; plot(gentwaECG_n); hold off; % plot
%clean and noisy ecg with TWA
