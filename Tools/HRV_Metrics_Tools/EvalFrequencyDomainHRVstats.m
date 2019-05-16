function out = EvalFrequencyDomainHRVstats(NN, tNN, sqi, HRVparams, tWin)
%
% out = EvalFrequencyDomainHRVstats (NN, tNN, , sqi, settings)
%   
%   OVERVIEW:   This function returns frequency domain HRV metrics 
%               calculated on input NN intervals.
%
%   INPUT:      MANDATORY:
%               NN          : a single row of NN (normal normal) interval
%                             data in seconds
%               tNN         : a single row of time of the rr interval 
%                             data (seconds)
%               sqi         : (Optional) Signal Quality Index; Requires a 
%                             matrix with at least two columns. Column 1 
%                             should be timestamps of each sqi measure, and 
%                             Column 2 should be SQI on a scale from 0 to 1.
%                             Additional columns can be included with
%                             additional sqi at the same timestamps
%               HRVparams   : struct of settings for hrv_toolbox analysis
%               tWin .      : vector containing the starting time of each
%                             windows (in seconds) 
%
%   OUTPUT:     out.ulf        : (ms^2) Power in the ultra low frequency range (default < 0.003 Hz)
%               out.vlf        : (ms^2) Power in very low frequency range (default 0.003 <= vlf < 0.04 Hz)
%               out.lf         : (ms^2) Power in low frequency range (default 0.04Hz  <= lf < 0.15 Hz)
%               out.hf         : (ms^2) Power in high frequency range (default 0.15 <= hf < 0.4 Hz)
%               out.lfhf       : Ratio LF [ms^2]/HF [ms^2]
%               out.ttlpwr     : (ms^2) Total spectral power (approximately <0.4 Hz)
%               out.fdflag     : 1 - Lomb Periodogram or other method failed
%                                2 - Not enough high SQI data
%                                3 - Not enough data in the window to analyze
%                                4 - Window is missing too much data
%                                5 - Success
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Gari Clifford HRV Tools
%                   G. Clifford 2001 gari@mit.edu, calc_lfhf.m 
%                   http://www.robots.ox.ac.uk/~gari/CODE/HRV/
%       ChengYu Lui
%       Adriana Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information  
%
%   10-28-2017 Modified by Giulia Da Poian: 
%   Initlialiezed output as NaN
%   Removed fft, Welch and Burg method without resampling. 
%   Added function for resampling when required by the method 
%   Acceppting only one method at time to make the function output more
%   clear
%   Removed option to write results on file directly from this function
%   Converted output to ms^2



% Verify input arguments

if nargin < 5
    error('Eval_FrequencydomainHRVstats: wrong number of input arguments!')
end

if isempty(sqi) 
     sqi(:,1) = tNN;
     sqi(:,2) = ones(length(tNN),1);
end

% Set Defaults
windowlength = HRVparams.windowlength;
fd_threshold1 = HRVparams.sqi.LowQualityThreshold;
fd_threshold2 = HRVparams.RejectionThreshold;
limits = HRVparams.freq.limits;
method = HRVparams.freq.method;
plot_on = HRVparams.freq.plot_on;
debug_sine = HRVparams.freq.debug_sine;   % if use the sine wave to debug
f_sine = HRVparams.freq.debug_freq;       % the frequency of the added sine wave
weight = HRVparams.freq.debug_weight;
sf = HRVparams.freq.resampling_freq;


% Preallocate arrays before entering the loop

out.ulf = nan(1,length(tWin));
out.vlf = nan(1,length(tWin));
out.lf = nan(1,length(tWin));
out.hf = nan(1,length(tWin));
out.lfhf = nan(1,length(tWin));
out.ttlpwr = nan(1,length(tWin));
out.fdflag = nan(1,length(tWin));
% Window by Window Analysis

% Loop through each window of RR data
for iWin = 1:length(tWin)
    % Check window for sufficient data
    if ~isnan(tWin(iWin))    
        % Isolate data in this window
        idx_NN_in_win = find(tNN >= tWin(iWin) & tNN < tWin(iWin) + windowlength);
        idx_sqi_win = find(sqi(:,1) >= tWin(iWin) & sqi(:,1) < tWin(iWin) + windowlength);

        sqi_win = sqi(idx_sqi_win,:);
        t_win = tNN(idx_NN_in_win);
        nn_win = NN(idx_NN_in_win);

        % Analysis of SQI for the window
        lowqual_idx = find(sqi_win(:,2) < fd_threshold1);

        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < fd_threshold2

            % Initialize variables
            % maxF=fs/2; % This calculation works for regularly sampled data
            N     =  length(nn_win);       % RR interval series length
            % m_fs  = 1/mean(nn_win);     % mean frequency of heart rate, i.e., the mean sample rate (fs) of RR sereis  
            % max_f = .5/(min(nn_win));   % max frequency of RR interval series
            % nfft  = 2^nextpow2(N);      % Next power of 2 from N
            nfft = 1024;
            %F = [1/nfft:1/nfft:m_fs];  % setting up frequency vector
            F = [1/nfft:1/nfft:.5];  % setting up frequency vector

            % add sine wave to RR signal
            if debug_sine
                s_sin  = weight*sin(2*pi*f_sine*t);
                nn_win = nn_win + s_sin;
            end

            % subtract mean of segment
            if HRVparams.freq.zero_mean
                rr_0 = nn_win - mean(nn_win); %rudimentary detrending
            else
                rr_0 = nn_win;
            end
            
            switch method
                % 1. Lomb-Scargle Periodogram (no resampling)
                case 'lomb'
                    try
                        [PSDlomb,Flomb] = CalcLomb(t_win,rr_0,F,nfft,HRVparams.freq.normalize_lomb);
                        % plomb equivalent to CalcLomb when normalized 
                        % lomb-scargle periodogram
                        %[PSDmatlablomb,fplombout] = plomb(vv,tt,F,'normalized');
                        [out.ulf(iWin), out.vlf(iWin), out.lf(iWin), out.hf(iWin), out.lfhf(iWin),...
                            out.ttlpwr(iWin)] = CalcLfHfParams(PSDlomb, Flomb, limits, plot_on);
                        
                        out.fdflag(iWin) = 5; %'sucess';

                    catch
                        out.fdflag(iWin) = 1; %'lomb_failed'
                    end

                % 2. Pwelch (after resampling)
                case 'welch'                    
                    % Resampling
                    rr_int = ApplyResamplingToRR(t_win,rr_0,HRVparams);
                    try
                        [PSDwelch,Fwelch] = pwelch(rr_int,[],[],2^nextpow2(length(rr_int)),sf);
                        [out.ulf(iWin), out.vlf(iWin), out.lf(iWin), out.hf(iWin), out.lfhf(iWin),...
                                out.ttlpwr(iWin)] = CalcLfHfParams(PSDwelch, Fwelch, limits, plot_on) ;  
                        out.fdflag(iWin) = 5; %'sucess';
                    catch
                        out.fdflag(iWin) = 1; %'Pwelch faild';
                    end
                    
                % 3. FFT (after resampling)
                case 'fft'
                    % Resampling
                    rr_int = ApplyResamplingToRR(t_win,rr_0,HRVparams);
                    try
                        Y_FFT2 = fft(rr_int);  % fft using the original sample length
                        PSDfft = Y_FFT2.*conj(Y_FFT2)/length(rr_int);
                        Ffft = sf*(0:floor(length(rr_int)/2)+1)/length(rr_int);
                        [out.ulf(iWin), out.vlf(iWin), out.lf(iWin), out.hf(iWin), out.lfhf(iWin),...
                                out.ttlpwr(iWin)] = CalcLfHfParams(PSDfft(1:length(Ffft)), Ffft, limits, plot_on);
                        out.fdflag(iWin) = 5; %'sucess';
                    catch
                        out.fdflag(iWin) = 1; %'fft faild';
                    end

                % 4. Burg (after resampling)
                case 'burg'
                    % Resampling
                    rr_int = ApplyResamplingToRR(t_win,rr_0,HRVparams);
                    try
                    % ChengYu's Burg Method with resampled data
                    pb = HRVparams.freq.resampled_burg_poles; % pole setting
                    [PSD_Burg1,f_Burg1] = pburg(rr_int,pb,2^nextpow2(length(rr_int)),sf);
                    [out.ulf(iWin), out.vlf(iWin), out.lf(iWin), out.hf(iWin), out.lfhf(iWin),...
                            out.ttlpwr(iWin)] = CalcLfHfParams(PSD_Burg1,f_Burg1, limits, plot_on);
                        out.fdflag(iWin) = 5; %'sucess';
                    catch
                        out.fdflag(iWin) = 1; %'Burg faild';
                    end
            end
        else
            out.fdflag(iWin) = 2; %'nt_enuf_hi_sqi_data'
        end % end conditional statements that run only if SQI is adequate

    else
        out.fdflag(iWin) = 3; %'nt_enuf_data_in_win';   
    end % end check for sufficient data
end % end of loop through windows


%   References:
%               5 minutes is the shortest period that HRV spectral metrics
%               should be calculated according to 184(Clifford Thesis). With a 5 min
%               window, the lowest frequency that can be theoretically resolved is
%               1/300 (.003 Hz).