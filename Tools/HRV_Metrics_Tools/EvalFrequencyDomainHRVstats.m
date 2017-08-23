function [ulfL, vlfL, lfL, hfL, lfhfL, ttlpwrL, methods, fdflag, windows_all] = EvalFrequencyDomainHRVstats(NN, tNN, sqi, HRVparams, windows_all)
%
% [ulfL, vlfL, lfL, hfL, lfhfL, ttlpwrL, methods, fdflag, windows_all] = ...
%         EvalFrequencyDomainHRVstats (NN, tNN, , sqi, settings, windows_all)
%   
%   OVERVIEW:   This function returns frequency domain HRV metrics 
%               calculated on input NN intervals.
%
%   INPUT:      MANDATORY:
%               NN          : a single row of NN (normal normal) interval
%                             data in seconds
%               
%               OPTIONAL:
%               tNN         : a single row of time indices of the rr interval 
%                             data (seconds)
%               sqi         : Signal Quality Index; Requires a matrix with
%                             at least two columns. Column 1 should be
%                             timestamps of each sqi measure, and Column 2
%                             should be SQI on a scale from 1 to 100.
%                             Additional columns can be included with
%                             additional sqi at the same timestamps
%               HRVparams   : struct of settings for hrv_toolbox analysis
%               windows_all :
%
%   OUTPUT:     ulfL        :
%               vlfL        :
%               lfl         :
%               hfL         :
%               lfhfL       :
%               ttlpwrL     :
%               methods     :
%               fdflag      : 1 - Lomb Periodogram Failed
%                             2 - Not enough high SQI data
%                             3 - Not enough data in the window to analyze
%                             4 - Window is missing too much data
%                             5 - Success
%               windows_all : 
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
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
%%

% Verify input arguments

if nargin< 1
    error('Eval_FrequencydomainHRVstats: wrong number of input arguments!')
end
if nargin<2 || isempty(tNN)
        tNN = cumsum(NN);
end
if nargin<3 || isempty(sqi) 
        sqi(:,1) = tNN;
        sqi(:,2) = 100* ones(length(tNN),1);
end
if nargin<4 || isempty(HRVparams) 
        HRVparams = initialize_HRVparams('demo');
end
if nargin<6 || isempty(windows_all)
    windows_all = CreateWindowRRintervals(tNN, NN, HRVparams);
end

% Set Defaults


windowlength = HRVparams.windowlength;
fd_threshold1 = HRVparams.freq.threshold1;
fd_threshold2 = HRVparams.freq.threshold2;
limits = HRVparams.freq.limits;
methods = HRVparams.freq.methods;
plot_on = HRVparams.freq.plot_on;
dataoutput = HRVparams.freq.dataoutput;
debug_sine = HRVparams.freq.debug_sine;   % if use the sine wave to debug
f_sine = HRVparams.freq.debug_freq;       % the frequency of the added sine wave
weight = HRVparams.freq.debug_weight;

flagLomb=false;
flagBurg=false; 
flagWelch=false;  
flagFft = false;

for m = 1:length(methods)
    if strcmpi(methods{m},'welch')
        flagWelch = true;
    elseif strcmpi(methods{m},'burg')
        flagBurg = true;
    elseif strcmpi(methods{m},'lomb')
        flagLomb = true;
    elseif strcmpi(methods{m},'fft')
        flagFft = true;
    end
end


% Preallocate arrays before entering the loop

ulfL = zeros(1,length(windows_all));
vlfL = zeros(1,length(windows_all));
lfL = zeros(1,length(windows_all));
hfL = zeros(1,length(windows_all));
lfhfL = zeros(1,length(windows_all));
ttlpwrL = zeros(1,length(windows_all));
fdflag = zeros(1,length(windows_all));

%% Window by Window Analysis

% Loop through each window of RR data
for i_win = 1:length(windows_all)
    % Check window for sufficient data
    if ~isnan(windows_all(i_win))    
        % Isolate data in this window
        idx_NN_in_win = find(tNN >= windows_all(i_win) & tNN < windows_all(i_win) + windowlength);
        idx_sqi_win = find(sqi(:,1) >= windows_all(i_win) & sqi(:,1) < windows_all(i_win) + windowlength);

        sqi_win = sqi(idx_sqi_win,:);
        t_win = tNN(idx_NN_in_win);
        nn_win = NN(idx_NN_in_win);

        % Analysis of SQI for the window
        lowqual_idx = find(sqi_win(:,2) < fd_threshold1);

        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < fd_threshold2

            % Initialize variables
            % maxF=fs/2; % This calculation works for regularly sampled data
            N     = length(nn_win);     % RR interval series length
            % m_fs  = 1/mean(nn_win);     % mean frequency of heart rate, i.e., the mean sample rate (fs) of RR sereis  
            % max_f = .5/(min(nn_win));   % max frequency of RR interval series
            % nfft  = 2^nextpow2(N);      % Next power of 2 from N
            nfft = 1024;
            %F = [1/nfft:1/nfft:m_fs];  % setting up frequency vector
            F = [1/nfft:1/nfft:.5];  % setting up frequency vector

            % add  sine wave to RR signal
            if debug_sine
                s_sin  = weight*sin(2*pi*f_sine*t);
                nn_win     = nn_win + s_sin;
            end

            % subtract mean of segment
            if HRVparams.freq.zero_mean
                rr_0 = nn_win - mean(nn_win); %rudimentary detrending
            else
                rr_0 = nn_win;
            end
            
            %% Initialize txt file output if needed
            if dataoutput == 1
                if i_win == 1
                    writing_dir = [HRVparams.writedata filesep]; 
                    % filesep returns the platform-specific file separator character
                    fileID = fopen(strcat(writing_dir,'freq_dom.txt'),'w');
                    fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s %6s %10s %10s %10s\n', ...
                        'segment#','start','end', 'ulf','vlf','lf','hf','lfhf','ttlpwr','method','resampled');
                end
            end

            %% 1. Lomb-Scargle Periodogram
            if flagLomb
                try
                    [PSDlomb,Flomb] = CalcLomb(t_win,rr_0,F,nfft,HRVparams.freq.normalize_lomb);
                    % plomb equivalent to CalcLomb when normalized 
                    % lomb-scargle periodogram
                    %[PSDmatlablomb,fplombout] = plomb(vv,tt,F,'normalized');
                    [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSDlomb, Flomb, limits, plot_on);
                    meth = 'lomb';
                    resamp = 0;
                    if dataoutput == 1
                        fprintf(fileID,'%5i %8.4f %8.4f %1.8f %1.8f %1.8f %1.8f %3.4f %1.8f %10s %1i\n',...
                                i_win,first,last,ulf,vlf,lf,hf,lfhf,ttlpwr,meth,resamp);
                    end
                    if plot_on
                        title('PSD - Lomb Method')
                    end

                    ulfL(i_win) = ulf;
                    vlfL(i_win) = vlf;
                    lfL(i_win) = lf; 
                    hfL(i_win) = hf;
                    lfhfL(i_win) = lfhf;
                    ttlpwrL(i_win) = ttlpwr;
                    methods{i_win} = 'lomb';
                    fdflag(i_win) = 5; %'sucess';

                catch
                    % warning('Lomb Failed. Using Burg After Resampling')
                    % sf = 7; % (Hz) resampling frequency
                    % ti = [t_win(1):1/sf:t_win(end)];% time values for interp.
                    % % rr_lin = interp1(t_win,rr_0,ti','linear')'; %linear interpolation
                    % rr_int = interp1(t_win,rr_0,ti','spline')'; % cubic spline interpolation
                    % pb = 100; % pole setting,
                    % [PSD_Burg1,f_Burg1] = pburg(rr_int,pb,2^nextpow2(length(rr_int)),sf);
                    % [ulf vlf lf hf lfhf ttlpwr] = CalcLfHfParams(PSD_Burg1,f_Burg1, limits, plot_on);
                    % meth = 'burg';
                    % resamp = 0;

                    ulfL(i_win) = NaN;
                    vlfL(i_win) = NaN;
                    lfL(i_win) = NaN; 
                    hfL(i_win) = NaN;
                    lfhfL(i_win) = NaN;
                    ttlpwrL(i_win) = NaN;
                    methods{i_win} = 'NA';
                    fdflag(i_win) = 1; %'lomb_failed'
                    if dataoutput == 1
                        fprintf(fileID,'%5i %8.4f %8.4f %1.8f %1.8f %1.8f %1.8f %3.4f %1.8f %10s %1i\n',...
                            i_win,first,last,ulf,vlf,lf,hf,lfhf,ttlpwr,meth,resamp);
                    end
                end
            end

            %% 2. Burg Method
            if flagBurg
                p = HRVparams.freq.burg_poles; % pole setting
                [PSDburgBRS,FburgBRS] = pburg(rr_0,p,F,m_fs);
                [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSDburgBRS, FburgBRS, limits, plot_on);
                meth = 'burg';
                resamp = 0;
                if dataoutput == 1
                    fprintf(fileID,'%5i %8.4f %8.4f %1.8f %1.8f %1.8f %1.8f %3.4f %1.8f %10s %1i\n',...
                        i_win,first,last,ulf,vlf,lf,hf,lfhf,ttlpwr,meth,resamp);
                end
                if plot_on
                    title('ChengYu Burg Method')
                end
            end
            %% 3. FFT (before resampling)
            if flagFft
                
                %Y = fft(rr_0); 
                %P2 = abs(Y); % The two sided spectrum
                %P1 = P2(1:N/2+1); % Single sided spectrum
                %P1(2:end-1) = 2*P1(2:end-1);
                %f = 7*(0:(N/2))/N;
                
                
                Y_FFT0       = fft(rr_0);  % fft using the original sample length
                PSD_FFT0     = Y_FFT0.*conj(Y_FFT0)/N; % double sided spectrum
                %m_fs0 = 1/mean(rr_0);
                f_FFT0       = m_fs*(0:N/2+1)/N; 
                % temporary!
                f_FFT0      = 7*(0:(N/2))/N;
                [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSD_FFT0(1:length(f_FFT0)), f_FFT0, limits, plot_on);
                meth = 'fft';
                resamp = 0;
                if dataoutput == 1
                    fprintf(fileID,'%5i %8.4f %8.4f %1.8f %1.8f %1.8f %1.8f %3.4f %1.8f %10s %1i\n',...
                        i_win,first,last,ulf,vlf,lf,hf,lfhf,ttlpwr,meth,resamp);
                end
                if plot_on
                    title('Chengyus FFT Method 1 Before Resampling')
                end
            end
            %% 4. Pwelch (before resampling)
            if flagWelch
                % temporary!
                m_fs      = 7;
                [PSDwelchBRS,FwelchBRS] = pwelch(rr_0,[],[],2^nextpow2(length(rr_0)),m_fs);
                [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSDwelchBRS, FwelchBRS, limits, plot_on);
                meth = 'welch';
                resamp = 0;
                if dataoutput == 1
                    fprintf(fileID,'%5i %8.4f %8.4f %1.8f %1.8f %1.8f %1.8f %3.4f %1.8f %10s %1i\n',...
                        i_win,first,last,ulf,vlf,lf,hf,lfhf,ttlpwr,meth,resamp);
                end
                if plot_on
                    title('Pwelch method: without resampling')
                end
            end
            %% Resampling
            sf = HRVparams.freq.resampling_freq; % (Hz) resampling frequency
            ti = [t_win(1):1/sf:t_win(end)];% time values for interp.
            interp_method = HRVparams.freq.resample_interp_method;

            if strcmp(interp_method,'cub')
                rr_int = interp1(t_win,rr_0,ti','spline')'; % cubic spline interpolation
            elseif strcmp(interp_method,'lin')
                rr_int = interp1(t_win,rr_0,ti','linear')'; %linear interpolation
            else
                warning('using cubic spline method')
                rr_int = interp1(t_win,rr_0,ti','spline')'; % cubic spline interpolation
            end

            %% 5. Welch 
            if flagWelch
                [PSDwelch,Fwelch] = pwelch(rr_int,[],[],2^nextpow2(length(rr_int)),sf);
                [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSDwelch, Fwelch, limits, plot_on) ;  
                meth = 'welch';
                resamp = 1;
                if dataoutput == 1
                    fprintf(fileID,'%5i %8.4f %8.4f %1.8f %1.8f %1.8f %1.8f %3.4f %1.8f %10s %1i\n',...
                        i_win,first,last,ulf,vlf,lf,hf,lfhf,ttlpwr,meth,resamp);
                end
                if plot_on
                    title(strcat('Pwelch method: cubic spline interpolation resampling (',num2str(sf),' Hz)'));
                end
            end

            %% 6. FFT
            if flagFft
                Y_FFT2 = fft(rr_int);  % fft using the original sample length
                PSDfft = Y_FFT2.*conj(Y_FFT2)/length(rr_int);
                Ffft = sf*(0:floor(length(rr_int)/2)+1)/length(rr_int);
                [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSDfft(1:length(Ffft)), Ffft, limits, plot_on);
                meth = 'fft';
                resamp = 1;
                if dataoutput == 1      
                    fprintf(fileID,'%5i %8.4f %8.4f %1.8f %1.8f %1.8f %1.8f %3.4f %1.8f %10s %1i\n',...
                        i_win,first,last,ulf,vlf,lf,hf,lfhf,ttlpwr,meth,resamp);
                end
                if plot_on
                    title(strcat('Chengyus FFT Method with resampling (cubic spline ',num2str(sf),' Hz)'));
                end
            end
            %% 7. Burg
            if flagBurg
                % ChengYu's Burg Method with resampled data
                pb = HRVparams.freq.resampled_burg_poles; % pole setting
                [PSD_Burg1,f_Burg1] = pburg(rr_int,pb,2^nextpow2(length(rr_int)),sf);
                [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSD_Burg1,f_Burg1, limits, plot_on);
                meth = 'burg';
                resamp = 1;
                if dataoutput == 1
                    fprintf(fileID,'%5i %8.4f %8.4f %1.8f %1.8f %1.8f %1.8f %3.4f %1.8f %10s %1i\n',...
                        i_win,first,last,ulf,vlf,lf,hf,lfhf,ttlpwr,meth,resamp);
                end
                if plot_on
                    title(strcat('Chengyus Burg Method with resampling (cubic spline ',num2str(sf),' Hz)'));
                end
            end
        else
            ulfL(i_win) = NaN;
            vlfL(i_win) = NaN;
            lfL(i_win) = NaN;
            hfL(i_win) = NaN;
            lfhfL(i_win) = NaN;
            ttlpwrL(i_win) = NaN;
            methods{i_win} = 'NA';
            fdflag(i_win) = 2; %'nt_enuf_hi_sqi_data'
        end % end conditional statements that run only if SQI is adequate

    else
        ulfL(i_win) = NaN;
        vlfL(i_win) = NaN;
        lfL(i_win) = NaN;
        hfL(i_win) = NaN;
        lfhfL(i_win) = NaN;
        ttlpwrL(i_win) = NaN;
        methods{i_win} = 'NA';
        fdflag(i_win) = 3; %'nt_enuf_data_in_win';
        
    end % end check for sufficient data
    

    
    if dataoutput == 1
        if i_win == i_win(end)
            fclose(fileID);    
        end
    end
end % end of loop through windows
end % end function


%   References:
%               5 minutes is the shortest period that HRV spectral metrics
%               should be calculated according to 184(Clifford Thesis). With a 5 min
%               window, the lowest frequency that can be theoretically resolved is
%               1/300 (.003 Hz).