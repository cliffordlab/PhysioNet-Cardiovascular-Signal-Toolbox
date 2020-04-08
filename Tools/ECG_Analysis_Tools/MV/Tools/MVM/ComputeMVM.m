function [energyinband_array,sqi_array,hr_array] = ComputeMVM(ecg,ann_ecg,Fs,segment_size,normalize)

% [energyinband_array,sqi_array] = ComputeMVM(ecg,ann,Fs,segment_size,normalize)
%   OVERVIEW:   This function returns the QRS MVM measured in non overlapping
%   windows of size 'segment_size' in minutes. It also returns the
%   signal quality and mean heart rate for each window.
%
%   INPUT:      MANDATORY:
%               ecg             : a single row of ECG data in samples. Must be sampled at 1 KHz.
%
%               ann             : annotation for the ecg passed as input. The ann variable
%                                 is a struct with Q,R,S fiduical points in sample number as fields
%                                 for the struct. i.e. ann.R, ann.Q, ...
%
%               Fs              : sampling frequency for the ecg signal (Hz)
%
%               segment_size        : the length of non overlapping windows
%                                 being analyzed in minutes. Default size
%                                 is 5 minutes
%
%               normalize       : boolean value (0 or 1). When set to 1
%                               will normalize each analysis window by the median R
%                               amplitude. Default is 1.
%
%
%   OUTPUT:
%               energyinband_array     : the energy for QRS complex morphological variability measured in
%                                        the every 2-7 beats region for each analysis window.
%
%               sqi_array        : the signal quality index measured for each analysis window.
%
%               hr_array          : the median heart rate measured for analysis
%                                 window
%   
%   REF:
%   The script is based on the algorithm given in the following paper,
%   Liu Y, Syed Z, Scirica BM, Morrow DA, Guttag JV, Stultz CM. 
%   ECG morphological variability in beat space for risk stratification after acute coronary syndrome. 
%   J Am Heart Assoc. 2014;3(3):e000981. Published 2014 Jun 24. doi:10.1161/JAHA.114.000981
%	REPO:
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

% divide into 5 minute segments
% initialize variables
segment_size_samples = Fs*60*segment_size; % size of analysis window in samples
Nsegments = floor(length(ecg)/segment_size_samples); % number of segments analyzed
energyinband_array = NaN(1, Nsegments); % array for storing the MVM energy measured in each analysis window
sqi_array = NaN(1, Nsegments); % array of storing the signal quality for each window analyzed
hr_array = NaN(1,Nsegments); % array for storing the heart rate estimate for each window analyzed
lengthflag = 0;

sqi_threshold = 0.7;

% convert ecg to row vector
if (size(ecg,1) > size(ecg,2))
   ecg = ecg'; 
end

if (Nsegments == 0)   % length of data < segment_size_samples, analyze available data
    Nsegments = 1; lengthflag = 1;
end

for segmentidx = 1:Nsegments
    
    try
        
        
        % isolate segment
        if(lengthflag)
            % if the length of signal is < segment_size in samples,
            % determine annotations till the end of the record.
            fiveminsegment = ecg(segmentidx*segment_size_samples+1-segment_size_samples:end);
            ann.R = ann_ecg.R(ann_ecg.R < length(ecg));
            ann.R = ann.R(ann.R > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;% + 1;
            ann.Q = ann_ecg.Q(ann_ecg.Q < length(ecg));
            ann.Q = ann.Q(ann.Q > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;% + 1;
            ann.S = ann_ecg.S(ann_ecg.S < length(ecg));
            ann.S = ann.S(ann.S > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;% + 1;
            ann.QRSon = ann_ecg.QRSon(ann_ecg.QRSon < length(ecg));
            ann.QRSon = ann.QRSon(ann.QRSon > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;
            ann.QRSoff = ann_ecg.QRSoff(ann_ecg.QRSoff < length(ecg));
            ann.QRSoff = ann.QRSoff(ann.QRSoff > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;
            lengthflag = 0;
        else
            % get ecg signal and corresponding annotations in the current
            % window
            fiveminsegment = ecg(segmentidx*segment_size_samples+1-segment_size_samples:segmentidx*segment_size_samples);
            ann = ann_ecg;
            ann.R = ann_ecg.R(ann_ecg.R < segmentidx*segment_size_samples);
            ann.R = ann.R(ann.R > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;% + 1;
            ann.Q = ann_ecg.Q(ann_ecg.Q < segmentidx*segment_size_samples);
            ann.Q = ann.Q(ann.Q > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;% + 1;
            ann.S = ann_ecg.S(ann_ecg.S < segmentidx*segment_size_samples);
            ann.S = ann.S(ann.S > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;% + 1;
            ann.QRSon = ann_ecg.QRSon(ann_ecg.QRSon < segmentidx*segment_size_samples);
            ann.QRSon = ann.QRSon(ann.QRSon > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;
            ann.QRSoff = ann_ecg.QRSoff(ann_ecg.QRSoff < segmentidx*segment_size_samples);
            ann.QRSoff = ann.QRSoff(ann.QRSoff > segmentidx*segment_size_samples+1-segment_size_samples) - (segmentidx-1)*segment_size_samples;
            %figure(2); plot(fiveminsegment); hold on; % Check annotations after loading 
            %scatter(ann.R(~isnan(ann.R)), fiveminsegment(ann.R(~isnan(ann.R)))); hold off;
        end
        
        ecg_mv = fiveminsegment; clear fiveminsegment;
        
        
        
        % remove any nan values from the annotations
        R = round(ann.R);
        sum(isnan(R)); R = R(~isnan(R)); sum(isnan(R));
        Q = round(ann.Q);
        sum(isnan(Q)); Q = Q(~isnan(Q)); sum(isnan(Q));
        S = round(ann.S);
        sum(isnan(S)); S = S(~isnan(S)); sum(isnan(S));
        
        % compute the median Q,S and R amplitude, used to determine
        % normalization factor
        median_Q_amp = nanmedian(ecg_mv(1,Q));
        median_R_amp = nanmedian(ecg_mv(1,R));
        median_S_amp = nanmedian(ecg_mv(1,S));
        
        % compute normalization factor
        if (median_S_amp <= median_Q_amp || isnan(median_Q_amp))
            norm_factor = median_R_amp - median_S_amp;
        else
            if (~isempty(median_Q_amp))
                norm_factor = median_R_amp - median_Q_amp;
            else
                % if unable to determine, set to 1
                disp('unable to compute norm factor, setting to 1')
                norm_factor = 1;
            end
        end
        
        % normalize signal
        if (normalize)
            ecg_mv = ecg_mv / norm_factor;
        end
        
        % generate second set of annotations
        refqrs = ann.R;
        % mean impute missing val
        ecg_mv(isnan(ecg_mv)) = nanmean(ecg_mv);
        testqrs = wqrsm(ecg_mv, Fs);
        thres = 0.1; margin = 0; windowlen = 60*segment_size;
        % Determine signal quality
        refqrs = refqrs./Fs; testqrs = testqrs./Fs; % convert to time (s)
        [current_win_sqi,Se,PPV,Nb] = run_sqi(refqrs,testqrs,thres,margin,windowlen,Fs);
        sqi_array(segmentidx) = current_win_sqi;
        if (current_win_sqi > sqi_threshold)   % may need to be set higher
            
            % perform morphological variability on QRS
            
            % determine qrs onset and qrs offset annotations over 5 minute segment
            qrson = ann.QRSon;
            qrsoff = ann.QRSoff;
            if (qrsoff(1) < qrson(1))
                startidx = find(qrsoff < qrson(1));
                qrsoff = qrsoff(startidx+1:end);
                qrson = qrson(1:length(qrsoff));
            else
                qrson = qrson(1:length(qrsoff));
            end
            
            % make sure sufficient annotations to estimate median S
            % amplitude
            if (-1*median_S_amp > median_R_amp && (length(S) > 200))
                alignmentpoint = 'S';
            else
                alignmentpoint = 'R';
            end
            
            % re-compute R peak
            radj = zeros(1,length(qrson)); % r-adjusted, used in ectopic beat removal
            % setting max value between qrson and qrsoff as r-pk
            for qrsonidx = 1:length(qrson)
                %[pks,locs] = findpeaks(ecg(qrson(qrsonidx):qrsoff(qrsonidx)));
                [pks,locs] = findpeaks(ecg_mv(qrson(qrsonidx):qrsoff(qrsonidx)));
                maxpk = max(pks); maxloc = locs((maxpk == pks));
                if (isempty(maxloc))
                    radj(qrsonidx) = qrson(qrsonidx);
                else
                    radj(qrsonidx) = qrson(qrsonidx)-1+maxloc(1);   % 1 in case of duplicates
                end
            end
            % compute alignemnt point
            % find optimal alignment points max r-pk or s-valley
            alignpoint = zeros(1,length(qrson)); % used in qrs complex isolation
            for qrsonidx = 1:length(qrson)
                
                if (strcmp(alignmentpoint, 'R'))
                    [pks,locs] = findpeaks((ecg_mv(qrson(qrsonidx):qrsoff(qrsonidx))));
                else
                    [pks,locs] = findpeaks(-1*(ecg_mv(qrson(qrsonidx):qrsoff(qrsonidx))));
                end
                
                if (isempty(pks))   % if no positive peaks define mid-point between qrson and qrsoff as alignment point
                    locs = round(mean(qrson(qrsonidx), qrsoff(qrsonidx)));
                    pks = ecg_mv(locs);
                    alignpoint(qrsonidx) = locs; continue;
                end
                maxpk = max(pks); maxloc = locs((maxpk == pks));
                alignpoint(qrsonidx) = qrson(qrsonidx)-1+maxloc(1);   % 1 in case of duplicates
            end
            
            % detect pvc beats
            signal = ecg_mv; th = 0.1;
            [qrs_time pvc_output] = detectpvc2(signal',Fs,th);
            qrs_time = qrs_time(logical(pvc_output));
                        
            % only use 5 minute signals with complete annotations to maintane validity
            % of fourier analysis (may want to replace with interpolation techhnique)
            if (sum(isnan(qrson)) || sum(isnan(qrsoff)))
                disp(['annotations incomplete ' num2str(fileidx)]);
                continue;
            end
            
            % ectopic beat removal
            qrsonectopicfree = qrson(1);
            rectopicfree = radj(1);
            qrsoffectopicfree = qrsoff(1);
            alignpointsectopicfree = alignpoint(1);
            for ridx = 2:length(radj)-1   % start from idx 2 and goto length(r)-1
                if (ridx <= length(radj)-40)
                    %meanRR = mean(diff(radj(ridx-1:ridx+38)));   % compute mean over 40 beat segment
                    meanRR = median(diff(radj(ridx-1:ridx+38)));   % compute median over 40 beat segment
                else
                    %meanRR = mean(diff(radj(end-39:end)));   % compute mean over last 40 beat segment
                    meanRR = median(diff(radj(end-39:end)));   % compute median over last 40 beat segment
                end
                preRR = radj(ridx) - radj(ridx-1);    % pre-RR interval
                postRR = radj(ridx+1) - radj(ridx);   % post-RR interval
                % if the pre-RR and post-RR vary from the mean RR by more than
                % 20% consider as ectopic
                if ((((meanRR-preRR) > 0.2*meanRR) && ((postRR-meanRR) > 0.2*meanRR)) || (((meanRR-preRR) > 0.2*meanRR) && ((meanRR-postRR) > 0.2*meanRR))) % && ((meanRR-postRR) > 0.2*meanRR)))  % 326, 367
                    continue;
                else % include the QRS complex, may also want to include t-wave onset/offset/amplitude/QT-interval
                    qrsonectopicfree = [qrsonectopicfree qrson(ridx)];
                    rectopicfree = [rectopicfree radj(ridx)];
                    qrsoffectopicfree = [qrsoffectopicfree qrsoff(ridx)];
                    alignpointsectopicfree = [alignpointsectopicfree alignpoint(ridx)];
                end
            end
            
            % remove ectopic beats
            [C,ia,ib] = intersect(rectopicfree, qrs_time);
            rectopicfree(ia) = [];
            qrsonectopicfree(ia) = [];
            qrsoffectopicfree(ia) = [];
            alignpointsectopicfree(ia) = [];
            
            % calculate NTWDseries and TWDseries
            % Isolate complexes
            [Complexes] = IsolateQRSComplexes(ecg_mv, qrsoffectopicfree, alignpointsectopicfree, qrsonectopicfree);
            
            % NTWD and TWD
            [NTWDseries TWDseries] = SquaredDiffComputation_beatbybeat(Complexes);
            
            % compute spectrum
            nfft = 256; % standardize the length of the fft to 256
            pxx = pwelch(NTWDseries,[],[],nfft);
            %             pxx = pwelch(TWDseries,[],[],nfft);
            N = nfft;
            stepsize = 0.5/length(pxx);
            
            dl = 2; du = 7;
            indexeverytwobeats = round(N/dl, 0);    % corresponding to every 2 beats
            indexeverysevenbeats = floor(N/du); % corresponding to every 7 beats
            
            % compute energy in every 2-7 beats region
            energyinband = trapz((indexeverysevenbeats:indexeverytwobeats) * stepsize, pxx(indexeverysevenbeats:indexeverytwobeats));
            energyinband_array(segmentidx) = energyinband; disp(['MVM for current window is ' num2str(energyinband)]);
            
            % hr est
            hr_array(segmentidx) = 60/(nanmedian(diff(alignpointsectopicfree))/Fs);
        else
            disp('unclean segment')
        end
        
    catch
        disp('error in segment')
        energyinband_array(segmentidx) = -inf;
        sqi_array(segmentidx) = -inf;
        hr_array(segmentidx) = -inf;
    end
    
end

end

