function [ energyinbandcollection ] = MVQRSRealData( ecgdata, annotations, fs )
%MVQRSRealData, calculate morphological variability in QRS complexes for real
%data. To compensate for inaccurate annotations, The variability is
%calculated from median onset interval to r-pk and r-pk to median offset
%interval between consecutive QRS complexes.
%   Inputs:     ecgdata - input ECG.
%               annotations - annotations structure (as obtained from Dr.
%               Qiao).
%   Outputs:    energyinbandcollection - array containing energy in every
%               2-7 beats band for each 5 minute segment

% Last modified by Ismail Sadiq 7/3/2017
% Dependencies: medianfilter_is, wqrsm, run_sqi, MDSeriesCalc,
%               DiffSeriesCalc, DiffSeriesComplexesCalc

bw = medianfilter_is(ecgdata, fs);
ecgdatafiltered = ecgdata - bw;
%     figure; plot(ecgdatafiltered(1:5000));
%     hold on;
%     stem(annotations.QRSon(1:5), ones(1, 5));
%     stem(annotations.Q(1:5), ones(1, 5));
%     stem(annotations.R(1:5), ones(1, 5));
%     stem(annotations.S(1:5), ones(1, 5));
%     stem(annotations.QRSoff(1:5), ones(1, 5));
%     hold off;
figure; plot(ecgdatafiltered ./ nanmean(annotations.R));
hold on;
stem(annotations.QRSon, 1e-4*ones(1, length(annotations.QRSon)));
stem(annotations.Q, 1e-4*ones(1, length(annotations.Q)));
stem(annotations.R, 1e-4*ones(1, length(annotations.R)));
stem(annotations.S, 1e-4*ones(1, length(annotations.S)));
stem(annotations.QRSoff, 1e-4*ones(1, length(annotations.QRSoff)));
hold off;
%% divide into 5 minute segments
N5 = fs*60*5;
numberof5minutesegments = floor(length(ecgdatafiltered)/N5);
energyinbandcollection = NaN(1, numberof5minutesegments);
for segmentidx = 1:numberof5minutesegments
    segmentidx
    % isolate segment
    fiveminsegment = ecgdatafiltered(segmentidx*N5+1-N5:segmentidx*N5);
    % generate second set of annotations
    refqrs = annotations.QRSon(annotations.QRSon < segmentidx*N5);
    refqrs = refqrs(refqrs > segmentidx*N5+1-N5) - (segmentidx-1)*N5;
    testqrs = wqrsm(fiveminsegment, fs);
    thres = 0.1; margin = 0; windowlen = 10;
    if ((length(testqrs) - length(refqrs)) == 1)
        testqrs = testqrs(2:end);
    else if ((length(refqrs) - length(testqrs)) == 1)
            refqrs = refqrs(1:end-1);
        end
    end
    % Determine signal quality
    refqrs = refqrs./fs; testqrs = testqrs./fs; % convert to time (s)
    [F1,Se,PPV,Nb] = run_sqi(refqrs,testqrs(2:end-1),thres,margin,windowlen,fs);
    if (F1 > 0.95)   % may need to be set higher
        %% perform morphological variability on QRS
        % normalize signal
        q = annotations.Q(annotations.Q < segmentidx*N5);
        q = q(q > segmentidx*N5+1-N5);
        r = annotations.R(annotations.R < segmentidx*N5);
        r = r(r > segmentidx*N5+1-N5);
        s = annotations.S(annotations.S < segmentidx*N5);
        s = s(s > segmentidx*N5+1-N5);
        t = annotations.T(annotations.T < segmentidx*N5);
        t = t(t > segmentidx*N5+1-N5);
        % onset/offset
        qrson = annotations.QRSon(annotations.QRSon < segmentidx*N5);
        qrson = qrson(qrson > segmentidx*N5+1-N5);
        qrsoff = annotations.QRSoff(annotations.QRSoff < segmentidx*N5);
        qrsoff = qrsoff(qrsoff > segmentidx*N5+1-N5);
        % adjust qrsoff so that it is greater than qrson(1)
        if (qrsoff(1) < qrson(1))
            startidx = find(qrsoff < qrson(1));
            qrsoff = qrsoff(startidx+1:end);
            qrson = qrson(1:length(qrsoff));
        else
            qrson = qrson(1:length(qrsoff));
        end
        % adjust r pks (may want to optimize r pk detection)
        %             if (r(end) < qrson(end))    % temporary solution incase r peaks inconsistant with qrson
        %                 continue;
        %             end
        %             radj = NaN(size(qrson));
        %             for qrsonidx = 1:length(qrson)
        %                 ridx = find(r < qrson(qrsonidx));
        %                 if(isempty(ridx))
        %                     ridx = 0;
        %                 end
        %                 if (r(ridx(end) + 1) < qrsoff(qrsonidx))
        %                     radj(qrsonidx) = r(ridx(end) + 1);
        %                 end
        %             end
        radj = zeros(1,length(qrson));
        % compute means points to determine alignment point
        medianQ = nanmedian(ecgdatafiltered(q));
        medianR = nanmedian(ecgdatafiltered(r));
        medianS = nanmedian(ecgdatafiltered(s));
        meanQ = medianQ;
        meanR = medianR;
        meanS = medianS;
        % make sure sufficient annotations (may need changing)
        if (-1*meanS > meanR && (length(s) > 200))
            alignmentpoint = 'S';
        else
            alignmentpoint = 'R';
        end
        % setting max value between qrson and qrsoff as r-pk
        for qrsonidx = 1:length(qrson)
            [pks,locs] = findpeaks(ecgdatafiltered(qrson(qrsonidx):qrsoff(qrsonidx)));
            maxpk = max(pks); maxloc = locs((maxpk == pks));
            if (isempty(maxloc))
                radj(qrsonidx) = qrson(qrsonidx);
            else
                radj(qrsonidx) = qrson(qrsonidx)-1+maxloc(1);   % 1 in case of duplicates
            end
        end
        % find optimal alignment points max r-pk or s-valley
        alignpoint = zeros(1,length(qrson));
        for qrsonidx = 1:length(qrson)
            if (strcmp(alignmentpoint, 'R'))
                [pks,locs] = findpeaks((ecgdatafiltered(qrson(qrsonidx):qrsoff(qrsonidx))));
            else
                [pks,locs] = findpeaks(-1*(ecgdatafiltered(qrson(qrsonidx):qrsoff(qrsonidx))));
            end
            if (isempty(pks))   % if no positive peaks search for s valley
                %                     [pks,locs] = findpeaks(abs(ecgdatafiltered(qrson(qrsonidx):qrsoff(qrsonidx))));
                locs = round(mean(qrson(qrsonidx), qrsoff(qrsonidx)));
                pks = ecgdatafiltered(locs);
                alignpoint(qrsonidx) = locs; continue;
            end
            maxpk = max(pks); maxloc = locs((maxpk == pks));
            alignpoint(qrsonidx) = qrson(qrsonidx)-1+maxloc(1);   % 1 in case of duplicates
        end
        % end editting qrson
        ton = annotations.Ton(annotations.Ton < segmentidx*N5);
        ton = ton(ton > segmentidx*N5+1-N5);
        toff = annotations.Toff(annotations.Toff < segmentidx*N5);
        toff = toff(toff > segmentidx*N5+1-N5);
        % first offset needs to be greater than first onset
        %             % aligning remaining annotations
        %             q = q(q > qrson(1));
        %             r = r(r > q(1));
        %             s = s(s > r(1));
        %             qrsoff = qrsoff(qrsoff > qrson(1));
        %             ton = ton(ton > qrsoff(1));
        %             t = t(t > ton(1));
        %             toff = toff(toff > t(1));
        % compute mean points
        % meanQ = mean(ecgdatafiltered(q));
        % meanR = mean(ecgdatafiltered(r));
        % meanS = mean(ecgdatafiltered(s));
        lowerbound = meanS;
        if (meanQ < meanS || isnan(meanS))
            lowerbound = meanQ;
        end
        normalizationfactor = meanR - lowerbound;
        ecgdatafilterednormalized = ecgdatafiltered ./ normalizationfactor;
        % plot segment of normalized data
        figure; plot(ecgdatafilterednormalized(r(1):r(10))); title(['segidx' num2str(segmentidx) 'sqi ' num2str(F1)])
        %             saveas(gcf,[resultpath '\ECGSegmentidx' num2str(segmentidx) '.png']);
        % only use signals with complete annotations to maintane validity
        % of fourier analysis (may want to replace with interpolation techhnique)
        if (sum(isnan(qrson)) || sum(isnan(qrsoff)))
            disp(['annotations incomplete ' num2str(fileidx)]);
            continue;
        end
        % may want to remove ectopic beats
        qrsonectopicfree = qrson(1);
        rectopicfree = radj(1);
        qrsoffectopicfree = qrsoff(1);
        alignpointsectopicfree = alignpoint(1);
        for ridx = 2:length(radj)-1   % start from idx 2 and goto length(r)-1
            if (ridx <= length(radj)-40)
                meanRR = mean(diff(radj(ridx-1:ridx+38)));   % compute mean over 40 beat segment
            else
                meanRR = mean(diff(radj(end-39:end)));   % compute mean over last 40 beat segment
            end
            preRR = radj(ridx) - radj(ridx-1);    % pre-RR interval
            postRR = radj(ridx+1) - radj(ridx);   % post-RR interval
            % if the pre-RR and post-RR vary from the mean RR by more than
            % 20% consider as ectopic
            if (abs(meanRR-preRR) > 0.15*meanRR || abs(meanRR-postRR) > 0.15*meanRR)
                continue;
            else % include the QRS complex, may also want to include t-wave onset/offset/amplitude/QT-interval
                qrsonectopicfree = [qrsonectopicfree qrson(ridx)];
                rectopicfree = [rectopicfree radj(ridx)];
                qrsoffectopicfree = [qrsoffectopicfree qrsoff(ridx)];
                alignpointsectopicfree = [alignpointsectopicfree alignpoint(ridx)];
            end
        end
        
        % plot ectopic free qrs
        % figure; hold on;
        % for qrsidx = 1:length(qrsonectopicfree)
        %     plot(ecgdatafilterednormalized(qrsonectopicfree(qrsidx):qrsoffectopicfree(qrsidx)));
        % end
        
        hold off; title(['QRS ' num2str(length(qrsonectopicfree)) ' segidx ' num2str(segmentidx)]);
        %             saveas(gcf,[resultpath '\QRSSegmentidx' num2str(segmentidx) ' sqi ' num2str(F1) '.png']);
        
        % calculate MD series
        MDSeries = MDSeriesCalc( ecgdatafilterednormalized, qrsoffectopicfree, qrsonectopicfree);
        MDSeriesLength = length(MDSeries);
        % calculate NTWDseries and TWDseries
        % Isolate complexes (may want to replace this function)
        [NTWDseriesold TWDseriesold TD Complexes] = DiffSeriesCalc( ecgdatafilterednormalized, qrsoffectopicfree, alignpointsectopicfree, qrsonectopicfree);
        % plot complexes
        figure; hold on;
        for qrsidx = 1:size(Complexes, 1)
            plot(Complexes(qrsidx, :));
        end
        title(['aligned qrs' num2str(size(Complexes, 1))]);
        % perform autoccorrelation to remove inconsistent complexes
        % r = xcorr(x,y)
        avgetemplate = mean(Complexes);
        correlation = zeros(1,size(Complexes, 1));
        for qrsidx = 1:size(Complexes,1)
            [r, lags] = xcorr(Complexes(qrsidx, :),avgetemplate);
            correlation(qrsidx) = max(r);
        end
        % plot correlation
        figure; hold on;
        stem(correlation);
        threshold = median(correlation)/2 * ones(size(correlation));
        plot(threshold, 'LineWidth', 2);
        hold off; title(['threshold' ]);
        
        % remove erroneous complexes
        % ComplexesAdjusted = NaN(size(Complexes));
        for qrsidx = 1:size(Complexes,1)
            if (correlation(qrsidx) < threshold)
                Complexes(qrsidx,:) = NaN;  % may cause error, double check
            end
        end
        % remove nans
        originalbeats = size(Complexes,1);
        Complexes(~any(~isnan(Complexes), 2),:)=[];
        cleanbeats = size(Complexes,1);
        % fraction of clean beats
        cleanfraction = (cleanbeats/originalbeats);
        figure(gcf); title(['threshold' num2str(cleanfraction)]);
        %             saveas(gcf,[resultpath '\threshold QRSSegmentalignidx' num2str(segmentidx) ' sqi ' num2str(F1) '.png']);
        % plot clean beats
        figure; hold on;
        for qrsidx = 1:size(Complexes, 1)
            plot(Complexes(qrsidx, :));
        end
        hold off;
        title(['aligned QRS ' num2str(size(Complexes, 1)) ' segidx ' num2str(segmentidx) ' cleanfrac ' num2str(cleanfraction)]);
        %             saveas(gcf,[resultpath '\QRSSegmentalignidx' num2str(segmentidx) ' sqi ' num2str(F1) '.png']);
        % Recompute NTWD and TWD
        [ NTWDseries TWDseries ] = DiffSeriesComplexesCalc( Complexes )
        % plot MD
        figure; stem(MDSeries); title(['MD ' 'sqi ' num2str(F1)]);
        %             saveas(gcf,[resultpath '\MDSegmentidx' num2str(segmentidx) '.png']);
        % plot NTWDseries
        figure; stem(NTWDseries); title(['NTWD ' 'sqi ' num2str(F1)]);
        %             saveas(gcf,[resultpath '\NTWDSegmentidx' num2str(segmentidx) '.png']);
        % plot TWDseries
        figure; stem(TWDseries); title(['TWD ' 'sqi ' num2str(F1)]);
        %             saveas(gcf,[resultpath '\TWDSegmentidx' num2str(segmentidx) '.png']);
        nfft = 256; % standardize the length of the fft to 256
        %             pxx = pwelch(MDSeries,[],[],nfft);
        pxx = pwelch(NTWDseries,[],[],nfft);
        %             pxx = pwelch(TWDseries,[],[],nfft);
        N = nfft;
        stepsize = 0.5/length(pxx);
        
        % bquency = (1:129)*stepsize;
        
        dl = 2; du = 7;
        indexeverytwobeats = round(N/dl, 0);    % corresponding to every 2 beats
        indexeverysevenbeats = floor(N/du); % corresponding to every 7 beats
        figure; plot((1:129)*stepsize, log10(pxx)); title(['pwelch 5 min' 'sqi ' num2str(F1)]); hold on;
        scatter(indexeverysevenbeats*stepsize, log10(pxx(indexeverysevenbeats)));
        scatter(indexeverytwobeats*stepsize, log10(pxx(indexeverytwobeats)));
        hold off;
        % compute energy in every 2-7 beats region
        energyinband = trapz((indexeverysevenbeats:indexeverytwobeats) * stepsize, pxx(indexeverysevenbeats:indexeverytwobeats));
        energyinbandcollection(segmentidx) = energyinband;
        xlabel(['energy = ' num2str(energyinband)]); hold off;
        %             saveas(gcf,[resultpath '\SpectrumSegmentidx' num2str(segmentidx) '.png']);
        close all;
    end
end
end