function [ NTWDseries TWDseries TD complexes] = DiffSeriesCalc( signalnormalized, QRSoff, R, QRSon)
%MDSeriesCalc: Calculates the morphological variability between QRS
%complexes for an ECG sequence.
%
%   Inputs:
%   signalnormalized - ecg signal normalized by the mean R amplitude
%   QRSoff - indices corresponding to QRS-offset
%   R - indices corresponding to R-points
%   QRSon - indices corresponding to QRS-onset
%
%   Outputs:
%   NTWDseries - the sum of squared differences between morphology calculated
%   between succesive qrs complexes in ecg aligned at the r points
%   TWDseries - the sum of squared differences between morphology calculated
%   between time warped succesive qrs complexes in ecg aligned at the r points


noofconsecutivecomplexes = length(QRSoff);
%     maxinterval = max(diff(QRSoffprocessed(1:noofconsecutivecomplexes) - QRSonprocessed(1:noofconsecutivecomplexes)));
onsets = QRSon;
offsets = QRSoff;
rpoints = R;

medianonset = round(median(R - QRSon));
medianoffset = round(median(QRSoff - R));



%     QRScomplexes(1:noofconsecutivecomplexes-1, 1:maxinterval) = 0;
%     QRScomplexlengths(1:noofconsecutivecomplexes-1) = 0;
%     complex(1:noofconsecutivecomplexes-1, 1:maxinterval) = 0;
%     complexlengths(1:noofconsecutivecomplexes-1) = 0;

%     %for index = 1:(noofconsecutivecomplexes-1)
%     for index = 1:(noofconsecutivecomplexes-1)
%         currentcomplexlength = QRSoffprocessed(index) - QRSonprocessed(index);
%         complex(index, 1:currentcomplexlength) = signalnormalized(QRSonprocessed(index):QRSoffprocessed(index)-1);
%         complexlengths(index) = currentcomplexlength;
%     end

%     figure; hold on;
%     for index = 1:(noofconsecutivecomplexes-1)
%         plot(complex(index, 1:complexlengths(index)));
%     end
%     hold off;

%% DWT and MD series
noofcomplexes = length(QRSoff);
NTWDseries(1:(noofcomplexes-1)) = 0;
complex1 = []; complex2 = [];
complexes = zeros(1,medianonset+medianoffset+1);
for index = 1:(noofcomplexes-1-1)
    
    index
    
    if (index == 334)
        wait = 1;
    end
    %         precomplex1 = complex(index, 1:complexlengths(index));
    %         precomplex2 = complex(index+1, 1:complexlengths(index+1));
    %         figure; hold on;
    %         plot(precomplex1); plot(precomplex2);
    %         hold off; title('Before processing');
    
    complex1 = signalnormalized(rpoints(index) - medianonset:rpoints(index) + medianoffset);
    complex2 = signalnormalized(rpoints(index+1) - medianonset:rpoints(index+1) + medianoffset);
    complexes(index,:) = complex1;
    
    SD = sum((complex1 - complex2).^2);
    SD = SD ./ length(complex1);
    NTWDseries(index) = SD;
    
    % testing
    %         figure; title('NTWDseries');
    %         plot(complex1); hold on;
    %         plot(complex2); hold off;
    
    % dtw
    [dist,ix,iy] = dtw(complex1, complex2);
    templt1 = complex1(ix);
    templt2 = complex2(iy);
    TD(index) = dist;
    
    % testing
    %         figure; title('TWDseries');
    %         plot(templt1); hold on;
    %         plot(templt2); hold off;
    
    MD = sum((templt1 - templt2).^2);
    MD = MD ./ length(templt1);
    
    TWDseries(index) = MD;
    
    %         close all;
end

end

