function [ MDseries ] = MDSeriesCalc( signalnormalized, QRSoffprocessed, QRSonprocessed)
%MDSeriesCalc: Calculates the morphological variability between QRS
%complexes for an ECG sequence.
%
%   Inputs:
%   signalnormalized - ecg signal normalized by the mean R amplitude
%   QRSoffprocessed - indices corresponding to QRS-offset
%   QRSonprocessed - indices corresponding to QRS-onset
%
%   Outputs:
%   MDseries - the sum of squared differences between morphology calculated
%   between succesive qrs complexes in ecg

noofconsecutivecomplexes = length(QRSoffprocessed);
maxinterval = max(diff(QRSoffprocessed(1:noofconsecutivecomplexes) - QRSonprocessed(1:noofconsecutivecomplexes)));

% need to correct
QRScomplexes(1:noofconsecutivecomplexes-1, 1:maxinterval) = 0;
QRScomplexlengths(1:noofconsecutivecomplexes-1) = 0;
%     complex(1:noofconsecutivecomplexes-1, 1:maxinterval) = 0;
%     complexlengths(1:noofconsecutivecomplexes-1) = 0;

%for index = 1:(noofconsecutivecomplexes-1)
for index = 1:(noofconsecutivecomplexes-1)
    currentcomplexlength = QRSoffprocessed(index) - QRSonprocessed(index);
    complex(index, 1:currentcomplexlength) = signalnormalized(QRSonprocessed(index):QRSoffprocessed(index)-1);
    complexlengths(index) = currentcomplexlength;
end

%     figure; hold on;
%     for index = 1:(noofconsecutivecomplexes-1)
%         plot(complex(index, 1:complexlengths(index)));
%     end
%     hold off;

%% DWT and MD series

noofcomplexes = noofconsecutivecomplexes;
MDseries(1:(noofcomplexes-1)) = 0;
noofcomplexes = length(QRSoffprocessed);
complex1 = []; complex2 = [];

% MDseries(1:noofcomplexes-1) = 0;
%     figure; hold on;
for index = 1:(noofcomplexes-1-1)
    
    %         if (find([10 31 63 155 195] == index) > 0)
    if (index == 150)
        wait = 1;
    end
    
    % index
    complex1 = complex(index, 1:complexlengths(index));
    complex2 = complex(index+1, 1:complexlengths(index+1));
    
    % %         % testing
    %         figure;
    %         plot(complex1); hold on;
    %         plot(complex2); hold off;
    
    % dtw
    [dist,ix,iy] = dtw(complex1, complex2);
    templt1 = complex1(ix);
    templt2 = complex2(iy);
    
    %         plot(templt1);
    %         plot(templt2);
    
    
    %         %         % testing
    %                 figure;
    %                 plot(templt1); hold on;
    %                 plot(templt2); hold off;
    
    MD = sum((templt1 - templt2).^2);
    MD = MD ./ length(templt1);
    
    MDseries(index) = MD;
    
    %         if (MD > 6e-4)
    %            wait = 0;
    %         end
    
    %                 close all;
end
%     hold off;

end

