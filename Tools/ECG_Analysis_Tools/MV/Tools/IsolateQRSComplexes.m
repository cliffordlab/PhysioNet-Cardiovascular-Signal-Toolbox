function [Complexes] = IsolateQRSComplexes(ecg, QRSoff, R, QRSon)

% [Complexes] = IsolateQRSComplexes(ecg, QRSoff, R, QRSon)
%   OVERVIEW:   Using the annotations in QRSon, R, QRSoff the QRS complexes
%               are segmented over a median QRS onset to QRS offset
%               interval in the ECG.
%
%   INPUT:      MANDATORY:
%               ecg             : ecg signal as input
%               QRSoff          : QRS offset for each QRS complex in ecg
%                               
%               R               : R peak for each QRS complex in ecg 
%                                 
%                                
%               QRSon           : QRS onset for each QRS complex in ecg
%                                                
%                
%   OUTPUT:     
%               Complexes       : 2D array with dimension (number of complexes by median QRS complex interval).
%                               Stores the QRS complexes in the ecg.  
%                              
%
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

medianonset = round(median(R - QRSon));
medianoffset = round(median(QRSoff - R));

noofcomplexes = length(QRSoff);
Complexes = zeros(noofcomplexes, medianonset+medianoffset+1);
for index = 1:(noofcomplexes)
    Complexes(index,:) = ecg(R(index) - medianonset:R(index) + medianoffset);
end

end

