function [TWAresult] = ComputeTWA(ecg,ann,fs)
%OVERVIEW, This is main function which calls the Analyze_TWA_segment function
% which measures TWA using the statistical reshuffling method.
%
% INPUT:        MANDATORY:      DESCRIPTION:
%               ecg             N by M (row by col) array of ecg data with each of M columns corresponding to a single channel of ecg.
%                               Each channel has N datapoints.
%
%               ann             the annotations structure contains Q,R,S,T
%                               fiducial points.
%
%   REF:
%   This script is based on the algorithm presented in the following paper,
%   Nemati S, Abdala O, Monasterio V, Yim-Yeh S, Malhotra A, Clifford GD.
%   A nonparametric surrogate-based test of significance for T-wave alternans detection.
%   IEEE Trans Biomed Eng. 2011;58(5):1356–1364. doi:10.1109/TBME.2010.2047859
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

if (fs ~= 1000)
    disp('The ECG needs to be sampled at 1000 Hz in order to perform analysis ... exiting analysis')
    return
else
    [TWAresult] = Analyze_TWA_segment(ecg, 1000, ann);
end

end

