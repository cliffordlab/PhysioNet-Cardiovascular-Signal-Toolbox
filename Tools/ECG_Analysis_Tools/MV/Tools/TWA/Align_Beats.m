function [Align HR]= Align_Beats(beats, q, r, s, STLen, ecg, freq, Param)
%OVERVIEW, This function aligns beats with a QRS template to ensure the
% ecg data is clean enough for analysis.
%
% INPUT         MANDATORY           DESCRIPTION
%
%               ecg                 N by M array of ecg data. Contains M channels of ecg data with N datapoints each.
%
%               beats               Number of analysis beats.
%
%               q                   Array of location of q points.
%
%               r                   Array of location of r points.
% 
%               s                   Array of location of s points.
% 
%               STlen               ST interval segment length.
% 
%               freq                Sampling frequency.
%
%               Param               Structure contains parameter
%                                   definitions used for computing TWAs.
%
% OUTPUT:
% 
%               Align                       Structure contains alignment and correaltion information of beats in current analysis window. 
%                                           The structure containing the fields:
%               fid                         (1 x num_of_beats): fiducial base marker after aligning to maximize correlation with S-T segment. 
%               amp                         (num_of_beats x num_of_leads): average of max and min amplitude for each qrs complex.
%               q2f (scaler):               length of Q-R segment
%               f2s (scaler):               length of R-S segment
%               valid (1 x num_of_beats):   1 if the beat is valid and 0 otherwise
%               st :                        estimate of median S-T offset segment for analysis window.                         
%               lead :                      lead used to estimate q-s s-t intervals for aligning.
%               orientation :               1 indicates the QRS complex has a maximum positive deflection.
%               fidBase :                   R point locations in the analysis window.
%               template :                  flag, set to 1 after the fid and fidQRS markers 
%                                           have been adjusted to maximize correlation between qrs, qrs
%                                           template and s-t, s-t template.
%               fidQRS :                    fiducial base points adjusted after maximizing correlation between qrs and ars template.
%               QRScorr :                   Correlation between each QRS complex and the QRS template in analysis window.
%               Tcorr :                     Correlation between each S-T offset segment and the S-T offset template.
%               qs_templ :                  QRS template for analysis window.
%               st_templ :                  S-T offset template for analysis window.
%               validleads :                list of noise free leads valid for analysis.
%
%               HR                          Heart rate of analysis window
%
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Shamim Nemati, 
%       editted by Ismail Sadiq on 10/26/2019. 
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.   

HR=NaN;
Align.fid = [];
Align.valid = [];
Align.st = STLen;

% Get alignment information
Align = Find_Fiducials(q, r, s, Align, ecg, Param);
if (isempty(Align.fid))
    disp('no lead is suitable for alignment, finishing..');
    return;
end

% Determine number of invalid or noisy beats in the window 
for lead = 1:size(ecg, 2)
    NInvalid = length(find(1 - Align.valid(:, lead)));
    if NInvalid
        s = sprintf('AlignBeats: %d invalid beats in lead %d', NInvalid, lead);
        disp(s)
    end
end

% If number of valid beats greater than analysis window length, limit to analysis window length. 
if (length(Align.fid) > beats)
    Align.fid = Align.fid(1:beats);
end

% compute HR for window
if (nargout>1)
    HR = median(60./medfilt1(diff(Align.fid/freq),5)) ;
end

end
