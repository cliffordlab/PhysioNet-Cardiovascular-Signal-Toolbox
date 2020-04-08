function [fidQRS, QRScorr, fidST, Tcorr] = Adjust_Fiducials(ecg, fidBase, q2f, f2s, qs_templ, st, st_templ, Param)
% OVERVIEW, The function adjusts the fiducial base marker for each
% beat so that the correlation between the QRS template and ST template is
% maximium.
%
% INPUTS        MANDATORY       DESCRIPTION
%               ecg             N by M array of ecg data. Contains M
%                               channels of ecg with N data points each.
%
%               fidBase         The location of R peaks for an upright QRS complex.
%
%               q2f             The number of samples between q point and fiducial base in the qrs complex template.
%
%               f2s             The number of samples between fiducial base
%                               and the s point in the qrs complex template.
%
%               qs_templ        qrs complex template.
%
%               st              s-t segment length estimate.
%
%               st_templ        s-t segment template.
%
%               Param           Structure containing parameters using for
%                               computing TWAs.
%
% OUTPUTS
%               fidQRS          The adjusted fiducial base so that the
%                               correlation between the QRS for each beat and template
%                               QRS is maximum.
%
%               QRScorr         The maximum correlation between QRS complex
%                               for each beat and template QRS.
%
%               fidST           The adjusted fiducial base so that the
%                               correlation between the S-T segment for each beat and template S-T segment is maximum.
%
%               Tcorr           The maximum correlation between S-T segment
%                               for each beat and the template S-T segment.
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


NInapr = 0;

qoff = -q2f - Param.stAdjIntv - 1;
soff = f2s - Param.stAdjIntv - 1;
toff = f2s + st - Param.stAdjIntv - 1;
Tcorr=[];
warning off MATLAB:divideByZero     % to supress corrcoef warning when data is constant, in that case corrcoef is NaN and max index is any

for i = 1:length(fidBase)
    cc = zeros(1, 4.5 * Param.stAdjIntv + 1);
    ccc = zeros(1, 4.5 * Param.stAdjIntv + 1);
    for j = 1:(4.5 * Param.stAdjIntv + 1)
        if((fidBase(i) + j + qoff)<1), continue; end
        qs_ecg = ecg((fidBase(i) + j + qoff):(fidBase(i) + j + soff));
        if(length(qs_templ)~=length(qs_ecg)), continue;
        else
            a = corrcoef(qs_templ, qs_ecg);
            cc(j) = a(1, 2);
        end
        
        if strcmp(Param.Alignment, 'st')
            if((fidBase(i) + j + soff)<1), continue;end
            st_ecg = ecg((fidBase(i) + j + soff):(fidBase(i) + j + toff));
            if (length(st_templ)~=length(st_ecg)), continue;
            else
                b = corrcoef(st_templ, st_ecg);
                ccc(j) = b(1, 2);
            end
        end;
    end;
    [m ind] = max(cc);
    indQRS = ind;
    if (strcmp(Param.Alignment, 'st'))
        [mm ind] = max(ccc);
    end;
    
    if (m < Param.corrQRS || (strcmp(Param.Alignment, 'st') && mm < Param.corrT))
        NInapr = NInapr + 1;
    end;
    
    if (NInapr > 0.1 * length(fidBase))
        fidST = [];
        fidQRS = [];
        QRScorr = [];
        Tcorr = [];
        break;
    else
        fidST(i) = fidBase(i) + ind - Param.stAdjIntv - 1; % with respect to ST alignment
        fidQRS(i) = fidBase(i) + indQRS - Param.stAdjIntv - 1; % with respect to QRS alignment
        QRScorr(i) = m;
        if (strcmp(Param.Alignment, 'st'))
            Tcorr(i) = mm;
        end;
        
    end;
end;

return;
