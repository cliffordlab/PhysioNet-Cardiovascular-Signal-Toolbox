function Align = CheckLead_Correlations(Align, ecg, Param)
% OVERVIEW, This function determines which beats in all leads being analyzed, are valid for analysis based 
%on correlation with the QRS template and S-T offset segment template.
% 
% INPUTS        MANDATORY           DESCRIPTION
%               Align               Alignment information for current
%                                   analysis window.
%
%               ecg                 N by M array of ecg data. Contains M
%                                   channels of ecg with N data points.
%
%               Param               Structure which contains parameters
%                                   used for computing TWAs.
% OUTPUTS
%               Align               Alignment information for current
%                                   analysis window. Contains information
%                                   on whether beat valid for analysis.
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


for lead = 1:size(ecg, 2)
    if lead == Align.lead
        continue;
    end;
    % checking correlations vs template beat makes us very vulnerable to a
    % noise within this particular beat (which is possible in a lead other than alignment
    % lead), thus we better check correlations vs average
    qs_templ = zeros(1, Align.q2f + Align.f2s + 1)';
    if strcmp(Param.Alignment, 'st')
        st_templ = zeros(1, Align.st + 1)';
    end;
    for i = 1:length(Align.fidQRS)
        qs_templ = qs_templ + (ecg((Align.fidQRS(i) - Align.q2f):(Align.fidQRS(i) + Align.f2s), lead) - Align.amp(i, lead));
        if strcmp(Param.Alignment, 'st')
            st_templ = st_templ + (ecg((Align.fid(i) + Align.f2s):(Align.fid(i) + Align.f2s + Align.st), lead) - Align.amp(i, lead));
        end;
    end;
    Align.qs_templ(:, lead) = qs_templ / length(Align.fidQRS);
    Align.st_templ(:, lead) = st_templ / length(Align.fidQRS);
    
    for i = 1:length(Align.fid)
        qs_ecg = ecg((Align.fidQRS(i) - Align.q2f):(Align.fidQRS(i) + Align.f2s), lead);
        a = corrcoef(qs_templ, qs_ecg);
        Align.QRScorr(i, lead) = a(1, 2);
        if (a(1, 2) < Param.corrQRS)
            Align.valid(i, lead) = false;
            continue;
        end;
    
        if strcmp(Param.Alignment, 'st')
            st_ecg = ecg((Align.fid(i) + Align.f2s):(Align.fid(i) + Align.f2s + Align.st), lead);
            a = corrcoef(st_templ, st_ecg);
            Align.Tcorr(i, lead) = a(1, 2);
            if (a(1, 2) < Param.corrT)
                Align.valid(i, lead) = false;
                continue;
            end;
        end;    
        Align.valid(i, lead) = true;        
    end;
    
    Align.validleads(lead) = (length(find(Align.valid(:, lead))) >= 0.9 * length(Align.fid));
    
end;

return;