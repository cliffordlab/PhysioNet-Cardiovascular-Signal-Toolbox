function Align = Find_Fiducials(q, r, s, Align, ecg, Param)
%OVERVIEW, Finds the lead with maximum QRS amplitude, takes a beat as QRS template and finds fiducial points such that
% aligned beats have maximum cross-correlation (over Q to S or S to T interval) in the chosen
% lead against the template.
%
% INPUTS            MANDATORY           OUTPUTS
%                   q                   q point locations
% 
%                   r                   r point locations
%
%                   s                   s point locations
% 
%                   Align               Alignment information for beats in
%                                       analysis window.
%
%                   ecg                 N by M array of ecg data. Contains M channels of ecg with N data points each.
%
%                   Param               Structure containing the parameters
%                                       used for evaluating TWAs.
%
% OUTPUTS           
%                   Align               Alignment information used to
%                                       determine if current window is clean enough for
%                                       analysis.
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

q = q(1:min(length(q), length(s)));
s = s(1:min(length(q), length(s)));

% Sort leads by amplitude.
leadSorted = Sort_Leads(q, s, Align.st, ecg, Param);

for ilead = 1:length(leadSorted)
    Align.lead = leadSorted(ilead);

    % Get the orientation of the QRS complex, max positive or max negative deflection.
    Align.orientation = MaxFrontOrientation(ecg(:, Align.lead), q, s);

    Align.fidBase = r;

    Align.f2s = floor(median(s - Align.fidBase) + Param.stAdjIntv);   % reserve to avoid including part of QRS in ST interval computations
    Align.q2f = min([floor(median(Align.fidBase - q)),Align.f2s-Param.stAdjIntv]);

    Align.fid = [];
    NFailed = 0;

    Align.template = 0;
    while isempty(Align.fid)
        Align.template = Align.template + 1;
        % Get alignment information
        Align = Adjust_e(Align, ecg, 1000, Param);
        if isempty(Align.fid)

            NFailed = NFailed + 1;
            if (NFailed > 0.1 * length(s))
                break
            end
        end
    end
    if (~isempty(Align.fid))
        for lead=1:size(ecg, 2)
            for i = 1:length(Align.fid)
                Align.amp(i, lead) = (max(ecg(q(i) : s(i), lead))+min(ecg(q(i) : s(i), lead)))/2;  % align to zero mean on QS segment
            end
        end
        % get correlation for each beat with template over Q-S and S-T
        % interval
        Align = CheckLead_Correlations(Align, ecg, Param);
        disp(['FindFiducials: Alignment in lead ' num2str(Align.lead)]);
        break;
    end;
    disp(['FindFiducials: lead ' num2str(Align.lead) ' is too noisy or ectopic..']);
end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function leads = Sort_Leads(q, s, st, ecg, Param)

NLead = size(ecg, 2);

%   QRS or T amplitudes in each lead
ptp = zeros(min(length(q), 8),NLead);
for i = 1:min(length(q), 8)
    for l = 1:NLead
        if strcmp(Param.Alignment, 'st')
            ptp(i, l) = max(ecg(s(i):s(i) + st, l)) - min(ecg(s(i):s(i) + st, l));
        else
            ptp(i, l) = max(ecg(q(i):s(i), l)) - min(ecg(q(i):s(i), l));
        end;
    end;
end;
ptpmed = zeros(1,NLead);
for i = 1:NLead
    ptpmed(i) = median(ptp(:, i));
end;

[sorted indices] = sort(ptpmed);
leads = indices(length(indices):-1:1);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Align = Adjust_e(Align, ecg, freq, Param)
% OVERVIEW, This function determines which beats in analysis window for
% the lead given by Align.lead are valid for analysis, based on correlation
% with a QRS and S-T offset segment template.
Align.fidBase(Align.template);

len_ecg = size(ecg,1);
qs_templ = ecg(max(Align.fidBase(Align.template) - Align.q2f,1):min(Align.fidBase(Align.template) + Align.f2s,len_ecg), Align.lead);
st_templ = ecg(max(Align.fidBase(Align.template) + Align.f2s,1):min(Align.fidBase(Align.template) + Align.f2s + Align.st,len_ecg), Align.lead);

[a b c d] = Adjust_Fiducials(ecg, Align.fidBase, Align.q2f, Align.f2s, qs_templ, Align.st, st_templ, Param);
% fid a -> with respect to QRS alignment
% fid c -> with respect to ST alignment

if ~isempty(c)
    Align.fidQRS = a;
    Align.QRScorr(:, Align.lead) = b;
    Align.fid = c;
    Align.Tcorr(:, Align.lead) = d;
    amp = (max(qs_templ)+min(qs_templ))/2;
    Align.qs_templ(:, Align.lead) = qs_templ - amp;
    Align.st_templ(:, Align.lead) = st_templ - amp;

    for i = 1:length(Align.fid)
        Align.valid(i, Align.lead) = (Align.QRScorr(i, Align.lead) >= Param.corrQRS) && (Align.Tcorr(i, Align.lead) >= Param.corrT);
    end;
    Align.validleads(Align.lead) = true;
end;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function orientation = MaxFrontOrientation(ecg, q, s)

NPos = 0;

for i = 1:min(20, length(s))
    [ma maI] = max(ecg(q(i):s(i)));
    [mi miI] = min(ecg(q(i):s(i)));

    NPos = NPos + (maI > miI);
end;

orientation = (NPos > 10);

return;
