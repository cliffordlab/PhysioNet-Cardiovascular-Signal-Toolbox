function [res,q,r,s,qs,stend,ecg] = get_ann(ecg, beat, Param, freq)%, exp)
% OVERVIEW, This function reads the ecg data and inverts it in case their
% is an inverted QRS complex. The annotations are read and stored in the q,
% r,s,qs and stend variables. If annotations are missing they are replaced
% using estimates from the other annotations.
%
% INPUT         MANDATORY       DESCRIPTION
%
%               ecg (uV)        N by M array of ecg data. Contains M
%                               channels of data with N data points each.
%
%               beat            Structure contains the Q,R,S,T-offset fiducal
%                               point locations.
%
%               Param           Structure contains the parameters used in
%                               computing TWA.
%
%               freq (Hz)       Sampling frequncy of the ecg, needs to be
%                               1000 Hz.
%
% OUTPUT        res             flag, if set to 1 indicates the annotations
%                               have been successfully assigned to the q,r,s variables and related intervals have been computed.
%
%               q               array of q fiducial points for ecg.
%
%               r               array of r fiducial points for ecg.
%
%               s               array of s fiducial points for ecg.
%
%               qs              array of qs intervals computed from q and s fiducial
%                               points.
%
%               stend           array of s-t end intervals computed from the s and t
%                               off set fiducal points.
%
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


res = 0;

len = size(ecg,1);min(ecg(:,1));max(ecg(:,1));

% if the QRS is inverted multiply with -1
for lead=1:size(ecg,2)
    [lu_q1] = quantile(ecg(ceil(0.2*len:.3*len),lead),[0.05 , 0.95]);
    [lu_q2] = quantile(ecg(ceil(0.8*len:.9*len),lead),[0.05 , 0.95]);
    if (abs(lu_q1(1))>lu_q1(2) && abs(lu_q2(1))>lu_q2(2))
        ecg(:,lead) = -1*ecg(:,lead);
    end
end
%
try
    % Ectopic detection in R wave
    if (numel(find(isnan(beat.R))) > 0.2*length(beat.R))
        return
    end
    [tm,ids,LOG]=incidencias(beat.R/freq);
    [normales,itm,itorg]=intersect(tm,beat.R/freq);
catch
    warning('Function get_ann failed ...')
    return;
end

% assign q,r,s,toff fiducial points.
if (isempty(itorg)) return;
else
    q = nan(size(beat.Q)); r = q; s = q; Toff = q;
end
q(itorg)=beat.Q(itorg);
r(itorg)=beat.R(itorg);
s(itorg)=beat.S(itorg);
Toff(itorg) = beat.Toff(itorg);

% estimate R peaks for missing points
ind_valid = find(~isnan(r));
r(ind_valid) = improv_fid(ecg(:,1),r(ind_valid),freq,'R',[],Param.stAdjIntv*2);

auxS=size(find(~isnan(s)),2)/size(find(~isnan(r)),2);
auxQ=size(find(~isnan(q)),2)/size(find(~isnan(r)),2);
% S wave detection based on wavelet detector (aux>0.95) or minimum value of
% the ECG in a time window, of length 0.1 seconds, beginning at the peak of
% the R-wave (aux<0.95)[Mason's thesis]
if auxS<0.95
    valid_beat=find(~isnan(r));
    SS = calculo_S(ecg,r(valid_beat),freq); %If S-detection was poor, use R peaks to redetect S points
    s(valid_beat) = SS;
end
if auxQ<0.95
    valid_beat=find(~isnan(r));
    QQ = calculo_Q(ecg,r(valid_beat),freq);
    q(valid_beat) = QQ;
end

% estimate qs and stend interval
qs = s-q;
stend = Toff-s;

ind_invalid = find(isnan(q) | isnan(r) | isnan(s));
i=1;
while (i <= length(ind_invalid))
    %if 2 consecutive invalid indices remove both
    if (i~=length(ind_invalid) & ind_invalid(i+1) == ind_invalid(i)+1)
        q(ind_invalid(i:i+1))=[];
        r(ind_invalid(i:i+1))=[];
        s(ind_invalid(i:i+1))=[];
        Toff(ind_invalid(i:i+1))=[];
        qs(ind_invalid(i:i+1))=[];
        stend(ind_invalid(i:i+1))=[];
        ind_invalid=ind_invalid-2;
        i=i+2;
        %if 1 non-consecutive invalid index take the average of previous and next fiducials
    else
        if (ind_invalid(i)==1 | ind_invalid(i)==length(q))
            q(ind_invalid(i))=[]; r(ind_invalid(i))=[]; s(ind_invalid(i))=[];
            Toff(ind_invalid(i))=[]; qs(ind_invalid(i))=[]; stend(ind_invalid(i))=[];
            ind_invalid=ind_invalid-1;
            i=i+1;
        else
            r(ind_invalid(i)) = round(0.5*(r(ind_invalid(i)-1)+r(ind_invalid(i)+1)));
            r(ind_invalid(i)) = improv_fid(ecg(:,1),r(ind_invalid(i)),freq,'R',[],Param.stAdjIntv*2);
            q(ind_invalid(i)) = r(ind_invalid(i))-2*Param.stAdjIntv;
            s(ind_invalid(i)) = r(ind_invalid(i))+2*Param.stAdjIntv;
            Toff(ind_invalid(i)) = round(0.5*(Toff(ind_invalid(i)-1)+Toff(ind_invalid(i)+1)));
            qs(ind_invalid(i)) = round(0.5*(qs(ind_invalid(i)-1)+qs(ind_invalid(i)+1)));
            stend(ind_invalid(i)) = round(0.5*(stend(ind_invalid(i)-1)+stend(ind_invalid(i)+1)));
            i=i+1;
        end
    end
end

% improve fiducial estimates of the fiducial points.
r = improv_fid(ecg(:,1),r,freq,'R',[],Param.stAdjIntv*2);
s = improv_fid(ecg(:,1),s,freq,'S',r,Param.stAdjIntv/2);
q = improv_fid(ecg(:,1),q,freq,'Q',r,Param.stAdjIntv);

res=1;
%plot(ecg(:,1)),hold on,plot(r,ecg(r,1),'r+'),hold on,plot(s,ecg(s,1),'co'),hold on,plot(q,ecg(q,1),'kx'); % check accuracy of fiducial points by plotting.
end




