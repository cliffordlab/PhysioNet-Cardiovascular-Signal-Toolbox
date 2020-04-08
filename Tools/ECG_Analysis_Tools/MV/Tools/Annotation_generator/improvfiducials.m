function [detection_flg, beato] = improvfiducials(beat,freq,ecg)

% OVERVIEW, This function improves the fiducial point detection for the
% Q,R,S,Toff points in the beat variable. The improved fiducial points are
% stored in the beato variable.
%
% INPUTS    MANDATORY       DESCRIPTION
%
%           beat            Structure variable which must contain the following
%                           fields.
%           .QRSon          Location for QRS onset for each beat.
%           .Q              Location for Q for each beat.
%           .R              Location for R for each beat.
%           .S              Location for S for each beat.
%           .QRSoff         Location for QRS offset for each beat.
%           .Toff           Location for T offset for each beat.
%
%           freq            Sampling frequncy of the ecg variable.
%           ecg             Single row of ecg data. Fiducial points for
%                           this data are provided in the beat variable.
% OUTPUTS
%           beato           Structure variable which contains the following
%                           fields.
%           .QRSon          Improved detection for the QRS onset.
%           .Q              Improved detection for the Q point.
%           .R              Improved detection for the R point.
%           .S              Improved detection for the S point.
%           ,QRSoff         Improved detection for the QRS offset.
%           .Toff           Improved detection for the T offset.
%
% NOTE: Other fiducial points not listed are detected but are not required for MVM
% and TWA computation.
% 
% REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Shamim Nemati
%	editted by Ismail Sadiq on 10/26/2019.
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  


%s=[];r=[]; q=[]; qs = []; stend = []; Toff=[]; qlatest = -1;

% if the QRS is inverted multiply with -1
% for lead=1:size(ecg,2)
%    [lu_q1] = quantile(ecg(ceil(0.2*len:.3*len),lead),[0.05 , 0.95]);
%    [lu_q2] = quantile(ecg(ceil(0.8*len:.9*len),lead),[0.05 , 0.95]);
%    if (abs(lu_q1(1))>lu_q1(2) && abs(lu_q2(1))>lu_q2(2))
%        ecg(:,lead) = -1*ecg(:,lead);
%    end
% end
%
detection_flg = 0;
% Ectopic detection in R wave
if (numel(find(isnan(beat.R))) > 0.2*length(beat.R))
    return
end

q = beat.Q;
r = beat.R;
s = beat.S;
Toff = beat.Toff;
qrson = beat.QRSon;
qrsoff = beat.QRSoff;

ind_valid = find(~isnan(r));
r(ind_valid) = improvfid(ecg',r(ind_valid),freq,'R',[],floor(40 * freq / 1000)*2);

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


ind_invalid = find(isnan(q) | isnan(r) | isnan(s) | isnan(qrson) | isnan(qrsoff));
i=1;
while (i <= length(ind_invalid))
    %if 2 consecutive invalid indices remove both
    if (i~=length(ind_invalid) & ind_invalid(i+1) == ind_invalid(i)+1)
        q(ind_invalid(i:i+1))=[];
        r(ind_invalid(i:i+1))=[];
        s(ind_invalid(i:i+1))=[];
        Toff(ind_invalid(i:i+1))=[];
        qrson(ind_invalid(i:i+1))=[];
        qrsoff(ind_invalid(i:i+1))=[];
        ind_invalid=ind_invalid-2;
        i=i+2;
        %if 1 non-consecutive invalid index take the average of previous and next fiducials
    else
        if (ind_invalid(i)==1 | ind_invalid(i)==length(q))
            q(ind_invalid(i))=[]; r(ind_invalid(i))=[]; s(ind_invalid(i))=[];
            Toff(ind_invalid(i))=[]; qrson(ind_invalid(i))=[]; qrsoff(ind_invalid(i))=[];
            ind_invalid=ind_invalid-1;
            i=i+1;
        else
            r(ind_invalid(i)) = round(0.5*(r(ind_invalid(i)-1)+r(ind_invalid(i)+1)));
            r(ind_invalid(i)) = improvfid(ecg',r(ind_invalid(i)),freq,'R',[],floor(40 * freq / 1000)*2);%(:,1),r(ind_invalid(i)),freq,'R',[],Param.stAdjIntv*2);
            q(ind_invalid(i)) = r(ind_invalid(i))-2*floor(40 * freq / 1000);
            s(ind_invalid(i)) = r(ind_invalid(i))+2*floor(40 * freq / 1000);
            Toff(ind_invalid(i)) = round(0.5*(Toff(ind_invalid(i)-1)+Toff(ind_invalid(i)+1)));
            qrson(ind_invalid(i)) = round(0.5*(qrson(ind_invalid(i)-1)+qrson(ind_invalid(i)+1)));
            qrsoff(ind_invalid(i)) = round(0.5*(qrsoff(ind_invalid(i)-1)+qrsoff(ind_invalid(i)+1)));
            i=i+1;
        end
    end
end

r = improvfid(ecg',r,freq,'R',[],floor(40 * freq / 1000)*2);
s = improvfid(ecg',s,freq,'S',r,floor(40 * freq / 1000)/2);
q = improvfid(ecg',q,freq,'Q',r,floor(40 * freq / 1000));
% improvfid t offset
t_nan = find(isnan(Toff));
median_qrsoff_toff_intv = nanmedian(Toff - qrsoff); % est median qrsoff toff interval 
for ii = 1:length(t_nan)
    if (qrsoff(t_nan(ii)) + median_qrsoff_toff_intv <= length(ecg))
        Toff(t_nan(ii)) = qrsoff(t_nan(ii)) + median_qrsoff_toff_intv;
    end
end

beato.R = r; beato.Q = q; beato.S = s; beato.QRSon = qrson; beato.QRSoff = qrsoff; beato.Toff = Toff;
detection_flg = 1;
%plot(ecg),hold on,plot(r,ecg(r),'r+'),hold on,plot(s,ecg(s),'co'),hold on,plot(q,ecg(q),'kx'), plot(Toff,ecg(Toff),'g*'); plot(beato.QRSon,ecg(beato.QRSon),'yx'); plot(beato.QRSoff,ecg(beato.QRSoff),'yx'); hold off
return
