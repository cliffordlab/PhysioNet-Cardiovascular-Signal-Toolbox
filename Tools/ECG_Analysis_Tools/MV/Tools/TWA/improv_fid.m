function [fidout] = improv_fid(ecg,fidin,fs,fidtype,fidhelp,N2)
% OVERVIEW, This function improves the fiducal point detections for Q, R, S
% and T points.
% 
% INPUTS        MANDATORY       DESCRIPTION
%
%               ecg             N by M array of ecg data, with M channels
%                               of ecg and N data points in each channel.
%
%               fidin           Fiducial points for improving. e.g. S, Q,..
%
%               fs (Hz)         Sampling frequency of ecg, needs to be 1000
%                               Hz.
%
%               fidtype         char type specifying type of fiducial
%                               point. E.g. 'S', 'Q'.
%
%               fidhelp         R fiducial point used to locate or improve
%                               S or Q detection.
%
%               N2              Window size in samples used for searching for or improving Q or S detection.
%
% OUTPUTS
%
%               fidout          Improved detection for the fiducial points
%                               specified in fidin.
%   REPO:
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

fidout=zeros(1,size(fidin,2));

if (nargin<5 || isempty(fidhelp))
    Rp=[]; Sp=[];
elseif (strcmp(fidtype,'S'))
    Rp=fidhelp; Sp=[];
elseif (strcmp(fidtype,'T'))
    Sp = fidhelp; Rp=[];
else
    Rp=[]; Sp=[];
end

if (nargin<4)
    direc=mean(ecg(fidin(:,1)))-mean(ecg);
else
    if (strcmp(fidtype,'R'))
        direc = 1;
    elseif (strcmp(fidtype,'T') & ~isempty(fidhelp))
        direc = median(ecg(fidin))-median(ecg(Sp));
    else
        direc=-1;
    end
end

if (nargin<6)
    N2 = round(0.04*fs);
end

for i=1:size(fidout,2)
    if (~isempty(Rp) & fidin(i)>Rp(i))%S
        upper = min(fidin(i)+N2,Rp(i)+2*N2+10);
        v=ecg(Rp(i):upper)';% if there's a R right before use it
    else
        v=ecg(fidin(i)-N2:fidin(i)+N2)';
    end
    vi=v;
    if direc>0
        [val,ind]=max(vi);
        fidout(i)=(fidin(i)-N2-1)+ind;
    else
        [val,ind]=min(vi);
        if (~isempty(Rp) & fidin(i)>Rp(i))%S
            fidout(i)=Rp(i)-1+ind;
        else
            fidout(i)=(fidin(i)-N2-1)+ind;
        end
    end
end

