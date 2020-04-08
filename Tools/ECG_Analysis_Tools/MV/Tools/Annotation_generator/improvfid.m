function [fidout] = improvfid(ecg,fidin,fs,fidtype,fidhelp,N2)
%
% OVERVIEW, This function improves the fiducial point detection for the
% Q,R,S,Toff points in the beat variable. The improved fiducial points are
% stored in the beato variable.
%
% INPUTS    MANDATORY       DESCRIPTION
%
%           ecg             input ecg data in num of samples by 1
%                           dimensions
%           fidin           fiducial point input in samples number.
%           fs              sampling frequency of the input ecg data.
%           fidtype         type of fiducial point input. i.e. Q, R, S or T.              
%           fidhelp         R peak annotations used for improving S and Q point detections.
%           N2              window size of 2*N2 is used to find improved fiducial point detection for fiducial point specified as fidtype.
%
% OUTPUTS
%           fidout          Adjusted fiducial point location after
%                           improving fiducial point detection.
%
% 
% REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Shamim Nemati
%       editted by Ismail Sadiq on 10/26/2019.
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
        %upper=fidin(i)+N2;
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

