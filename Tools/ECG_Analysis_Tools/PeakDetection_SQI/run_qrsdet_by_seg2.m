function QRS = run_qrsdet_by_seg2(ecg,fs,window,thres,ecgType)
% this function is used to run the QRS detector for each window window (non overlapping)
% this is needed because in the case of big artefacts at the begining of
% the record (i.e. where the P&T threshold is computed) it can make the detection fail.
%
% inputs
%   ecg:        ecg signal
%   fs:         sampling frequency
%   window:     size of the window onto which to perform QRS detection (in seconds)
%   thres:      threshold to be used
%   ecgType:    'FECG' or 'MECG'
% 
% output
%   QRS:        QRS location in nb samples (ms)
%
%
% FECG-ESN toolbox, version 1.0
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 28-01-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.


% == managing inputs
if nargin<1; error('run_qrsdet_by_seg: wrong number of input arguments \n'); end;
if nargin<2; fs=1000; end;
if nargin<3; window=15; end;
if nargin<4; thres=0.6; end;
if nargin<5; ecgType='FECG'; end;

% == general
segsizeSamp = window*fs; % convert window into nb of samples
NbSeg = floor(length(ecg)/segsizeSamp); % nb of segments
QRS = [];
start = 1;
stop  = segsizeSamp;
signForce = 0; % if we want to force the sign of the peak we are looking for

% == core function
try
    for ch=1:NbSeg
        % for each segment perform QRS detection
        QRStemp = 0;

        % take +/-1sec around selected subsegment exept for the borders. This
        % is in case there is a QRS in between segments -> allows to locate
        % them well.
        if ch==1
            % first subsegment
            dTplus  = fs;
            dTminus = 0;
        elseif ch==NbSeg
            % last subsegment
            dTplus  = 0;
            dTminus = fs;    
        else
            % any other subsegment
            dTplus  = fs;
            dTminus = fs;  
        end

        % lowering the threshold in case not enough beats are
        % detected. Also changed the refractory period to be different between
        % mother and foetus. sign of peaks is determined by the sign on the
        % first window and then is forced for the following windows.
        if strcmp(ecgType,'FECG')
            thresTrans = thres;
            while length(QRStemp)<20 && thresTrans>0.1
                [QRStemp,signForce] = qrs_detect2(ecg(start-dTminus:stop+dTplus),0.150,thresTrans,fs,[],signForce);
                thresTrans = thresTrans-0.1;
            end
        else
            [QRStemp,signForce] = qrs_detect2(ecg(start-dTminus:stop+dTplus),0.250,thres,fs,[],signForce);
        end

        NewQRS = (start-1)-dTminus+QRStemp;
        NewQRS(NewQRS>stop) = [];
        NewQRS(NewQRS<start) = [];

        if ~isempty(QRS) && ~isempty(NewQRS) && NewQRS(1)-QRS(end)<0.25*fs
            % this is needed to avoid double detection at the transition point
            % between two windows
            NewQRS(1) = [];
        end
        QRS = [QRS NewQRS];

        start = start+segsizeSamp;
        stop = stop+segsizeSamp;
    end

catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    QRS = [1000 2000];
    rethrow(ME);
end

end