function [qrs,jpoints] = wqrsm_fast(data,fs,PWfreq,TmDEF,jflag)
% The swqrsm detector rewrited from swqrsm.c
% rewrited by Giulia Da Poian, August 7, 2018
% use struct swqrsm instead of global varaibles
% giulia.da.poian@emory.com
%
% input:
%   data:   ECG data, the swqrsm.c used the data in raw adus, if the input 
%           data is in physical units, when the function found the peak-peak 
%           value < 10, it will multiply the value by WFDB_DEFGAIN (200)
%   fs:     sampling frequency (default: 125)
%   PWfreq: power line (mains) frequency, in Hz (default: 60)
%   TmDEF:  minimum threshold value (default: 100)
%   jflag:  annotate J-points (ends of QRS complexes) (default: 0)
% output:
%   qrs:    QRS fiducial mark in samples
%   jpoints:J-points annotation, if jflag==1
%
% Matlab code based on:
%  file: swqrsm.c		Wei Zong      23 October 1998
% 			Last revised:  9 April 2010 (by G. Moody)
% -----------------------------------------------------------------------------
% swqrsm: Single-lead QRS detector based on length transform
% Copyright (C) 1998-2010 Wei Zong
% 
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
% PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 59 Temple
% Place - Suite 330, Boston, MA 02111-1307, USA.
% 
% You may contact the author by e-mail (wzong@mit.edu) or postal mail
% (MIT Room E25-505, Cambridge, MA 02139, USA).  For updates to this software,
% please visit PhysioNet (http://www.physionet.org/).
% ------------------------------------------------------------------------------
% 
% This program analyzes an ECG signal, detecting QRS onsets and J-points, using
% a nonlinearly-scaled ECG curve length feature.  This version has been optimized
% for ECGs sampled at 125 Hz, but it can analyze ECGs sampled at any frequency
% using on-the-fly resampling provided by the WFDB library.
% 
% `swqrsm' can process records containing any number of signals, but it uses only
% one signal for QRS detection (signal 0 by default;  this can be changed using
% the `-s' option, see below).  'swqrsm' has been optimized for adult human ECGs.
% For other ECGs, it may be necessary to experiment with the input sampling
% frequency and the time constants.
% ----------------------------------------------------------------------
% This matlab version has updates to deal with a variable sampling 
% frequency but has not been exhanstively tested. 
% ----- Add notes here:
%
%

if nargin<5
jflag = 0;
end
if nargin<4
TmDEF = 100;
end
if nargin<3
PWfreq = 60;
end
if nargin<2
fs = 125;
end

swqrsm.BUFLN = 16384; %	/* must be a power of 2, see ltsamp() */
EYE_CLS = 0.25; %   /* eye-closing period is set to 0.25 sec (250 ms) */ 
MaxQRSw = 0.13; %    /* maximum QRS width (130ms) */                        
NDP	= 2.5; %    /* adjust threshold if no QRS found in NDP seconds */
PWFreqDEF = PWfreq; %    /* power line (mains) frequency, in Hz (default) */
% TmDEF = 100; %	/* minimum threshold value (default) */
WFDB_DEFGAIN = 200.0; %  /* default value for gain (adu/physical unit) */

timer_d=0;
gain=WFDB_DEFGAIN;
swqrsm.lfsc = fix(1.25*gain*gain/fs);	% /* length function scale constant */
PWFreq = PWFreqDEF;	% /* power line (mains) frequency, in Hz */

% test data is physical units (mV) or raw adus units
datatest=data(1:fix(length(data)/fs)*fs);
if length(datatest)>fs
datatest=reshape(datatest,fs,[]);
end
test_ap=median(max(datatest)-min(datatest));
if test_ap<10 % peak-peak < 10 mV, may physical units
data=data*gain;
end

swqrsm.data = data;
swqrsm.lbuf = zeros(swqrsm.BUFLN,1);
swqrsm.ebuf = zeros(swqrsm.BUFLN,1);
swqrsm.ebuf(1:end)=fix(sqrt(swqrsm.lfsc));
swqrsm.lt_tt = 0;
swqrsm.aet = 0;
swqrsm.Yn = 0;
swqrsm.Yn1 = 0;
swqrsm.Yn2 = 0;

qrs = [];
jpoints = [];

Tm = fix(TmDEF / 5.0);
%     spm = 60 * fs;
%     next_minute = from + spm;
swqrsm.LPn = fix(fs/PWFreq); %		/* The LP filter will have a notch at the power line (mains) frequency */
if (swqrsm.LPn > 8)  
swqrsm.LPn = 8;	% /* avoid filtering too agressively */
end
swqrsm.LP2n = 2 * swqrsm.LPn;
EyeClosing = fix(fs * EYE_CLS); % /* set eye-closing period */
ExpectPeriod = fix(fs * NDP);   % /* maximum expected RR interval */
swqrsm.LTwindow = fix(fs * MaxQRSw);   % /* length transform window size */

%     for i=1:2000
%         ltdata(i)=ltsamp(i);
%     end
%     /* Average the first 8 seconds of the length-transformed samples
%        to determine the initial thresholds Ta and T0. The number of samples
%        in the average is limited to half of the ltsamp buffer if the sampling
%        frequency exceeds about 2 KHz. */

t1 = fs*8;
if t1> fix(swqrsm.BUFLN*0.9)
    t1=swqrsm.BUFLN/2;
end

T0=0;
for t=1:t1
    swqrsm = ltsamp(t,swqrsm);
    T0 = T0 + swqrsm.lt_data;
end

T0 = T0 / t1;
Ta = 3 * T0;

% /* Main loop */
t=1;
learning=1;
while t<length(data)
    if (learning==1)
        if (t > t1) 
            learning = 0;
            T1 = T0;
            t = 1;	% /* start over */
        else
            T1 = 2*T0;
        end
    end

    % /* Compare a length-transformed sample against T1. */
    swqrsm = ltsamp(t,swqrsm);
    if ( swqrsm.lt_data > T1) %	/* found a possible QRS near t */
        timer_d = 0; % /* used for counting the time after previous QRS */
        maxd = swqrsm.lt_data;
        mind = maxd;
        for tt = t+1:t + fix(EyeClosing/2)
            swqrsm = ltsamp(tt,swqrsm);
            if (swqrsm.lt_data > maxd) 
                maxd = swqrsm.lt_data;
            end
        end
        for tt = t-1:-1:t - fix(EyeClosing/2)
            swqrsm = ltsamp(tt,swqrsm);
            if (swqrsm.lt_data < mind) 
                mind = swqrsm.lt_data;
            end
        end
        if (maxd > mind+10) % /* There is a QRS near tt */
        % /* Find the QRS onset (PQ junction) */
            onset = fix(maxd/100) + 2;
            tpq = t - 5;
            for tt = t:-1:t - fix(EyeClosing/2)
                swqrsm = ltsamp(tt,swqrsm);
                swqrsm_1 = ltsamp(tt-1,swqrsm);
                swqrsm_2 = ltsamp(tt-2,swqrsm_1);
                swqrsm_3 = ltsamp(tt-3,swqrsm_2);
                swqrsm_4 = ltsamp(tt-4,swqrsm_3);
                if (swqrsm.lt_data - swqrsm_1.lt_data < onset && swqrsm_1.lt_data - swqrsm_2.lt_data < onset && swqrsm_2.lt_data - swqrsm_3.lt_data < onset && swqrsm_3.lt_data - swqrsm_4.lt_data < onset) 
                    swqrsm = swqrsm_4;
                    tpq = tt - swqrsm.LP2n; %	/* account for phase shift */
                    break;
                end
            end

            if (learning~=1) 
                % /* Check that we haven't reached the end of the record. */
                if tpq>length(data)
                    break;
                end
                % /* Record an annotation at the QRS onset */
                qrs = [qrs tpq];

                % J-points processing
                if (jflag)
                    tj = t+5;
                    for tt=t:t + fix(EyeClosing/2)
                        swqrsm = ltsamp(tt, swqrsm);
                        if swqrsm.lt_data > maxd - fix(maxd/10)
                            tj = tt;
                            break;
                        end
                    end
                    if tj>length(data)
                        break;
                    end
                    % Record an annotation at the J-point
                    jpoints = [jpoints tj];
                end
            end
            % /* Adjust thresholds */
            Ta = Ta + (maxd - Ta)/10;
            T1 = Ta / 3;

            % /* Lock out further detections during the eye-closing period */
            t = t + EyeClosing;
        end

    elseif (learning~=1) 
    % 	    /* Once past the learning period, decrease threshold if no QRS
    % 	       was detected recently. */
        timer_d = timer_d + 1;
        if timer_d > ExpectPeriod && Ta > Tm
            Ta = Ta -1;
            T1 = Ta / 3;
        end
    end
    t = t+1;
end

end

function swqrsm = ltsamp(t,swqrsm)
    % /* ltsamp() returns a sample of the length transform of the input at time t.
    %    Since this program analyzes only one signal, ltsamp() does not have an
    %    input argument for specifying a signal number; rather, it always filters
    %    and returns samples from the signal designated by the global variable
    %    'sig'.  The caller must never "rewind" by more than swqrsm.BUFLN samples (the
    %    length of ltsamp()'s buffers). */
    while (t > swqrsm.lt_tt) 
        swqrsm.Yn2 = swqrsm.Yn1;
        swqrsm.Yn1 = swqrsm.Yn;
        v0=swqrsm.data(1);
        v1=swqrsm.data(1);
        v2=swqrsm.data(1);
        if swqrsm.lt_tt>0 && swqrsm.lt_tt<=length(swqrsm.data)
            v0 = swqrsm.data(swqrsm.lt_tt);
        end
        if swqrsm.lt_tt-swqrsm.LPn>0 && (swqrsm.lt_tt-swqrsm.LPn)<=length(swqrsm.data)
            v1 = swqrsm.data(swqrsm.lt_tt-swqrsm.LPn);
        end        
        if swqrsm.lt_tt-swqrsm.LP2n>0 && (swqrsm.lt_tt-swqrsm.LP2n)<=length(swqrsm.data)
            v2 = swqrsm.data(swqrsm.lt_tt-swqrsm.LP2n);
        end
        if v0~=-32768 && v1~=-32768 && v2~=-32768
            swqrsm.Yn = 2*swqrsm.Yn1 - swqrsm.Yn2 + v0 - 2*v1 + v2;
        end
        dy = fix((swqrsm.Yn - swqrsm.Yn1) / swqrsm.LP2n);	%	/* lowpass derivative of input */
        swqrsm.lt_tt=swqrsm.lt_tt+1;
        swqrsm.et = fix(sqrt(swqrsm.lfsc +dy*dy)); % /* length transform */
        id = mod(swqrsm.lt_tt,swqrsm.BUFLN);
        if id == 0
            id = swqrsm.BUFLN;
        end
        swqrsm.ebuf(id) = swqrsm.et;
        id2 = mod(swqrsm.lt_tt-swqrsm.LTwindow,swqrsm.BUFLN);
        if id2 == 0
            id2 = swqrsm.BUFLN;
        end
        swqrsm.aet = swqrsm.aet + (swqrsm.et - swqrsm.ebuf(id2));
        swqrsm.lbuf(id) = swqrsm.aet;
    % 	/* swqrsm.lbuf contains the average of the length-transformed samples over
    % 	   the interval from tt-swqrsm.LTwindow+1 to tt */
    end
    
    id3 = mod(t,swqrsm.BUFLN);
    if id3 == 0
        id3 = swqrsm.BUFLN;
    end
    swqrsm.lt_data = swqrsm.lbuf(id3);
end