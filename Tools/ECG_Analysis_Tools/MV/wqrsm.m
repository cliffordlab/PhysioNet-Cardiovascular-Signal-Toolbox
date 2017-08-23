function [qrs]=wqrsm(data,fs,PWfreq,TmDEF)
% The WQRS detector rewrited from wqrs.c
% rewrite by Qiao Li, March 6, 2017
% qiaolibme@gmail.com
%
% input:
%   data:   ECG data, the wqrs.c used the data in raw adus, if the input 
%           data is in physical units, when the function found the peak-peak 
%           value < 10, it will multiply the value by WFDB_DEFGAIN (200)
%   fs:     sampling frequency (default: 125)
%   PWfreq: power line (mains) frequency, in Hz (default: 60)
%   TmDEF:  minimum threshold value (default: 100)
% output:
%   qrs:    QRS fiducial mark in samples
%
% /* file: wqrs.c		Wei Zong      23 October 1998
% 			Last revised:  9 April 2010 (by G. Moody)
% -----------------------------------------------------------------------------
% wqrs: Single-lead QRS detector based on length transform
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
% `wqrs' can process records containing any number of signals, but it uses only
% one signal for QRS detection (signal 0 by default;  this can be changed using
% the `-s' option, see below).  'wqrs' has been optimized for adult human ECGs.
% For other ECGs, it may be necessary to experiment with the input sampling
% frequency and the time constants indicated below.
% 
% To compile this program under GNU/Linux, MacOS/X, MS-Windows, or Unix, use gcc:
% 	gcc -o wqrs wqrs.c -lwfdb -lm
% You must have installed the WFDB library, available at	
% 	http://www.physionet.org/physiotools/wfdb.shtml
% gcc is standard with GNU/Linux and is available for other platforms from:
% 	http://www.gnu.org/		(sources and Unix binaries)
% 	http://fink.sourceforge.net	(Mac OS/X only)
% 	http://www.cygwin.com/   	(MS-Windows only)
% 	
% For a usage summary, see the help text at the end of this file.  The input
% record may be in any of the formats readable by the WFDB library, and it may
% be anywhere in the WFDB path (in a local directory or on a remote web or ftp
% server).  The output of 'wqrs' is an annotation file named RECORD.wqrs (where
% RECORD is replaced by the name of the input record).  Within the output
% annotation file, the time of each NORMAL annotation marks a QRS onset;  if
% the '-j' option is used, additional JPT annotations mark the J points (the
% ends of the QRS complexes).
% 
% For example, to mark QRS complexes in record 100 beginning 5 minutes from the
% start, ending 10 minutes and 35 seconds from the start, and using signal 1, use
% the command:
%     wqrs -r 100 -f 5:0 -t 10:35 -s 1
% The output may be read using (for example):
%     rdann -a wqrs -r 100
% To evaluate the performance of this program, run it on the entire record, by:
%     wqrs -r 100
% and then compare its output with the reference annotations by:
%     bxb -r 100 -a atr wqrs
% */

if nargin<4
    TmDEF = 100;
end
if nargin<3
    PWfreq = 60;
end
if nargin<2
    fs = 125;
end

clear global Yn Yn1 Yn2 lt_tt lbuf ebuf aet et LPn LP2n lfsc wqrs_data BUFLN LTwindow

global Yn Yn1 Yn2;
global lt_tt lbuf ebuf;
global aet et;
global LPn LP2n lfsc;
global wqrs_data;
global BUFLN;
global LTwindow;

BUFLN = 16384; %	/* must be a power of 2, see ltsamp() */
EYE_CLS = 0.25; %   /* eye-closing period is set to 0.25 sec (250 ms) */ 
MaxQRSw = 0.13; %    /* maximum QRS width (130ms) */                        
NDP	= 2.5; %    /* adjust threshold if no QRS found in NDP seconds */
PWFreqDEF = PWfreq; %    /* power line (mains) frequency, in Hz (default) */
% TmDEF = 100; %	/* minimum threshold value (default) */
WFDB_DEFGAIN = 200.0; %  /* default value for gain (adu/physical unit) */

timer_d=0;
gain=WFDB_DEFGAIN;
lfsc = fix(1.25*gain*gain/fs);	% /* length function scale constant */
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

wqrs_data = data;

lbuf = zeros(BUFLN,1);
ebuf = zeros(BUFLN,1);
ebuf(1:end)=fix(sqrt(lfsc));

lt_tt = 0;
aet = 0;
Yn = 0;
Yn1 = 0;
Yn2 = 0;
    
qrs = [];

    Tm = fix(TmDEF / 5.0);
%     spm = 60 * fs;
%     next_minute = from + spm;
    LPn = fix(fs/PWFreq); %		/* The LP filter will have a notch at the power line (mains) frequency */
    if (LPn > 8)  
        LPn = 8;	% /* avoid filtering too agressively */
    end
    LP2n = 2 * LPn;
    EyeClosing = fix(fs * EYE_CLS); % /* set eye-closing period */
    ExpectPeriod = fix(fs * NDP);   % /* maximum expected RR interval */
    LTwindow = fix(fs * MaxQRSw);   % /* length transform window size */

%     for i=1:2000
%         ltdata(i)=ltsamp(i);
%     end
    
%     /* Average the first 8 seconds of the length-transformed samples
%        to determine the initial thresholds Ta and T0. The number of samples
%        in the average is limited to half of the ltsamp buffer if the sampling
%        frequency exceeds about 2 KHz. */
    t1 = fs*8;
    if t1> fix(BUFLN*0.9)
        t1=BUFLN/2;
    end
    
    T0=0;
    for t=1:t1
        T0 = T0 + ltsamp(t);
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
	if (ltsamp(t) > T1) %	/* found a possible QRS near t */
	    timer_d = 0; % /* used for counting the time after previous QRS */
	    maxd = ltsamp(t);
        mind = maxd;
	    for tt = t+1:t + fix(EyeClosing/2)
            if (ltsamp(tt) > maxd) 
                maxd = ltsamp(tt);
            end
        end
	    for tt = t-1:-1:t - fix(EyeClosing/2)
    		if (ltsamp(tt) < mind) 
                mind = ltsamp(tt);
            end
        end
	    if (maxd > mind+10) % /* There is a QRS near tt */
		% /* Find the QRS onset (PQ junction) */
            onset = fix(maxd/100) + 2;

            tpq = t - 5;
            for tt = t:-1:t - fix(EyeClosing/2)
                if (ltsamp(tt)   - ltsamp(tt-1) < onset && ltsamp(tt-1) - ltsamp(tt-2) < onset && ltsamp(tt-2) - ltsamp(tt-3) < onset && ltsamp(tt-3) - ltsamp(tt-4) < onset) 
                    tpq = tt - LP2n; %	/* account for phase shift */
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

function lt_data = ltsamp(t)
% /* ltsamp() returns a sample of the length transform of the input at time t.
%    Since this program analyzes only one signal, ltsamp() does not have an
%    input argument for specifying a signal number; rather, it always filters
%    and returns samples from the signal designated by the global variable
%    'sig'.  The caller must never "rewind" by more than BUFLN samples (the
%    length of ltsamp()'s buffers). */

global Yn Yn1 Yn2;
global lt_tt lbuf ebuf;
global aet et;
global LPn LP2n lfsc;
global wqrs_data;
global BUFLN;
global LTwindow;

while (t > lt_tt) 
	Yn2 = Yn1;
	Yn1 = Yn;
    v0=wqrs_data(1);
    v1=wqrs_data(1);
    v2=wqrs_data(1);
    if lt_tt>0 && lt_tt<=length(wqrs_data)
        v0 = wqrs_data(lt_tt);
    end
    if lt_tt-LPn>0 && (lt_tt-LPn)<=length(wqrs_data)
        v1 = wqrs_data(lt_tt-LPn);
    end        
    if lt_tt-LP2n>0 && (lt_tt-LP2n)<=length(wqrs_data)
        v2 = wqrs_data(lt_tt-LP2n);
    end
    if v0~=-32768 && v1~=-32768 && v2~=-32768
        Yn = 2*Yn1 - Yn2 + v0 - 2*v1 + v2;
    end
	dy = fix((Yn - Yn1) / LP2n);	%	/* lowpass derivative of input */
    lt_tt=lt_tt+1;
	et = fix(sqrt(lfsc +dy*dy)); % /* length transform */
    id = mod(lt_tt,BUFLN);
    if id == 0
        id = BUFLN;
    end
    ebuf(id) = et;
    id2 = mod(lt_tt-LTwindow,BUFLN);
    if id2 == 0
        id2 = BUFLN;
    end
    aet = aet + (et - ebuf(id2));
    lbuf(id) = aet;
% 	/* lbuf contains the average of the length-transformed samples over
% 	   the interval from tt-LTwindow+1 to tt */
end
    id3 = mod(t,BUFLN);
    if id3 == 0
        id3 = BUFLN;
    end
    lt_data = lbuf(id3);
end
