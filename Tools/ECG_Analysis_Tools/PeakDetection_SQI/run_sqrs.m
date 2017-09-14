function out = run_sqrs(ecg,HRVparams,rs)
%out = run_sqrs(ecg,s,rs)
%	OVERVIEW:
%       QRS onset detector
%   INPUT:
%       ecg = single channel of ECG in digital values
%       s   = settings struct
%       rs  = resampling option; 1 = resample original signal, 0 = don't
%   OUTPUT:
%       out = sample points of onset of each QRS complex
%   DEPENDENCIES & LIBRARIES:
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Converted from sqrs.c at:
%       http://www.physionet.org/physiotools/wfdb/app/sqrs.c
%       to  sqrs.m (Matlab)     by J. Perry,    July        2006
%           sqrs.c	(C)         by	G.B. Moody  27 October  1990
%
%	Last revised by John:  25 February 2006
%	Last revised by Gari:  07 February 2007
%   
%   03-06-2017
%   Edited by Adriana Vest
%   Now requires settings struct HRVparams, which defines debug mode and 
%   sampling frequency, and resampling option rs.
%   LICENSE AND COPYRIGHT:    
%       See Below
% 
% if debug>0, then plot final data (default; debug=0)
% 
% This algorithm is sensitive to high frequency noise, so you 
% might want to use an FIR band-pass filter first.
% 
% Avaialble under the GPL (see below)
%-------------------------------------------------------------------------------
% sqrs: Single-channel QRS detector
% Copyright (C) 1990-2006 George B. Moody
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
% 02111-1307, USA.
%
% You may contact the author by e-mail (george@mit.edu) or postal
% mail (MIT Room E25-505A, Cambridge, MA 02139 USA).  For updates to
% this software, please visit PhysioNet (http://www.physionet.org/).
% _________________________________________________________________
%
% The detector algorithm is based on example 10 in the WFDB
% Programmer's Guide, which in turn is based on a Pascal program
% written by W.A.H. Engelse and C. Zeelenberg, "A single scan
% algorithm for QRS-detection and feature extraction", Computers in
% Cardiology 6:37-42 (1979).  `sqrs' does not include the feature
% extraction capability of the Pascal program.  The output of `sqrs'
% is an annotation file (with annotator name `qrs') in which all
% detected beats are labelled normal;  the annotation file may also
% contain `artifact' annotations at locations which `sqrs' believes
% are noise-corrupted.
%
% 'sqrs' has been optimized for adult human ECGs. For other ECGs,
% it may be necessary to experiment with the input sampling
% frequency and the time constants indicated below. 
%
% This program is provided as an example only, and is not intended
% for any clinical application.  At the time the algorithm was
% originally published, its performance was typical of
% state-of-the-art QRS detectors.  Recent designs, particularly
% those that can analyze two or more input signals, may exhibit
% significantly better performance.
%%

if nargin<2
    error('Incorrect number of input arguments: must provide ECG signal and HRV paramters struct')
end
if nargin < 3
    rs = 1;
end

%debug = s.debug;
debug = 0;
fs = HRVparams.Fs;

if ~rs
    %% 1. Use Fs of Input
    time = 0;
    now = 10;

    %freq = 256;
    freq = fs;
    ms160 = ceil(0.16*freq);
    ms200 = ceil(0.2*freq);
    s2 = ceil(2*freq);
    scmin = 500;
    % number of ADC units corresponding to scmin microvolts
    % scmin = muvadu(signal, scmin); 
    scmax = 10 * scmin;
    slopecrit = 10 * scmin;
    maxslope = 0;
    nslope = 0;
    out = [];

    while (now < length(ecg)) %   && (to == 0)
        filter = [1 4 6 4 1 -1 -4 -6 -4 -1] * ecg((now-9):now);
        if (mod(time, s2) == 0)
	    % Adjust slope 
            if (nslope == 0)
	      slopecrit = max(slopecrit - slopecrit/16, scmin);
            elseif (nslope >= 5)
	      slopecrit = min(slopecrit + slopecrit/16, scmax);
            end
        end
        if (nslope == 0  && abs(filter) > slopecrit)
            nslope = nslope + 1; 
	    maxtime = ms160;
	    if (filter > 0) 
	      sign = 1;
	    else
	      sign = -1;
	    end
            qtime = time;
        end
        if (nslope ~= 0)
            if (filter * sign < -slopecrit)
                sign = -sign;
                nslope = nslope + 1;
                if (nslope > 4) 
                  maxtime = ms200;
                else
                  maxtime = ms160;
                end
            elseif (filter * sign > slopecrit &&  abs(filter) > maxslope)
                maxslope = abs(filter);
            end
            
            if (maxtime < 0)
                if (2 <= nslope && nslope <= 4)
                    slopecrit = slopecrit + ((maxslope/4) - slopecrit)/8;
                    if (slopecrit < scmin)
                        slopecrit = scmin;
                    elseif (slopecrit > scmax) 
                        slopecrit = scmax;
                    end
                    out = [out; now - (time - qtime) - 4];
                    %annot.anntyp = NORMAL; 
                    time = 0;
                elseif (nslope >= 5)
                    out = [out; now - (time - qtime) - 4];
                    %annot.anntyp = ARFCT; 
                end
                nslope = 0;
            end
            maxtime = maxtime - 1;
        end
	time = time + 1;
	now = now + 1;
    end

    out=out-1; % adjust for 1 sample offset problem.

    if debug > 0
       plot(ecg,'b');
       hold on;
       plot(out,ecg(out),'m*'); 
    end
elseif rs
    %% 2. Resample at 256 Hz
    % These time constants may need adjustment for pediatric or
    % small mammal ECGs.
    ecg_rs = resample(ecg,256,fs);
    
    time = 0;
    now = 10;

    freq = 256;
    ms160 = ceil(0.16*freq);
    ms200 = ceil(0.2*freq);
    s2 = ceil(2*freq);
    %scmin = 500;
    % ANV changed this to work with real units
    scmin = 500;
    % number of ADC units corresponding to scmin microvolts
    % scmin = muvadu(ecg_rs, scmin); 
    scmax = 10 * scmin;
    slopecrit = 10 * scmin;
    maxslope = 0;
    nslope = 0;
    out = [];

    while (now < length(ecg_rs)) %   && (to == 0)
        filter = [1 4 6 4 1 -1 -4 -6 -4 -1] * ecg_rs((now-9):now);
        if (mod(time, s2) == 0)
            % Adjust slope 
            if (nslope == 0)
                slopecrit = max(slopecrit - slopecrit/16, scmin);
            elseif (nslope >= 5)
                slopecrit = min(slopecrit + slopecrit/16, scmax);
            end
        end
        if (nslope == 0  && abs(filter) > slopecrit)
            nslope = nslope + 1; 
            maxtime = ms160;
            if (filter > 0) 
                sign = 1;
            else
                sign = -1;
            end
            qtime = time;
        end
        if (nslope ~= 0)
            if (filter * sign < -slopecrit)
                sign = -sign;
                nslope = nslope + 1;
                if (nslope > 4) 
                    maxtime = ms200;
                else
                  maxtime = ms160;
                end
            elseif (filter * sign > slopecrit &&  abs(filter) > maxslope)
                maxslope = abs(filter);
            end
            if (maxtime < 0)
                if (2 <= nslope && nslope <= 4)
                    slopecrit = slopecrit + ((maxslope/4) - slopecrit)/8;
                    if (slopecrit < scmin)
                        slopecrit = scmin;
                    elseif (slopecrit > scmax) 
                        slopecrit = scmax;
                    end
                    out = [out; now - (time - qtime) - 4];
                    %annot.anntyp = NORMAL; 
                    time = 0;
                elseif (nslope >= 5)
                    out = [out; now - (time - qtime) - 4];
                    %annot.anntyp = ARFCT; 
                end
                nslope = 0;
            end
            maxtime = maxtime - 1;
        end
	time = time + 1;
	now = now + 1;
    end

    out=out-1; % adjust for 1 sample offset problem.

    if debug > 0
       plot(ecg_rs,'b');
       hold on;
       plot(out,ecg_rs(out),'m*'); 
    end

end