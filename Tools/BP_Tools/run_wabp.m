function r = run_wabp(abp)
% WABP  ABP waveform onset detector.
%   r = run_wabp(abp) obtains the onset time (in samples) 
%       of each beat in the ABP waveform.
%
%   In:   ABP (125Hz sampled)
%   Out:  Onset sample time
% 
%   Usage:
%   - ABP waveform must have units of mmHg
%
%   Written by James Sun (xinsun@mit.edu) on Nov 19, 2005.  This ABP onset
%   detector is adapted from Dr. Wei Zong's wabp.c.
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

%% Input checks
if nargin ~=1
    error('exactly 1 argment needed');
end

if size(abp,2)~=1
    error('Input must be a <nx1> vector');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scale physiologic ABP
offset   = 1600;
scale    = 20;
Araw     = abp*scale-offset;

% LPF
A = filter([1 0 0 0 0 -2 0 0 0 0 1],[1 -2 1],Araw)/24+30;
A = (A(4:end)+offset)/scale; % Takes care of 4 sample group delay

% Slope-sum function
dypos          = diff(A);
dypos(dypos<0) = 0;
ssf            = [0; 0; conv(ones(16,1),dypos)];

% Decision rule
avg0       = sum(ssf(1:1000))/1000;   % average of 1st 8 seconds (1000 samples) of SSF
Threshold0 = 3*avg0;                  % initial decision threshold

% ignoring "learning period" for now
lockout    = 0;    % lockout >0 means we are in refractory
timer      = 0;
z          = zeros(100000,1);
counter    = 0;

for t = 50:length(ssf)-17
    lockout = lockout - 1;
    timer   = timer   + 1;      % Timer used for counting time after previous ABP pulse

    if (lockout<1) & (ssf(t)>avg0+5)  % Not in refractory and SSF has exceeded threshold here
        timer = 0;
        maxSSF = max(ssf(t:t+16));  % Find local max of SSF
        minSSF = min(ssf(t-16:t));  % Find local min of SSF
        if maxSSF > (minSSF + 10)
            onset = 0.01*maxSSF ;  % Onset is at the time in which local SSF just exceeds 0.01*maxSSF

            tt       = t-16:t;
            dssf     = ssf(tt) - ssf(tt-1);
            BeatTime = find(dssf<onset,1,'last')+t-17;
            counter  = counter+1;

            if isempty(BeatTime)
                counter = counter-1;
            else
                z(counter) = BeatTime;
            end
            Threshold0 = Threshold0 + 0.1*(maxSSF - Threshold0);  % adjust threshold
            avg0 = Threshold0 / 3;        % adjust avg

            lockout = 32;   % lock so prevent sensing right after detection (refractory period)
        end
    end

    if timer > 312  % Lower threshold if no pulse detection for a while
        Threshold0 = Threshold0 - 1;
        avg0       = Threshold0/3;
    end
end
r = z(find(z))-2;