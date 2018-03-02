function r = abpfeature(abp,OnsetTimes, Fs)
% ABPFEATURE  ABP waveform feature extractor.
%   r = ABPFEATURE(ABP,ONSETTIMES) extracts features from ABP waveform such
%   as systolic pressure, mean pressure, etc.
%
%   In:     ABP = signal (default 125Hz sampled)
%           ONSETTIMES = times of onset (in samples)
%   Out:    Beat-to-beat ABP features
%           Col 1:  Time of systole   [samples]
%               2:  Systolic BP       [mmHg]
%               3:  Time of diastole  [samples]
%               4:  Diastolic BP      [mmHg]
%               5:  Pulse pressure    [mmHg]
%               6:  Mean pressure     [mmHg]
%               7:  Beat Period       [samples]
%               8:  mean_dyneg
%               9:  End of systole time  0.3*sqrt(RR)  method
%              10:  Area under systole   0.3*sqrt(RR)  method
%              11:  End of systole time  1st min-slope method
%              12:  Area under systole   1st min-slope method
%              13:  Pulse             [samples]
% 
%   Usage:
%   - OnsetTimes must be obtained using wabp.m
%
%   Written by James Sun (xinsun@mit.edu) on Nov 19, 2005.
%   Updated by Alistair Johnson, 2014.
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information




if length(OnsetTimes)<30 % don't process anything if too little onsets
    r = [];
    return
end

%% P_sys, P_dias
if nargin<3
Fs        = 125; % Sampling Frequency
end
Window    = ceil(0.32*Fs);
OT        = OnsetTimes(1:end-1);
BeatQty   = length(OT);

% [MinDomain,MaxDomain] = init(zeros(BeatQty,Window));

for i=1:Window % Vectorized version
    MinDomain(:,i) = OT-i+1;
	MaxDomain(:,i) = OT+i-1;
end

%MinDomain = bsxfun(@minus,OT,(0:Window-1)); 
%MaxDomain = bsxfun(@plus,OT,(0:Window-1)); 
MinDomain(MinDomain<1) = 1;  % Error protection
MaxDomain(MaxDomain<1) = 1;
MinDomain(MinDomain>length(abp)) = length(abp);
MaxDomain(MaxDomain>length(abp)) = length(abp);

[P_dias Dindex]  = min(abp(MinDomain),[],2);% Get lowest value across 'Window' samples before beat onset
[P_sys  Sindex]  = max(abp(MaxDomain),[],2);% Get highest value across 'Window' samples after beat onset

DiasTime         = MinDomain(sub2ind(size(MinDomain),(1:BeatQty)',Dindex)); % Map offset indices Dindex to original indices
SysTime          = MaxDomain(sub2ind(size(MaxDomain),(1:BeatQty)',Sindex)); % Map offset indices Sindex to original indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulse Pressure [mmHg]
PP         = P_sys - P_dias;

% Beat Period [samples]
BeatPeriod = diff(OnsetTimes);

% Mean,StdDev, Median Deriv- (noise detector)
dyneg          = diff(abp);
dyneg(dyneg>0) = 0;

[MAP,stddev,mean_dyneg] = init(zeros(BeatQty,1));
if OnsetTimes(end)==numel(abp)
    OnsetTimes(end) = numel(abp) - 1;
end

for i=1:BeatQty
    interval       = abp(OnsetTimes(i):OnsetTimes(i+1));
    MAP(i)         = mean(interval);
    stddev(i)      = std(interval);

    dyneg_interval = dyneg(OnsetTimes(i):OnsetTimes(i+1));
    dyneg_interval(dyneg_interval==0) = [];
    if min(size(dyneg_interval))==0
        dyneg_interval = 0;
    end
    mean_dyneg(i)  = mean(dyneg_interval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Systolic Area calculation using 0.3*sqrt(RR)
RR                  = BeatPeriod/Fs;  % RR time in seconds
sys_duration        = 0.3*sqrt(RR);
EndOfSys1           = round(OT + sys_duration*Fs);
SysArea1            = localfun_area(abp,OT,EndOfSys1',P_dias);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Systolic Area calculation using 'first minimum slope' method
SlopeWindow         = ceil(0.28*Fs);
ST                  = EndOfSys1; % start 12 samples after P_sys

if ST(end) > (length(abp)-SlopeWindow)   % error protection
    ST(end) = length(abp)-SlopeWindow;
end

SlopeDomain = zeros(BeatQty,SlopeWindow);
for i=1:SlopeWindow
    SlopeDomain(:,i) = ST+i-1;
end
% y[n] = x[n]-x[n-1]
Slope              = diff(abp(SlopeDomain),1,2);
Slope(Slope>0)     = 0; % Set positive slopes to zero

[MinSlope index]   = min(abs(Slope),[],2); % Find first positive slope
EndOfSys2          = SlopeDomain(sub2ind(size(SlopeDomain),(1:BeatQty)',index));
SysArea2           = localfun_area(abp,OT,EndOfSys2,P_dias);
Pulse              = 60./BeatPeriod;


% OUTPUT:
% Ensure that there is no concatenation error by using
% (:) indexing to force column vectors
r = [SysTime(:),P_sys(:),DiasTime(:),P_dias(:),PP(:),...
    MAP(:),BeatPeriod(:),mean_dyneg(:),EndOfSys1(:),SysArea1(:),...
    EndOfSys2(:),SysArea2(:),Pulse(:)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function:

function SysArea = localfun_area(abp,onset,EndSys,P_dias)
%% Input:  abp,
%%         onset, P_dias, end of systole, beat duration in unit of samples, same length
%% Output: systolic area, warner correction factor

BeatQty = length(onset);
SysArea = init(zeros(BeatQty,1));
for i=1:BeatQty
    SysArea(i) = sum(abp(onset(i):EndSys(i))); % faster than trapz below
    %b(i) = trapz(abp(onset(i):EndSys(i)));
end
EndSys = EndSys(:); onset = onset(:); % force col vectors
SysPeriod  = EndSys     - onset;

% Time scale and subtract the diastolic area under each systolic interval
SysArea = (SysArea - P_dias.*SysPeriod)/125; % Area [mmHg*sec]
end

% for initializing a variable
function varargout = init(z)
for i=1:nargout
    varargout(i) = {z};
end
end