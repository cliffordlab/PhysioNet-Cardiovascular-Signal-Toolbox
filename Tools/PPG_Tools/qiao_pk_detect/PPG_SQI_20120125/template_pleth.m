% [template t2 valid] = template_pleth(wave,anntime,temp_ahead,samp_freq)
% 
% template_pleth.m - PPG waveform template creation.
% by Qiao Li 21 Feb 2011
% 
% input: 
%     wave:       PPG data; 
%     anntime:    PPG annotation time (sample), read from ple annot file,
%                 the time is a offset based on wave(1)
%     temp_ahead: N samples before the beginning of PPG waveform mark,
%                 default is 0
%     samp_freq:  sampling frequency, default is 125Hz
%     
% output:
%     template:   PPG waveform template based on normal-length beats
%     t2:         PPG waveform template2 based on the beats which have 
%                 good correlation with the template
%     valid:      1 for valid template, 2 for valid t2, 3 for both valid,
%                 0 for invalid template and t2
%	
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

function [template t2 valid] = template_pleth(wave,anntime,temp_ahead,samp_freq)

if nargin < 4
    samp_freq = 125;
end

if nargin < 3
    temp_ahead = 0;
end

% according to heart rate max(300bpm) and min(20bpm) to get max and min
% beat-by-beat interval
hr_max=300;
bb_interval_min=samp_freq*60/hr_max;
hr_min=20;
bb_interval_max=samp_freq*60/hr_min;

normal_beat_length_min=0.7;
normal_beat_lentth_max=1.5;
normal_beat_percent_threshold=0.5;

% using xcorr to get the basic period of the PPG as the length of template
y=xcorr(detrend(wave));

len=length(wave);
lena=length(anntime);
i=len+1;
n=1;
peaki=[];
while i<len*2-1
    if y(i)>y(i-1) && y(i)>y(i+1)
        peaki(n)=y(i);
        peakp(n)=i;
        n=n+1;
    end
    i=i+1;
end
if (length(peaki)<1) || (lena<1)
    template=[];
    t2=[];
    valid=0;
    return;
end
[c i]=max(peaki);
peakp=peakp-len;
i=peakp(i);

cycle=samp_freq;
if i<len-1
    cycle=i;
end


% cumulate the beats with reasonable length to get template
p0=1;
i=anntime(p0);
while i-temp_ahead <1 
    p0=p0+1;
    if (p0>lena)
        template=wave;
        valid=0;
        return;
    end
    i=anntime(p0);
end

if p0+1>=lena
    template=[];
    t2=[];
    valid=0;
    return;
end

beat_interval=diff(anntime(p0:length(anntime)));
median_bi=median(beat_interval);
if ~isnan(median_bi)
    temp_peak=abs(peakp-median_bi);
    [m i]=min(temp_peak);
    cycle=peakp(i);
else
    template=[];
    t2=[];
    valid=0;
    return;
end
   
% the length of template valid detection
valid=1;
if cycle > bb_interval_max || cycle < bb_interval_min
    valid=0;
    template=zeros(1,cycle);
    t2=zeros(1,cycle);
    return;
end
    
n=0;
d1=0;
invalidn=0;
currentbeatlength=anntime(p0+1)-anntime(p0);
if currentbeatlength>0 %cycle*normal_beat_length_min %% && currentbeatlength < cycle*normal_beat_lentth_max
    d1=wave(i-temp_ahead:i+cycle-1);%-temp_ahead);
    n=1;
else
    invalidn=invalidn+1;
    d1=zeros(cycle+temp_ahead,1);
end
    
p0=p0+1;
if p0<lena-1
    i=anntime(p0);
    n=1;
    invalidn=0;
    while i<len-cycle && p0<lena-1
        currentbeatlength=anntime(p0+1)-anntime(p0);
        if currentbeatlength>0 %cycle*normal_beat_length_min %% && currentbeatlength < cycle*normal_beat_lentth_max
            d1=d1+wave(i-temp_ahead:i+cycle-1);%-temp_ahead);
            n=n+1;
        else
            invalidn=invalidn+1;
        end
        p0=p0+1;
        i=anntime(p0);
    end
    d1=d1./n;
    % normal beat is less than the reasonable percentage of all beats
    if (n/(n+invalidn))<normal_beat_percent_threshold
        valid=0;
    end
else
    valid=0;
end

% Compare each beat to the template to get template 2
d2=0;
if (valid)
    p0=2;
    i=anntime(p0);
    n=0;
    while i<len-cycle && p0<lena-1
        cc=corrcoef(d1,wave(i-temp_ahead:i+cycle-1));
        if cc(1,2)>0.8
            d2=d2+wave(i-temp_ahead:i+cycle-1);
            n=n+1;
        end
        p0=p0+1;
        i=anntime(p0);
    end
    d2=d2./n;
    valid = uint8(valid);
    if n>length(anntime)*normal_beat_percent_threshold
        valid = bitor(valid, uint8(2));
    end
end
template=d1;
t2=d2;
