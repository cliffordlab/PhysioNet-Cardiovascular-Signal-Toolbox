function ecgout = Filter_TWA(ecgin,mains,fc)
%OVERVIEW, This function filters the ecg for powerline interference at a
%frequency of mains Hz if needed. It also removes baseline wander and
% low pass filters the ecg with cut off frequency fc.
%
% INPUTS        MANDATORY
%
%               ecgin (uV)      N by M array of ecg data. Contains M
%                               channels of ecg with N datapoints.
%
%               mains (Hz)      Frequency for mains electricity.
%
%               fc (Hz)         cut off frequency for low pass filter.
%
% OUTPUTS       ecgout (uV)     filtered ecg
%
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Ismail Sadiq
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  


out_len = ceil(size(ecgin,1));
ecgout = zeros(out_len,size(ecgin,2));
% to remove maines noise
if (fc > mains-10)
    [bm,am] = butter(5,[(mains-0.5)/(Fs_in/2) (mains+0.5)/(Fs_in/2)],'stop');
end
% low-pass filter (removes anything above fc)
if (fc < 1000/2)
    [bl,al] = butter(max([min([ceil(1000/100),11]),5]),40/(1000/2));
end

% apply filters
for i = 1:size(ecgin, 2)
    tmp = detrend(ecgin(:,i));
    tmp = medianfilter(tmp,1000); % remove baseline
        if (fc < 1000/2)
            tmp = filtfilt(bl,al,tmp); % lowpass noise
        end
    if (fc > mains-10)
        tmp = filtfilt(bm,am,tmp); % mains noise
    end
    ecgout(:,i)= tmp;
end

end

% if (0)
%     t = linspace(0,length(ecgin(:,1))/Fs_in,length(ecgin(:,1)));
%     tt = linspace(0,length(ecgin(:,1))/Fs_in,length(ecgin(:,1))/Fs_in*Fs_out);
%     plot(t,ecgin(:,1))
%     hold on
%     plot(tt,ecgout(:,1),'r')
% end