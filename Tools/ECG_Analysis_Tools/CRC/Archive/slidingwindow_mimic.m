% This script will compute the CRC Index for the mimic data provided you
% already have the RR Interval (saved in variable "rr"), R ampltude (save
% in variable "ramp") and RR interval time stamp (saved in variable "t").
% These variables were obtained using your matlab version of Pan and
% Tompkins algo.


close all

% initialize
fs = 125;
win = 180;
slide = 120;

MAX_C = [];
MAX_F = [];
T = [];

N = floor((t(end) - win)/slide) + 1; % Compute the number of windows

    % compute CRC Index
    for k = 1:N
        st = (k-1) * slide;
        T = [T; st];    % Add to time vector
        fn = st + win;
        index_rr = find(t> st & t< fn);
        rr_seg = rr(index_rr); t_seg = t(index_rr);
        ramp_seg = ramp(index_rr);
        [t_seg,rr_seg, ramp_seg] = ectopic_rejection(t_seg, rr_seg, ramp_seg, .15);
        if length(rr_seg) > 80
            [t_interp,rr_interp] =  interp_RR(t_seg, rr_seg, 4, 3);
            [t_interp,ramp_interp] =  interp_RR(t_seg, ramp_seg, 4, 3);
            
            [c,f] = crc_index(rr_interp,ramp_interp, 2^9, 4, hanning(400), 200);
            MAX_C = [MAX_C; c];
            MAX_F = [MAX_F; f];
        else 
            MAX_C = [MAX_C; 0];
            MAX_F = [MAX_F; 0];
        end
        
    end
    
%plot

figure
subplot(211)
plot(T,MAX_C)
axis([0 T(end) 0 1.3])
xlabel('Time (s)')
ylabel('Cardiorespiratory Coupling Index')

subplot(212)
plot(T,MAX_F)
axis([0 T(end) 0 .5])
ylabel('frequency at maximum CRC between 1.5 and 5 Hz')
xlabel('Time (s)')