function [acm dcm AC DC] = prsa(rr,rri,fs,win_length,thresh_per,testmode) 
% [acm dcm AC DC] = prsa(rr,rri,fs,win_length,thresh_per,testmode); 
%
% calculates acceleration and deceleration matrix  
% rr - rr intervals, rri - index of rr intervals, fs - sampling freq, 
% if fs=[] then rri are assumed to be in seconds
% win_length-length of window taken around anchors
% Thresh_per - percent that beat can differ from previous beat - default =5
% Testmode = 1 creates a figure 
% empty arrays are permissable except for first arguemnt:
% ... e.g.:
% [acm dcm AC DC] = prsa(rr,[], [],[],[],1);
%
% License - GNU GPL
% Author: L. Campana
% Last change - G. Clifford Sept 2009
%

% See following URLs for more info:
% http://adsabs.harvard.edu/abs/2007Chaos..17a5112K
% http://h-r-t.org/prsa/en/

% defaults
tp=20;
wl=30;
fsd=512;
time_stamp=[];

rr=rr(:);

if nargin < 6
   testmode = 0;
end

if nargin <5
       thresh_per = tp;
end

if nargin < 4
   win_length = wl; % seconds
end

if(isempty(thresh_per))
       thresh_per = tp;
end

if(isempty(win_length))
   win_length = wl; % seconds
end

threshlow = 1-thresh_per/100-.0001;
threshhigh = 1+thresh_per/100;

if nargin<3
    fs=512; %Hz
end

if nargin<3
    fs=fsd; %Hz
end

if(isempty(fs)==1)
  time_stamp = rri;
  fs=fsd;
end 

if nargin < 3
 time_stamp = rri.*fs;
end


% make time stamps
if nargin < 2
 rri = cumsum(rr).*fs;
end

% or
if(isempty(rri)==1)
 rri = cumsum(rr).*fs;
end

if(isempty(time_stamp)==1)
  time_stamp = rri*fs;
end 
  


% intialize
j =0;
k = 0; 
drrpi = [];
drrni = [];
acm = [];
dcm = [];
drr = rr(2:end)-rr(1:end-1);
drr_per = rr(2:end)./rr(1:end-1);

ac_anchor = (drr_per > threshlow) & (drr_per <= .9999); %defines AC achors, no changes greater then 5%
dc_anchor = (drr_per > 1) & (drr_per <= threshhigh);
ac_anchor = [0; ac_anchor];
dc_anchor = [0; dc_anchor];
ac_anchor(1:win_length) = 0;
ac_anchor(length(ac_anchor)-win_length+1:length(ac_anchor)) = 0;
dc_anchor(1:win_length) = 0;
dc_anchor(length(dc_anchor)-win_length+1:length(dc_anchor)) = 0;
ac_ind = find(ac_anchor);
dc_ind = find(dc_anchor);

for i = 1:length(dc_ind);
    dcm(i,:) = 1000*rr(dc_ind(i)-win_length:dc_ind(i)+win_length); %%convert from sec to msec
end

for i = 1:length(ac_ind);
    acm(i,:) = 1000*rr(ac_ind(i)-win_length:ac_ind(i)+win_length);
end

m_acm = mean(acm);
m_dcm = mean(dcm);

DC = (m_dcm(win_length+1)+m_dcm(win_length+2)-m_dcm(win_length)-m_dcm(win_length-1))/4;
AC = (m_acm(win_length+1)+m_acm(win_length+2)-m_acm(win_length)-m_acm(win_length-1))/4;




% do we plot?
if testmode == 1
    figure(1);
    plot(time_stamp,rr,'k+');
    hold on;
    figure(2)
    subplot(2,1,2)
    title('Deceleration');
    plot(dcm,'k--');
    hold on;

    figure(2)
    subplot(2,1,1)
    title('Acceleration');
    plot(acm,'k--');
    hold on;

figure
plot(time_stamp,rr,'k-+');
hold on;
plot(time_stamp(ac_ind),rr(ac_ind),'g+');
plot(time_stamp(dc_ind),rr(dc_ind),'r+');
legend('non-anchor','ac anchor','dc anchor');
title('RR anchors');

figure(3);
subplot(2,1,1)
plot([-2:2*win_length-2],dcm,'k--')
title('DC matrix')
hold on
plot([-2:2*win_length-2],m_dcm,'r');
subplot(2,1,2)
plot([-2:2*win_length-2],acm,'k--')
hold on
plot([-2:2*win_length-2],m_acm,'r');
title('AC matrix')

end