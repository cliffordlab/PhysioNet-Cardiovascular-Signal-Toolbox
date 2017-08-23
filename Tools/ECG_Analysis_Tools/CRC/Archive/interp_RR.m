function [t_out,x_out] =  interp_RR(t,x, samp_freq, interp_type)
 
%[t_out,x_out] =  interp_RR(t,x, samp_freq,interp_type);
% performs the following interpolations:
% interp_type == 1 - 'nearest' - nearest neighbor interpolation
% interp_type == 2 - 'linear'  - linear interpolation
% interp_type == 3 - 'spline'  - cubic spline interpolation
% interp_type == 4 - 'cubic'   - cubic interpolation
% interp_type == 5 - 'smooth cubic' - smooth cubic interpolation
% Defaults are: samp_freq=7Hz, interp_type = 1 (linear);
%
% G Clifford, Oxford University 2001. 
% Please send error reports to gari@mit.edu


format long e
 
% check if defaults are changed
if nargin < 4
   interp_type=1
end
 
if nargin < 3
   samp_freq=7
end
 
if nargin < 2
   error('You need at least a time vector and a signal')
end
 
if length(t)~=length(x)
     error('x and t must be the same length')
end

%find out starting times and finishing times
start_time = t(1);
end_time = t(length(t));
t=t';

% adjust so that we calculate the starting and finishing times 
% of the new time scale (if we run over the ends we may get garbage)

scale_start = ceil(samp_freq*start_time);
scale_end   = fix(samp_freq*end_time);

start  = scale_start/samp_freq;
finish = scale_end/samp_freq;

t_out = start:1/samp_freq:finish;

if (interp_type == 1)
  x_out = interp1(t,x,t_out,'linear'); 
end
if (interp_type == 2)
  x_out = interp1(t,x,t_out,'nearest'); 
end
if (interp_type == 3)
  x_out = interp1(t,x,t_out,'spline'); 
end
if (interp_type == 4)
  x_out = interp1(t,x,t_out,'cubic'); 
end
if (interp_type == 5)
  x_out = csaps(t,x,0.5,t_out); 
end











