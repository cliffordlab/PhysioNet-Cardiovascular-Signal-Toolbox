function output = pnna(rr, alpha)
%
%   output = pnna(rr,alpha)
%
%	OVERVIEW:   finds the fraction of consecutive beats that differ by
%               more than a specified time.
%               Setting <alpha> to 50 msec will make this equivalent to pnn50
%               50 msec is the default value of <alpha>
%               When no arguments are given, the program documents itself%
%
%   INPUT:     
%               rr    - is a vector containing RR intervals in seconds
%               alpha - is a number which specifies the limit of successive 
%                       differences in milliseconds (50 ms by default) 
%   OUTPUT:     
%               pnna - Probability of intervals greater alpha ms or 
%                      smaller -alpha ms (e.g., alpha = 50 ms)
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   SOURCE:     Adapted from Kaplan Tools
%               See: http://www.robots.ox.ac.uk/~gari/CODE/HRV/
%	AUTHORS:    Adriana Vest <adriana.vest@gmail.com>
%	COPYRIGHT (C) 2016 AUTHORS (see above)
%   LICENSE:    This code (and all others in this repository) are under 
%               the GNU General Public License v3

% set default parameters
if nargin < 2
    % assume alpha is in millisecs
    alpha = 50;
end

dif = diff(rr); 

% pnnalpha is the fraction of differences larger than alpha

output = sum( abs(dif) >= alpha/1000 )/length(dif);


end