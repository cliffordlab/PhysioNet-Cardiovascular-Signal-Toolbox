function [SD1, SD2, SD1_SD2_ratio] = PoincareMetrics(RRints,flag)
%   [SD1, SD2, SD12] = PoincareMetrics(rr)
%   OVERVIEW:
%   The Poincaré plot analysis is a geometrical and nonlinear method to 
%   assess the dynamics of HRV. 
%   Poincaré plot is generated via plotting each RR interval (RR[n]) 
%   against the subsequent RR interval (RR[n+1]).
%   This function allows the user to calculate a number of Poincaré
%   features, which are detailed below.
%
%   INPUTS:    
%       rr  : Row vector of NN-intervals in seconds
%     flag  : (Optional) 1: show Poincaré plot, 
%                        0: do not show Poincaré plot
%   OUTPUTS:
%       SD1           : (ms) standard  deviation  of  projection  of  the  
%                       PP  on  the line perpendicular to the line of 
%                       identity (y=-x)
%       SD2           : (ms) standard deviation of the projection of the PP  
%                       on the line of identity (y=x)
%       SD1_SD2_ratio : SD1/SD2 ratio
% 
% 
%   Written by: Giulia Da Poian <giulia.dap@gmail.com>
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%

if nargin<2
    flag = 0; % default do not show plot
end

SDSD = std(diff(RRints));
SDRR = std(RRints);
SD1 = (1 / sqrt(2)) * SDSD; % measures the width of poincare cloud
SD2 = sqrt((2 * SDRR^2) - (0.5 * SDSD^2)); % measures the length of the poincare cloud

SD1_SD2_ratio = SD1/SD2;

% Poincare plot  
if flag
    ax1 = RRints(1:end-1);
    ax2 = RRints(2:end);
    scatter(ax1, ax2)
    xlabel('RR_n (s)')
    ylabel('RR_n+1 (s)')
end

% Convert to ms
SD1 = SD1 * 1000;
SD2 = SD2 * 1000;

end


% Brennan, Michael, Marimuthu Palaniswami, and Peter Kamen. 
% "Do existing measures of Poincare plot geometry reflect nonlinear 
% features of heart rate variability?." IEEE transactions
% on biomedical engineering 48.11 (2001): 1342-1347.