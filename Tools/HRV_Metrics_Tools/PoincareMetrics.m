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
%     flag  : ...
%   OUTPUTS:
%       SD1           : standard  deviation  of  projection  of  the  PP  on  
%                       the line perpendicular to the line of identity (y=-x)
%       SD2           : standard deviation of the projection of the PP on 
%                       the line of identity (y=x)
%       SD1_SD2_ratio : SD1/SD2 ratio
% 
% 
%   Written by: Giulia Da Poian <giulia.dap@gmail.com>
%   REPO:       
%       https://github.com/cliffordlab/HRVToolbox1.0  
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%


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

end


% Brennan, Michael, Marimuthu Palaniswami, and Peter Kamen. 
% "Do existing measures of Poincare plot geometry reflect nonlinear 
% features of heart rate variability?." IEEE transactions
% on biomedical engineering 48.11 (2001): 1342-1347.