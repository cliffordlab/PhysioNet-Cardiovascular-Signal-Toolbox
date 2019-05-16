function [ac_results, dc_results, prsa_ac, prsa_dc] = prsa(rr, t_rr, HRVparams, sqi, tWin)
%   [ac_results, dc_results] = prsa(rr, t_rr, HRVparams, swi, tWin )
%
%   OVERVIEW:
%       Calculates acceleration and deceleration capacity values  
%   INPUT
%       rr           - (seconds) rr intervals
%       t_rr         - (seconds) time of the rr intervals
%       HRVparams    - struct of settings for hrv_toolbox analysis
%       sqi          - Signal Quality Index; Requires a matrix with
%                      at least two columns. Column 1 should be
%                      timestamps of each sqi measure, and Column 2
%                      should be SQI on a scale from 0 to 1.
%       tWin         - Starting time of each windows to analyze 
%
%   OUTPUT
%       ac_results   - Acceleration Capacity value
%       dc_results   - Deceleration Capacity value
%       prsa_ac      - PRSA signal for AC anchor points
%       prsa_dc      - PRSA signal for DC anchor points
%
%   AUTHORS
%       L. Campana
%       G. Clifford
%   REFERENCE   
%       http://adsabs.harvard.edu/abs/2007Chaos..17a5112K
%       http://h-r-t.org/prsa/en/
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   DEPENDENCIES & LIBRARIES
%	COPYRIGHT (C) 2016 AUTHORS (see above)
%       This code (and all others in this repository) are under the GNU General Public License v3
%       See LICENSE in this repo for details.

% Feb 13, 2017
% Code altered by Adriana N Vest for the HRV Toolbox. Most of code
% unchanged except for HRVparams and Initializations Section and the
% inclusion of windowing options.
% Oct 13, 2017
% Code altered by Giulia Da Poian to compute AC and DC value at different
% scale s 
%
% These default parameters were a part of the original PRSA code developed
% by the authors above. 
% tp = 20; thresh_per
% wl = 30; prsa_win_length
% fsd = 512; fs
%%
% Make vector a column
rr = rr(:);

if nargin < 5
    error('no data provided')
end
if isempty(sqi)
    sqi(:,1) = t_rr;
    sqi(:,2) = ones(length(t_rr),1);
end

prsaWinLength = HRVparams.prsa.win_length;
s = HRVparams.prsa.scale;
PRSA_th = HRVparams.prsa.thresh_per; 
plot_results = HRVparams.prsa.plot_results;
windowlength = HRVparams.windowlength;
SQI_th = HRVparams.sqi.LowQualityThreshold;        % SQI threshold
WinQuality_th = HRVparams.RejectionThreshold; % Low quality windows threshold

Anchor_Low_th = 1-PRSA_th/100-.0001; % The lower limit for the AC anchor selection
Anchor_High_th = 1+PRSA_th/100;      % The upper limit for the DC anchor selection

% Preallocation

ac_results = nan(length(tWin),1);
dc_results = nan(length(tWin),1);
prsa_ac = [];
prsa_dc = [];

% Run PRSA by Window
% Loop through each window of RR data
for i_win = 1:length(tWin)
	% Check window for sufficient data
    if isnan(tWin(i_win))

    end
    if ~isnan(tWin(i_win))
        % Isolate data in this window
        sqi_win = sqi( sqi(:,1) >= tWin(i_win) & sqi(:,1) < tWin(i_win) + windowlength,:);
        nn_win = rr( t_rr >= tWin(i_win) & t_rr < tWin(i_win) + windowlength );
         
        % Analysis of SQI for the window
        lowqual_idx = find(sqi_win(:,2) < SQI_th);

        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < WinQuality_th

            % intialize
            acm = [];
            dcm = [];
            drr_per = nn_win(2:end)./nn_win(1:end-1);

            ac_anchor = (drr_per > Anchor_Low_th) & (drr_per <= .9999); % defines ac anchors, no changes greater then 5%
            dc_anchor = (drr_per > 1) & (drr_per <= Anchor_High_th);
            ac_anchor = [0; ac_anchor(:)];
            dc_anchor = [0; dc_anchor(:)];
            ac_anchor(1:prsaWinLength) = 0;                                        % sets first window to 0
            ac_anchor(length(ac_anchor)-prsaWinLength+1:length(ac_anchor)) = 0;    % sets last window to 0
            dc_anchor(1:prsaWinLength) = 0;                                        % sets first window to 0
            dc_anchor(length(dc_anchor)-prsaWinLength+1:length(dc_anchor)) = 0;    % sets last window to 0
            ac_ind = find(ac_anchor);
            dc_ind = find(dc_anchor);

            for i = 1:length(dc_ind)
                dcm(i,:) = 1000*nn_win(dc_ind(i)-prsaWinLength:dc_ind(i)+prsaWinLength-1); % convert from sec to msec
            end
            
            for i = 1:length(ac_ind)
                acm(i,:) = 1000*nn_win(ac_ind(i)-prsaWinLength:ac_ind(i)+prsaWinLength-1); % convert from sec to msec
            end

            prsa_ac = mean(acm,1);
            prsa_dc = mean(dcm,1);
            
            % Edited by Adriana Vest: Added the following if statements to
            % deal with scenarios when ac or dc is not computable for a
            % particular window. An example scenario is when the RR
            % intervals are flat.
            if (~isnan(sum(prsa_dc)) && numel(dc_ind)>= HRVparams.prsa.min_anch) % Edited by Giulia Da Poian
                dc = (sum(prsa_dc(prsaWinLength+1:prsaWinLength+s)) - sum(prsa_dc(prsaWinLength-(s-1):prsaWinLength))) ./ (2*s);
                dc_results(i_win,1) = dc; % assign output of window
            end
            if ~isnan(sum(prsa_ac)) && numel(ac_ind)>= HRVparams.prsa.min_anch % Edited by Giulia Da Poian
                ac = (sum(prsa_ac(prsaWinLength+1:prsaWinLength+s)) - sum(prsa_ac(prsaWinLength-(s-1):prsaWinLength))) ./ (2*s);
                ac_results(i_win,1) = ac; % assign output of window
            end

            % Load custom colors
            custom_colors.red = [72 11 11] / 100;
            custom_colors.green = [30 69 31] / 100;

            % Plot results
            if plot_results == 1
                figure(1);
                plot(t_rr,nn_win,'k+');
                hold on;

                figure(2)
                subplot(2,1,2)
                plot(dcm,'k--');
                legend('Deceleration');
                hold on;

                figure(2)
                subplot(2,1,1)
                plot(acm,'k--');
                legend('Acceleration');    
                hold on;

                figure(3)
                plot(t_rr,nn_win,'k-+');
                hold on;
                plot(t_rr(ac_ind),nn_win(ac_ind), 'color', custom_colors.green, 'marker', '+');
                plot(t_rr(dc_ind),nn_win(dc_ind), 'color', custom_colors.red, 'marker', '+');
                legend('non-anchor','ac anchor','dc anchor');
                title('RR anchors');

                figure(4);
                subplot(2,1,1)
                plot([-2:2*prsaWinLength-2],dcm,'k--')
                title('dc matrix')
                hold on
                plot([-2:2*prsaWinLength-2],prsa_dc,'r');
                subplot(2,1,2)
                plot([-2:2*prsaWinLength-2],acm,'k--')
                hold on
                plot([-2:2*prsaWinLength-2],prsa_ac,'r');
                title('ac matrix')
            end % end of plotting condition
        else % else, if SQI is not adequate
        end % end of conditional statements run when SQI is adequate
    else % else, if window is NaN
    end % end of check for sufficient data
end % end of loop through windows


end % end of function
