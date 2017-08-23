function [ac_results, dc_results, windows_all] = prsa(rr, rri, sqi, windows_all, HRVparams)
%   [ac_results, dc_results, windows_all] = prsa(rr, rri, sqi, windows_all, HRVparams )
%
%   OVERVIEW:
%       Calculates acceleration and deceleration matrix  
%   INPUT
%       rr           - rr intervals
%       rri          - index of rr intervals
%       sqi          - 
%       windows_all  -
%       HRVparams    - struct of settings for hrv_toolbox analysis
%
%       empty arrays are permissable except for first argument, e.g.
%       [acm dcm ac dc] = prsa(rr, [], [], [], [], 1);
%   OUTPUT
%       ac_results   - ?
%       dc_results   - ?
%       windows_all  - ?
%
%   AUTHORS
%       L. Campana
%       G. Clifford
%   REFERENCE   
%       http://adsabs.harvard.edu/abs/2007Chaos..17a5112K
%       http://h-r-t.org/prsa/en/
%	REPO
%       https://github.com/cliffordlab/hrv_toolbox
%   DEPENDENCIES & LIBRARIES
%	COPYRIGHT (C) 2016 AUTHORS (see above)
%       This code (and all others in this repository) are under the GNU General Public License v3
%       See LICENSE in this repo for details.

% Feb 13, 2017
% Code altered by Adriana N Vest for the HRV Toolbox. Most of code
% unchanged except for HRVparams and Initializations Section and the
% inclusion of windowing options.

% These default parameters were a part of the original PRSA code developed
% by the authors above. 
% tp = 20; thresh_per
% wl = 30; prsa_win_length
% fsd = 512; fs
%%
% Make vector a column
rr = rr(:);

%% HRVparams and Initializations


if nargin < 1
    error('no data provided')
end
if nargin < 2 || isempty(rri)
    rri = cumsum(rr);
end
if nargin <3 || isempty(sqi)
    sqi(:,1) = rri;
    sqi(:,2) = 100* ones(length(rri),1);
end
if nargin<4 || isempty(windows_all)
    windows_all = CreateWindowRRintervals(rri, NN, HRVparams);
end
if nargin<5 || isempty(HRVparams)
    HRVparams = initialize_HRVparams;
end

time_stamp = rri;

prsa_win_length = HRVparams.prsa.win_length;
thresh_per = HRVparams.prsa.thresh_per;
plot_results = HRVparams.prsa.plot_results;
windowlength = HRVparams.windowlength;
prsa_threshold1 = HRVparams.prsa.threshold1;
prsa_threshold2 = HRVparams.prsa.threshold2;

threshlow = 1-thresh_per/100-.0001;
threshhigh = 1+thresh_per/100;

%% Run PRSA by Window
% Loop through each window of RR data
for i_win = 1:length(windows_all)
	% Check window for sufficient data
    if isnan(windows_all(i_win))

    end
    if ~isnan(windows_all(i_win))
        % Isolate data in this window
        idx_NN_in_win = find(rri >= windows_all(i_win) & rri < windows_all(i_win) + windowlength);
        idx_sqi_win = find(sqi(:,1) >= windows_all(i_win) & sqi(:,1) < windows_all(i_win) + windowlength);
        sqi_win = sqi(idx_sqi_win,:);
        nn_win = rr(idx_NN_in_win);
         
        % Analysis of SQI for the window
        lowqual_idx = find(sqi_win(:,2) < prsa_threshold1);

        % If enough data has an adequate SQI, perform the calculations
        if numel(lowqual_idx)/length(sqi_win(:,2)) < prsa_threshold2

            % intialize
            j =0;
            k = 0; 
            drrpi = [];
            drrni = [];
            acm = [];
            dcm = [];
            drr = nn_win(2:end)-nn_win(1:end-1);
            drr_per = nn_win(2:end)./nn_win(1:end-1);

            ac_anchor = (drr_per > threshlow) & (drr_per <= .9999); % defines ac anchors, no changes greater then 5%
            dc_anchor = (drr_per > 1) & (drr_per <= threshhigh);
            ac_anchor = [0; ac_anchor(:)];
            dc_anchor = [0; dc_anchor(:)];
            ac_anchor(1:prsa_win_length) = 0;                                        % sets first window to 0
            ac_anchor(length(ac_anchor)-prsa_win_length+1:length(ac_anchor)) = 0;    % sets last window to 0
            dc_anchor(1:prsa_win_length) = 0;                                        % sets first window to 0
            dc_anchor(length(dc_anchor)-prsa_win_length+1:length(dc_anchor)) = 0;    % sets last window to 0
            ac_ind = find(ac_anchor);
            dc_ind = find(dc_anchor);

            for i = 1:length(dc_ind)
                dcm(i,:) = 1000*nn_win(dc_ind(i)-prsa_win_length:dc_ind(i)+prsa_win_length); % convert from sec to msec
            end

            for i = 1:length(ac_ind)
                acm(i,:) = 1000*nn_win(ac_ind(i)-prsa_win_length:ac_ind(i)+prsa_win_length); % convert from sec to msec
            end

            m_acm = mean(acm);
            m_dcm = mean(dcm);
            
            % Edited by Adriana Vest: Added the following if statements to
            % deal with scenarios when ac or dc is not computable for a
            % particular window. An example scenario is when the RR
            % intervals are flat.
            if ~isnan(m_dcm)
                dc = (m_dcm(prsa_win_length+1)+m_dcm(prsa_win_length+2)-m_dcm(prsa_win_length)-m_dcm(prsa_win_length-1))/4;
                dc_results(i_win,1) = dc; % assign output of window
            else
                dc_results(i_win,1) = NaN;
            end
            if ~isnan(m_acm)
                ac = (m_acm(prsa_win_length+1)+m_acm(prsa_win_length+2)-m_acm(prsa_win_length)-m_acm(prsa_win_length-1))/4;
                ac_results(i_win,1) = ac; % assign output of window
            else
                ac_results(i_win,1) = NaN;
            end

            % Load custom colors
            custom_colors.red = [72 11 11] / 100;
            custom_colors.green = [30 69 31] / 100;

            % Plot results
            if plot_results == 1
                figure(1);
                plot(time_stamp,nn_win,'k+');
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
                plot(time_stamp,nn_win,'k-+');
                hold on;
                plot(time_stamp(ac_ind),nn_win(ac_ind), 'color', custom_colors.green, 'marker', '+');
                plot(time_stamp(dc_ind),nn_win(dc_ind), 'color', custom_colors.red, 'marker', '+');
                legend('non-anchor','ac anchor','dc anchor');
                title('RR anchors');

                figure(4);
                subplot(2,1,1)
                plot([-2:2*prsa_win_length-2],dcm,'k--')
                title('dc matrix')
                hold on
                plot([-2:2*prsa_win_length-2],m_dcm,'r');
                subplot(2,1,2)
                plot([-2:2*prsa_win_length-2],acm,'k--')
                hold on
                plot([-2:2*prsa_win_length-2],m_acm,'r');
                title('ac matrix')
            end % end of plotting condition
        else % else, if SQI is not adequate
            dc_results(i_win,1) = NaN;
            ac_results(i_win,1) = NaN;
        end % end of conditional statements run when SQI is adequate
    else % else, if window is NaN
        dc_results(i_win,1) = NaN;
        ac_results(i_win,1) = NaN;
    end % end of check for sufficient data
end % end of loop through windows

%prsaout(1).acm = acm;
%prsaout(1).dcm = dcm;


end % end of function
