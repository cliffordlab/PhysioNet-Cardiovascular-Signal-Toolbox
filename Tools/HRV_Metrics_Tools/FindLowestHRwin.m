function windows_output = FindLowestHRwin(t_win, t, rr, num_segs, HRVparams)
%
%	windows_output = FindLowestHRwin(t_win, t, rr, num_segs, HRVparams)
%
%	OVERVIEW
%       Given time series data, calculate and return median HR within window
%       satisfying condition (e.g. lowest or highest median HR)
%
%	INPUT
%       t_win            - 
%       t                - time stamps for 'hr'; in units of seconds
%       rr          	 - input time series of RR intervals from QRS file
%       num_segs         - number of windows to find;
%                          i.e. if this == 3 and 'condition' == 'lowest', then return
%                          windows with three lowest median HR in order of
%                          lowest to highest
%       HRVparams        - struct of settings for hrv_toolbox analysis
%
%	OUTPUT
%          windows_output.t  - vector of start times for windows
%          windows_output.hr - vector of median HR values for windows
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%	AUTHORS
%		Erik Reinertsen <erikrtn@gmail.com>
%       (Adapted for use in HRV Toolbox)
%	COPYRIGHT (C) 2016 AUTHORS (see above)
%		This code (and all others in this repository) are under the GNU General Public License v3
%		See LICENSE in this repo for details.

% Convert RR intervals into HR
hr = 60 ./ rr;

% Set value of increment window (sec)
increment = HRVparams.increment;
windowlength = HRVparams.windowlength;

% Initialize struct to hold t and hr of each window
windows = [];

% Initialize start time
t_window_start = 0;

% Initialize index counter for storing results
ii = 1;

% Loop through time series until we reach the end
% NOTE: vectorizing is non-trivial because the data is not evenly
% sampled. Otherwise I'd just stagger the matrix and calculate
% medians down each column containing each window of data.
while t_window_start <= t(end) - windowlength + increment
    
    % Find indices of HR values in this window
    idx_window = find(t > t_window_start & t < t_window_start + windowlength);
    
    % Store the time
    windows(ii).t = t_window_start;
    
    % Calculate the median of these HR values
    windows(ii).hr = median(hr(idx_window));
    
    % Increment time by sliding window length (sec)
    t_window_start = t_window_start + increment;
    
    % Increment struct index
    ii = ii + 1;
end

% Find indices for top non-overlapping windows that meet conditions
[~, idx_sorted] = sort([windows.hr]');

% Remove indexes where the data is disqualified from analysis, ie t_win =
% NaN
idxN = find(isnan(t_win));
idx_sorted(idxN) = [];

% The first non-NaN sorted index is the window with the lowest median HR,
% so we definitely return that one, as long as it's corresponding t_win is
% not a NaN value
windows_output(1).t = [];

while isempty(windows_output(1).t)
    
    windows_output(1).t = windows(idx_sorted(1)).t;
    windows_output(1).hr = windows(idx_sorted(1)).hr;
    
    % Erase the first entry from 'idx_sorted'
    idx_sorted(1) = [];
end

% Until we've found all of the lowest median HR windows,
% go to the next ordered window, check if overlaps with any other
% windows; if not, add. If yes, erase it.
while length(windows_output) < num_segs

    % Check dt between all best windows and the next sorted window
   dt = [windows_output.t] - windows(idx_sorted(1)).t;

   % Check if any dt's are less than the window length;
   % any that are suggests new window overlaps with a best window
   bool_overlaps = abs(dt) < windowlength;

   % If all bools are zero, the candidate window doesn't overlap
   if sum(bool_overlaps) == 0
       % Append t and hr
       windows_output(end + 1).t = windows(idx_sorted(1)).t;
       windows_output(end).hr = windows(idx_sorted(1)).hr; % note 'end' is +1 compared to the line before
   end

   % If there exist any bools == 1, the candidate window does
   % overlap, so don't append t and hr

   % Erase candidate window from 'idx_sorted'
   idx_sorted(1) = [];
end

end % end function