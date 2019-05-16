function windowRRintervals = CreateWindowRRintervals(tNN, NN, HRVparams,option)
%
% windowRRintervals = CreateWindowRRintervals(NN, tNN, settings, options)
%   
%   OVERVIEW:   This function returns the starting time (in seconds) of 
%               each window to be analyzed.
%
%   INPUT:      tNN       : (seconds) a single row of time of the rr interval 
%               NN        : (seconds) a single row of NN (normal normal) interval
%                           data
%               HRVparams : struct of settings for hrv_toolbox analysis
%               option    : 'normal', 'af', 'sqi', 'mse', 'dfa'
%               
%   OUTPUT:     windowRRintervals : array containing the starting time 
%                                   (in seconds) of each window to be
%                                   analyzed
%                                   
%
%   DEPENDENCIES & LIBRARIES:
%       PhysioNet Cardiovascular Signal Toolbox
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   REFERENCE: 
%   Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%   Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Main script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)   
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

% Verify input arguments

if nargin< 1
    error('Need to supply time to create windows')
end
if nargin<4 || isempty(HRVparams) 
     option = 'normal';
end

% Set Defaults

increment = HRVparams.increment;
windowlength = HRVparams.windowlength;
win_tol = HRVparams.MissingDataThreshold;

switch option
    case 'af'
        increment = HRVparams.af.increment;
        windowlength = HRVparams.af.windowlength;
    case 'mse'
        increment = HRVparams.MSE.increment * 3600;
        windowlength = HRVparams.MSE.windowlength * 3600;
        if isempty(increment)
            windowRRintervals = 0;
            return % no need to crate windows , use entair signal
        end 
    case 'dfa'
        increment = HRVparams.DFA.increment * 3600;
        windowlength = HRVparams.DFA.windowlength * 3600;
        if isempty(increment)
            windowRRintervals = 0;
            return  % no need to crate windows , use entair signal
        end       
    case 'sqi'
        increment = HRVparams.sqi.increment;
        windowlength = HRVparams.sqi.windowlength;
    case 'HRT'
        increment = HRVparams.HRT.increment;
        windowlength = HRVparams.HRT.windowlength * 3600;
        if windowlength > tNN(end) 
            windowRRintervals = 0;
            return;
        end
            
end

nx = floor(tNN(end));               % length of sequence
overlap = windowlength-increment;   % number of overlapping elements
Nwinds = fix((nx-overlap)/(windowlength-overlap));    % number of sliding windows
% Initialize output matrix
windowRRintervals = (0:(Nwinds-1))*(windowlength-overlap);  % starting index of each windows

% Initialize loop variables
t_window_start = 0;     % Window Start Time
i = 1;                       % Counter

% for j = 1:(floor(tNN(end))-(windowlength))/increment
%     indicies = [];
%     try
%         indicies = find(tNN >= (j-1)*increment & tNN < windowlength+(j-1)*increment);
%         % windows_all: first coloumn is start time of sliding window
%         % and second column is end time of sliding window
%         windows_all(j,1) = tNN(indicies(1));
%         windows_all(j,2) = tNN(indicies(end));
%     catch
%         windows_all(j,1) = NaN;
%         windows_all(j,2) = NaN;
%     end
% end

if ~strcmp(option,'af') && ~strcmp(option,'sqi')
    
    while t_window_start <= tNN(end) - windowlength + increment

        % Find indices of time values in this segment
        t_win = tNN(tNN >= t_window_start & tNN < t_window_start + windowlength);

        % if NN intervals are supplied, assign them to the current window
        % if not, put in place a vector of ones as a place holder
        if ~isempty(NN)
            nn_win = NN(tNN >= t_window_start & tNN < t_window_start + windowlength);
        else
            nn_win = (windowlength/length(t_win))* ones(length(t_win),1);
        end

        % Store the begin time of window
        windowRRintervals(i) = t_window_start;

        % Increment time by sliding segment length (sec)
        t_window_start = t_window_start + increment;

        % Check Actual Window Length and mark windows that do not meet the
        % crieria for Adequate Window Length
        % First remove unphysiologic beats from candidates for this
        % measurement:
        idxhi = find(nn_win > HRVparams.preprocess.upperphysiolim);
        idxlo = find(nn_win < HRVparams.preprocess.lowerphysiolim);
        comb = [idxhi(:) ; idxlo(:)];
        nn_win(comb) = [];
        % Now query the true length of the window by adding up all of the NN
        % intervals
        truelength = sum(nn_win(:));
        if truelength < (windowlength * (1-win_tol))
                windowRRintervals(i) = NaN; 
        end

        % Increment loop index
        i = i + 1;
    end
end

end % end of function