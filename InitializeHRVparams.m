function HRVparams = InitializeHRVparams(project_name)
%
%   settings = InitializeHRVparams('project_name')
%
%   OVERVIEW:   
%       This file stores settings and should be configured before
%       each use of the HRV toolbox:
%       1.  Project Specific Input/Output Data type and Folders
%       2.  How much does the user trust the data
%       3.  Global Settings (Window Size) for signal segmentation
%       4.  Quality Threshold Settings
%       5.  Debug Settings
%       6.  SQI Settings
%       7.  Output Settings 
%       8.  Time of Process and Filename to Save Data
%       9.  Preprocess Settings
%       10. AF Detection Settings
%       11. Time Domain Analysis Settings
%       12. Frequency Domain Analysis Settings
%       13. SDANN and SDNNI Analysis Settings
%       14. PRSA Analysis Settings
%       15. Peak Detection Settings
%       16. Entropy  and Multiscale Entropy - MSE - Settings
%       17. Detrended Fluctuation Analysis - DFA - Settings
%       18. Poincaré plot - Settings
%
%   INPUT:      
%       project_name = a string with the name of the project - this
%       will determine the naming convention of file folders 
%
%   OUTPUT:
%       HRVparams - struct of various settings for the hrv_toolbox analysis
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       This script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%


% Set up input options for a specific project
switch project_name   
    % Define new project name and parameters
    case 'MyProjectName'                   % Update with your project name
        HRVparams.Fs = Nan;                % Spacify sampling frequency
        HRVparams.readdata = '';           % (Optional) Specify name for data input folder
        HRVparams.writedata = '';          % (Optional) Specify name for data output folder
        HRVparams.datatype = '';           % (Optional) Spacify Data type of input
        HRVparams.ext = '';                % (Optional) Spacify file extension of input (e.g., 'mat','qrs')
        
    % Existing demo projects
    case 'demo_NSR'      % Parameters for demo using MIT nsr data
        HRVparams.readdata = strcat('TestData', filesep, 'Physionet_nsr2db');
        HRVparams.writedata = strcat('OutputData', filesep, 'ResultsNSR');
        HRVparams.ext = 'ecg';
        HRVparams.Fs = 125;
    case 'demo_RRgen' 
        HRVparams.writedata = strcat('OutputData', filesep, 'ResultsRRgen');
        HRVparams.Fs = 125;
        HRVparams.demo.length = 5*60*60;   % Length of demo RR intervals in seconds
        HRVparams.demo.pe = 0.0003;        % Probability of ectopy ~ 1 per hr 
        HRVparams.demo.pn = 0.0048;        % Probability of noise ~ 16 per hr 
        HRVparams.demo.seed = 1;           % Seed for RRGEN
    case 'demoICU'      % Parameters for ICU demo 
        HRVparams.readdata = strcat('TestData');
        HRVparams.writedata = strcat('OutputData', filesep, 'ResultsICU');
        HRVparams.ext = 'mat';
        HRVparams.Fs = 128;
    case 'demoAF'      % Parameters for AF demo 
        HRVparams.readdata = strcat('TestData');
        HRVparams.writedata = strcat('OutputData', filesep, 'ResultsAFData');
        HRVparams.ext = 'mat';
        HRVparams.Fs = 128;
           
end

if  isempty(HRVparams.writedata)    
    % Default data OUTPUT folder name based on project name
    HRVparams.writedata = strcat(project_name,'_Results');  
    fprintf('New OUTPUT folder: "%s"\n', HRVparams.writedata)
    mkdir(HRVparams.writedata);          % Create output folder and 
elseif ~exist([pwd filesep HRVparams.writedata], 'dir')
    fprintf('New OUTPUT folder: "%s"\n',HRVparams.writedata)
    mkdir(HRVparams.writedata);          % Create output folder and 
end
addpath(genpath(HRVparams.writedata));   % Add folder to search path

%% 2. How much does the user trust the data:
% This setting determines how stringently filtered the data is prior to
% calculating HRV metrics. Raw ECG or Pulse data that has been labeled
% using a peak detector (like jqrs) would require the most stringent
% filtering, whereas RR interval data that has been reviewed by a human
% technician would require the least amount of filtering. 

%   - qrs detection no beats labeled - most stringent filters
%   - automatic beat detection beat typed - moderately stringent filtered
%   - hand corrected - least filtered, maybe no filter

% EXPLAIN THE THRESHOLDS BETTER
% ADD THIS TO A DEMO

HRVparams.data_confidence_level = 1;
% 1 - raw data with automatic beat detection
% 2 - raw data with automatic beat detection, but beat typed (ie N, SV,etc)
% 3 - technician reviewed data

% * ^^^^ NOT YET IN USE ^^^^  *

%% 3. Global Settings (Window Size)

HRVparams.windowlength = 300;	% Default: 300, seconds
HRVparams.increment = 30;       % Default: 30, seconds increment
HRVparams.numsegs = 5;          % Default: 5, number of segments to collect with lowest HR

%% 4. Quality Threshold Settings
HRVparams.threshold1 = 0.9;       % Default: 0.9, Threshold for which SQI represents good data
HRVparams.threshold2 = .20;       % Default: 0.2, amount (%) of data that can be rejected before a
                                  % window is considered too low quality for analysis
HRVparams.win_tol = .15;          % Default: 0.15, maximum percentage of data allowable to be missing
                                  % from a window .15 = 15%
%% 5. Debug Settings

HRVparams.rawsig = 0;           % Load raw signal if it is available for debugging
HRVparams.debug = 0;

%% 6. SQI Settings

HRVparams.sqi.windowlength = 10; % Default: 10, seconds
HRVparams.sqi.increment = 1;     % Default: 1, seconds
HRVparams.sqi.threshold = 0.1;   % Default: 0.1, seconds
HRVparams.sqi.margin = 2;        % Default: 2, seconds

%% 7. Output Settings

HRVparams.gen_figs = 0;             % Generate figures
HRVparams.save_figs = 0;            % Save generated figures
if HRVparams.save_figs == 1
    HRVparams.gen_figs = 1;
end

% Format settings for HRV Outputs
HRVparams.output.format = 'csv';        % 'csv' - creates csv file for output
                                        % 'mat' - creates .mat file for output
HRVparams.output.separate = 1;          % Default : 1 = separate files for each subject
                                        %           0 = all results in one file
HRVparams.output.num_win = [];          % Specify number of lowest hr windows returned
                                        % leave blank if all windows should be returned
                                        % Format settings for annotations generated
HRVparams.output.ann_format = 'binary'; % 'binary'  = binary annotation file generated
                                        % 'csv'     = ASCII CSV file generated
                            

%% 8. Filename to Save Data
HRVparams.time = datestr(now, 'yyyymmdd');              % Setup time for filename of output
HRVparams.filename = [HRVparams.time '_' project_name];

%% 9. Preprocess Settings

HRVparams.preprocess.figures = HRVparams.gen_figs;  % Figures on = 1, Figures off = 0
HRVparams.preprocess.gaplimit = 4;                  % Default: 2, seconds; maximum believable gap in rr intervals
HRVparams.preprocess.per_limit = 0.2;               % Default: 0.2, Percent limit of change from one interval to the next
HRVparams.preprocess.forward_gap = 3;	            % Default: 3, Maximum tolerable gap at beginning of timeseries in seconds
HRVparams.preprocess.method_outliers = 'rem';       % Default: 'rem', Method of dealing with outliers
                                                    % 'cub' = replace outlier points with cubic spline method
                                                    % 'rem' = remove outlier points
                                                    % 'pchip' = replace with pchip method
HRVparams.preprocess.lowerphysiolim = 60/160;       % Default: 60/160
HRVparams.preprocess.upperphysiolim = 60/30;        % Default: 60/30
HRVparams.preprocess.method_unphysio = 'rem';       % Default: 'rem', Method of dealing with unphysiologically low beats
                                                    % 'cub' = replace outlier points with cubic spline method
                                                    % 'rem' = remove outlier points
                                                    % 'pchip' = replace with pchip method

% The following settings do not yet have any functional effect on 
% the output of preprocess.m:                             
HRVparams.preprocess.threshold1 = 0.9 ;	        % Default: 0.9, Threshold for which SQI represents good data
HRVparams.preprocess.minlength = 30;            % Default: 30, The minimum length of a good data segment in seconds
                                
%% 10. AF Detection Settings

HRVparams.af.on = 1;              % Default: 1, AF Detection On or Off
HRVparams.af.windowlength = 30;   % Default: 30, Set to include ~30 beats in each window
HRVparams.af.increment = 30;      % Default: 30, No overlap necessary in AF feat calc

%% 11. Time Domain Analysis Settings

HRVparams.timedomain.on = 1;                             % Default: 1, Time Domain Analysis 1=On or 0=Off

HRVparams.timedomain.threshold1 = 0.9;                   % Default: 0.9, Threshold for which SQI represents good data
HRVparams.timedomain.threshold2 = 0.20;                  % Default: 0.20, Amount (%) of data that can be rejected before a
                                                         % window is considered too low quality for analysis
HRVparams.timedomain.dataoutput = 0;                     % 1 = Print results to .txt file
                                                         % Anything else = utputs to return variables only
                                                         % returned variables
HRVparams.timedomain.alpha = 0.050;                      % Default: 0.05 ,In seconds
HRVparams.timedomain.win_tol = .15;                      % Default: .15, Maximum percentage of data allowable 
                                                         % to be missing from a window

%% 12. Frequency Domain Analysis Settings

HRVparams.freq.on = 1;              % Default: 1, Frequency Domain Analysis 1=On or 0=Off

ULF = [0 .0033];                    % 
VLF = [0.0033 .04];                 % Requires at least 300 s window
LF = [.04 .15];                     % Requires at least 25 s window
HF = [0.15 0.4];                    % Requires at least 7 s window

HRVparams.freq.limits = [ULF; VLF; LF; HF];
HRVparams.freq.threshold1 = 0.9;	        % Default: 0.9, Threshold for which SQI represents good data
                                            % (Used only when SQI is provided)
HRVparams.freq.threshold2 = 0.05;           % Default: 0.05, Amount (%) of data that can be rejected before a
                                            % window is considered too low quality for analysis
                                            % (Used only when SQI is provided)
HRVparams.freq.zero_mean = 1;               % Default: 1, Option for subtracting the mean from the input data
HRVparams.freq.method = 'lomb';             % Default: 'lomb' 
                                            % Options: 'lomb', 'burg', 'fft', 'welch' 
HRVparams.freq.plot_on = 0;

% The following settings are for debugging spectral analysis methods
HRVparams.freq.debug_sine = 0;              % Default: 0, Adds sine wave to tachogram for debugging
HRVparams.freq.debug_freq = 0.15;           % Default: 0.15
HRVparams.freq.debug_weight = .03;          % Default: 0.03

% Lomb:
HRVparams.freq.normalize_lomb = 0;	        % Default: 0
                                            % 1 = Normalizes Lomb Periodogram, 
                                            % 0 = Doesn't normalize 

% Burg: (not recommended)
HRVparams.freq.burg_poles = 15;    % Default: 15, Number of coefficients 
                                   % for spectral estimation using the Burg
                                   % method (not recommended)

% The following settings are only used when the user specifies spectral
% estimation methods that use resampling : 'welch','fft', 'burg'
HRVparams.freq.resampling_freq = 7;             % Default: 7, Hz 
HRVparams.freq.resample_interp_method = 'cub';  % Default: 'cub'
                                                % 'cub' = cublic spline method
                                                % 'lin' = linear spline method
HRVparams.freq.resampled_burg_poles = 100;      % Default: 100

%% 13. SDANN and SDNNI Analysis Settings
HRVparams.sd.on = 1;                        % Default: 1, SD analysis 1=On or 0=Off
HRVparams.sd.segmentlength = 300;           % Default: 300, windows length in seconds

%% 14. PRSA Analysis Settings

HRVparams.prsa.on = 1;             % Default: 1, PRSA Analysis 1=On or 0=Off
HRVparams.prsa.win_length = 30;    % Default: 30, The length of the PRSA signal 
                                   % before and after the anchor points
                                   % (the resulting PRSA has length 2*L)
HRVparams.prsa.thresh_per = 20;    % Default: 20%, Percent difference that one beat can 
                                   % differ from the next in the prsa code
HRVparams.prsa.plot_results = 0;   % Default: 0                            
HRVparams.prsa.scale = 2;          % Default: 2, scale parameter for wavelet analysis (to compute AC and DC)

%% 15. Peak Detection Settings

% The following settings are for jqrs.m

HRVparams.REF_PERIOD = 0.250;   % Default: 0.25
HRVparams.THRES = .6;           % Default: 0.6
HRVparams.fid_vec = [];         % Default: []
HRVparams.SIGN_FORCE = [];      % Default: []
HRVparams.debug = 0;            % Default: 0

%% 16. Entropy Settings
% Multiscale Entropy
HRVparams.MSE.on = 1;                      % Default: 1, MSE Analysis 1=On or 0=Off
HRVparams.MSE.windowlength = 60*60*0.5;           % Default: [], perform DFA on the entair signal
HRVparams.MSE.increment = 60*60*0.5;              % Default: [];
HRVparams.MSE.RadiusOfSimilarity = 0.15;   % Default: 0.15, Radius of similarity (% of std)
HRVparams.MSE.patternLength = 2;           % Default: 2, 
HRVparams.MSE.maxCoarseGrainings = 20;     % Default: 7, Maximum number of coarse-grainings
% SampEn an ApEn 
HRVparams.Entropy.on = 1;                     % Default: 1, MSE Analysis 1=On or 0=Off
HRVparams.Entropy.RadiusOfSimilarity = 0.15;  % Default: 0.15, Radius of similarity (% of std)
HRVparams.Entropy.patternLength = 2;          % Default: 2, 

%% 17. DFA Settings

HRVparams.DFA.on = 1;               % Default: 1, DFA Analysis 1=On or 0=Off
HRVparams.DFA.windowlength = 60*60*0.5;    % Default [], perform DFA on the entair signal
HRVparams.DFA.increment = 60*60*0.5;       % Default: [];
HRVparams.DFA.minBoxSize = 4 ;      % Default: 4, Smallest box width
HRVparams.DFA.maxBoxSize = [];      % Largest box width (default in DFA code: signal length/4) 
%

%% 18. Poincaré plot

HRVparams.poincare.on = 1;     % Default: 1, Poincare Analysis 1=On or 0=Off

%% Export Parameter as Latex Table
% Note that if you change the order of the parameters or add parameters 
% this might not work

ExportHRVparams2Latex(HRVparams);


end