function scaledData = coarsegrain(data,tau,cg_moment,cg_method)
%   scaledData = coarsegrain(data,tau)
%
%   OVERVIEW:  generate the consecutive coarse-grained time series
%   
%   INPUTS:    
%     data      : Vector containing the original seires
%     tau       : scale factor (if tau = 1 returs the original series)
%     cg_moment : [optional] - moment used to coarse-grain the time series
%                 when using the traditional method,by default use the
%                  mean 'mean'. Options: 'mean', 'varaince'
%     cg_method : [optional] - method use to generate the coarse-grain 
%                 time series. By default use the original method proposed
%                 by Costa et al. ('fir' with mean).
%                 Options: 'fir' [Costa et al.], 'butter' [Valencia et al.]
%
%   OUTPUTS:
%     scaledData : the coarse-grained time series at the scale factor tau
%
%   Written by: Giulia Da Poian <giulia.dap@gmail.com>
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox 
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

if nargin < 2
    error('Wrong number of input arguments');
end
if nargin < 3 || isempty(cg_moment)
    cg_moment = 'mean';    
end
if nargin < 4 || isempty(cg_method)
    cg_method = 'fir';   
end


L = length(data);
J = fix(L/tau);




if tau < 1
    error('The scale factor must be a positive integer number');
elseif tau == 1
    switch cg_method
        case 'fir'
            scaledData = data; % return orignal series
        case 'butter'
            undec_scaledData = data;
            scaledData= data;        
    end
    return
else
    switch cg_method
        case 'fir' 
            for idx = 1:J
                switch cg_moment
                    case 'mean'
                        scaledData(idx) = mean(data((idx-1)*tau+1:idx*tau));
                    case 'variance'
                        scaledData(idx) = var(data((idx-1)*tau+1:idx*tau));
                    otherwise % still use mean
                        scaledData(idx) = mean(data((idx-1)*tau+1:idx*tau));                
                end
            end
            
        case 'butter'
            fs = 1;        % sample frequency 1Hz
            fn = fs/2;     % Nyquist frequency
            fc = fn/tau;   % cutoff frequency in function of the scale
            [b,a] = butter(6,fc/fn); % Butterworth low pass filter order 6,
            undec_scaledData = filter(b,a,data);     %series filtered and without decimate
            idxd=tau:tau:length(data);                       %index to do decimate
            scaledData = undec_scaledData(idxd);    %series filtered and with decimate
            clear a b idxd fc
        

    end
end


% REFERENCES

% Costa et al. "Multiscale entropy analysis of complex physiologic time
% series." Physical review letters 89.6 (2002): 068102

% Valencia et al. "Refined multiscale entropy: Application to 24-h holter 
% recordings of heart period variability in healthy and aortic stenosis
% subjects." IEEE Transactions on Biomedical Engineering 56, no. 9 (2009).

