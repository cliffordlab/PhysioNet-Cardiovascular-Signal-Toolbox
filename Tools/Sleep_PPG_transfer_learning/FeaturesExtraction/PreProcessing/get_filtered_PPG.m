function ppg_filt = get_filtered_PPG(ppg,params,saveFlag,varargin)

% ppg_filt = get_filtered_PPG(ppg,params,saveFlag,varargin)
%
% INPUTS: 
%     ppg       : ppg signal
%     parmas    : struct with parameters used for processing cardiovascular
%                 signals
%     saveFlag  : 1 - save filtered signal
%                 0 - do not save filtered signal
%     vararging : (required when saveFlag = 1) signal name to save filterd
%                 signal
%                
%
% OUTPUTS : 
%     ppg_filt  : filtered ppg signal
%
% Written by Giulia Da Poian, 07-Dec-2018

if (saveFlag == 1 && nargin<4)
    saveFlag =0 ;
end

if ~isempty(varargin)
    SaveTo = params.writedata;
    SigName = varargin{1};
end
    

% remove baseline wander
ppg_filt = bw_filter(ppg,params.Fs);
% search better up-slope edge
domi_factor = dominant_slope(ppg_filt,params.Fs);
ppg_filt=ppg_filt*domi_factor;
ppg_filt(isnan(ppg_filt))=0;
% lowpass FIR filter to remove high frequency
ppg_filt=doFilter_FIR_lowpass(ppg_filt,params.Fs);

% save the pre-processed PPG waveform
if saveFlag == 1
    if ~exist([SaveTo filesep 'FiltSigs'],'dir')
        mkdir([SaveTo filesep 'FiltSigs']);
    end
    save([SaveTo filesep 'FiltSigs' filesep SigName '_filt.mat'],'ppg_filt');
end
