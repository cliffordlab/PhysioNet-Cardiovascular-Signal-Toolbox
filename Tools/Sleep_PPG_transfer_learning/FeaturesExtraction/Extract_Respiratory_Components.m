function  RIAV = Extract_Respiratory_Components(ppg_sig,onsets,Fs)

%
% Use time-series of peaks and onsets to derive a new time-series which 
% represents information related to respiration, in particular respiratory 
% induced  amplitude  variation
%
% Inputs :
%    ppg_sig : 
%     onsets :
%         Fs :
%
% Outputs :
%       RIAV : respiratory-induced  amplitude  variation
%
% Written by Giulia Da Poian, 07-Dec-2018
%

obsWin = round(Fs*0.75); % find PPG peak within 75ms from onset

IBI = diff(onsets)./Fs; % Intra beat intervals 

for iOns = 1:length(onsets)
    % look for max peak after onset
    try
        [~, pos] = max(ppg_sig(onsets(iOns): onsets(iOns)+obsWin));
        ppgPeak(iOns) = pos + onsets(iOns)-1;
        RIAV(iOns) = ppg_sig(ppgPeak(iOns)) - ppg_sig(onsets(iOns)); 
    catch
        ppgPeak(iOns) = NaN;
        RIAV(iOns) = NaN;
    end 
end

