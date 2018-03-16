function [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSD, F, limits,plot_on)
% [ulf, vlf, lf, hf, lfhf, ttlpwr] = CalcLfHfParams(PSD, F, limits,plot_on)
%
%   OVERVIEW: Compute the frequency domain features for a given PSD and
%             frequency bans limits
%         
%   INPUT:      
%        PSD     - power spectral density 
%        F       - frequency vector
%        limits  - frequency domain analysis limits
%        plot_on - 
%
%   OUTPUT:     
%	- ulf     : (ms^2) Power in the ultra low frequency range (default < 0.003 Hz)
%	- vlf     : (ms^2) Power in very low frequency range (default 0.003 <= vlf < 0.04 Hz)
%	- lf      : (ms^2) Power in low frequency range (default 0.04Hz  <= lf < 0.15 Hz)
%	- hf      : (ms^2) Power in high frequency range (default 0.15 <= hf < 0.4 Hz)
%	- lfhf    : Ratio LF [ms^2]/HF [ms^2]
%	- ttlpwr  : (ms^2) Total spectral power (approximately <0.4 Hz)
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
%%
if nargin <3
    ULF = [0 .003];
    VLF = [0.003 .04];
    LF = [.04 .15];
    HF = [0.15 0.4];
    limits = [ULF; VLF; LF; HF];
end
if nargin < 4
    plot_on =1;
end

Indx_ULF = find( (limits(1,1) <= F) & (F <= limits(1,2)) );
Indx_VLF = find( (limits(2,1) <= F) & (F <= limits(2,2)) );
Indx_LF = find( (limits(3,1) <= F) & (F <= limits(3,2)) );
Indx_HF = find( (limits(4,1) <= F) & (F <= limits(4,2)) );
space = F(2)-F(1);

ulf = sum(PSD(Indx_ULF)*space) * 1e6; % convert to ms^2
vlf = sum(PSD(Indx_VLF)*space) * 1e6; % convert to ms^2
lf = sum(PSD(Indx_LF)*space) * 1e6;   % convert to ms^2
hf = sum(PSD(Indx_HF)*space) * 1e6;   % convert to ms^2

ttlpwr = sum([ulf vlf lf hf]);

lf_n = lf/ttlpwr; % normalized
hf_n = hf/ttlpwr;
lfhf = round(lf_n/hf_n*100)/100; % lf/hf ratio

if plot_on
    figure
    % plot PSD
    plot(F,10*log10(PSD),'b','linewidth',2)
    hold on
    % plot limits on graph for lf and hf
    plot([F(Indx_LF(1)) F(Indx_LF(1))],[-80 40],'k:')
    hold on
    plot([F(Indx_LF(end)) F(Indx_LF(end))],[-80 40],'k:')
    hold on
    plot([F(Indx_HF(end)) F(Indx_HF(end))],[-80 40],'k:')

    % labelsc
    text(0.07,30,'LF','Fontname','Times New Roman','Fontsize',10)
    text(0.25,30,'HF','Fontname','Times New Roman','Fontsize',10)
    %text(0.15, 35, 'Power Spectral Density','Fontname','Times New Roman','Fontsize',10)
    text(0.3, -60, strcat('LF/HF=',num2str(lfhf)),'Fontname','Times New Roman','Fontsize',10)
    ylabel('Normalized PSD (db/Hz)','Fontname','Times New Roman','fontsize',10)
    xlabel('Frequency (Hz)','Fontname','Times New Roman','fontsize',10)
    axis([0 .45 -80 40]);
    box off
end % end plot

end % end function
