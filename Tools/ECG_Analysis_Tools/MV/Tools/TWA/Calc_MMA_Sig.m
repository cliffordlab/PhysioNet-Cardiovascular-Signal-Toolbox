function [VAlt, VAlt_Sig, noise_median, noise_95, VAltPt]= Calc_MMA_Sig(beat_matrix, Param, Align)
%OVERVIEW, This function computes the TWAs from the beats stored in
% beat_matrix. Performs the reshuffling to estimate a gamma distribution to
% TWA ampitudes from the reshuffled sequences. It estimates a noise
% threshold from the noise distribution using the value in Param.threshold.
% If the TWA estimate in the original series is greater than this threshold
% the TWA amplitude is considered statisitcally significant.
%
% INPUTS        MANDATORY       DESCRIPTION
%               beat_matrix     () array with beats for each stored.
%
%               Param           Structure contains parameters used for
%                               computing the TWAs.
%
%               Align           Structure contains the alignment
%                               information for beats in the current 
%                               analysis window.
%
% OUTPUTS
%               VAlt            TWA amplitude estimates for the analysis
%                               window.
%
%               VAlt_Sig        TWA amplitude estimates which are statistically
%                               significant compared to the noise threshold
%                               estimate.
%
%               noise_median    median of the gamma distribution used to model the noise.
%
%               noise_95        95th percentile of the gamma distribution
%                               used to model the noise.
%
%               VAltPt          The location of the point with maximum
%                               difference between the average even and odd beats.
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Shamim Nemati, 
%       editted by Ismail Sadiq on 10/26/2019. 
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  

[~ , num_beats] = size(beat_matrix);
num_beats = num_beats-mod(num_beats,2);

TOTAL_ATMPT=5; % Total # of attempts to find a significant difference between even and odd beats
Flg=true; atmpt=0; VAlt = zeros(1,TOTAL_ATMPT); VAltPt=zeros(1,TOTAL_ATMPT);

while (Flg && atmpt<TOTAL_ATMPT)
    atmpt=atmpt+1;

    [VAlt(atmpt),VAltPt(atmpt)] = Calc_MMA(beat_matrix, Param, Align);
    rand_VAlt = zeros(Param.NumRuns,1);
    rand_ind=zeros(1,num_beats);
    for trial = 1:Param.NumRuns
        even_beats = 2:2:num_beats;
        odd_beats = 1:2:num_beats;
        inde=randperm(length(even_beats));
        indo=randperm(length(odd_beats));
        rand_ind=[even_beats(inde) odd_beats(indo)];
        rand_beat_matrix = beat_matrix(VAltPt(atmpt),rand_ind);

        rand_VAlt(trial) = Calc_MMA(rand_beat_matrix, Param, Align);

    end
    [parmhat,parmci] = gamfit(rand_VAlt(:));
    P = gamcdf(VAlt(atmpt),parmhat(1),parmhat(2),parmci);
    if (P < 0.01*str2double(Param.threshold)) % 0.05 threshold for significent % 0.01?
       beat_matrix(VAltPt(atmpt),:)=0; % ignore, maybe noise
    else
       Flg=0;
    end
end

if (P < 0.01*str2double(Param.threshold)) % 0.05 threshold for significent
    VAlt_Sig=NaN;
    VAlt = VAlt(1);
    VAltPt = VAltPt(1);
else
    % If p value greater than threshold, consider significant alternan.
    VAlt_Sig = VAlt(atmpt);
    VAlt = VAlt(atmpt);
    VAltPt = VAltPt(atmpt);
end
qtl= quantile(rand_VAlt,[0.5 0.95]);
noise_median = qtl(1);
noise_95 = qtl(2);
return;
