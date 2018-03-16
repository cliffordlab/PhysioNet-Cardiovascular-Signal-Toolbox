function [ RR_gapFilled, t_gapFilled] = RR_Preprocessing_for_MSE_DFA(rr, t )
%
%    [ RR_gapFilled, t_gapFilled] = RR_Preprocessing_for_MSE_DFA(rr, t )
%
%   OVERVIEW:   This function deals with missing data in rr time series. 
%               When rr intervals is 'too' long, i.e., grather that 2 times
%               the median value of the 10 rr inetrvals before and after, 
%               this function will fill the gap with median values
%
%   INPUT:      rr            - a single row of rr interval data in seconds
%               t             - the time indices of the rr interval data 
%                                 (seconds)
%   OUTPUT:     RR_gapFilled   - RR intervals with gaps filled with median
%                                values
%               t_gapFilled    - time indices of RR_gapFilled
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
%   Written by Giulia Da Poian (giulia.dap@gmail.com) and Pradyumna Suresha
%   (pradyumna.suresha@gmail.com)
%
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

t_gapFilled = [];
RR_gapFilled = [];


NmbOfBeatsForMedian = 10;

tmp_RR_gapFilled = abs(diff(t));

ii = 0;
for idx = 1:length(tmp_RR_gapFilled)
    try
        medianValue = median(rr(idx-NmbOfBeatsForMedian: idx+NmbOfBeatsForMedian));
    catch 
        try
            medianValue = median(rr(idx-NmbOfBeatsForMedian: idx));
        catch
            try
                medianValue = median(rr(idx:idx+NmbOfBeatsForMedian));
            catch
                medianValue = median(rr(1:idx));
            end
        end
    end
    
    thresh = 2 * medianValue;
       
            
    if tmp_RR_gapFilled(idx) > thresh
       
       MissingBeats = floor(tmp_RR_gapFilled(idx)/medianValue); 
              
       RR_gapFilled = [RR_gapFilled ; rr(idx)];
       t_gapFilled = [t_gapFilled; t(idx)];
       ii=ii+1;
       % adding 'samps' points
       t_gapFilled = [t_gapFilled; ones(MissingBeats,1)];
       RR_gapFilled = [RR_gapFilled; ones(MissingBeats,1)];
       ii=ii+MissingBeats;
       t_gapFilled(end) = abs(t(idx+1) - rr(idx+1));
       RR_gapFilled(end) = medianValue;
       
       for k = 1: MissingBeats-1
       
           t_gapFilled(end-k) = abs(t_gapFilled(end-k+1) - RR_gapFilled(end-k+1));
           RR_gapFilled(end-k) = medianValue;
       
       end
       
       % chack if fisrt value in the gap is correct
       
       if abs(t_gapFilled(ii-MissingBeats+1) - t(idx)) < 0.5*medianValue
           t_gapFilled(ii-MissingBeats+1) = [];
           RR_gapFilled(ii-MissingBeats+1) = [];
           ii=ii-1;
           MissingBeats = MissingBeats -1;
       end   
       RR_gapFilled(ii-MissingBeats+1) = abs(t_gapFilled(ii-MissingBeats+1) - t(idx));

        

    else
        RR_gapFilled = [RR_gapFilled ; rr(idx)];
        t_gapFilled = [t_gapFilled; t(idx)];
        ii=ii+1;
    end

end

RR_gapFilled = [RR_gapFilled; rr(end)];
t_gapFilled = [t_gapFilled; t(end)];


