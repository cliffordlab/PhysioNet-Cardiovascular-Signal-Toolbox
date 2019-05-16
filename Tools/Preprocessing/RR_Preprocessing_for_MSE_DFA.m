function [ RR_gapFilled, t_gapFilled] = RR_Preprocessing_for_MSE_DFA(rr, tt_rr )
%
%    [ RR_gapFilled, t_gapFilled] = RR_Preprocessing_for_MSE_DFA(rr, t )
%
%   OVERVIEW:   This function deals with missing data in rr time series for 
%               long range HRV measure
%               
%               
%   INPUT:      rr            - a single row of rr interval data in seconds
%               t             - the time of the rr interval data 
%                               (seconds)
%   OUTPUT:     RR_gapFilled   - RR intervals with gaps filled with median
%                                values
%               t_gapFilled    - time indices of RR_gapFilled
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
%   Written by Giulia Da Poian (giulia.dap@gmail.com) and Pradyumna Suresha
%   (pradyumna.suresha@gmail.com)
%
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

t_gapFilled = tt_rr;
RR_gapFilled = rr;

% NOTE: 26 April 2018 This function is under development and will be released soon
