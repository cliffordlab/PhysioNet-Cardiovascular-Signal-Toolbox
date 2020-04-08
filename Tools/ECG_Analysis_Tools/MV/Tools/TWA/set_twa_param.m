function [res Param] = set_twa_param()
%OVERVIEW, This function sets the params, like analysis window size, needed
% to perform th twa analysis.
%
% INPUTS        MANDATORY           DESCRIPTION
%               none
%
% OUTPUTS       
%
%               res                 flag which indicates if script evaluated
%                                   successfully, 1 indicates parameters were succesfully set.
%
%               Param               A structure which contains the different paramters set, used in measuring TWAs.
%                                   Contains the following fields,     
%               .NumRuns            Perform 100 reshuffles to estimate gamma distribution for noise in T wave amplitude series.
%               .MethodForEctopy    Missing or abnormal beats are to be replaced with the average even or odd beat.
%               .Alignment          Specify alignment is to be performed over S-T offset interval.
%               .corrQRS            Correlation threshold with QRS template to determine clean QRS complex in beat
%               .corrT              Correlation threshold with S-T offset template to determine clean S-T segment in beat.
%               .Interval           Beats, window of analysis for the MMA.
%               .simple             If simple == true, odd and even average beats are just that - standard average beats.
%               .threshold          Confidence interval for significance.
%
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Ismail Sadiq
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  

res=0;

Param.NumRuns =100; % Perform 100 reshuffles to estimate gamma distribution for noise in T wave amplitude series.
Param.MethodForEctopy = 'replace';  % Missing or abnormal beats are to be replaced with the average even or odd beat.
Param.Alignment = 'st'; % Find max amplitude in S-T offset segment and use to sort leads by this amplitude. Specify alignment is to be performed over S-T offset interval.
Param.corrQRS = 0.90; % Correlation threshold with QRS template to determine clean QRS complex in beat.
Param.corrT = 0.9;  % Correlation threshold with S-T offset template to determine clean S-T segment in beat.
Param.Interval = 60;    % Beats, window of analysis for the MMA.
Param.simple = true;   % If simple == true, odd and even average beats are just that - standard average beats.
Param.threshold = '95'; % Confidence interval for significance.

res = 1;


