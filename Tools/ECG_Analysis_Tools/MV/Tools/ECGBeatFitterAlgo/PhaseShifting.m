function phase = PhaseShifting(phasein,teta)
%
% phase = PhaseShifting(phasein,teta),
% Phase shifter.
%
% inputs:
% phasein: calculated ECG phase
% teta: desired phase shift. teta>0 and teta<0 corresponds with phase leads
% and phase lags, respectively.
%
% output:
% phase: the shifted phase.
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

phase = phasein + teta;
phase = mod(phase + pi, 2*pi) - pi;
