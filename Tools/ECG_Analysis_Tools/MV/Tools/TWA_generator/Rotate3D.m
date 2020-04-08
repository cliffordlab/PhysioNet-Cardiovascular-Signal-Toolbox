function R = Rotate3D(tetax,tetay,tetaz);
%
% R = Rotate3D(tetax,tetay,tetaz)
% 3D rotation matrix
%
% inputs:
% tetax: rotation angle around the x axis (in rads.)
% tetay: rotation angle around the y axis (in rads.)
% tetaz: rotation angle around the z axis (in rads.)
%
% output:
% R: The rotation matrix
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
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.

Rx = [1 0 0 ; 0 cos(tetax) sin(tetax) ; 0 -sin(tetax) cos(tetax)];
Ry = [cos(tetay) 0 -sin(tetay) ; 0 1 0 ; sin(tetay) 0 cos(tetay)];
Rz = [cos(tetaz) sin(tetaz) 0; -sin(tetaz) cos(tetaz) 0 ; 0 0 1];

R = Rx*Ry*Rz;
