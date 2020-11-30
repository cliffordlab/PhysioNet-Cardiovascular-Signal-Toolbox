function E = ECGModelError(X,ECGmn,Phasemn,flag)
%
% Synthetic ECG model error
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

L = (length(X)/3);

alphai = X(1:L);
bi = X(L+1:2*L);
tetai = X(2*L+1:3*L);

Z = zeros(size(ECGmn));
for j = 1:length(alphai),
    dtetai = rem(Phasemn - tetai(j) + pi,2*pi)-pi;
    Z = Z + alphai(j) .* exp(-dtetai .^2 ./ (2*bi(j) .^ 2));
end

if(flag==0)
    E = (Z-ECGmn);
elseif(flag==1)
    E = Z;
end