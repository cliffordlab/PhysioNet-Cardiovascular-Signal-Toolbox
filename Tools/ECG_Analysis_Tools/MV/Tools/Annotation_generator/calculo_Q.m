function S = calculo_Q(ecg,R,fsecg)
% Q wave detection based on R detection (minimum value of the ECG in a time window, of length 0.1 seconds, ending at the peak of the R-wave).
%
%   In:   ecg signal
%         R peaks (in samples)
%         fsecg ecg sampling frequency
%
%   Out:  Q peaks (in samples)
% 
%
%   Created by Eduardo Gil (edugilh@unizar.es) on Sept 3, 2008.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

%% Input checks
if nargin < 2
    error('Provide ECG signal & R peaks');
end

if nargin < 3
    fsecg = 500;
end

w_size=round(0.1*fsecg);

for i=1:size(R,2)
    [value,index]=min(ecg(R(i)-w_size:R(i)));
    S(i)=R(i)-w_size+index-1;
end
