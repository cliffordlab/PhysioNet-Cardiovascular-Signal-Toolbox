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
