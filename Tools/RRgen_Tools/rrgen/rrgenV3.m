function data = rrgen(tmax, Pe, Pn, seed, filename);
% data = rrgen(tmax, Pe, Pn, seed, 'filename');
% generates artificial rr time series over a time tmax seconds
% and dumps the output to 'filename'. If filename is not specified,
% the output is returned into the data array. If it is specified, then the
% data array is zero-valued and of length(tmax). Pe and Pn are the 
% probability of ectopy and noise.
% Defaults : tmax=300; NoEctopic=1; NoNoise=1; seed=1; 
% P(ectopy), Pe = 0.0003;  ...  Probability of ectopy ~ 1 per hr 
% P(noise),  Pn = 0.0048;  ...  Probability of noise ~ 16 per hr 
% P(noise) in sleep is Pn/12, P(ectopy) in sleep is unchanged. 
%
% Requires rrgen binary - compilation of rrgen.c on your system:
% gcc -Wall rrgen.c -lm -o rrgen
%
% (C) P.E. McSharry & G.D. Clifford 2002 under the 
% GNU Public License. Contact P. McSharry (patrick@mcsharry.net) 
% or G. Clifford (gari@mit.edu)
% Also see - McSharry P.E., Clifford G.D., Tarassenko L., 
% Smith L.: Method for generating an artificial RR tachogram of a typical 
% healthy human over 24-hours, Computers in Cardiology, 29:225-228, IEEE 
% Computer Society Press, September 2002. 
 
if nargin < 5
  filename = 'tmp.dat';
end
if nargin < 4
  seed=1;    
end
if nargin < 3
    Pn=0.0048;
end
if nargin < 2
    Pe=0.0003;
end
if nargin < 1
  tmax=300;
end

if tmax<115
    fprintf('Warning, tmax minimum value is 115, using this as default\n')    
end

tmax=round(tmax);
eval(['!rrgen ' int2str(seed) ' ' int2str(tmax) ' ' int2str(Pe) ' ' int2str(Pn) ' > ' filename]);
fid = fopen(filename,'r');
data = zeros(1,tmax);
for(i=1:tmax)
 fscanf(fid,'%f\n',data(i)); 
end
fclose(fid);

% if not file name read in from tmp.dat and then delete it.
if nargin < 5
 load tmp.dat;
 data = tmp;
 if(ispc==1)
  system('del tmp.dat');
 else
  system('rm tmp.dat');
 end
end
