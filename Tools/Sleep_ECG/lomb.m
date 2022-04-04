function varargout=lomb(varargin)
%
% [Pxx,F]=lomb(x,dcOffset,smooth)
%
%    Wrapper to WFDB LOMB:
%         http://www.physionet.org/physiotools/wag/lomb-1.htm
%
% Transforms a real-valued time series 'x' into a power spectrum 'X', using a 
% technique known as the Lomb periodogram. The input is a Nx2 matrix containing 
% a sampled time series, presented as two columns of numbers (the sample times 
% in the first column and the sample values in the second). The intervals between 
% consecutive samples need not be uniform.  
%
%Input Parameters:
% x    
%       Nx2 vector of doubles. First column is sample time index (in
%       seconds), and second column is the sample value of the signal at
%       that time.
%
% dcOffset (Optional)
%       Booelan. If present add constant to input samples ( x(:,2) ), such that the mean
%       values of the time series is zero (default=1).
%
% smooth   (Optional)
%       Boolean String specifying the if the output should be smoothed (default =1).
% 
%
%Output Parameters:
%
%Pxx 
%       Mx1 Double. Estimated power spectrum.
%
%F 
%       Mx1 Double. Frequency of the estimated power spectrum (Hz).
%
%
% CITING CREDIT: To credit this function, please cite the following paper at your work:
%
%Moody, G.B.
%    Spectral analysis of heart rate without resampling. Computers in Cardiology 1993, pp. 715-718 (IEEE Computer Society Press, 1993). http://www.physionet.org/physiotools/lomb/lomb.html . 
%
%
%Additional References:
%Lomb, N.R.
%    Least-squares frequency analysis of unequally spaced data. Astrophysics and Space Science 39:447-462 (1976). 
%Press, W.H, and Rybicki, G.B.
%    Fast algorithm for spectral analysis of unevenly sampled data. Astrophysical J. 338:277-280 (1989). 
%Press, W.H. Teukolsky, S.A., Vetterling, W.T., and Flannery, B.P.
%    Numerical Recipes in C: the Art of Scientific Computing, pp. 575-584 (Cambridge Univ. Press, 1992). 
%Moody, G.B.
%    Spectral analysis of heart rate without resampling. Computers in Cardiology 1993, pp. 715-718 (IEEE Computer Society Press, 1993). http://www.physionet.org/physiotools/lomb/lomb.html . 
%
%
%
% MATLAB wrapper written by Ikaro Silva, 2013
% Last Modified: -
% Version 1.0
% Since 0.9.0 
%
%
% %Example: Heart Rate Spectral Analysis:
% [tm, signal]=rdsamp('mitdb/100',1);
% [ann]=rdann('mitdb/100','atr');
% [Pxx,F]=lomb([tm(ann) signal(ann)]);
% plot(F,Pxx);grid on;hold on
%
% See also RDANN, TACH, SQRS, WQRS

%endOfHelp

persistent javaWfdbExec
if(isempty(javaWfdbExec))
    javaWfdbExec=getWfdbClass('lomb');
end

%Set default pararamter values
%[Pxx,F]
inputs={'x','dcOffset','smooth'};
dcOffset=1;
smooth=1;
for n=1:nargin
    if(~isempty(varargin{n}))
        eval([inputs{n} '=varargin{n};'])
    end
end

wfdb_argument={'-P'};

if(dcOffset)
     wfdb_argument{end+1}='-z';
end
if(smooth)
     wfdb_argument{end+1}='-s';
end

wfdb_argument{end+1}='-';
del=repmat([' '],size(x(:,1)));
data=[num2str(x(:,1)) del num2str(x(:,2))];
javaWfdbExec.setArguments(wfdb_argument);
pxx=char(javaWfdbExec.execWithStandardInput(cellstr(data)));
pxx=sscanf(pxx(2:end-1), '%f %f,');

varargout{1}=pxx(2:2:end);
varargout{2}=pxx(1:2:end);

