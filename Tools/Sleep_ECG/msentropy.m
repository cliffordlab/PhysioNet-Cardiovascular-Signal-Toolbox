function varargout=msentropy(varargin)
%
% [y,scale,info]=msentropy(x,dn,dm,dr,N,N0,minM,maxM,maxScale,minR,maxR)
%
%    Wrapper to the Multiscale Entropy C code written by Madalena Costa (mcosta@fas.harvard.edu):
%         http://physionet.org/physiotools/mse/mse-1.htm
%
% Calculates the multi scale entropy of a signal 'x'. A tutorial on Mulsticale
% entropy is available at:
% http://www.physionet.org/physiotools/mse/tutorial/
%
%
% Please cite these publications when referencing this material:
%     Costa M., Goldberger A.L., Peng C.-K. Multiscale entropy analysis of biological signals. Phys Rev E 2005;71:021906.
%     Costa M., Goldberger A.L., Peng C.-K. Multiscale entropy analysis of physiologic time series. Phys Rev Lett 2002; 89:062102.
%
% Also include the standard citation for PhysioNet:
%     Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
%     Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank, PhysioToolkit, and PhysioNet: components of a new research resource for complex physiologic signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages; http://circ.ahajournals.org/cgi/content/full/101/23/e215]; 2000 (June 13)
%
% Readers of may also wish to read:
%     Costa M, Peng C-K, Goldberger AL, Hausdorff JM. Multiscale entropy analysis of human gait dynamics. Physica A 2003;330:53-60.
%
% Required Parameters:
%
% x
%       Nx1 vector of doubles in which to caculate the multiscale entropy.
%
% Optional Parameters are:
% dn
%       1x1 double. Sets the scale increment to dn (1-40; default: 1).
% dm
%       1x1 double. Sets the m increment to dm (1-10; default: 1).
% dr
%       1x1 double. Sets the scale increment to dr (>0; default: 0.05).
% N
%       1x1 integer. Stop the analysis with row N.
%       By default, analysis ends at row 39999, or at the end of the data set if there are fewer rows.
% N0
%       1x1 integer. Begin the analysis with row N0.
%       By default, analysis begins with row 1.
% minM
%       1x1 integer betwee 1-10. Set the minimum pattern length for SampEn to minN (1-10; default: 2).
% maxM
%        1x1 integer betwee 1-10. Set the maximum m to maxM (1-10; default: 2).
% maxScale
%        1x1 integer betwee 1-40. Set the maximum scale for coarse-graining to maxScale (1-40; default: 20).
% minR
%        1x1 double >0. Set the minimum similarity criterion for SampEn to minR (>0; default: 0.15).
% maxR
%        1x1 double > 0. Set the maximum m to maxR (>0; default: 0.15).
%
%
% Outputs:
% y
%       A 1xM vector of doubles corresponding to estimated sample entropies at each scale.
% scale
%       A 1xM vector of integers specifying the scales in which 'y' was
%       estimated.
%
% info
%       An optional 3x1 cell array of strings providing loggin and verbose information from
%       the calculation.
%
% Wrapper written by Ikaro Silva, 2013
% Last Modified: March 20, 2014
% Version 0.0.1
%
% Since 0.9.5
%
% %Example
% N=30000;
% noise=randn(N,1);
% maxScale=10;
%[entropyNoise,scale1]=msentropy(noise,[],[],[],[],[],[],[],maxScale);
% %Simulate determistic system with noise-like 2nd order statistics
% nlinear=zeros(N,1);nlinear(1)=0.2;u=4;
% for n=2:N;nlinear(n)=u*nlinear(n-1)*(1-nlinear(n-1));end
%[entropyDeterm,scale2]=msentropy(nlinear,[],[],[],[],[],[],[],maxScale);
%subplot(2,1,1);
%plot(noise(1:1000));hold on;grid on;plot(nlinear(1:1000),'r');legend('Stochastic','Deterministic')
%subplot(2,1,2);
%plot(scale1,entropyNoise);hold on;grid on;plot(scale2,entropyDeterm,'r');legend('Stochastic','Deterministic')
%
%
% See also SURROGATE, DFA, WFDBDESC, PHYSIONETDB, RDANN, ANN2RR, MAPRECORD

%endOfHelp
persistent javaWfdbExec config
if(isempty(javaWfdbExec))
    [javaWfdbExec,config]=getWfdbClass('mse');
end

%Set default pararamter values
inputs={'x','dn','dm','dr','N','N0','minM','maxM','maxScale','minR','maxR'};
outputs={'y','scale','info'};
dn=[];
dm=[];
dr=[];
N=[];
N0=[];
minM=[];
maxM=[];
maxScale=[];
minR=[];
maxR=[];
wfdb_argument={};
info=[];
scale=[];
y=[];
x=[];
for n=1:nargin
    if(~isempty(varargin{n}))
        eval([inputs{n} '=varargin{n};'])
    end
end
if(~isempty(dn))
    wfdb_argument{end+1}='-a';
    wfdb_argument{end+1}=[num2str(dn)];
end
if(~isempty(dm))
    wfdb_argument{end+1}='-b';
    wfdb_argument{end+1}=[num2str(dm)];
end
if(~isempty(dr))
    wfdb_argument{end+1}='-c';
    wfdb_argument{end+1}=[num2str(dr)];
end
if(~isempty(N0))
    wfdb_argument{end+1}='-i';
    wfdb_argument{end+1}=[num2str(N0-1)];
end
if(~isempty(N))
    wfdb_argument{end+1}='-I';
    wfdb_argument{end+1}=[num2str(N-1)];
end
if(~isempty(minM))
    wfdb_argument{end+1}='-m';
    wfdb_argument{end+1}=[num2str(minM)];
end
if(~isempty(maxM))
    wfdb_argument{end+1}='-M';
    wfdb_argument{end+1}=[num2str(maxM)];
end
if(~isempty(maxScale))
    wfdb_argument{end+1}='-n';
    wfdb_argument{end+1}=[num2str(maxScale)];
end
if(~isempty(minR))
    wfdb_argument{end+1}='-r';
    wfdb_argument{end+1}=[num2str(minR)];
end
if(~isempty(maxR))
    wfdb_argument{end+1}='-R';
    wfdb_argument{end+1}=[num2str(maxR)];
end
javaWfdbExec.setArguments(wfdb_argument);

if(config.inOctave)
    x=cellstr(num2str(x));
    x=java2mat(javaWfdbExec.execWithStandardInput(x));
    Nx=x.size;
    out=cell(Nx,1);
    for n=1:Nx
        out{n}=x.get(n-1);
    end
else
    out=cell(javaWfdbExec.execWithStandardInput(x).toArray);
end
M=length(out);
if(M<4)
    error(['Error calculating MSE:' out{:}])
end
info=out(1:3);
out(1:4)=[];
M=M-4;
scale=zeros(M,1)+NaN;
y=zeros(M,1)+NaN;
for m=1:M
    str=out{m};
    sep=regexp(str,'\s');
    scale(m)=str2num(str(1:sep));
    y(m)=str2num(str(sep(1):sep(2)));
end

for n=1:nargout
    eval(['varargout{n}=' outputs{n} ';'])
end

