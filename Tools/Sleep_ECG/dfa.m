function [ln,lf]=dfa(varargin)
%
% [ln,lf]=dfa(x,p,integrateFlag,minBoxSize,maxBoxSize,slideWindowFlag)
%
%
% Wrapper to the DFA Algorithm in:
%    http://www.physionet.org/physiotools/dfa/
%
% References: 
% Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL. Mosaic organization of DNA nucleotides. Phys Rev E 1994;49:1685-1689.
% Peng C-K, Havlin S, Stanley HE, Goldberger AL. Quantification of scaling exponents and crossover phenomena in nonstationary heartbeat time series. Chaos 1995;5:82-87.
%
% Please cite at least one of the above publications when referencing this
% material.
%
% Required Input Options are:
%
% x
%       A Nx1 vector of doubles. 
%
% Optional Input Options are:
%
% p
%       Detrend using a polynomial of degree p (default: p=1, linear
%       detrending).
%
% integrateFlag
%       Input series is already integrated ( default= false ).
%
% minBoxSize
%       Smallest box width (default: 2p+2)
%
% maxBoxSize
%       Largest box width (default: N/4)
%
% slideWindowFlag
%       Sliding window DFA (default =false);
%
%
% The Output variables are:
%
% ln
%       A (MaxBoxSize -MinBoxSize) x 1 vector of log(boxsize)
%
% lf
%       A (MaxBoxSize -MinBoxSize) x 1 vector of the log of the root
%       mean square fluctuation for the given boxsize. 
%
%
% Written by Ikaro Silva, 2014
% Last Modified: November 21, 2014
% Version 1.0
%
% Since 0.9.8
%
% %Example:
%
%  gqrs('mitdb/117');
%  [rr]=ann2rr('mitdb/117','qrs');
%  [ln,lf]=dfa(rr);
%  plot(ln,lf)

%
%
% See also MSENTROPY, SURROGATE

%endOfHelp


persistent javaWfdbExec config
if(isempty(javaWfdbExec))
    [javaWfdbExec,config]=getWfdbClass('dfa');
end

%Set default pararamter values
inputs={'x','p','integrateFlag','minBoxSize','maxBoxSize','slideWindowFlag'};
p=[];
integrateFlag=false;
minBoxSize=[];
maxBoxSize=[];
slideWindowFlag=false;
wfdb_argument={};
for n=1:nargin
    if(~isempty(varargin{n}))
        eval([inputs{n} '=varargin{n};'])
    end
end
if(~isempty(p))
    wfdb_argument{end+1}='-d';
    wfdb_argument{end+1}=[num2str(p)];
end
if(integrateFlag)
    wfdb_argument{end+1}='-i';
end
if(~isempty(minBoxSize))
    wfdb_argument{end+1}='-l';
    wfdb_argument{end+1}=[num2str(minBoxSize)];
end
if(~isempty(maxBoxSize))
    wfdb_argument{end+1}='-u';
    wfdb_argument{end+1}=[num2str(maxBoxSize)];
end
if(slideWindowFlag)
    wfdb_argument{end+1}='-s';
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
ln=zeros(M,1)+NaN;
lf=zeros(M,1)+NaN;
for m=1:M
    str=out{m};
    sep=regexp(str,'\s');
    ln(m)=str2num(str(1:sep));
    lf(m)=str2num(str(sep(1):end));
end

