function varargout = dfaScalingExponent(x, minBoxSize, maxBoxSize, midBoxSize, pflag)
%
% varargout = dfaScalingExponent(xminBoxSize, midBoxSize, maxBoxSize, pflag) 
% calculates the detrended fluctuation analysis estimate of the scaling 
% exponent alpha. 
%
% INPUTS
%         x          : A Nx1 vector containing the series to be analyzed
%         minBoxSize : Smallest box width (default: 4)
%         maxBoxSize : Largest box width (default: N/4)
%         midBoxSize : Medium time scale  (default: empty)
%         pflag      : (Optional) pflag=1 plot,  pflag=0 
% OUTPUTS     
%         alpha      : estimate of scaling exponent, +
%                      minBoxSize <= n <= maxBoxSize
%         alpha1     : estimate of scaling exponent for short-term 
%                      fluctuations, minBoxSize <= n < midBoxSize
%         alpha2     : estimate of scaling exponent for long-term 
%                      fluctuations, midBoxSize <= n <= maxBoxSize
%
% The raw time series x(i) is first integrated to give y(i); i=1,...,N. 
% For each length scale, n, y(i) is divided into segments of equal length, n.
% In each segment, the data is detrended by subtracting the local linear least 
% squares fit, yn(k).  The root-mean-square fluctuation of this integrated 
% and detrended time series is given by 
% F(n) = sqrt( (1/N) sum_{k=1}^N [y(k) - yn(k)]^2 )
% We calculate the average fluctuation F(n) for each segment n. 
% If the scaling approximately given by F(n) = c n^alpha, 
% we can estimate alpha by calculating the slope of log F(n) versus log n.
% Such a linear relationship on a log-log plot indicates the presence of 
% power law (fractal) scaling. 
% A log-log plot of F(n) against n is provided when pflag=1.  Default: plag=0.
% Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL. 
% Mosaic organization of DNA nucleotides. Phys Rev E 1994;49:1685-1689.
%
%
% 09-20-2017 Modified by Giulia Da Poian (GDP) to be included in the VOSIM 
%            HRV toolbox. (Original function name: dfa)
%
% Copyright (c) 2005 Patrick E. McSharry (patrick@mcsharry.net)
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

if nargin < 2 || isempty(maxBoxSize)
    minBoxSize = 4;
end
if nargin < 3 || isempty(maxBoxSize)
    maxBoxSize = length(x)/4;
end
if nargin < 4 
    midBoxSize = [];
end
if nargin < 5
   pflag = 0;
end

if size(x,1)<size(x,2)
    x=x';
end

N = length(x);     
y = cumsum(x);

n1 = round(log2(minBoxSize)); % modified GDP, was 3
n2 = round(log2(maxBoxSize)); % modified GDP, was n2 = round(log2(N/2))
ns = (2.^(n1:n2))';           % modified GDP, was ns =[2.^[n1:n2] N]' 

nn = length(ns);
F = zeros(nn,1);
for n=1:nn
   t = trend(y, ns(n));
   z = y - t;
   F(n) = sqrt(mean(z.^2));
 
end

lns = log10(ns);
lF = log10(F);
A = ones(nn,2);
A(:,2) = lns;
a = pinv(A)*lF;
alpha = a(2);  
lFpred = A*a;

% output
varargout{1} = alpha;

% Introduced by GDP to include the possibility of computing alpha1 and alpha2 
if ~isempty(midBoxSize)   
    nMid = round(log2(midBoxSize-minBoxSize));
    % alpha1 from minBoxSize to midBoxSize (excluded)
    lns1 = log10(ns(1:nMid-1));
    lF1 = log10(F(1:nMid-1));
    A1 = ones(length(ns(1:nMid-1)),2);
    A1(:,2) = lns1;
    a1 = pinv(A1)*lF1;
    alpha1 = a1(2); % alpha1  
    % alpha2 from midBoxSize (included) to maxBoxSize
    lns2 = log10(ns(nMid:end));
    lF2 = log10(F(nMid:end));
    A2 = ones(length(ns(nMid:end)),2);
    A2(:,2) = lns2;
    a2 = pinv(A2)*lF2;
    alpha2 = a2(2); % alpha2

    % output
    varargout{2} = alpha1;
    varargout{3} = alpha2;       
end
    

if pflag == 1
    figure;
    loglog(10.^lns, 10.^lF,'b.-','MarkerSize',16);
    hold;
    loglog(10.^[lns(1) lns(nn)], 10.^[lFpred(1) lFpred(nn)],'k');
    xlabel('n');
    ylabel('F(n)');
    title(['F(n) ~ n^{\alpha} with \alpha = ' num2str(a(2)) ]);
end

end % dfaScalingExponent function



function t = trend(y, n)
    N = length(y);
    t = zeros(N,1);
    r = floor(N/n);
    for i=1:r 
       v = y((i-1)*n+1:i*n);
       t((i-1)*n+1:i*n) = linfit(v); 
    end 
    v = y(r*n+1:N);
    t(r*n+1:N) = linfit(v);
end % trend function
   
function up = linfit(v)
    k = length(v);
    A = ones(k,2);
    u = [1:k]';
    A(:,2) = u;
    a = pinv(A)*v;
    up = A*a;
end % linfit function