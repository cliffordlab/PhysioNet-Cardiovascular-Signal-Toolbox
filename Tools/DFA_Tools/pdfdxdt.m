function P = pdfdxdt(x, tau)
% P = pdfdxdt(x, tau) calculates the probability density function (PDF)
% of changes in x with respect to changes in time.  The PDF is estimated at 
% time intervals of tau. The technique uses a Gaussian kernel 
% with bandwidth  given by 1.06*std(x)^(1/5) [Silverman, 1986].
% Copyright (c) 2007 Patrick E. McSharry (patrick@mcsharry.net)

n = length(x);
if nargin < 2
   tau = [1 2 3 4];
end
if nargin < 3
   k = 100;
end

for i=1:length(tau)
   taustr{i} = ['\tau=' int2str(tau(i))]; 
end

% Select k evenly spaced points
taumax = max(tau);
d = x(1+taumax:n) - x(1:n-taumax);
dmax = max(abs(d));
dc = linspace(-dmax, dmax, k);

for i=1:length(tau)   
   d = x(1+tau(i):n) - x(1:n-tau(i));

   % calculate the probability of r at k points in rc
   p = pdfkernel(d, dc);
   
   % set zero values to eps in order to take the log 
   izero = find(p<=eps);
   p(izero) = eps;
   P(:,i) = p;
end

plot(dc,P,'LineWidth',2)
xlabel('\Deltax','FontSize',12)
ylabel('p(\Deltax)','FontSize',12);
set(gca,'FontSize',12);
legend(taustr);

% function for kernel estimation of the PDF
function p = pdfkernel(x, xc)
n = length(x);
k = length(xc);
p = zeros(k,1);
sigma = 1.06*std(x)*(n.^(-1/5));
sigma2 = sigma*sigma;
norm = sqrt(2*pi*sigma2);
d = x*ones(1,k) - ones(n,1)*xc;
d2 = d.*d;
g = exp(-0.5*d2/sigma2)/norm;
p = mean(g);
