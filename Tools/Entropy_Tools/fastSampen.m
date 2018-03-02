function s = fastSampen(y,m,r)
% INPUTS
%  y     num_var x num_samples
%  m     template length
%  r     radius
% OUTPUTS
%  s     Sample Entropy
%
% Written by Shamim Nemati <shamim.nemati@alum.mit.edu>
%
% 01-31-2018 : modified by Giulia Da Poian <giulia.dap@gmail.com>, see
%              comments in the code 
%              NOTE : Sample entropy quantifies the likelihood that a 
%              sequence of m consecutive data points that matches another
%              sequence of the same length (match within a tolerance of r) 
%              will still match the other sequence when their length is 
%              increased of one sample (sequences of length m + 1); 
%              References : 
%              1. Richman JS, Moorman JR. Physiological time-series 
%              analysis using approximate entropy and sample entropy. 
%              American Journal of Physiology-Heart and Circulatory 
%              Physiology. 2000 Jun 1;278(6):H2039-49.
%              2. Humeau-Heurtier A. The multiscale entropy 
%              algorithm and its variants: A review. Entropy. 2015 May 
%              12;17(5):3110-23.
% LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information

% Ensure y is a row vector
if size(y, 1) > size(y,2)
    y = y';
end

xx = convert_to_lagged_form(y, m)';   % Giulia: replaced m-1 with m to match SampEn definition
Dxx = pdist(xx,'chebychev');

yy = convert_to_lagged_form(y, m+1)'; % Giulia: replaced m with m+1 to match SampEn definition
Dyy = pdist(yy,'chebychev');

A = mean( Dxx < r ) ;
B = mean( Dyy < r );

s = -log(B/A);


function yy = convert_to_lagged_form(y, k)
% Create an observation vector yy(:,t) containing the last k values of y, newest first
% e.g., k=2, y = (a1 a2 a3)     yy  = a2 a3
%                (b1 b2 b3)           b2 b2
%                                     a1 a2
%                                     b1 b2
[s, T] = size(y);
bs = s*ones(1,k);
yy = zeros(k*s, T-k+1);
for i=1:k,     yy(block(i,bs), :) = y(:, k-i+1:end-i+1); end

function sub = block(blocks, block_sizes)
% BLOCK Return a vector of subscripts corresponding to the specified blocks.
% sub = block(blocks, block_sizes)
%
% e.g., block([2 5], [2 1 2 1 2]) = [3 7 8].
blocks = blocks(:)';
block_sizes = block_sizes(:)';
skip = [0 cumsum(block_sizes)];
start = skip(blocks)+1;
fin = start + block_sizes(blocks) - 1;
sub = [];
for j=1:length(blocks)
    sub = [sub start(j):fin(j)];
end
