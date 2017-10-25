function s = fastSampen(y,m,r)
% INPUTS
%  y     num_var x num_samples
%  m     max template length
%  r     radius
% OUTPUTS
%  s     Sample Entropy
%
% Written by Shamim Nemati <shamim.nemati@alum.mit.edu>

% Ensure y is a row vector
if size(y, 1) > size(y,2)
    y = y';
end

xx = convert_to_lagged_form(y, m-1)';
Dxx = pdist(xx,'chebychev');

yy = convert_to_lagged_form(y, m)';
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
