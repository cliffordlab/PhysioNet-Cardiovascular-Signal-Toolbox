function [fe_output, fe_ci] = fuzzyEntropy(data, m, r)
 
%   [fe_output, fe_ci] = fuzzyEntropy(data, m, r)
%
%	AUTHORS: Chengyu Liu <bestlcy@sdu.edu.cn>
%
%	OVERVIEW:
%       Calculates the fuzzy measure entropy for a physiological signal time sequence. 
%       Using fuzzy entropy aims to reduce the poor statistical stability of traditional entropy (ApEn and SampEn). 
%       The algorithm can refer to [Ref1]
%       The recommendatory parameter settings of m=2 and r=0.1 can refer to [Ref2]
%
%	INPUTS & VARIABLES:
%       data - time series
%       m - template length
%       r - similarity criteria
%   	x is signal sequence; m is embedding dimension (usually m=2)
%       r is threshold value (usually r = 0.1-0.25), in here:
%       local threshold r_l = 0.2, global threshold r_g = 0.2, 
%       local weight of sequence segments' similarity n_l = 3
%       global weight of sequence segments' similarity n_g = 2
%
%	OUTPUTS:
%		fe_output - fuzzy entropy output; single value > 0, usually ~ 0 to 2.5
%       fe_ci - fuzzy entropy 95% confidence interval
%
%	REFERENCES:
%		[Ref1]. Liu, C Y. et al. Analysis of heart rate variability using fuzzy measure entropy. Comput. Biol. Med. 43, 100???108 (2013).
%               [Ref2]. Zhao???L N. et al. Determination of sample entropy and fuzzy measure entropy parameters for distinguishing congestive heart failure from normal sinus rhythm subjects. Entropy, 2015, 17(9): 6270-6288.
%	
%	LAST REVISED
%		by Erik Reinertsen <ereiner@emory.edu>
%		on 2 Apr 2016
%		for complete revision history of this code, refer to https://github.com/cliffordlab/icelandcode
%	
%    COPYRIGHT (C) 2016 AUTHORS (see above)
%		This program is free software: you can redistribute it and/or modify
%		it under the terms of the GNU General Public License as published by
%		the Free Software Foundation, either version 3 of the License, or
%		(at your option) any later version.
%
%		This program is distributed in the hope that it will be useful,
%		but WITHOUT ANY WARRANTY; without even the implied warranty of
%		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%		GNU General Public License for more details.
%
%		You should have received a copy of the GNU General Public License
%		along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Check on number of input arguments and set default values
if nargin < 3
	r = 0.1;
end
if nargin < 2
	m = 2;
end

% Convert data from double to single to improve memory efficiency
data = single(data);

% If row vector, transpose to column vector for readability
if size(data,1) < size(data,2)
    data = data';
end

N = length(data);
r_l = r; % r for local similarity
r_g = r; % r for global similarity
n_l = 3; % weight for local similarity
n_g = 2; % weight for global similarity

% Compute indices
indb = hankel(1:N-m, N-m:N-1); % m --> vector b
inda = hankel(1:N-m, N-m:N);   % m+1 --> vector a

% Global vectors
yb_g   = data(indb); % b global
if m  == 1
    yb_g = yb_g(:);
end
ya_g = data(inda); % a global

% local vectors
yb_mean = mean(data(indb)');
yb_mean = repmat(yb_mean',1,size(yb_g,2));
yb_l = yb_g - yb_mean; % local b
ya_mean = mean(data(inda)');
ya_mean = repmat(ya_mean',1,size(ya_g,2));
ya_l = ya_g - ya_mean; % local a

clearvars data yb_mean ya_mean

% using pdist for saving time but need a large memory

% for global b
db_g = pdist(yb_g, 'chebychev'); % maximum coordinate difference
Db_g = sum(exp(-(db_g.^n_g/r_g)));
B_g = Db_g*2/(size(yb_g,1)*(size(yb_g,1)-1));

% for global a
da_g = pdist(ya_g, 'chebychev');
Da_g = sum(exp(-(da_g.^n_g/r_g)));
A_g = Da_g*2/(size(ya_g,1)*(size(ya_g,1)-1));

% for local b
db_l = pdist(yb_l, 'chebychev'); % maximum coordinate difference
Db_l = sum(exp(-(db_l.^n_l/r_l)));
B_l = Db_l*2/(size(yb_l,1)*(size(yb_l,1)-1));

% for local a
da_l = pdist(ya_l, 'chebychev');
Da_l = sum(exp(-(da_l.^n_l/r_l)));
A_l = Da_l*2/(size(ya_l,1)*(size(ya_l,1)-1));

% calculate local and global fuzzy measure entropy
fuzzylmen = -log(A_l/B_l);
fuzzygmen = -log(A_g/B_g);

fe_output = fuzzylmen + fuzzygmen;

fe_ci = NaN;

end