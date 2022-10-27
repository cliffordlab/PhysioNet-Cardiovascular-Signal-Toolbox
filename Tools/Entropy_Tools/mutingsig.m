function [I, Is, p] = mutingsig(x,y,Ns,g_flag)
% [I, Is, p] = mutingsig(x,y,Ns,g_flag);
% Calculates the mutual information I between vectors x and y which 
% take on continuous values, but is quantised by partitioning data to 
% produce an equal number of data poits in each bin.
% Ns permutated surrogates are used to calculate values Is and the 
% significance probability level p that Is < I.
% Use Ns = 19 for 95% and Ns = 99 for 99% (one-sided hypothesis test)
% graphics off (g_flag=0) by default, else plot on figure(ceil(g_flag))

% BSD 3-Clause License
%
% Copyright (c) 2001, P. McSharry
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



if nargin < 4
   g_flag=0;
end
if(g_flag>0) % ensure it is an integer for plotting
     g_flag=ceil(g_flag);
end;
if nargin < 3
Ns = 99;
end
if nargin < 2
   y=x;% calculate auto mutingsig
   fprintf('calculating mutingsig between vector and itself\n');
end



N = length(x);

I = muting(x,y);
 
for i=1:Ns
   sx = x(randperm(N));
   sy = y(randperm(N));
   Is(i) = muting(sx,sy);
end

% find number of Is < I
Nf = length(find(Is < I));
I;
p = Nf/(Ns+1);

if(g_flag>0)
figure(g_flag);
boundRegion=0.01;
subplot(3,1,1)
plot(x,'b')
axis([1 length(x) min(x)-(boundRegion*abs(min(x))) max(x)+(boundRegion*abs(max(x)))])
subplot(3,1,2)
plot(y,'c')
axis([1 length(y) min(y)-(boundRegion*abs(min(y))) max(y)+(boundRegion*abs(max(y)))])
subplot(3,2,5)
plot(x,y,'*g')
axis([min(x)-(boundRegion*abs(min(x))) max(x)+(boundRegion*abs(max(x))) min(y)-(boundRegion*abs(min(y))) max(y)+(boundRegion*abs(max(y)))])
hold on
%lsline
plot(x,y,'*m')
subplot(3,2,6)
plot(Is,'*r')
hold on
plot([1,Ns],[I I],'k')
axis([1 Ns 0 I+0.2])
end
