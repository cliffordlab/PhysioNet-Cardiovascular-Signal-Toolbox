function [entropy, conf95] = sampenMaxim(data, m, r)
% SAMPEN Calculate SampEn value
%
% Description:
%   The function takes a column vector of data and calculates Sample
%   Entropy value of order m and similarity r.
%
% Arguments:
%   data - column vector with data
%   m - pattern length
%   r - similarity criteria (absolute value)
%
% Results:
%   entropy - Sample Entropy value
%   conf95 - 95% confidence interval
%
% Copyright (C) 2011-2013, Maxim Osipov
%
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%
%  - Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%  - Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%  - Neither the name of the University of Oxford nor the names of its
%    contributors may be used to endorse or promote products derived from this
%    software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
% OF THE POSSIBILITY OF SUCH DAMAGE.
%

    len = length(data);
    % create an array of all m-templates
    tpl_len = len-m+1;
    tpl_data = zeros(tpl_len, m);
    i_tpl = 1;      % index of first template
    n_tpl = 0;      % length of first template series
    for i_data = 1:m
        n_data = floor((len-i_data+1)/m)*m;
        i_tpl = i_tpl + n_tpl;
        n_tpl = floor(n_data/m);
        tpl_data(i_tpl:(i_tpl+n_tpl-1), :) =...
            reshape(data(i_data:(i_data+n_data-1)), m, n_tpl)';
    end
    % count template matches, excluding self-matches
    BB = zeros(1, tpl_len);
    tpl_data = repmat(tpl_data, 2, 1);
    % shifting the vector until floor((tpl_len/2)+1 to avoid
    % double count of the same match
    for i_tpl = 2:floor((tpl_len/2)+1)
        tpl_diff = abs(tpl_data(1:tpl_len, :) -...
            tpl_data(i_tpl:i_tpl+tpl_len-1, :));
        BB(i_tpl-1) = sum(sum((tpl_diff < r), 2) == m);
    end
    B = sum(BB)/(len-m);
    % create an array of all m+1 templates
    mp = m+1;
    tpl_len = len-mp+1;
    tpl_data = zeros(tpl_len, mp);
    i_tpl = 1;      % index of first template
    n_tpl = 0;      % length of first template series
    for i_data = 1:mp,
        n_data = floor((len-i_data+1)/mp)*mp;
        i_tpl = i_tpl + n_tpl;
        n_tpl = floor(n_data/mp);
        tpl_data(i_tpl:(i_tpl+n_tpl-1), :) =...
            reshape(data(i_data:(i_data+n_data-1)), mp, n_tpl)';
    end
    % count template matches, excluding self-matches
    AA = zeros(1, tpl_len);
    tpl_data = repmat(tpl_data, 2, 1);
    % shifting the vector until floor((tpl_len/2)+1) to avoid
    % double count of the same match
    for i_tpl = 2:floor((tpl_len/2)+1),
        tpl_diff = abs(tpl_data(1:tpl_len, :) -...
            tpl_data(i_tpl:i_tpl+tpl_len-1, :));
        AA(i_tpl-1) = sum(sum((tpl_diff < r), 2) == mp);
    end
    A = sum(AA)/(len-mp);
    
    % calculate value of sample entropy
    if (A<=0) || (B<=0),
        entropy = -log(1/(m*mp));
    else
        entropy = -log((A/(len-mp))/(B/(len-m)));
    end

    % calculate 95% confidence interval
    EE = AA(1:len-mp)./BB(1:len-mp);
    conf95 = std(EE)*tinv(0.95,length(EE))/sqrt(length(EE));
end