% //This software is licensed under the BSD 3 Clause license: http://opensource.org/licenses/BSD-3-Clause 
% 
% 
% //Copyright (c) 2013, University of Oxford
% //All rights reserved.
% 
% //Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 
% //Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% //Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% //Neither the name of the University of Oxford nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% //THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%   The method implemented in this file has been patented by their original
%   authors. Commercial use of this code is thus strongly not
%   recocomended.
%
% //Authors: 	Gari D Clifford - 
% //            Roberta Colloca -
% //			Julien Oster	-

function [ MAD ] = comput_MAD(data)

%comput_MAD computes the MAD value for the vector  'data'

%Input
%data RR interval series (column vector)

%Output
%MAD median of median value for heart rate deviation


len=length(data);
J=3; %number of segments
N=len/J;  % number of RR intervals per segment

%heart rates segments: J-1=seg(;,1) J=seg(;,2) J+1=seg(;,3)
data=reshape(data,N,J);
seg=1./data ;

%remove mean from each segment: J-1, J, J+1
seg_m=detrend(seg,'constant');

%remove linear trend
seg_m_lt=detrend(seg_m,'linear');

%determine absolute deviation from mean, absdev_i for i=1..N
%within segments J-1, J, J+1
abs_dev=abs(detrend(seg_m_lt,'constant'));

%determine the Median of absolute deviation from mean for each RR interval, mABSdev_j within segments J-1, J, J+1
 mABSDEV= median(abs_dev);
 
 %determine median absolute deviation of segment J, compared to adjacent
 %segments , J-1 and J+1
  MAD= median(mABSDEV);
end

