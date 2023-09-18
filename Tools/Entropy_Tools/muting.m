function xcor=muting(data1, data2)

% xcor=muting(data1, data2)
% computing mutual information between
% columns of the data matrix "data"

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


data = [data1' data2']; 
[N  dim]=size(data);
[xdat  idat]=sort(data);
clear xdat
ydat=idat;
for d=1:dim
      ydat(idat(:,d),d)=(1:N)';%'
end
ddim=2^dim; dim2=2*dim;
clear idat
xcor=0; npar=1; poc(1)=1; kon(1)=N; poradi=1:N;
NN=zeros(1,ddim); marg=zeros(8*ddim,dim2);
marg(1,:)=[ones(1,dim)   N*ones(1,dim)];
Imm=[0; 1];
for d=2:dim
      Imm=[zeros(length(Imm),1)  Imm; ones(length(Imm),1)  Imm];
end
chi2=[0   7.81   13.9   25.0   42.0];
run=0;
while npar>0
    run=run+1;
    apoc=poc(npar); akon=kon(npar);
    apor=poradi(apoc:akon); Nex=length(apor);
    ave=floor((marg(npar,1:dim)+marg(npar,dim+1:dim2))/2);
    J=(ydat(apor,:)<=ones(Nex,1)*ave);
    I=zeros(Nex,ddim);
    amarg=ones(ddim,1)*marg(npar,:);
    for d=1:ddim
          I(:,d)=ones(Nex,1);
          for k=1:dim
                if Imm(d,k)
                    I(:,d)=I(:,d)&~J(:,k);
                    amarg(d,k)=ave(k)+1;
                else
                    I(:,d)=I(:,d)&J(:,k);
                    amarg(d,k+dim)=ave(k);
                end
          end
    end
    NN=sum(I);
    tst=ddim*sum((NN-Nex/ddim*ones(1,ddim)).^2)/Nex;
    if (tst>chi2(dim))|(run==1)
          npar=npar-1;
          for ind=1:ddim
                if NN(ind) > ddim
                    npar=npar+1;
                    akon=apoc+NN(ind)-1;
                    poc(npar)=apoc; kon(npar)=akon;
                    marg(npar,:)=amarg(ind,:);
                    poradi(apoc:akon)=apor(find(I(:,ind)));
                    apoc=akon+1;
               else
                    if NN(ind) > 0

Nxx=prod(amarg(ind,dim+1:dim2)-amarg(ind,1:dim)+ones(1,dim));
                        xcor=xcor+NN(ind)*log(NN(ind)/Nxx);
                    end
               end
         end
     else
         Nxx=prod(marg(npar,dim+1:dim2)-marg(npar,1:dim)+ones(1,dim));
         xcor=xcor+Nex*log(Nex/Nxx);
         npar=npar-1;
     end
end
xcor=xcor/N+(dim-1)*log(N);
