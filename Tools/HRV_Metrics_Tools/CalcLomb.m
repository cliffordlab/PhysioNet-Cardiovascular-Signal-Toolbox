function [PSD,F] = CalcLomb(tNN,NN,F,nfft,norm)
%
% [PSD,F] = CalcLomb(tNN,NN,F,nfft,norm);
% TODO: Change Function Name
%   OVERVIEW:   Calculates the Lomb-Scargle normalized periodogram values 
%               "Px" as a function of the supplied vector of frequencies 
%               "f" for input vectors "tNN" (time) and "NN" (observations).
%               Also returns the probability "Prob" that the null hypothesis
%               is valid (same length as Px and freq). Time stamps, tNN and 
%               amplitudes "NN" must be the same length. 
%
%               Equivalent to the following native Matlab function:
%               [PSD,f] = plomb(y,t,F,'normalized');
%
%   INPUT:      tNN  : the time indices of the rr interval data (seconds)
%               NN   : a single row of NN (normal normal) interval data in seconds
%               F    : frequencies at which PSD should be estimated
%                      default will be determined based on the maximal
%                      frequency detected in RR intervals (1/minimal interval)
%               nfft : (optional) Only needed for normalization
%               norm : normalize output; 1 for normalize, 0 for not    
%
%   OUTPUT:     PSD  :
%               F    :
   
%   REFERENCE:  See Scargle J.D.:"Studies in astronomical time series analysis. II. 
%               Statistical aspects of spectral analysis of unevenly spaced data,"
%               Astrophysical Journal, vol 263, pp. 835-853, 1982.  ... and
%               Lomb N.R: "Least-squares frequency analysis of unequally spaced data",
%               Astrophysical and Spcae Science, vol 39, pp. 447-462, 1976.
%
%               Sources 43, 162, 198 in Clifford Thesis
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL
%   SOURCE:     Gari Clifford HRV
%               http://www.robots.ox.ac.uk/~gari/CODE/HRV/
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%
 
% Lomb Scargle Method
% maxfreq = 2/(min(diff(t)));
% nfft = 1024;
% F = [1/nfft:1/nfft:maxfreq]; % setting up frequency vector

var = std(NN);          % calc var for normalization
PSD  = zeros(size(F));  

for i=1:length(F)
	%pcnt = 100*i/length(F);
	%fprintf('Percent Complete: %2.2f\n', pcnt);
	w = 2*pi*F(i);
	if w > 0 
        twt = 2*w*tNN;
        tau = atan2(sum(sin(twt)),sum(cos(twt)))/2/w;
        wtmt = w*(tNN - tau);
        PSD(i) = (sum(NN.*cos(wtmt)).^2)/sum(cos(wtmt).^2) + ...
        (sum(NN.*sin(wtmt)).^2)/sum(sin(wtmt).^2);
    else
        PSD(i) = (sum(NN.*tNN).^2)/sum(tNN.^2);
    end
    
end

if norm
    for i=1:length(PSD)
        if var~=0            % check for divide by zero
            PSD(i)=PSD(i)/2/var.^2;
            Prob(i) = 1-(1-exp(-PSD(i)))^nfft;
        else
            PSD(i)=inf; 
            Prob(i)=1;
        end
        if Prob(i) < .001  % allow for possible roundoff error
            Prob(i) = nfft*exp(-PSD(i));
        end
    end
end


end