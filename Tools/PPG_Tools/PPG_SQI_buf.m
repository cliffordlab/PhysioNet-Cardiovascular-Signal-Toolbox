% [annot sqimatrix template valid] = PPG_SQI_buf(wave,anntime,template_ahead,windowlen,Fs)
% 
% PPG_SQI.m - PPG SQI based on beat template correlation.
% (as an advice, the algorithm get 30 beats at each call and run in loop)
% by Qiao Li 30 Mar 2011
% 
% input: 
%     wave:       PPG data; 
%     anntime:    PPG annotation time (samples), read from ple annot file,
%                 But the ann-time is the OFFSET based on wave(1)
%     template:   Last PPG beat template 
%     windowlen:  length of window to calculate template(default: 30s)
%     Fs       :  sampling frequency (default Fs=125Hz)
% output:
%     annot:      ppg sqi annotation
%                     E - excellent beat; 
%                     A - acceptable beat; 
%                     Q - unacceptable beat
%     sqimatrix:  ppg sqi matrix   
%                     [N,1]: SQI based on Direct compare
%                     [N,2]: SQI based on Linear resampling
%                     [N,3]: SQI based on Dynamic time warping
%                     [N,4]: SQI based on Clipping detection
%     template:   Current PPG beat template
%     valid:      1 or greater for valid template, 
%                 0 for invalid template
%	
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%
% 03-03-2017
% Edits by Adriana Vest
% - Changed output variable annot from numeric to cell to preserve
%   characters
% - Style changes to align loops and conditional statements
%
% 12-01-2017 Modified by Giulia Da Poian: sampling frequency as input
% parameter instead of fixed fs = 125
% 
% 12-19-2017 Modified by Giulia Da Poian: replaced dp_dtw with dpfast 


function [annot sqimatrix template valid] = PPG_SQI_buf(wave,anntime,template,windowlen,Fs)

    if nargin < 5
        Fs = 125; 
    end
    if nargin < 4 || isempty(windowlen)
        windowlen=30*Fs;
    end
    if nargin < 3 || isempty(template)
        template=[];
    end
    if nargin < 2
        sprintf('Error: must provide wave, anntime');
        annot=[];
        sqimatrix=[];
        template=[];
        valid=0;
        return;
    end

    annot=[];
    sqimatrix=[];
    
    % get PPG template
    [t t2 v]=template_pleth(wave(1:min(windowlen,length(wave))), anntime(find(anntime<min(windowlen,length(wave)))),0, Fs);

    if v<1 && length(template)<1 % Current template invalid && no previous template available
        for j=1:length(anntime)
            annot{j}='Q'; % Unacceptable
        end
        t=[];
    else
        % Using previous template as the template
        if v<1
            t=template;
        end
        % Using t2 if available
        if v>1
            t=t2;
        end
        
        % Calculate the PLA of template for dynamic time warping
        d1=t;
        d1=(d1-min(d1))/(max(d1)-min(d1)).*100;
        [y1 pla1]=PLA(d1,1,1);

        % Main Loop
        for j=1:length(anntime)-1
            %% SQI1: Direct compare
            % Calculate correlation coefficients based on the template
            % length
            beatbegin=anntime(j);
            beatend=anntime(j+1);
            % 07/11/2011 ADD max beat length <= 3s detection 
            if beatend-beatbegin>3*Fs
                beatend=beatbegin+3*Fs;
            end
            templatelength=length(t);
            if beatbegin+templatelength-1 > length(wave) || beatend > length(wave) || beatbegin < 1
                continue;
            end
            currentb = j;
            cc = corrcoef(t,wave(beatbegin:beatbegin+templatelength-1));
            c1(j) = cc(1,2);
            if (c1(j)<0)
                c1(j)=0;
            end
            sqimatrix(currentb,1)=int8(c1(j)*100);

            %% SQI2: Linear resampling
            % Calculate correlation coefficients based on the 
            % linear resampling (interp1)
            
            y=interp1(1:beatend-beatbegin, wave(beatbegin:beatend-1),1:(beatend-beatbegin-1)/(templatelength-1):(beatend-beatbegin),'spline');
            y(isnan(y))=0;
            cc=corrcoef(t,y);
            c2(j)=cc(1,2);
            if (c2(j)<0)
                c2(j)=0;
            end
            sqimatrix(currentb,2)=int8(c2(j)*100);

            %% SQI3: Dynamic Time Warping                
            % Calculate correlation coefficients based on the dynamic time
            % warping
            d2=wave(beatbegin:beatend-1);
            
            % if beat too long, set SQI = 0;
            if (length(d2)>length(d1)*10)
                c3(j)=0;
            else
                d2=(d2-min(d2))/(max(d2)-min(d2)).*100;
               [y2 pla2]=PLA(d2,1,1);

               [w ta tb] = simmx_dtw(y1,pla1,y2,pla2);
               try % try to use the fast version if possible
                   [p,q,Dm] = dpfast(w);
               catch
                   [p,q,Dm] = dp_dtw(w); 
               end
               [ym1, ym2, yout1] = draw_dtw(y1,pla1,p,y2,pla2,q); 
                cc=corrcoef(y1,ym2);
                c3(j)=cc(1,2);
                if (c3(j)<0)
                    c3(j)=0;
                end
            end
            sqimatrix(currentb,3)=int8(c3(j)*100);
                
            %% SQI4: Clipping detection   
            d2=wave(beatbegin:beatend-1);
            y=diff(d2);
            clipthreshold=0.5;
            c4(j)=int8(length(find(abs(y)>clipthreshold))/length(y)*100);
            sqimatrix(currentb,4)=c4(j);

            %% SQI: Combined

            sqibuf=[sqimatrix(currentb,1) sqimatrix(currentb,2) sqimatrix(currentb,3) sqimatrix(currentb,4)];
            if min(sqibuf)>=90
                annot{currentb}='E'; % Excellent
            else
                if length(find(sqibuf>=90))>=3 || (median(sqibuf(1:3))>=80 && sqibuf(1)>=50 && sqibuf(4)>=70) || min(sqibuf)>=70
                    annot{currentb}='A'; % Acceptable
                else
                    annot{currentb}='Q'; % Unacceptable
                end
            end
                
        end
    end
    
    template=t;
valid=v;