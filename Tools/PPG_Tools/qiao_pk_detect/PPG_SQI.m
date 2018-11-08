% [template t2 valid] = PPG_SQI(wave,anntime,annot,template_ahead,windowlen,Fs)
% 
% PPG_SQI.m - PPG SQI based on beat template correlation.
% (as an advice, the algorithm get 30 beats at each call and run in loop)
% by Qiao Li 30 Mar 2011
% 
% PPG sampling frequency is Fs
%
% input: 
%     wave:       PPG data; 
%     anntime:    PPG annotation time (samples), read from ple annot file,
%                 But the ann-time is the OFFSET based on wave(1)
%     annot:      Annotation of beats, read from from ple annot file
%                 directly
%     template:   Last PPG beat template 
%     windowlen:  length of window to calculate template(default: 30s)
%     Fs       :  sampling frequency (defatult: 125 to work with pervious code)
% output:
%     annot:      ppg sqi annotation
%                     annot.typeMnemonic: E - excellent beat; 
%                                         A - acceptable beat; 
%                                         Q - unacceptable beat
%                     annot.subtype: SQI based on Direct compare
%                     annot.chan:    SQI based on Linear resampling
%                     annot.num:     SQI based on Dynamic time warping
%                     annot.auxInfo: SQI based on Clipping detection
%     template:   Current PPG beat template
%     valid:      1 or greater for valid template, 
%                 0 for invalid template
%	
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%
% 12-01-2107 Modified by Giulia Da Poian: replace fixed samplig frequency
% (125) with a variable Fs

function [annot template valid] = PPG_SQI(wave,anntime,annot,template,windowlen,Fs)

    if nargin < 6
        Fs =125;
    end
    if nargin < 5
        windowlen=30*Fs;
    end
    if nargin < 4
        template=[];
    end
    if nargin < 3
        sprintf('Error: must provide wave, anntime and annot');
        annot=[];
        template=[];
        valid=0;
        return;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 19/09/2011 ADD baseline wander filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wave = PPGmedianfilter(wave, Fs,Fs);
    % get PPG template
    [t t2 v]=template_pleth(wave(1:min(windowlen,length(wave))), anntime(find(anntime<min(windowlen,length(wave)))), 0,Fs);

    if v<1 && length(template)<1 % Current template invalid && no previous template available
        for j=1:length(annot)
            annot(j).typeMnemonic='Q'; % Unacceptable
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
            for j=1:length(annot)-1

%             SQI1: Direct compare

%             Calculate correlation coefficients based on the template
%             length
                beatbegin=anntime(j);
                beatend=anntime(j+1);
% 07/11/2011 ADD max beat length <= 3s detection 
		if beatend-beatbegin>3*Fs
			beatend=beatbegin+3*Fs;
		end
		templatelength=length(t);
		complength=min(templatelength,beatend-beatbegin-1);
                if beatbegin+complength-1 > length(wave) || beatend > length(wave) || beatbegin < 1
                    continue;
                end
                currentb=j;
                cc=corrcoef(t(1:complength),wave(beatbegin:beatbegin+complength-1));
                c1(j)=cc(1,2);
                if (c1(j)<0)
                    c1(j)=0;
                end
                annot(currentb).subtype=int8(c1(j)*100);


%             SQI2: Linear resampling

%             Calculate correlation coefficients based on the linear
%             resampling (interp1)
            
                y=interp1(1:beatend-beatbegin, wave(beatbegin:beatend-1),1:(beatend-beatbegin-1)/(templatelength-1):(beatend-beatbegin),'spline');
                y(isnan(y))=0;
                cc=corrcoef(t,y);
                c2(j)=cc(1,2);
                if (c2(j)<0)
                    c2(j)=0;
                end
                annot(currentb).chan=int8(c2(j)*100);

%             SQI3: Dynamic Time Warping                
                
%             Calculate correlation coefficients based on the dynamic time
%             warping
            
                d2=wave(beatbegin:beatend-1);
            
                % if beat too long, set SQI=0;
                if (length(d2)>length(d1)*10)
                    c3(j)=0;
                else
                    d2=(d2-min(d2))/(max(d2)-min(d2)).*100;
                   [y2 pla2]=PLA(d2,1,1);
            
                   [w ta tb] = simmx_dtw(y1,pla1,y2,pla2);
                   [p,q,Dm] = dp_dwt2(w,ta,tb);
                   [ym1 ym2 yout1]=draw_dtw(y1,pla1,p,y2,pla2,q);
                    cc=corrcoef(y1,ym2);
                    c3(j)=cc(1,2);
                    if (c3(j)<0)
                        c3(j)=0;
                    end
                end
                annot(currentb).num=int8(c3(j)*100);
                
%             SQI4: Clipping detection   

                d2=wave(beatbegin:beatend-1);
                y=diff(d2);
                clipthreshold=0.5;
                c4(j)=int8(length(find(abs(y)>clipthreshold))/length(y)*100);
                
%             SQI: Combined

                sqibuf=[annot(currentb).subtype annot(currentb).chan annot(currentb).num c4(j)];
                if min(sqibuf)>=90
                    annot(currentb).typeMnemonic='E'; % Excellent
                else
                    if length(find(sqibuf>=90))>=3 || (median(sqibuf(1:3))>=80 && sqibuf(1)>=50 && sqibuf(4)>=70) || min(sqibuf)>=70
                        annot(currentb).typeMnemonic='A'; % Acceptable
                    else
                        annot(currentb).typeMnemonic='Q'; % Unacceptable
                    end
                end
                annot(currentb).auxInfo=num2str(int8(c4(j)));
                
            end

    end
    template=t;
    valid=v;
