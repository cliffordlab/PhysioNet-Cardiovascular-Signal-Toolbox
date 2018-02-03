function [TO, TS, TOsingle, TSsingle] = HRT_Analysis(RRInts,Labels,BeatsBefore, BeatsAfter, GraphOn)

%   [TO, TS, TOsingle, TSsingle] = HRT_Analysis(RRInts,Labels,BeatsBefore, BeatsAfter, GraphOn)
%   OVERVIEW:
%       This function return TO and TS, i.e., the basic parameters of 
%       heart rate turbulence (HRT) used to quantify the return to 
%       equilibrium of heart rate  after a premature ventricular 
%       contraction (PVC).
%       TO is the Turbulence Onset, and it is defined as the percentage 
%       difference between average value of the first two normal RR
%       following the PVC and the last two normal intervals preciding the
%       PVC
%       TS is the Turbulence Slope, and it is calculated by constructing an
%       average ectopy-centered time series and determining the steepest
%       slope each possible five consecutive normal intervals in the
%       post-PVC 
%
%   INPUTS:
%       RRInts          : Vector containing RR intervals data (in seconds) 
%       Labels          : Vector containing annotations of the RR data at 
%                         each point indicating the type of the beat (see 
%                         https://www.physionet.org/physiobank/annotations.shtml)
%       BeatsBefore     : number of RR intervals to consider before the PVC
%                         (default 2)
%       BeatsBefore     : number of RR intervals to consider after the PVC
%                         (default 16)
%       GraphOn         : (optional) set GraphOn = 1 for a plot of the HRT analys 
%
%   OUTPUTS:
%       TO         : average turbulence onset (TO) 
%       TS         : turbulence slop (TS) of the average tachogram 
%       TOsingle   : a vector containing the turbulence onset (TO) of each PVC
%       TSsingle   : a vector containing the turbulence slop (TS) of each PVC
%
%	REPO:       
%       https://github.com/cliffordlab/Physionet-HRV-toolbox-for-MATLAB
%   ORIGINAL SOURCE AND AUTHORS:     
%       Giulia Da Poian    
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
  

if nargin<2 
    error('')
end
if nargin<3
    BeatsBefore = 2;
end
if nargin<4
    BeatsAfter = 16;
end
if nargin<5
    GraphOn = 0;
end

% Convert to ms
RRInts = RRInts *1000;

% Find 
PVCs = find(Labels=='V'); 


% Pre-process the original RR intervals to get tachograms time series from
% which unphysiological beats and not Normal 'N' beats are remove, also
% computes the average tachogram

[tachograms, avgTachogram] = HRT_preprocessing(RRInts,Labels, PVCs, BeatsBefore, BeatsAfter); 
% Analyze tachograms and avgTachogram to get HRT values
  
[TO, TS, TOsingle, TSsingle] = computeHRT(tachograms, avgTachogram,BeatsBefore);

% Plots
if GraphOn
    PlotHRT(tachograms, avgTachogram, BeatsBefore,TO,TS)
end

end % end of Eval_HRT function



function [CleanTachos, AvgTacho] = HRT_preprocessing(RRInts,Labels, PVCs, BeatsBefore, BeatsAfter)       
        
% Function that performs a complete HRT preprocessing steps that includes:
%     1.  Find all VPC-tachograms (pre_pvc_beats+PVC+CP+post_pvc_beats)
%     2.  Filter all VPC-tachograms; i.e. find all the availables VPC-tachograms
%         to be used on the HRT analysis
% REFERENCE: Clifford et al. 'Advanced Methods and Tools for ECG Data
% Analysis' (pp. 89-92)

    % Initialize variables for:
    Tachos = [];             % all possible tachograms
    TachoLabels = [];      % and their labels
    CleanTachos = [];
    
    idxPVC = BeatsBefore+1;
    idxCP = idxPVC+1;
    
    %  Guarantee last VPC position allows post_pvc_beats beats after it.        
    if PVCs(end) + BeatsAfter > length(RRInts)
        PVCs(end) = [];         
    end

    for m = 1: length(PVCs)
        idxTachos = PVCs(m)-BeatsBefore : PVCs(m)+1+BeatsAfter;
        if idxTachos(end) >= sum(RRInts) - BeatsAfter
            continue
        end
        Tachos = [Tachos RRInts(idxTachos)];            
        TachoLabels = [TachoLabels Labels(idxTachos)];          
    end

    % Filtering process of tachograms:        
        % Every possible tachogram is going to be "filtered" 
        % to determine if it is a valid tachogram to be used 
        % in HRT analysis

    % Indices used to check condition (excluding the PVC and the one after the PVC)
    idxSinus = [1:BeatsBefore idxCP+1:idxCP+BeatsAfter];

    for idx = 1:size(Tachos,2)

        tmpTacho = Tachos(:,idx);
        tmpLabs = TachoLabels(:,idx);
        
        % compute all the conditions that a VPC tachogram must fulfill
        % All beats, except for VPC, must be 'N'
        if ~CheckForNormalBeats(tmpLabs, BeatsBefore)
                continue;
        end          
           
        % Do not consider RR intervals < 300ms or > 2000ms
        if  ~isempty(find(tmpTacho(idxSinus) < 300,1)) 
            continue;
        end
        if ~isempty(find(tmpTacho(idxSinus) > 2000,1))
            continue;
        end

        % Do not consider RR intervals where || RR(n-1) - RR(n)|| > 200ms
        if ~isempty(find(abs(diff(tmpTacho(1:BeatsBefore)))>=200,1)) || ...
                    ~isempty(find(abs(diff(tmpTacho(BeatsBefore+3:end)))>=200,1))
            continue;
        end

        % Reference RR interval, mean of the N RR-intervals precedding the VPC
        refRRinterval = mean(tmpTacho(1:BeatsBefore));
        % All the RR-intervals should be lower than 1.2*refRRinterval
        if ~isempty(find(tmpTacho(idxSinus) >1.2*refRRinterval,1))
            continue;
        end

        % PVC should be at least 20% shorter than refRRinterval, i.e must be 
        % lower than 80% than refRRinterval
        if tmpTacho(BeatsBefore+1) > 0.8*refRRinterval
            continue;
        end
        % extrasystolic pause should be greater than 1.2*refRRinterval
        if tmpTacho(BeatsBefore+2) < 1.2*refRRinterval
            continue;
        end

        CleanTachos = [CleanTachos  tmpTacho];
    end

    if isempty(CleanTachos) 
        % there is no tachograms that fulfill all conditions
        AvgTacho = [];
    else
        AvgTacho = mean(CleanTachos,2);
    end
end


function validTach = CheckForNormalBeats(tachLab, rri_before)
       
    aux = tachLab([1:rri_before rri_before+2:end]);
    cond1 = tachLab(rri_before+1) == 'V';
    OtherBeats = find(aux ~= 'N',1);
    if isempty(OtherBeats) && cond1
        validTach =  1;
    else
        validTach = 0;
    end
    
end


function [TS, iTS, p]= TurbulenceSlope(Tach,BeatsBefore) 

        % make sure Tach is an Nx1 vector
        if size(Tach,1)<size(Tach,2)
            Tach = Tach';
        end
        
        %Computes the Turbulence Slope on the tachogram given by parameter
        SegLeng = 5;
        posPC = BeatsBefore+1;
        posFin = length(Tach) - (SegLeng -1) ;

        slopes = zeros(1,length(Tach));        
        
        idx = 1;
        for m = posPC+2:posFin
            seg = Tach(m:m+SegLeng-1);
            p(idx,:) = polyfit(1:5, seg', 1);
            slopes(idx) = p(idx,1);
            idx =idx+1;
        end
        
        [TS, iTS] = max(slopes);            
end
        
        
function [TO, num, den] = TurbulenceOnset(Tach,PVCidx)
        
        % Computes the Turbulence Onset on the tachogram given by parameter
        num = (Tach(PVCidx+2)+Tach(PVCidx+3))/2-(Tach(PVCidx-1)+Tach(PVCidx-2))/2;
        den = Tach(PVCidx-1)+Tach(PVCidx-2)/2;
        TO = (num/den)*100; 

end
 

function [TO, TS, TOsingle, TSsingle ] = computeHRT(CleanTachos, AvgTacho, BeatsBefore)
    
    %Computes the TS and TO on all conditions tachograms (tachograms_ok) and on the mean tachogram
    TSsingle = nan(1,size(CleanTachos,2));
    TOsingle = nan(1,size(CleanTachos,2));
    if sum(isnan(AvgTacho)) > 0
        % if there is no tachograms that fulfill all conditions, the value of the mean tachogram is nan,
        % so the TS and TO values will be nan too
        TS = NaN;
        TO = NaN;
    else            
        % compute TS and TO for each VPC tachogram
        for tac = 1 : size(CleanTachos,2)
            TSsingle(tac) = TurbulenceSlope(CleanTachos(:,tac),BeatsBefore);
            TOsingle(tac) = TurbulenceOnset(CleanTachos(:,tac),BeatsBefore+1);       
        end
        % compute the TS and TO from the average VPC tachogram    
        TS = TurbulenceSlope(AvgTacho,BeatsBefore);
        TO = mean(TOsingle);
    end
end
         

function PlotHRT(tachograms, avgTachogram, pre_pvc_beats,TO,TS)

    N = size(tachograms,2);

    figure1 = figure('Position',[150 300 1400 600]);

    subplot1 = subplot(1,2,1,'Parent',figure1);
    plot(tachograms,'Parent',subplot1,'Marker','x');
    xlabel('# of RR Interval');ylabel('RR Intervals [ms]');
    title({['RR interval profiles of ' num2str(N) ' ventricular premature'], [' complexes (VPC) aligned by VPC position']})
    set(subplot1,'FontSize',18);
    
    % Create subplot
    subplot2 = subplot(1,2,2,'Parent',figure1);
    hold(subplot2,'on');

    plot(avgTachogram,'Marker','o','LineWidth',1.5)
    xlabel('# of RR Interval');ylabel('RR Intervals [ms]');
    set(subplot2,'FontSize',18);

    % Plot TO
    l1 = (avgTachogram(pre_pvc_beats-1)+ avgTachogram(pre_pvc_beats))/2;
    x1 = pre_pvc_beats-2 : pre_pvc_beats+3;
    line(x1, l1*ones(length(x1)),'LineWidth',2,'LineStyle','--',...
    'Color',[0.929411768913269 0.694117665290833 0.125490203499794]);
    l2 = (avgTachogram(pre_pvc_beats+3)+ avgTachogram(pre_pvc_beats+4))/2;
    x2 = pre_pvc_beats : pre_pvc_beats+5;
    line(x2, l2*ones(length(x2)),'LineWidth',2,'LineStyle','--',...
    'Color',[0.929411768913269 0.694117665290833 0.125490203499794]);
    
    % Plot TS
    [TS, iTS,p] = TurbulenceSlope(avgTachogram,pre_pvc_beats);
    f = polyval(p(iTS,:),1:5);
    plot(pre_pvc_beats+1+iTS:pre_pvc_beats+1+iTS+4,f,'LineWidth',3,...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625])
    title('Average Ectopy-Centered Time Series')

    % Add TS and TO values
    descr = {['Turbulence Onset TO:' num2str(TO) '%'] ;
             ' ';
             ['Turbulence Slope TS:' num2str(TS)]};
    % Create textbox
    annotation(figure1,'textbox',...
    [0.74 0.13 0.13 0.2],'Color',[1 1 1],...
    'String',descr,'LineWidth',2,'FitBoxToText','off','EdgeColor',[0.08 0.17 0.56],...
    'BackgroundColor',[0 0.45 0.75],'FontSize',18);
         
          



end

