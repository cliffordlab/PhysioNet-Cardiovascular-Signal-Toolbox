function Plot_SignalDetection_SQI(time, signal, ann, sqi,signalType)

% Plot_SignalDetection_SQI(time, ecg, qrs, sqi)
%
% INPUTS
%       time       : vector time
%       signal     : vector signal 
%       ann        : vector containing signal evenu locations, e.g., QRS 
%                    locations, PPG onsets
%       sqi        : vector containing Signal Quality Index (sqi) values 
%       signalType : string containing the type of signal, e.g., 'ECG', 
%                    'ABP', 'PPG'
%	REPO:       
%       https://github.com/cliffordlab/Physionet-HRV-Toolbox-for-Matlab
%   ORIGINAL SOURCE AND AUTHORS:     
%       Giulia Da Poian
%
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% New figure
figure1 = figure('position', [100 200 1200 500]);

% First supblot
axes1 = axes('Parent',figure1, 'Position',[0.13 0.39259 0.775 0.45]);
hold(axes1,'on');

% Plot ECG
plot(time,signal,'Parent',axes1,'DisplayName', signalType); 
% Plot QRS locations over ECG
plot(time(ann), signal(ann),'Parent',axes1,'DisplayName','Onsets','Marker','o','LineStyle','none');

% Set axes
set(axes1,'FontSize',18,'XTick',zeros(1,0));
if signalType == 'ECG'
    title('ECG - Detected QRS complexes and SQI')
    ylabel('Amplitude (mV)');
else
    ylabel('Amplitude');
    title([signalType '- Detected beat onsets and SQI'])
end

set(axes1,'FontSize',18);
legend1 = legend(axes1,'show');
set(legend1,'Location','southeast');

% Second subplot for SQI
axes2 = axes('Parent',figure1, 'Position',[0.13 0.28 0.775 0.08]);
clims = [0 1];
imagesc([0.5 291],0,sqi,'Parent',axes2,'CDataMapping','scaled',clims)
xlabel('Time (s)');
ylh=ylabel('SQI');
pos1=get(ylh,'Position'); pos1(1,2)=pos1(1,2)-0.05;
set(ylh,'Position',pos1)
axis(axes2,'ij');
% Set the remaining axes properties
% set(axes2,'DataAspectRatio',[10 1 1],'FontSize',18,'Layer','top','YTick',zeros(1,0));
set(axes2,'FontSize',18,'Layer','top','YTick',zeros(1,0));
% Create colorbar
colorbar('peer',axes2,'Position',[0.914 0.28 0.0139 0.0829]);

linkaxes([axes1,axes2],'x')