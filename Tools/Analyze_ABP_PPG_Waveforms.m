function [rr,t,sqi] = Analyze_ABP_PPG_Waveforms(Waveform,Type,HRVparams,detectedQRS,subjectID)
%
%   Analyze_ABP_PPG_Waveforms(Waveform,Type,HRVparams,detectedQRS,subjectID)
%	OVERVIEW:
%       Analyze ABP or PPG waveform 
%
%   INPUT:
%       Waveform    - matrix containing the a raw signal in each column
%       Type        - array containing the signal type of waveforms in each column: 
%                     'APB' for ABP waveform 
%                     'PPG' for PPG waveform
%       HRVparams   - struct of settings for HRV analysis  
%       subjectID   - string to identify current subject
%
%   OUTPUT
%       rr          - (seconds) Vector containing RR interval
%       t           - (seconds) Time of the rr interval data 
%       sqi         - (optional) Signal Quality Index; Requires a 
%                     matrix with at least two columns. Column 1 
%                     should be timestamps of each sqi measure, and 
%                     Column 2 should be SQI on a scale from 0 to 1.
%   
%   DEPENDENCIES & LIBRARIES:
%       PhysioNet Cardiovascular Signal Toolbox
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%
%   REFERENCE: 
%   Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%   Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Giulia Da Poian (giulia.dap@gmail.com) on Sep 6, 2017.
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


NmbOfSigs = size(Waveform,2);
AnnotationFolder = strcat(HRVparams.writedata, filesep, 'Annotation', filesep);
if ~exist(AnnotationFolder, 'dir')
   mkdir(AnnotationFolder)
end
addpath(AnnotationFolder)

rr = [];
t = [];
sqi = [];

for i = 1:NmbOfSigs
    
    current_type = Type{i};
    switch current_type
        case 'PPG'
            % PPG Detection - qppg
            [PPGann] = qppg(Waveform(:,i),HRVparams.Fs);
            % PPG SQI 
            [ppgsqi_numeric, ~, ppgsqi]= calculate_ppgsqi(PPGann,Waveform(:,i),HRVparams.Fs);
            % Write PPG  annotations
            write_ann(strcat(AnnotationFolder, subjectID),HRVparams,'ppg',PPGann);
            write_ann(strcat(AnnotationFolder, subjectID),HRVparams,'sqippg',PPGann(1:length(ppgsqi_numeric)),char(ppgsqi),ppgsqi_numeric);
            
            rr = diff(PPGann)./HRVparams.Fs;
            t = PPGann(2:end)./HRVparams.Fs;
            sqi = [PPGann'./HRVparams.Fs, ppgsqi_numeric'./100];

        case 'ABP'
            % ABP
            ABPann = run_wabp(Waveform(:,i));
            % ABP SQI
            ABPfeatures =  abpfeature(Waveform(:,i), ABPann, HRVparams.Fs);
            [BeatQ, ~] = jSQI(ABPfeatures, ABPann, Waveform(:,i));
            
            if ~isempty(detectedQRS)
                % Pulse Transit Time
                ptt = pulsetransit(detectedQRS, ABPann);
                % Plot BP vs PTT
                syst = ABPfeatures(:,2);
                if HRVparams.gen_figs
                    figure;
                    plot(syst,ptt(:,3)./HRVparams.Fs,'o');
                    xlabel('BP (mmHg)'); ylabel('PTT (s)');
                    title('Pulse Transit Time - BP vs PTT (ABP - QRS)')
                end
            end
            % Write ABP annotations
            write_ann(strcat(AnnotationFolder, subjectID),HRVparams,'abpm',ABPann);
            write_ann(strcat(AnnotationFolder, subjectID),HRVparams,'sqiabp',ABPann,'N',int8(BeatQ(:,1)));%write_ann(strcat(AnnotationFolder, subjectID),HRVparams,'sqiabp',BeatQ(:,1));
                      
    end
end
