function filename = SaveHRVoutput(sub_id,windows_all,results,titles,type,HRVparams, tNN, NN)
%   
%   SaveHRVoutput(sub_id,windows_all,results,titles,type,HRVparams, tNN, NN)
%
%   OVERVIEW:   Generates output based on the selections in struct HRVparams.
%
%   INPUT:      sub_id       :
%               windows_all  :
%               results      :
%               titles       :
%               type         : 'AF' results of AF detection 
%                              'MSE' results of Multiscale Entropy
%                              'SQI' average SQI index for the record
%                               [] otherwise
%               HRVparams    : struct of settings for hrv_toolbox analysis
%               tNN          : the time indices of the rr interval data (seconds)
%               NN           : a single row of NN (normal normal) interval
%                              data in seconds
%
%   OUTPUT:    filename      : current name of the file 
%              Outputs csv files or mat files based on the user's HRVparams   
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%
%   09-18-2017 - Modified by Giulia Da Poian
%   Convert to table with results in columns and then save csv, it makes  
%   easier to extract the results to perform statistical analysis from the
%   csv file (columns have the name of the variable, e.g.,SDNN, AC, DC...)

% What class of data is sub_id
if ischar(sub_id)
end
if isnumeric(sub_id)
    sub_id = num2str(sub_id);
end
if iscell(sub_id)
    sub_id = char(sub_id);
end



% Establish Filename Based on Type of Output
if strcmp(type,'AF') || strcmp(type,'MSE') || strcmp(type,'SQI')|| strcmp(type, 'DFA') || strcmp(type, 'HRT')
    filename = strcat( sub_id, '_', type, '_results_', HRVparams.time);
    if strcmp(HRVparams.output.format,'csv')
        % Add .csv extension to filename and directory
        fullfilename = strcat(HRVparams.writedata, filesep, filename, '.csv'); 
        % Write AF results to a table 
        T =  array2table(results,'VariableNames',titles);
        % Use writetable to geberate csv file with the results   
        writetable(T,fullfilename);        
    elseif strcmp(HRVparams.output.format,'mat')
        % Add .mat extension to filename and directory
        fullfilename = strcat(HRVparams.writedata, filesep, filename, '.mat');
        %Save results
        save(fullfilename, 'results', 'titles');
    end
      
else % HRV results
   
    % All windows or Lowest HR windows results 
    if ~isempty(HRVparams.output.num_win) 
        fileNameWind = strcat('HRV_results_', num2str(HRVparams.output.num_win), 'LowestHRwin_', HRVparams.time);
    else
        fileNameWind = strcat('HRV_results_allwindows_', HRVparams.time);
    end
    % All patients or Separate
    if HRVparams.output.separate
        % Generate a new file for each output
        filename = strcat(sub_id, '_', fileNameWind);
    else
        filename = strcat('AllPatients_', fileNameWind);
    end
    
    % Generate csv file
    if strcmp(HRVparams.output.format,'csv')  
        % Add .csv extension to filename
        fullfilename = strcat(HRVparams.writedata, filesep, filename, '.csv');    
        
        if ~isempty(HRVparams.output.num_win) 
            % Returns results based on the number of windows set by the HRVparams file
            x = size(results);
            idx = find(length(titles) == x);
            num_results = x(3-idx);
            
            % Find num_win windows with the lowest HR
            if num_results > 1
                windows_output = FindLowestHRwin(windows_all,tNN, NN, HRVparams.output.num_win,HRVparams);
                for i = 1:HRVparams.output.num_win
                    windx(i) = find(windows_output(i).t == windows_all);  
                end
                variables_names = ['patID' titles]; % Add colum with patient ID
                results = results(windx,:);
                patid_array  = string(repmat({sub_id},size(results,1),1));
                variables_vals = [patid_array results]; % Add colum with patient ID
            else
                variables_names = ['patID' titles]; % Add colum with patient ID
                patid_array  = string(repmat({sub_id},size(results,1),1));
                variables_vals = [patid_array results]; % Add colum with patient ID
            end

            % Write HRV results to a table 
            try
                T = readtable(fullfilename);
            catch
                T = [];
            end
            T = [T ; array2table(variables_vals,'VariableNames',variables_names)];
            % Use writetable to geberate csv file with the results   
            writetable(T,fullfilename);
                
        elseif isempty(HRVparams.output.num_win) 
        
            % Print out all the window values for all variables
            variables_names = ['patID' titles]; % Add colum with patient ID
            patid_array  = string(repmat({sub_id},size(results,1),1));
            variables_vals = [patid_array results]; % Add colum with patient ID

            % Write HRV results to a table
            try
                T = readtable(fullfilename);
            catch
                T = [];
            end
            T = [T ; array2table(variables_vals,'VariableNames',variables_names)];
            % Use writetable to geberate csv file with the results   
            writetable(T,fullfilename);

        else
            % Do nothing.
        end % End decision based on number of windows needed to be returned

    elseif strcmp(HRVparams.output.format,'mat')
        % Add .mat extension to filename
        fullfilename = strcat(HRVparams.writedata, filesep, filename, '.mat');
    
        if ~isempty(HRVparams.output.num_win)
            windows_output = FindLowestHRwin(tNN, NN,HRVparams.output.num_win);

            for i = 1:HRVparams.output.num_win
                windx(i) = find(windows_output(i).t ==  windows_all);  
            end

            output = results(wind,:);
            save(fullfilename, 'output','titles');

        elseif isempty(HRVparams.output.num_win) 

            save(fullfilename, 'results', 'titles');
        else
        end
    else
        % Return nothing.
    end
end



end
