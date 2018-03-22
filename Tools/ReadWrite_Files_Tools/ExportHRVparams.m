function [T, Tab] = ExportHRVparams(HRVparams)
%   ExportHRVparams2Latex(HRVparams)
%
%   INPUT:      
%       HRVparams - struct of various settings for the hrv_toolbox analysis
%
%   OUTPUT:
%       Latex file  
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Script written by Giulia Da Poian
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

% Note that if you change the order of the parameters or add parameters 
% this might not work

% Windows analysis parametrs 

GeneralParams = struct2table(HRVparams,'AsArray',true);

NameOfParameters = GeneralParams.Properties.VariableNames;
ValuesOfParameters = table2cell(GeneralParams);
T = [];
for idx=1:length(ValuesOfParameters)
    
    if isstruct(GeneralParams{1,idx})
    % extract structure and concatenate to Table
      
        temp = struct2table(GeneralParams{1,idx},'AsArray',true);
        tempName = strcat(NameOfParameters(idx),'_',temp.Properties.VariableNames);
        tempVal = table2cell(temp);
        % Check for matrix instead of a single vector 
        % This works only for this particular program for extracting the
        % frequency bands
             for kk =1:length(tempVal)
                 if ~isvector(tempVal{kk}) && ~isempty(tempVal{kk})
                     extractedArray = tempVal{kk};
                     freqLimits = cellstr(strcat(num2str(extractedArray(:,1)),'-', num2str(extractedArray(:,2))));
                     T = [T array2table(freqLimits','VariableNames',{'ulf','vlf','lf','hf'})];
                     tempVal{kk} = "[]";
                 elseif ~isnumeric(tempVal{kk})
                     tempVal{kk} = string(tempVal{kk});
                 elseif isempty(tempVal{kk})
                     tempVal{kk} = "[]";
                 end
             end
               
        try
           T = [T table( tempVal{:},'VariableNames',{tempName{:}})];
        catch
        end
    else
        T = [T GeneralParams(:,idx)];
    end
end

TableNames = T.Properties.VariableNames;
% Change Names to be in latex format
TableNames = strrep(TableNames,'_',' ');

TableValues = table2cell(T);
% Change Names to be in latex format
for idx = 1: length(TableValues)
    if iscellstr(TableValues(idx))
        TableValues(idx) = strrep(TableValues(idx),{'\'},{'$\$'});
        TableValues(idx) = strrep(TableValues(idx),{'_'},{' '});
    elseif ~ischar(TableValues{1,idx}) && ~isstring(TableValues{1,idx})
        TableValues{idx} = char(num2str(TableValues{1,idx}));
    end
end

Tab = [TableNames' TableValues'];

writetable(array2table(Tab),[HRVparams.writedata filesep 'ParametersTable' HRVparams.filename '.csv']);
try
    matrix2latex(Tab, [HRVparams.writedata filesep 'ParametersTable' HRVparams.filename '.tex']);
catch
end