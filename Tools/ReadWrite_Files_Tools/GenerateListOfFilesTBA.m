function [ids,fullfilenames] = GenerateListOfFilesTBA(ext,datadir,recursive)
%
%[ids,fullfilenames] = GenerateListOfFilesTBA(ext, datadir, recursive)
%
%   OBJECTIVE:
%       This function takes the extension ext and finds all relevant 
%       files within the working directory and the subfolders of the 
%       working directory. Works with records that are stored with a
%       filename in the format '####.ext'
%
%   INPUT:  ext       - string variable of the extension of filetype
%           	        example: 'mat' 'qrs'
%           datadir   - string variable representing the name of the
%           	        directory of the data
%           recursive -
%
%   OUTPUT: ids          - vector of identifiers (doubles) listing all records
%           fulfilenames - full file name for each record(strings)
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Main script written by Adriana N. Vest
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

% Check for folders inside the data directory so that a decision can be
% made about whether to include subfolders or not


if strcmp(computer,'GLNXA64')
    if recursive
    	[~, sysout] = system(['ls ' datadir filesep '*.',ext]);
    else
        [~, sysout] = system(['ls ' datadir filesep '*.',ext]);
        % ANV modified this line to work in linux for MIT NSR files on box
        % May need to revisit
    end
    allfiles = strsplit(sysout);
elseif strcmp(computer,'MACI64') || strcmp(computer,'PCWIN64')
    allfiles1 = [];
    if recursive
        allfiles1 = dir(fullfile(datadir, '*', filesep, ['*.' ext]));
    end
	allfiles2 = dir(fullfile(datadir, ['*.' ext]));
    allfiles = [allfiles1;allfiles2];
else
    error('no files or records found');
end

% Initialize arrays to store data
ids = cell(length(allfiles), 1); % Create ids vector to store subject identifiers
fullfilenames = cell(length(ids),1); % Create cell array to store filenames

for i = 1:length(fullfilenames)
    try
        if strcmp(computer,'GLNXA64')
            filename = dir(allfiles{i});
            fullfilenames{i} = [allfiles{i}];
            
            filename = filename.name;
            filename = strrep(filename, ['.' ext], '');
            ids{i} = filename;
        elseif strcmp(computer,'MACI64')
            fullfilenames{i} = [allfiles(i).folder filesep allfiles(i).name];
            filename = allfiles(i).name;
            filename = strrep(filename, ['.' ext], '');
            ids{i} = filename;
        elseif strcmp(computer,'PCWIN64')
            fullfilenames{i} = [datadir filesep allfiles(i).name];
            filename = allfiles(i).name;
            filename = strrep(filename, ['.' ext], '');
            ids{i} = filename;
        end
    catch
        
    end
end

if strcmp(computer,'GLNXA64')
	%idx_rem = find((isempty(ids))); %% This line isn't working on Linux
    emptyCells = cellfun(@isempty,ids);
    idx_rem = find(emptyCells);
    ids(idx_rem) = [];
	fullfilenames(idx_rem) = [];
end

if isnumeric(ids{1})
    for k = 1:length(ids)
        ids{k} = num2str(ids{k});
    end
end

