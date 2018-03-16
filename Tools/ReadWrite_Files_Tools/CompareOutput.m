function test = CompareOutput(file1,file2)
    
%   test = CompareOutput(file1,file2)
%	OVERVIEW: This function returns test=1 if the two files are equal 
%             test=0 otherwise
%   INPUT:
%       file1    - string with the reference file name 
%       file2    - string with the test file name 
%   OUTPUT 
%       test      - return 0 if the two files are different, 1 if they are
%                   the same
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       This script written by Giulia Da Poian
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
% 
if isempty(file1) || isempty(file2)
    test = 0;
    return;
else
    currentFile = javaObject('java.io.File', file1 );
    referenceFile = javaObject('java.io.File', file2);
    test = javaMethod('contentEquals','org.apache.commons.io.FileUtils', ...
        referenceFile, currentFile);
end

% Windows system return test = 0 when the two files are equal so we need to
% invert the results for cross-.platform compatibility
if strcmp(computer, 'PCWIN64')
    test = not(test);
end   
