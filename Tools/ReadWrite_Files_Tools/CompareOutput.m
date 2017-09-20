function test = CompareOutput(file1,file2)
    
% This function returns test=1 if the two files are equal test=0 otherwise
    
referenceFile = javaObject('java.io.File', file1 );
currentFile = javaObject('java.io.File', file2);
test = javaMethod('contentEquals','org.apache.commons.io.FileUtils', ...
    referenceFile, currentFile);


% Windows system return test = 0 when the two files are equal so we need to
% invert the results for cross-.platform compatibility
if strcmp(computer, 'PCWIN64')
    test = not(test);
end   
