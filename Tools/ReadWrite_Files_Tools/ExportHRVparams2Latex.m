function ExportHRVparams2Latex(HRVparams)
%   ExportHRVparams2Latex(HRVparams)
%
%   INPUT:      
%       HRVparams - struct of various settings for the hrv_toolbox analysis
%
%   OUTPUT:
%       Latex file  
%
%   DEPENDENCIES & LIBRARIES:
%       HRV_toolbox https://github.com/cliffordlab/hrv_toolbox
%       WFDB Matlab toolbox https://github.com/ikarosilva/wfdb-app-toolbox
%       WFDB Toolbox https://physionet.org/physiotools/wfdb.shtml
%   REFERENCE: 
%	REPO:       
%       https://github.com/cliffordlab/hrv_toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Script written by Giulia Da Poian
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%



% Table for frequency parameter
freqParams = struct2table(HRVparams.freq,'AsArray',true);
freqParams = freqParams(:,2:end); %  The first column contains a struct