# hrv_toolbox

------- Quick Start: ---------
1)  Review InitializeHRVparams.m and optimize the parameters for your 
    data. 
2)  The toolbox does not assume any format of data except that the input 
    of the Main_VOSIM.m fucntion are a two equal length vectors: RR interval
    and time in units of seconds or the ?raw? ECG signal (physical units,mV) 
    and time. 
    Additionaly, blood pressure waveform and photoplethysmographic/pulsatile
    data can be analyzed and they should be in the standard physical units 
    (mmHg or normalized units respectively). 
3)  Results will be stored in a subdirectory of the project's data file
    in a folder called as indicated in the InitializeHRVparams.m. 
    If the folder does not exist, it will be created.

------ Guide to Output: ------
The following metrics are output from the HRV Toolbox:]
	NNmean
	NNmode
	NNmedian
	NNskew
	NNkurt
	NNiqr
	t_win
	SDNN
	RMSSD
	pnn50
	fdflag      1 = lomb failed
        	    2 = not enough high SQI data in the window to process
                	(amount of data above threshold1 is greater than threshold2)
            	    3 = not enough data in the window 
                    4 = window is missing too much data
                    5 = success
	ulfL
	vlfL
	lfL
	hfL
	lfhfL
	ttlpwrL
	SDANN
	SDNNI

----- Full Instructions: -----
I. Introduction
The HRV Toolbox is a a cardiovascular dynamics analysis package, designed 
to meet the need in the clinical and scientific community for a validated, 
standardized, well-documented open-source toolkit to evaluate the 
relationships between physiological signals and disease. The package not 
only includes standard HRV tools to generate time and frequency domain 
metrics from ECG or pulsatile waveforms (like the blood pressure or 
photoplethysmographic waveforms), but more recent metrics such as 
acceleration and deceleration capacity and pulse transit time. The package
is designed to accommodate a variety of input data, from raw unprocessed 
and unannotated waveforms, to fully annotated tachogram data. In general, 
abnormal beat and noise removal and methods for dealing with the missing 
data are poorly described and highly variant in most of the literature. 
Therefore, we have included signal processing methods that include state 
of the art peak detectors, signal quality processing units, and beat/rhythm 
phenotyping.  The package can also analyze the interactions between 
multiple physiological signals.

II. Getting Started
System requirements:
- Matlab and License    https://www.mathworks.com/
- WFDB Toolbox          https://physionet.org/physiotools/wfdb.shtml
- WFDB Matlab Toolbox   https://www.physionet.org/physiotools/matlab/wfdb-app-matlab/

1)  Download and install Matlab, the WFDB Toolbox, and the WFDB Toolbox
    for Matlab.

2) Requires rrgen binary - compilation of rrgenV3.c on your system:
        1.  Compile rrgen
            Navigate to rrgen in HRV Toolbox & Compile using gcc
            gcc -Wall rrgenV3.c -lm -o rrgen
                or
            gcc -Wall -o rrgenV3 rrgenV3.c 
        2.  Ensure executable is on the system path, or move executable to
            usr/local/bin or similar location on the path
        3.  Ensure executable is on Matlab's path using the addpath fn

III. Starting Analysis

IV. FAQ
