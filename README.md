![alt text](https://zenodo.org/badge/DOI/10.5281/zenodo.1243112.svg)

# PhysioNet Cardiovascular Signal Toolbox

1. [Introduction](intro)
2. [Full Instructions](fullinst)
3. [Guide to Output](output) 
4. [Contributing to this project](cont)
5. [FAQ](faq)

If you are using this software, please cite:
```
 Vest A, Da Poian G, Li Q, Liu C, Nemati S, Shah A, Clifford GD, 
 "An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", 
  Physiological measurement 39, no. 10 (2018): 105004. DOI:10.5281/zenodo.1243111; 2018. 
```       
<a name="intro"></a>
## Introduction
The **PhysioNet Cardiovascular Signal Toolbox** is a a cardiovascular dynamics analysis package, designed 
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
phenotyping. The package can also analyze the interactions between 
multiple physiological signals.

<a name="fullinst"></a>
## Full Instructions: 

### I. Getting Started
System requirements:

- Matlab and License    https://www.mathworks.com/

1)  Download and install Matlab 2017b (v9.3) (required Matlab Toolboxes: 
    Signal Processing Toolbox, and Statistics and Machine Learning Toolbox, 
    Neural Network Toolbox)

2)  Add the PhysioNet Cardiovascular Signal Toolbox to your
    Matlab path: **run startup.m**

3)  (Optional) rrgen binary - compilation of rrgenV3.c on your system:

        1.  Compile rrgen
            Navigate to rrgen in HRV Toolbox & Compile using gcc
            gcc -Wall rrgenV3.c -lm -o rrgen
                or
            gcc -Wall -o rrgenV3 rrgenV3.c 
        2.  Ensure executable is on the system path, or move executable to
            usr/local/bin or similar location on the path
        3.  Ensure executable is on Matlab's path using the addpath fn

### II. Starting Analysis

#### Quick Start: 
1)  Review [InitializeHRVparams.m](https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox/blob/master/InitializeHRVparams.m) and optimize the parameters for your 
    data. 
2)  The toolbox does not assume any format of data except that the input 
    of the [Main_HRV_Analysis.m](https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox/blob/master/Main_HRV_Analysis.m) fucntion are a two equal length vectors: RR interval
    and time in units of seconds or the 'raw' ECG signal (physical units,mV) 
    and time. 
    Additionaly, blood pressure waveform and photoplethysmographic/pulsatile
    data can be analyzed and they should be in the standard physical units 
    (mmHg or normalized units respectively). 
3)  Results will be stored in folder called as indicated in the [InitializeHRVparams.m](https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox/blob/master/InitializeHRVparams.m)
    If the folder does not exist, it will be created.

<a name="output"></a>
## III. Guide to Output:
The following metrics are output from the HRV Toolbox:

    - t_start   : (s)  Start time of each window analyzed
    - t_end     : (s)  End time of each window analyzed

#### Time domain measures of HRV:

	- NNmean    : (ms) mean value of NN intervals
	- NNmode    : (ms) mode of NN intervals
	- NNmedian  : (ms) median value of NN intervals
	- NNskew    : skweness of NN intervals
	- NNkurt    : kurtosis of NN intervals
	- NNiqr     : interquartile range of NN intervals
	- SDNN      : (ms) Standard deviation of all NN intervals.
	- RMSSD     : (ms) The square root of the mean of the sum of the squares 
                       of differences between adjacent NN intervals.
	- pnn50     : (%) NN50 count divided by the total number of all NN intervals.
                  (Number of pairs of adjacent NN intervals differing by more than 50 ms )
	- tdflag    :   2 = not enough high SQI data in the window to process
                	(amount of data above threshold1 is greater than threshold2)
			3 = not enough data in the window 
			4 = window is missing too much data
			5 = success

#### Frequency domain measures of HRV (default using Lomb Periodogram method):

	- ulf         : (ms^2) Power in the ultra low frequency range (default < 0.003 Hz)
	- vlf         : (ms^2) Power in very low frequency range (default 0.003 <= vlf < 0.04 Hz)
	- lf          : (ms^2) Power in low frequency range (default 0.04Hz  <= lf < 0.15 Hz)
	- hf          : (ms^2) Power in high frequency range (default 0.15 <= hf < 0.4 Hz)
	- lfhf        : Ratio LF [ms^2]/HF [ms^2]
	- ttlpwr      : (ms^2) Total spectral power (approximately <0.4 Hz)
	- fdflag      : 1 = Lomb Periodogram or other method failed
	                2 = not enough high SQI data in the window to process
                	(amount of data above threshold1 is greater than threshold2)
			3 = not enough data in the window
			4 = window is missing too much data
			5 = success

#### Other HRV measures: 
    
    - PRSA - AC     : (ms) acceleration capacity
    - PRSA - DC     : (ms) deceleration capacity
    - SDANN         : (ms) Standard deviation of the average of NN intervals 
                       in all 5-minute segments of a long recording
    - SDNNI         : (ms) Mean of the standard deviation in all 5-minute 
                      segments of a long recording
		      
#### Entropy measures:
    - SampEn        : (a.u.) Sample entropy, which measures the regularity and complexity of a time series
    - ApEn          : (a.u.)Approximate entropy, which measures the regularity and complexity of a time series
    
#### Long range measures:
    
    - MSE           : First column contains the scale factors, and the second 
                      column provides the corresponding entropy values
     
    - DFA - alpha1  : Short range fractal scaling exponents (default 4<=n<16)
    - DFA - alpha2  : Long range fractal scaling exponents (default 16<=n<length/4)

#### Nonlinear HRV measures: 

    PoincarÃ© plot (PP)
     - SD1        : (ms) standard  deviation  of  projection  of  the  PP    
                    on the line perpendicular to the line of identity (y=-x)
     - SD2        : (ms) standard deviation of the projection of the PP on 
                    the line of identity (y=x)
     - SD2/SD1    : (ms) SD1/SD2 ratio

#### Heart Rate Turbulence HRT Analysis:

     - TO         : (%) turbulence onset
     - TS         : turbulence slope    

#### Detection Annotation Files 

Using [Main_HRV_Analysis.m](https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox/blob/master/Main_HRV_Analysis.m), [Analyze_ABP_PPG_Waveforms.m](https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox/blob/master/Tools/Analyze_ABP_PPG_Waveforms.m) to analyze the ECG, PPG and/or ABP the function 
will return an annotation file with the locations of detected QRS peaks or PPG/ABP onsets:

    ECG : *.jqrs (for jqrs detector)
          *.wqrs (for wqrs detector)
          *.sqrs (for sqrs detector)

    PPG : *.ppg (for PPG onset)

    ABP : *.abp (for ABP onset)

To read these files use the [read_ann.m] function included in the toolbox:

    QRS_locations = read_ann('fileName', 'jqrs')
    PPG_onsets = read_ann('fileName','ppg') 

Note that QRS locations and PPG/ABP onstets are in samples not in seconds

#### SQI Annotation Files

The SQI values are also saved as annotations files both for ECG and PPG/ABP

For ECG the SQI values are saved as a number from 0 to 100 in a file with extension:

    *.sqijw : comparison of jqrs wrt wqrs detection
    *.sqijs : comparison of jqrs wrt sqrs detection

read these files as follows

    [sqiTime,~,sqiValue] = read_ann('fileName' , 'sqijw')

For PPG and ABP two different values of SQI are seved in each annotation files
and they are related to a specific 'beat', one is a char value (E: excellent 
beat, A: acceptable beat, Q: unaceptable beat) and the other value is an integer 
in the range 0-100 given by the average of three SQI values (see PPG_SQI_buf.m)

read *.ppgsqi files as follows 

    [ppgAnn, ppgSQI, ppgSQInum] = read_ann('fileName', 'sqippg')
  

<a name="cont"></a>
## V. Contributing to this project

We are more than happy to accept contributions! If you like the project and find it useful, you can also start to improve the code or add new features yourself, it would be a great contribution to the community!

### Using the issue tracker

The issue tracker is the preferred channel for [bug reports](#bugs) but please **do not** use the issue tracker for personal support requests.


<a name="bugs"></a>
### Bug reports

A bug is a _demonstrable problem_ that is caused by the code in the repository. A good bug report is extremely important to solve the problem! Please, **check if the issue has already been reported** before opening a new issues. 

Please try to be as detailed as possible in your report. **Use the provided template!** and answer the necessary points: what is your environment? What steps will reproduce the issue? What Matlab verison and OS experience the problem? What would you expect to be the outcome? All these details will help people to fix any potential bugs.


<a name="pull-requests"></a>
### Pull requests


Follow this process if you'd like your work  - patches, improvements, new features - considered for inclusion in the
project.  

1. [Fork](http://help.github.com/fork-a-repo/) the project, clone your fork, and configure the remotes.

2. Create a new topic branch (off the main project development branch) to contain your feature, change, or fix.

3. Commit your changes in logical chunks. 

4. Locally merge (or rebase) the upstream development branch into your topic branch.

5. Push your topic branch up to your fork.

6. [Open a Pull Request](https://help.github.com/articles/using-pull-requests/) with a clear title and description.

**IMPORTANT**: By submitting a patch, you agree to allow the project owner to license your work under the same license as that used by the project.


<a name="faq"></a>
## IV. FAQ
