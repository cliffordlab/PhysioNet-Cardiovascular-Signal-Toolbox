# PhysioNet Cardiovascular Signal Toolbox: lombEC57.m
# August 2021, Advanced Algorithm Research Center, Philips Healthcare
# Zhang, Yu-He <yu-he.zhang@philips.com>

## Background
American National Standard ANSI/AAMI EC57: 2012 Sections 4.3.3 and A.4.3.3 provided guidelines to evaluate the performance of heart rate variability (HRV) algorithm. Specifically, Section A.4.3.3 described method to generate test patterns (TPs) to evaluate the HRV frequency domain (FD) algorithm with expected spectrum values in VLF, LF, and HF, in addition to the time domain (TD) stats, as shown in Table A.10. The frequency domain criteria is most significant since method of objective evaluation of FD performance with quantitative precision is rarely available. This version of the Lomb-Scargle algorithm in lombEC57.m was implemented to meet the aforementioned criteria. Compared to the original version in CalcLomb.m, lombEC57.m is much faster especially when running the TPs where the size of frequency vector is big, in addition to meeting the EC57 criteria.

## Using lombEC57.m
Replace CalcLomb.m with lombEC57.m in case 'lomb' in EvalFrequencyDomainHRVstats.m.

## EC57 Test Patterns (TPs)
The four test patterns (TPs) according to the method in EC57 Section 4.3.3 are available in ../Physionet_EC57TPs.

## Settings specific to EC57 TP testing
Since extraordinary long RR intervals are generated in TPs to evaluate VLF and LF, the following settings are required to survive the preprocessing stage of the Toolbox,

    HRVparams.preprocess.upperphysiolim = 60/15;
    HRVparams.preprocess.gaplimit = 4;

In addition, the window size of the HRV stats is the entire duration of the TP signal, instead of the default 5 minute window.

## Expected Results
As shown in EC57 Table A.10, the HRV algorithm is expected to generate 612.5 ms2 power in HF with TP2 (35 ms deviation, 0.25 Hz), 2450 ms2 power in LF with TP3 (70 ms deviation, 0.10 Hz), 39200 ms2 power in VLF with TP4 (280 ms deviation, 0.033333 Hz), and zero power in VLF, LF, and HF with TP5 (140 ms deviation, 0.000278 Hz).