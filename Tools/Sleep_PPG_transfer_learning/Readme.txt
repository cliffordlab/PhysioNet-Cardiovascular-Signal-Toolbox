# Sleep_PPG_transfer_learning
Sleep staging code from photoplethysmogram (PPG) signals, with additional use of actigraphy features when available

If you are using this software, please cite:
```
Li Q, Li Q, Cakmak AS, Da Poian G, Bliwise DL, Vaccarino V, Shah AJ, Clifford GD. 
Transfer learning from ECG to PPG for improved sleep staging from wrist-worn wearables. 
Physiol Meas. 2021 May 13;42(4):10.1088/1361-6579/abf1b0.
```     

Dependencies:
- Cardiovascular Signal Processing Toolbox : 
- VLFeat open source library : https://github.com/vlfeat/vlfeat
- LIBSVM : https://github.com/cjlin1/libsvm


## 1. Features Extraction

The first step consists in PPG onste detetction for generating the IBI tachogram,
 from which the **RSA**(Respiratory sinus arrhythmia, heart rate variability in synchrony with respiration)
 and **PDR** (PPG-Derived Respiration) time series are derived. 

List of PPG features
```
- CRC (Cross Spectral Coherence)
- AC and DC (Acceleration and Decelaration Capacity) 
- SampEn (Sample Entropy)
- SDNN (Standard Deviation of Normal-to-Normal intervals)
- LF/HF 
- SQI (Signal Quality Index)
- Energy Ratio LF and HF
- Ratio of the sum of two maximal choerent cross-power peaks in LF (0.01-0.1 Hz) and the sum of two maximal choerent cross-power peaks in HF (0.1-0.4 Hz) 
```
## 2. Code for Prediction given a trained ECG model and a transfer learning model

do_PPG_sleep_staging_prediction(ppg,HRVparams,class,trainingset,SavePath,SigName) 

## 3. Demo for PPG Sleep staging

do_PPG_sleep_staging_prediction_DEMO
