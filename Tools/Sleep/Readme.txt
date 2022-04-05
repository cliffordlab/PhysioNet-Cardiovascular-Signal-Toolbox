# ECG_Sleep_Staging_Prediction

If you are using this software, please cite:
```
Li Q, Li Q, Liu C, Shashikumar SP, Nemati S, Clifford GD. 
"Deep learning in the cross-time-frequency domain for sleep staging from a single lead electrocardiogram". 
Physiol. Meas. 2018 Dec 21;39(12):124005

Li Q, Li Q, Cakmak AS, Da Poian G, Bliwise DL, Vaccarino V, Shah AJ, Clifford GD. 
Transfer learning from ECG to PPG for improved sleep staging from wrist-worn wearables. 
Physiol Meas. 2021 May 13;42(4):10.1088/1361-6579/abf1b0.
```     


## Dependencies
- Physionet Cardiovascular Signal Toolbox (NOTE, an initialization struct with the patameters used for the specific project is provided. If not the same the default parameters 
have been used and can be generated using InitializeHRVparams([]) )

- WFDB library (for signal conversion)

## How to run the code

- Predict sleep staging from ECG using pre-trained model

  -- sleep_stage = ECG_sleep_staging_prediction(ECGdata,fs,class);

## DEMO

- prediction demo on slpdb: demo_ECG_sleep_staging_slpdb.m % predict sleep stages on slpdb from physionet.org, data were converted to Matlab format by wfdb2mat function before analysis.
