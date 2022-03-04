# ECG_Sleep_Staging_Prediction

## Dependencies
- Physionet Cardiovascular Signal Toolbox (NOTE, an initialization struct with the patameters used for the specific project is provided. If not the same the default parameters 
have been used and can be generated using InitializeHRVparams([]) )

- WFDB library (for signal conversion)

## How to run the code

- Predict sleep staging from ECG using pre-trained model

  -- sleep_stage = ECG_sleep_staging_prediction(ECGdata,fs,class);

## DEMO

- prediction demo on slpdb: demo_ECG_sleep_staging_slpdb.m % predict sleep stages on slpdb from physionet.org, data were converted to Matlab format by wfdb2mat function before analysis.
