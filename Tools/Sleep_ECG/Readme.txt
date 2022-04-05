# ECG_Sleep_Staging_Prediction_by_CNN_SVM_model

If you are using this software, please cite:
```
Li Q, Li Q, Liu C, Shashikumar SP, Nemati S, Clifford GD. 
"Deep learning in the cross-time-frequency domain for sleep staging from a single lead electrocardiogram". 
Physiol. Meas. 2018 Dec 21;39(12):124005
```     

## Dependencies
- Physionet Cardiovascular Signal Toolbox (NOTE, an initialization struct with the patameters used for the specific project is provided. If not the same the default parameters 
have been used and can be generated using InitializeHRVparams([]) )

- WFDB library (for signal conversion)

## How to run the code

- Predict sleep staging from ECG using pre-trained CNN and SVM model

  -- predict_label = do_sleep_staging_prediction(ECG,fs,class);

## DEMO

- prediction demo on slpdb: do_sleep_staging_prediction_DEMO % predict sleep stages on slpdb from physionet.org, data were converted to Matlab format by wfdb2mat function before analysis.
