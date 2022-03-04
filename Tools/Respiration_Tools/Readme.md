This code estimate the respiration rate from the ECG using multiple methods

The main demo function is rr_is.m. 

Currently it reads a sample signal using rdmat from the EHF_is.m function (This line can be commented and the ECG signal (ecg) and sampling frequency (Fs) passed in as parameters to the EHF_is function).


The original version only worked with a sample frequency of = 500 Hz. In case Fs is not equal to 500 Hz, the ecg signal is re-sampled to this frequency.

Then the algorithms runs Respirtion Rates (RR) algorithms on ECG and PPG signals using each possible combination of options, as specified in "setup_universal_params.m".

           
   Further Information:
       This version of the RRest is provided to facilitate reproduction of
       the analysis performed in:
           Charlton P.H. et al. Extraction of respiratory signals from the 
           electrocardiogram and photoplethysmogram: technical and physiological
           determinants, Physiological Measurement, 38(5), 2017
           DOI: https://doi.org/10.1088/1361-6579/aa670e
       Further information on this study can be obtained at:
           http://peterhcharlton.github.io/RRest/factors_assessment.html
       In addition, further information on RRest, including future
       versions, can be obtained at:
           http://peterhcharlton.github.io/RRest

   Comments, Questions, Criticisms, Feedback, Contributions:
       See: http://peterhcharlton.github.io/RRest/contributions.html

   Version:
       v.3 - published on 4th May 2017 by Peter Charlton

   Licence:
       please see the accompanying file named "LICENSE"





