This code estimate the respiration rate from the ECG using multiple methods

The main demo function is rr_is.m. 

Currently it reads a sample signal using rdmat from the EHF_is.m function (This line can be commented and the ECG signal (ecg) and sampling frequency (Fs) passed in as parameters to the EHF_is function).


The original version only worked with a sample frequency of = 500 Hz. In case Fs is not equal to 500 Hz, the ecg signal is re-sampled to this frequency.


Then the algorithms runs Respirtion Rates (RR) algorithms on ECG and PPG signals using each possible combination of options, as specified in "setup_universal_params.m".

               RRest('vortal_rest')

	Inputs:
		data            data files should be stored in the specified format
                       in the directory specified in
                       "setup_universal_params.m". Data can be downloaded
                       using the "
       period          this string specifies the dataset to be analysed.
                       Only the 'mimic' dataset has been used with this
                       version of the toolbox.

	Outputs:
       for each subject, N, the following files are made:
           N_int_respSigs      intermediate respiratory signals,
           N_respSigs          final respiratory signals
           N_rrEsts            RR estimates
           N_rrRef             Reference RR values
           N_sqi               Signal Quality Index values
       for the entire dataset, the following files are made:
           alg_names           Names of RR algorithms tested
           win_data            Data for every algorithm tested, every
                               window, and every subject.

   Context:    This is the main file used to run the algorithms. It calls
               lots of other functions, contained in separate files.
           
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





