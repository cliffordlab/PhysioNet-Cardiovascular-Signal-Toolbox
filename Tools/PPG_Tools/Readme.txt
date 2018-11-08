PPGSQI

ToDo:

Run call.m to test the algorithm. 

Note: Make sure you can access the MIMIC database when run the example (http://physionet.org/physiobank/database/mimicdb/)
It will read the 039/039 waveform data and 039.wpleth or 039.ple annotation
and create the PPGSQI annotation named 039.wplesqi
Make sure you have 'The WFDB Toolbox for Matlab' installed
see: http://www.physionet.org/physiotools/matlab/wfdb-swig-matlab/

Algorithm description:

    * Build Template:
          o Select 30s PPG data, do cross-correlation (xcorr) and let the length of the template (L) equals to the length between the two main peaks in the cross correlation sequence
          o Average the 30s beats, each beat beginning at the fiducial mark and ending at the length of the template (length L), get Template1 (T1). 
          o With T1, calculate the correlation coefficients (corrcoef) with each beat (length L) in the 30s window, C, drop any with C<0.8, average, get Template2 (T2). If more than half were dropped, un-believe T2.
          o If T2 valid, using T2; else but T1 valid, using T1; else, template invalid.

    * Get the PPG SQI based on the template matching and clipping detection:

          o Template matching:
                + Direct matching: select each beat (length L) and calculate C with template, -> SQI1
                + Linear resampling: interpolate each beat (length Li) to length L, calculate C, -> SQI2
                + Dynamic Time Warping:
                      # Transform template and each beat (length Li) to a sequence of lines using piecewise linear approximation(PLA)
                      # Acquire distance matrix based on the absolute difference between the slopes of lines
                      # DTW algorithm to find the minimal distance to wrap the beat  (length Li) to template (length L)
                      # Calculate C, -> SQI3
          o Clipping detection: detect the percentage (P) of saturation to maximum and minimum of each beat -> SQI4 = 1-P
          o Combine SQI1-SQI4, make decision

    * PPG SQI type: E: excellent beat
                    A: acceptable beat
                    Q: unacceptable beat

* You can use wabp_pleth to detect PPG beat and crate annotation file 039.wpleth by:
    ./wabp_pleth -r 039/039 -t 1:0:0


Or: run call_buf.m to test the algorithm without read and write WFDB annotation file.
