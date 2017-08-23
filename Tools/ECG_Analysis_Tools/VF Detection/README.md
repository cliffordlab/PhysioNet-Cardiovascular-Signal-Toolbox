# VF_detector

Current release: version 0.0.1

Last updated June 27, 2013

This program slides a 5 second non-overlapping window across an ECG and classifies a segment of an ECG as either ventricular fibrillation (and/or ventricular tachycardia) or absence of evidence for such a rhythm (1 or zero output respectively). Two features, called leakage and count2, are calculated from ECG waveform and a support vector machine (LIBSVM 3.11, http://www.csie.ntu.edu.tw/~cjlin/libsvm/) is used to classify each segment. The SVM was trained and tested to perform classification on expert-verified events in three annotated public domain ECG databases (the American Heart Association Database (https://www.ecri.org/components/Pages/AHA_ECG_DVD.aspx ), the Creighton University Ventricular Tachyarrhythmia Database, and the MIT-BIH Malignant Ventricular Arrhythmia Database (available from www.physionet.org)).

This program uses the WFDB Toolbox for MATLAB (http://physionet.org/physiotools/matlab/wfdb-app-matlab/) to read WFDB format data and the LIBSVM 3.11 toolbox (included, but also available from http://www.csie.ntu.edu.tw/~cjlin/libsvm/)

To run the program, call VF_Classification('recordname') in Matlab, e.g. 

VF_output = VF_Classification('cu01');

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Please cite this publication when referencing this program:

Q Li, C Rajagopalan, GD Clifford, Ventricular fibrillation and tachycardia classification using a machine learning approach, IEEE Transactions on Biomedical Engineering, 61 (6), 1607-1613, 2014.
