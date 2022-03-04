% Main function to estimate the respiration rate from ECG signals

% Set up initial parameters
addpath('./Algorithms')
up = setup_universal_params('RRSYNTH');

% Read the demo file and pre-process the data
[s_filt] = EHF_is();

% QRS wave detector
[ekg_RDtGC_EHF] = RDt_is(s_filt, up);

% Fiducial point detector
[ekg_FPt_RDtGC_EHF] = FPt_is(ekg_RDtGC_EHF,s_filt);

% Measures ECG features from peak and trough values such as QRS slopes and R-peaks
[feat_data] = FMe_is(ekg_FPt_RDtGC_EHF,s_filt,up);

%Resamples feat-based respiratory signals at a regular sampling rate.
[ekg_RSlinB_am, ekg_RSlinB_bw, ekg_RSlinB_fm] = RS_is(feat_data,up);

%ELF eliminates very low frequencies from resampled respiratory
[ekg_RSlinB_am_ELF, ekg_RSlinB_bw_ELF, ekg_RSlinB_fm_ELF] = ELF_is(ekg_RSlinB_am, ekg_RSlinB_bw, ekg_RSlinB_fm, up);

% Estimation of respiration rate (RR) using  using each possible set
% of RR estimation options and multiple methods:
% 	ekg_RSlinB_am_ELF_Cto.v respiration cycles per minute and ekg_RSlinB_am_ELF_Cto.t is the time index of the mid point. Estimation from amplitude
% 	ekg_RSlinB_bw_ELF_Cto.v change from baseline
%	ekg_RSlinB_fm_ELF_Cto.v frecuency modulation of RR intervals
[ekg_RSlinB_am_ELF_Cto, ekg_RSlinB_bw_ELF_Cto, ekg_RSlinB_fm_ELF_Cto] = estimate_rr_is(ekg_RSlinB_am_ELF, ekg_RSlinB_bw_ELF, ekg_RSlinB_fm_ELF, up);

% Fuse estimate RR if the am, bw, fm are similar it computes an average. 
[ekg_RSlinB__ELF_Cto_SFu] = fuse_rr_is(ekg_RSlinB_am_ELF_Cto, ekg_RSlinB_bw_ELF_Cto, ekg_RSlinB_fm_ELF_Cto, up);  


% Plot results
% figure(1); plot(s_filt.t, s_filt.v); hold on;
% stem(ekg_RDtGC_EHF.p.t, ekg_RDtGC_EHF.p.v); 
% stem(ekg_RDtGC_EHF.tr.t, ekg_RDtGC_EHF.tr.v); hold off
% 
% figure(1); plot(s_filt.t, s_filt.v); hold on;
% stem(ekg_FPt_RDtGC_EHF.det_p.t, ekg_FPt_RDtGC_EHF.det_p.v); 
% stem(ekg_FPt_RDtGC_EHF.det_tr.t, ekg_FPt_RDtGC_EHF.det_tr.v); hold off
% 
% figure(1); plot(feat_data.am.t, feat_data.am.v);
% figure(1); plot(ekg_RSlinB_am.t, ekg_RSlinB_am.v);
% plot(ekg_RSlinB_bw.t, ekg_RSlinB_bw.v);
% plot(ekg_RSlinB_fm.t, ekg_RSlinB_fm.v);
