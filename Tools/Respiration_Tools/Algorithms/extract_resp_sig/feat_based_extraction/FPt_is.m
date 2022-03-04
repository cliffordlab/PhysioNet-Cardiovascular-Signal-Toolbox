function [ekg_FPt_RDtGC_EHF] = FPt_is(rel_name1,rel_name2)
%FPT detects Fiducial Points from PPG peak and trough annotations
% as specified in PC's literature review.

eval(['rel_data.beats = ' getVarName(rel_name1) ';']);
eval(['rel_data.s = ' getVarName(rel_name2) ';']);
eval(['rel_data.fs = ' getVarName(rel_name1) '.fs;']);
eval(['rel_data.timings = ' getVarName(rel_name1) '.timings;']);

%% ECG Peaks
% Used to find peaks as max between detected onsets, but it now
% appears that this doesn't work if it detects the onsets as a
% mixture of Q and S waves.

temp.p_max = rel_data.beats.p;

%% ECG Troughs
% Find troughs as min within search range 0.1s before peaks
temp.tr_min.t = nan(length(temp.p_max.t),1);
temp.tr_min.v = nan(length(temp.p_max.t),1);
thresh = 0.1;
for beat_no = 1 : (length(temp.p_max.t))
    rel_range = find(rel_data.s.t >= (temp.p_max.t(beat_no) - thresh) & rel_data.s.t < temp.p_max.t(beat_no));
    % used to be:     rel_range = find(rel_data.s.t >= (temp.p_max.t(beat_no) - thresh) & rel_data.s.t < (temp.p_max.t(beat_no) + thresh));
    [~, rel_el] = min(rel_data.s.v(rel_range));
    if ~isempty(rel_el)   % it is empty if the peak is at the very first element of the signal
        temp.tr_min.t(beat_no) = rel_data.s.t(rel_range(rel_el));
        temp.tr_min.v(beat_no) = rel_data.s.v(rel_range(rel_el));
    end
end
% get rid of any nans (arise if the peak is at the very first element of the signal)
bad_els = isnan(temp.tr_min.t);
temp.tr_min.t = temp.tr_min.t(~bad_els);
temp.tr_min.v = temp.tr_min.v(~bad_els);

% very ocassionally it picks out the same trough or peak twice (if two consecutive peaks are ridiculously close together so share some of the same search range)
[~, rel_els, ~] = unique(temp.tr_min.t);
temp.tr_min.t = temp.tr_min.t(rel_els);
temp.tr_min.v = temp.tr_min.v(rel_els);

% Carry forward detected peaks and onsets
temp.det_p.t = rel_data.beats.p.t;
temp.det_p.v = rel_data.beats.p.v;
temp.det_tr.t = rel_data.beats.tr.t;
temp.det_tr.v = rel_data.beats.tr.v;

% carry forward fs and timings
temp.fs = rel_data.fs;
temp.timings = rel_data.timings;

ekg_FPt_RDtGC_EHF = temp;

end

function out = getVarName(var)
    out = inputname(1);
end
