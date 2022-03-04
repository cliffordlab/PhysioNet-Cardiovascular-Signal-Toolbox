function [wins] = identify_subj_wins_is(ekg_RSlinB_am_ELF, ekg_RSlinB_bw_ELF, ekg_RSlinB_fm_ELF, up)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%% Identify start and end times of each respiratory signal
respSigs.ekg_RSlinB_am_ELF = ekg_RSlinB_am_ELF; respSigs.ekg_RSlinB_bw_ELF = ekg_RSlinB_bw_ELF; respSigs.ekg_RSlinB_fm_ELF = ekg_RSlinB_fm_ELF;
respSig_names = fieldnames(respSigs);
[timings.start, timings.end] = deal(nan(length(respSig_names),1));
for sig_no = 1 : length(respSig_names)
    eval(['rel_data = respSigs.' respSig_names{sig_no} '.t;']);
    timings.start(sig_no) = rel_data(1);
    timings.end(sig_no) = rel_data(end);
end

latest_start = max(timings.start);
earliest_end = min(timings.end);

%% Define windows
duration_of_one_win = up.paramSet.winLeng;
no_of_secs_bet_con_wins = up.paramSet.winStep;
gap_between_win_starts = duration_of_one_win - no_of_secs_bet_con_wins;
wins.t_start = latest_start : gap_between_win_starts : (earliest_end - duration_of_one_win);
wins.t_end = wins.t_start + duration_of_one_win;

end

