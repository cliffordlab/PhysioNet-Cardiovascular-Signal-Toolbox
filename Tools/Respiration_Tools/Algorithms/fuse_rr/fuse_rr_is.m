function [ekg_RSlinB__ELF_Cto_SFu] = fuse_rr_is(ekg_RSlinB_am_ELF_Cto, ekg_RSlinB_bw_ELF_Cto, ekg_RSlinB_fm_ELF_Cto, up)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
rel_sub_comps = up.al.sub_components.fus_mod;
% Identify relevant data
feat_mods = {'bw', 'am', 'fm'};
for mod_no = 1:3
    mod = feat_mods(mod_no);
    eval(['rel_data.' mod{1,1} ' = ekg_RSlinB_', mod{1,1}, '_ELF_Cto' ';']);
    %eval(['rel_data.' mod{1,1} ' = rrEsts.' temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no} ';']);
    %eval(['rel_data.' mod{1,1} ' = rrEsts.' temp.bef{fusable_rr_est_no}, mod{1,1}, temp.aft{fusable_rr_est_no} ';']);
end

for sub_comp_no = 1 : length(rel_sub_comps)
    % Do fusion
    temp_rr = feval(rel_sub_comps{sub_comp_no}, rel_data, up); clear rel_data
    % store this series of rrs:
    save_name = ['ekg_RSlinB__ELF_Cto_SFu'];
    eval([save_name, ' = temp_rr;']);
end

end

