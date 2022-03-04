function [ekg_RSlinB_am_ELF, ekg_RSlinB_bw_ELF, ekg_RSlinB_fm_ELF] = ELF_is(ekg_RSlinB_am, ekg_RSlinB_bw, ekg_RSlinB_fm, up)
%ELF eliminates very low frequencies from resampled respiratory


for rel_var_name_no = 1 : length(up.al.options.FMe)
%     temp = strfind(rel_var_names{rel_var_name_no},'_');
%     start_el = temp(1); clear temp
%     eval(['save_name = ''' curr_sig, up.paths.filenames.elim_vlf2, rel_var_names{rel_var_name_no}(start_el:end) ''';']); clear start_el
%     exist_log = check_exists(savepath, save_name);
%     if exist_log
%         continue
%     end
    
    %% Load relevant data
    rel_name = ['ekg_RSlinB_' up.al.options.FMe{rel_var_name_no} ];
    %loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
    %load(loadpath, rel_name);
    eval(['old_data = ' rel_name ';']);
    
    %% Eliminate VLFs
    data.hpfilt.t = old_data.t;
    try
        data.hpfilt.v = elim_vlfs(old_data, up);
    catch
        % if there aren't enough points to use the filter, simply carry forward the previous data
        data.hpfilt.v = old_data.v;
    end
    data.hpfilt.fs = old_data.fs;
    %eval(['save_name = ' rel_name '_ELF']);
    save_name = [rel_name '_ELF'];
    eval([save_name ' = data.hpfilt;']);
    
    %% Save processed data
    %save_or_append_data
end

end

