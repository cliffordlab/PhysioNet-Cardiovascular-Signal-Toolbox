function [ekg_RSlinB_am_ELF_CtO, ekg_RSlinB_bw_ELF_CtO, ekg_RSlinB_fm_ELF_CtO] = estimate_rr_is(ekg_RSlinB_am_ELF, ekg_RSlinB_bw_ELF, ekg_RSlinB_fm_ELF, up)
% Estimation of respiration rate using multiple methods

subj = 1;
loaded_this_subj_respSigs = 0;
%% Make window timings if necessary
[wins] = identify_subj_wins_is(ekg_RSlinB_am_ELF, ekg_RSlinB_bw_ELF, ekg_RSlinB_fm_ELF, up);
%% Cycle through each resp signal
for respSig_no = 1:3    % 3 input signals
    rel_name = ['ekg_RSlinB_' up.al.options.FMe{respSig_no} '_ELF'];
    for option_no = 1 : length(up.al.options.estimate_rr)
        % Skip if this processing has been done previously
        if strcmp(up.al.options.estimate_rr{option_no}, 'GCE')
            %save_name = [ respSigs{respSig_no}(1:3) '_' up.al.options.estimate_rr{option_no} ];
            save_name = [ rel_name '_' up.al.options.estimate_rr{option_no} ];
        else
            %save_name = [ respSigs{respSig_no} '_' up.al.options.estimate_rr{option_no} ];
            save_name = [ rel_name '_' up.al.options.estimate_rr{option_no} ];
        end
%         savepath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.rrEsts, '.mat'];
%         exist_log = check_exists(savepath, save_name);
%         if exist_log
%             continue
%         end
        
%         % load data if it hasn't yet been loaded
%         if ~loaded_this_subj_respSigs
%             % Signals
%             load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.respSigs]);
%             % Window timings
%             load([up.paths.data_save_folder, num2str(subj), up.paths.filenames.win_timings, '.mat']);
%             loaded_this_subj_respSigs = 1;
%         end
        % Identify the relevant respSig data
        eval(['rel_data = ' rel_name ';']);
        % add current subject and resp sig for any methods that do not use a resp sig.
        rel_data.subj = subj;
        rel_data.respSig = rel_name;
        %eval(['rel_data.respSig = ' rel_name ';']); %respSigs{respSig_no};
        %% Calculate RR from this resp sig using each option for estimating RR
        if (length(rel_data.t) == 1 && isnan(rel_data.t)) || sum(isnan(rel_data.v))==length(rel_data.v)
            temp_rr.t = mean([wins.t_start(:)' ; wins.t_end(:)']); temp_rr.t = temp_rr.t(:);
            temp_rr.v = nan(length(temp_rr.t),1);
            [temp_rr.f, temp_rr.p] = deal(cell(length(temp_rr.t),1));
        else
            temp_rr = feval(up.al.options.estimate_rr{option_no}, rel_data, wins, up);
        end
        % store this series of rrs:
        eval([rel_name '_' up.al.options.estimate_rr{option_no} ' = temp_rr;']);
        
        clear temp_rr
%         %% Save RRs to file
%         save_or_append_data
    end
end

end

