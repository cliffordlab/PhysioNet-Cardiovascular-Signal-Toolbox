function [feat_data] = FMe_is(ekg_FPt_RDtGC_EHF,ekg_EHF,up)
% Measures ECG features from peak and trough values such as QRS slopes and R-peaks

sig_type = 'ekg';
curr_sig = 'ekg';

%% Skip if this processing has been done previously
% log_int_respSig = 1;
% subj = 1;
% iden_resp_sig_file_ending
% savepath = [up.paths.data_save_folder, num2str(subj), ending];
% filecontents = whos('-file', savepath);
var_names = {'ekg_EHF', 'ekg_FPt_RDtGC_EHf', 'ekg_RDtGC_EHF'};%extractfield(filecontents, 'name');
rel_log = zeros(size(var_names));
for s = 1 : length(var_names)
    if strfind(var_names{s}, [curr_sig, up.paths.filenames.fid_pts])
        rel_log(s) = 1;
    end
end
rel_var_names = var_names(logical(rel_log));
for rel_var_name_no = 1 : length(rel_var_names)
    for current_opt_no = 1 : length(up.al.options.FMe)
        % skip out ones that aren't relevant to the ppg
        %     if strcmp(sig_type, 'ppg') && ( strcmp(up.al.options.FMe{current_opt_no}, 'qrsW') ...
        %             || strcmp(up.al.options.FMe{current_opt_no}, 'qrsA') ...
        %             || strcmp(up.al.options.FMe{current_opt_no}, 'qrS') ...
        %             || strcmp(up.al.options.FMe{current_opt_no}, 'rsS') ...
        %             || strcmp(up.al.options.FMe{current_opt_no}, 'Rang'))
        %         continue
        %     end
        
        % skip out ones that aren't relevant to the ekg
        if strcmp(sig_type, 'ekg') && strcmp(up.al.options.FMe{current_opt_no}, 'pulW')
            continue
        end
        
            temp = strfind(rel_var_names{rel_var_name_no},'_');
            start_el = temp(1); clear temp
        eval(['save_name = ''' curr_sig, up.paths.filenames.feat_meas, up.al.options.FMe{current_opt_no}, rel_var_names{rel_var_name_no}(start_el:end) ''';']); clear start_el
        %     exist_log = check_exists(savepath, save_name);
        %     if exist_log
        %         continue
        %     end
        
        %% Load relevant data
        
        % only load signal data if not loaded already
        if ~exist('sig_data', 'var')
            %         loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
            %         rel_name = [curr_sig, up.paths.filenames.elim_vhf];
            %         load(loadpath, rel_name);
            %         eval(['sig_data = ' rel_name ';']);
            sig_data = ekg_EHF;
        end
        
        % only load fiducial point data if not loaded already
        if ~exist('rel_data', 'var')
            %loadpath = [up.paths.data_save_folder, num2str(subj), up.paths.filenames.int_respSigs];
            %rel_name = rel_var_names{rel_var_name_no};
            %load(loadpath, rel_name);
            %eval(['rel_data = ' rel_name ';']);
            rel_data = ekg_FPt_RDtGC_EHF;
            %% Settings
            % choose which peaks and onsets:
            rel_data.tr = rel_data.tr_min;
            rel_data.p = rel_data.p_max;
            
            if isempty(rel_data.tr.t) || isempty(rel_data.p.t)
                [peaks.t, peaks.v, onsets.t, onsets.v] = deal([]);
            else
                
                % Want an onset followed by a peak.
                if rel_data.tr.t(1) > rel_data.p.t(1)
                    rel_data.p.t = rel_data.p.t(2:end);
                    rel_data.p.v = rel_data.p.v(2:end);
                end
                % Want the same number of peaks and onsets
                diff_in_length = length(rel_data.p.t) - length(rel_data.tr.t);
                if diff_in_length > 0
                    rel_data.p.t = rel_data.p.t(1:(end-diff_in_length));
                    rel_data.p.v = rel_data.p.v(1:(end-diff_in_length));
                elseif diff_in_length < 0
                    rel_data.tr.t = rel_data.tr.t(1:(end-diff_in_length));
                    rel_data.tr.v = rel_data.tr.v(1:(end-diff_in_length));
                end
                % find onsets and peaks
                onsets.t = rel_data.tr.t; onsets.t = onsets.t(:);
                onsets.v = rel_data.tr.v; onsets.v = onsets.v(:);
                peaks.t = rel_data.p.t; peaks.t = peaks.t(:);
                peaks.v = rel_data.p.v; peaks.v = peaks.v(:);
                
                % exclude ectopics
                % following: Mateo, J. & Laguna, P., 2003. Analysis of heart rate variability in the presence of ectopic beats using the heart timing signal. IEEE transactions on bio-medical engineering, 50(3), pp.334–43. Available at: http://www.ncbi.nlm.nih.gov/pubmed/12669990.
                if 1%~isempty(strfind(up.paths.data_load_filename, 'vortal'))
                    tk_neg1 = peaks.t(1:(end-2));
                    tk = peaks.t(2:(end-1));
                    tk_pos1 = peaks.t(3:end);
                    r = 2*abs( (tk_neg1 - (2*tk) + tk_pos1)./ ...
                        ( (tk_neg1-tk).*(tk_neg1 - tk_pos1).*(tk-tk_pos1) ) );
                    thresh = min([4.3*std(r), 0.5]);
                    
                    %%%%%%%%%%%%%%%%%%%%%%
                    % additional rule inserted by PC:
                    % thresh = 0.5;   % so that artificial data with a very low variability doesn't trigger.
                    %%%%%%%%%%%%%%%%%%%%%%
                    temp = [0;r;0]; temp = logical(temp>thresh);
                    
                    tk_neg1 = onsets.t(1:(end-2));
                    tk = onsets.t(2:(end-1));
                    tk_pos1 = onsets.t(3:end);
                    r = 2*abs( (tk_neg1 - (2*tk) + tk_pos1)./ ...
                        ( (tk_neg1-tk).*(tk_neg1 - tk_pos1).*(tk-tk_pos1) ) );
                    thresh = min([4.3*std(r), 0.5]);
                    %%%%%%%%%%%%%%%%%%%%%%
                    % additional rule inserted by PC:
                    %thresh = 0.5;   % so that artificial data with a very low variability doesn't trigger.
                    %%%%%%%%%%%%%%%%%%%%%%
                    temp2 = [0;r;0]; temp2 = logical(temp2>thresh);
                    
                    peaks.t(temp | temp2) = nan; peaks.v(temp | temp2) = nan;
                    onsets.t(temp | temp2) = nan; onsets.v(temp | temp2) = nan;
                    %                     % to check:
                    %                     old_peaks  = rel_data.p; old_onsets = rel_data.tr;
                    %                     plot(diff(old_peaks.t), 'b'), hold on, plot(diff(peaks.t), 'r')
                    %                     close all
                    %                     plot(diff(old_onsets.t)), hold on, plot(diff(onsets.t))
                    %                     close all
                end
                clear temp temp2 r tk tk_neg1 tk_pos1 thresh
            end
        end
        
        %% Measure Features
        sig_data.wave_type = curr_sig;   % for PCA
        feat_data = feval(['calc_' up.al.options.FMe{current_opt_no}], peaks, onsets, sig_data.fs, sig_data, up);
        feat_data.timings = rel_data.timings;
        
        %% Save results
        %ekg_FMeam_FPt_RDtGC_EHF = feat_data;
        save_name = ['feat_data_' up.al.options.FMe{current_opt_no}];
        eval([save_name ' = feat_data;']);
        %save_or_append_data
    end

end

feat_data = [];
feat_data.am = feat_data_am;
feat_data.bw = feat_data_bw;
feat_data.fm = feat_data_fm;

end

function feat_data = calc_am(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

% Find am
am.t = mean([onsets.t, peaks.t], 2);
am.v = [peaks.v - onsets.v];

% Normalise
feat_data.v = am.v./nanmean(am.v);
feat_data.t = am.t;
feat_data.fs = fs;

end

function feat_data = calc_bw(peaks, onsets, fs, sig_data, up)

% eliminate any nans (which represent ectopics which have been removed)
peaks.t = peaks.t(~isnan(peaks.t));
peaks.v = peaks.v(~isnan(peaks.v));
onsets.t = onsets.t(~isnan(onsets.t));
onsets.v = onsets.v(~isnan(onsets.v));

% Find bw
bw.v = mean([onsets.v, peaks.v], 2);
bw.t = mean([onsets.t, peaks.t], 2);

% Find am
am.t = mean([onsets.t, peaks.t], 2);
am.v = [peaks.v - onsets.v];

% Normalise
feat_data.v = bw.v./nanmean(am.v);
feat_data.t = bw.t;
feat_data.fs = fs;

end

function feat_data = calc_fm(peaks, onsets, fs, sig_data, up)

% find fm
fm.v = [peaks.t(2:end) - peaks.t(1:(end-1))]/fs;
fm.t = mean([peaks.t(2:end), peaks.t(1:(end-1))], 2);

% eliminate any nans (which represent ectopics which have been removed)
fm.t = fm.t(~isnan(fm.t));
fm.v = fm.v(~isnan(fm.v));

% Normalise
feat_data.v = fm.v./nanmean(fm.v);
feat_data.t = fm.t;
feat_data.fs = fs;

end

