function [beat_matrix, avg_even, avg_odd] = CalcBeat_Matrix(fid, validleads, STLen, valid, TWARes, ecg)

% OVERVIEW This function calculates the alternans series for different points within st segment.
% If a given beat is marked as invalid it's replaced by the average even/odd beat
%
% INPUT:    MANDATORY       DESCRIPTION
%           ecg             (num_of_samples x num_of_leads)
%           fid             (1 x num_of_beats)
%           amp             (num_of_beats x num_of_leads)
%           STLen (scaler)  length of the S-Tend segment
%           valid           (1 x num_of_beats)
% OUTPUT:
%           beat_matrix     (num_of_timepoints x num_of_beats x num_of_leads): aligned ecg beats matrix
%           avg_even        (num_of_timepoints x num_of_leads) average even S-T offset beat
%           avg_odd         (num_of_timepoints x num_of_leads) average odd S-T offset beat.
%
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Shamim Nemati, 
%       editted by Ismail Sadiq on 10/26/2019. 
%	COPYRIGHT (C) 2019
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  

if(isfield(TWARes,'avg_even') & isfield(TWARes,'avg_odd'))
    avg_even = TWARes.avg_even;
    avg_odd = TWARes.avg_odd;
else
    avg_even=[];
    avg_odd=[];
end

beat_matrix = zeros(STLen,length(fid),size(ecg, 2));
for lead = 1:size(ecg, 2)
    if (validleads(lead))
        invalid_exist = ~isempty(find(valid(:,lead)==0));
        % if (invalid_exist)-> calculate the indeces of even and odd beats
        % minus the invalid beats. Then construct an average template for even and odd beats.
        % finally replace the invalid beat with the corresponding template,
        % ie. if even with an even template and if odd with an odd template
        vinds = find(valid(:,lead));
        odd = vinds(find(mod(vinds(:), 2)));  %   odd indices of valid QRS
        even = vinds(find((1 - mod(vinds(:), 2))));   %   even indices of valid QRS

        for timept = 1:STLen
            %   whenever there are invalid complexes compute odd and even
            %   averages to replace those
            
            avg_odd(timept, lead) = mean(ecg(fid(odd) + timept, lead)); % - amp(odd, lead));
            avg_even(timept, lead) = mean(ecg(fid(even) + timept, lead)); % - amp(even, lead));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for i = 1:(length(fid))

                if (valid(i)) %   beats that are followed e.g. by an extrasystole should be rejected as well
                    first = ecg(fid(i) + timept, lead);% - amp(i, lead);
                elseif mod(i, 2)
                    first = avg_odd(timept, lead);
                else
                    first = avg_even(timept, lead);
                end
                beat_matrix(timept, i, lead) = first;
            end
        end
    end
end
%


return;
