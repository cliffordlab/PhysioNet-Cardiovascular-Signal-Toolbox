% Compare the alpha scaling factor obtained by VOSIM function and Physionet
% one with respect to a given one.

% We test the two algorithm with realisation of noise processes with 
% different values of Beta 


InputDataType = 'NSR' ; % syntheticaData or NSR

switch InputDataType
    
    case 'syntheticaData'
        Beta = -2:6/24:4;
        SerieLength = 4096;

        for idx = 1: length(Beta)
            y = fftpowerlaw(Beta(idx), SerieLength);
            % Alpha using VOSIM
            alpha_VOSIM(idx) = dfaScalingExponent(y);
            % Compute Physionet DFA
            [lns,lF]=dfa(y);
            % Compute the alpha value based on lns and lF
            A = ones(length(lns),2);
            A(:,2) = lns;
            a = pinv(A)*lF;
            alpha_Physionet(idx) = a(2);
        end
        
    case 'NSR' % Normal Sinus Rhytm Data
        FileIDs = 1:54;
        for idx = 1: length(FileIDs)  % Number of signal in the Dataset
            if FileIDs(idx)<10
                thisFile = strcat('nsr2db/nsr00', num2str(FileIDs(idx)));
            else
                thisFile = strcat('nsr2db/nsr0', num2str(FileIDs(idx)));
            end
            
            [anntime,~] = rdann(thisFile,'ecg',[],[],[],'N');
            RRserie = diff(anntime);
            alpha_VOSIM(idx) = dfaScalingExponent(RRserie);
            
            % Physionet DFA
            [lns,lF]=dfa(RRserie);
            % Compute the alpha value based on lns and lF
            A = ones(length(lns),2);
            A(:,2) = lns;
            a = pinv(A)*lF;
            alpha_Physionet(idx) = a(2);
        end
 
end


    