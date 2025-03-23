clear

subs = {'KS' 'LS' 'S1' 'S2' 'S3' 'S4' 'S8' 'S6' 'S10' 'S14' 'S15' 'S13' 'S16' 'S21' 'S17'};

% select which MEG measure we're using for time series correlation across sensors
DVs = {'raw' ...                                     %  1
       'lowpass_5hz' 'lowpass_3hz' 'lowpass_1hz' ... %  2 -  4
       '3.33Hz_amp' '3.33Hz_phase' ...               %  5 -  6
       'delta_amp' 'delta_phase' ...                 %  7 -  8 
       'theta_amp' 'theta_phase' ...                 %  9 - 10
       'alpha_amp' 'alpha_phase' ...                 % 11 - 12
       'beta_amp' 'beta_phase' ...                   % 13 - 14
       'low_gamma_amp' 'low_gamma_phase' ...         % 15 - 16
       'mid_gamma_amp' 'mid_gamma_phase' ...         % 17 - 18
       'high_gamma_amp' 'high_gamma_phase' ...       % 19 - 20
       ...
       '3.33Hz_log_amp' ...                          % 21
       'delta_log_amp' 'theta_log_amp' ...           % 22 - 23
       'alpha_log_amp' 'beta_log_amp' ...            % 24 - 25
       'low_gamma_log_amp' 'mid_gamma_log_amp' ...   % 26 - 27
       'high_gamma_log_amp'};                        % 28
   
% subs = {'S21'};   

phase_option = 0;
for i_sub = 1:length(subs)
    i_sub
    for DV_ind = 1 %1:4
        DV_ind
        
        for betaID = 1:5
            stim_corr(subs{i_sub}, DV_ind, phase_option, betaID);
        end    
        
    end
end