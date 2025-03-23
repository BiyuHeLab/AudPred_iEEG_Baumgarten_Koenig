clear

subs = {'KS' 'LS' 'S1' 'S2' 'S3' 'S4' 'S8' 'S6' 'S10' 'S14' 'S15' 'S13' 'S16' 'S21' 'S17'};

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

saveplot = 1;

% subs = {'S21'};

for i_sub = 1:length(subs)
    for DV_ind = 1 %6 : 2 : 20
        for phase_option = 0 %1:1 %1:2
%             stim_corr_anova_topo(subs{i_sub}, DV_ind, saveplot, phase_option);

            for betaID = 1:5
                stim_corr_anova_topo(subs{i_sub}, DV_ind, saveplot, phase_option, betaID);
            end
        end
    end
end