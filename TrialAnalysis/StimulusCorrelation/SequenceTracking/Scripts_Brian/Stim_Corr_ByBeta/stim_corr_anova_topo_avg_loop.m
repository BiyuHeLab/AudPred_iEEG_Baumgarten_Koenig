clear

% subs = {'KS' 'LS' 'S1' 'S2' 'S3' 'S4' 'S8' 'S6' 'S10' 'S14' 'S15'};

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

pthresh = [.05 .005 .0005 1];
% pthresh = 1;

for DV_ind = 1 %1:4 %6 : 2 : 20
    
    for betaID = 1:5
    for i_p = 1:length(pthresh)
        stim_corr_anova_topo_avg('F', DVs{DV_ind}, pthresh(i_p), saveplot, betaID);
    end
    end
%     stim_corr_anova_topo_avg('p', DVs{DV_ind}, pthresh(1), saveplot);
end
