function stim_corr(sub, DV_ind, phase_option, input_betaID)

if ~exist('phase_option','var') || isempty(phase_option)
    phase_option = 0;
end

% phase_option values
% 0 --> do circular correlation on phase
% 1 --> unwrap phase, detrend, and do linear correlation
% 2 --> unwrap phase, subtract expected phase, and do linear correlation


%% initialize paths for fieldtrip and personal scripts

addpath('/data/gogodisk2/brian/analysis/');
brian_ft_path


%% initialize paths etc for analysis of this subject

% sub = 'KS';
addpath('/data/gogodisk2/brian/analysis/sfa_expt2_v2/');
subject_info_sfa_expt2_v2


%%

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


DV  = DVs{DV_ind};


isPhase = 0;
if length(DV) >= 5  &&  strcmp( DV(end-4:end), 'phase' )
    isPhase = 1;
end

if phase_option > 0 && isPhase == 0
    error('phase option selected for non-phase MEG DV')
end


isAmp = 0;
if length(DV) >= 3 && strcmp( DV(end-2:end), 'amp')
    isAmp = 1;
end

isLogAmp = 0;
if length(DV) >= 7 && strcmp( DV(end-6:end), 'log_amp')
    isAmp    = 0;
    isLogAmp = 1;
end




% define input paths and files for MEG
switch DV
    case {'raw' 'lowpass_5hz' 'lowpass_3hz' 'lowpass_1hz'}
        path_in_MEG  = si.path_ica;
        loadfile_MEG = [path_in_MEG  'data_clean.mat'];

% %     case {'3.33Hz_amp' '3.33Hz_phase'  '3.33Hz_log_amp'}
% %         path_in_MEG  = [path_sub 'analysis/hilbert/'];
% %         loadfile_MEG = [path_in_MEG  'hilbert_3.33Hz.mat'];        
% %         
% %     case {'delta_amp' 'delta_phase'  'delta_log_amp'}
% %         path_in_MEG  = [path_sub 'analysis/hilbert/'];
% %         loadfile_MEG = [path_in_MEG  'hilbert_delta.mat'];
% %         
% %     case {'theta_amp' 'theta_phase'  'theta_log_amp'}
% %         path_in_MEG  = [path_sub 'analysis/hilbert/'];
% %         loadfile_MEG = [path_in_MEG  'hilbert_theta.mat'];
% %         
% %     case {'alpha_amp' 'alpha_phase'  'alpha_log_amp'}
% %         path_in_MEG  = [path_sub 'analysis/hilbert/'];
% %         loadfile_MEG = [path_in_MEG  'hilbert_alpha.mat'];        
% %         
% %     case {'beta_amp' 'beta_phase'  'beta_log_amp'}
% %         path_in_MEG  = [path_sub 'analysis/hilbert/'];
% %         loadfile_MEG = [path_in_MEG  'hilbert_beta.mat'];
% %         
% %     case {'low_gamma_amp' 'low_gamma_phase'  'low_gamma_log_amp'}
% %         path_in_MEG  = [path_sub 'analysis/hilbert/'];
% %         loadfile_MEG = [path_in_MEG  'hilbert_low_gamma.mat'];
% %         
% %     case {'mid_gamma_amp' 'mid_gamma_phase'  'mid_gamma_log_amp'}
% %         path_in_MEG  = [path_sub 'analysis/hilbert/'];
% %         loadfile_MEG = [path_in_MEG  'hilbert_mid_gamma.mat'];        
% % 
% %     case {'high_gamma_amp' 'high_gamma_phase'  'high_gamma_log_amp'}
% %         path_in_MEG  = [path_sub 'analysis/hilbert/'];
% %         loadfile_MEG = [path_in_MEG  'hilbert_high_gamma.mat'];
        
end

% % path_in_behav  = [path_sub 'behavioral/'];
loadfile_behav = si.path_behav;


% define output paths and files
mkdir([si.path_analysis 'stimulus_tracking/stim_corr_by_beta/stim_corr_results/'], sub);
path_out = [si.path_analysis 'stimulus_tracking/stim_corr_by_beta/stim_corr_results/' sub '/'];

switch phase_option
    case 1,
        savefile = [path_out sub '_stim_corr_' DV '_unwrap_detrend.mat'];
        return

    case 2,
        savefile = [path_out sub '_stim_corr_' DV '_unwrap_subtr_exp.mat'];    
        return
        
    otherwise
        savefile = [path_out sub '_stim_corr_' DV '_betaID=' num2str(input_betaID) '.mat'];    
end

%% load preprocessed data

disp('loading...')
load(loadfile_MEG);
load(loadfile_behav);
disp('done.')



%% define the MEG data

switch DV
    case 'raw'
        data_MEG = data_clean;
        DV_label = 'timecourse with general preprocessing';
        
    % lowpass filtered
    case 'lowpass_5hz'
        
        disp('filtering...')
        data_lp = data_clean;
        for i_trial = 1:length(data_clean.trial)
            data_lp.trial{i_trial} = ft_preproc_lowpassfilter( ...
                data_clean.trial{i_trial}, data_clean.fsample, 5);
        end
        disp('done');
        
        data_MEG = data_lp;
        DV_label = 'timecourse lowpass filtered at 5 Hz';        
        
        
    case 'lowpass_3hz'
        
        disp('filtering...')
        data_lp = data_clean;
        for i_trial = 1:length(data_clean.trial)
            data_lp.trial{i_trial} = ft_preproc_lowpassfilter( ...
                data_clean.trial{i_trial}, data_clean.fsample, 3);
        end
        disp('done');
        
        data_MEG = data_lp;
        DV_label = 'timecourse lowpass filtered at 3 Hz';   
        
        
    case 'lowpass_1hz'
        
        disp('filtering...')
        data_lp = data_clean;
        for i_trial = 1:length(data_clean.trial)
            data_lp.trial{i_trial} = ft_preproc_lowpassfilter( ...
                data_clean.trial{i_trial}, data_clean.fsample, 1);
        end
        disp('done');
        
        data_MEG = data_lp;
        DV_label = 'timecourse lowpass filtered at 1 Hz';
% %         
% % 
% %     % 3.33Hz
% %     case {'3.33Hz_amp'  '3.33Hz_log_amp'}
% %         data_MEG = data_hilbert_amp;
% %         
% %         if isLogAmp
% %             DV_label = '3.33 Hz oscillation log amplitude';
% %             
% %             for i_trial = 1:length(data_MEG.trial)
% %                 data_MEG.trial{i_trial} = log( data_MEG.trial{i_trial} );
% %             end      
% %             
% %         else
% %             DV_label = '3.33 Hz oscillation amplitude';
% %         end
% %         
% %         
% %     case '3.33Hz_phase'
% %         data_MEG = data_hilbert_phase;
% %         DV_label = '3.33 Hz oscillation phase';
% %         freq_phase = 1/.3;
% %         
% %         
% %     % delta
% %     case {'delta_amp'  'delta_log_amp'} 
% %         data_MEG = data_hilbert_amp;
% %         
% %         if isLogAmp
% %             DV_label = 'delta band log amplitude (1-4 Hz)';
% %             
% %             for i_trial = 1:length(data_MEG.trial)
% %                 data_MEG.trial{i_trial} = log( data_MEG.trial{i_trial} );
% %             end            
% %             
% %         else
% %             DV_label = 'delta band amplitude (1-4 Hz)';
% %         end
% %         
% %        
% %     case 'delta_phase'
% %         data_MEG = data_hilbert_phase;
% %         DV_label = 'delta band phase (1-4 Hz)';
% %         freq_phase = 2.5;
% % 
% %         
% %     % theta
% %     case {'theta_amp'  'theta_log_amp'} 
% %         data_MEG = data_hilbert_amp;
% %         
% %         if isLogAmp
% %             DV_label = 'theta band log amplitude (4-8 Hz)';
% %             
% %             for i_trial = 1:length(data_MEG.trial)
% %                 data_MEG.trial{i_trial} = log( data_MEG.trial{i_trial} );
% %             end
% %             
% %         else
% %             DV_label = 'theta band amplitude (4-8 Hz)';
% %         end
% %         
% %     case 'theta_phase'
% %         data_MEG = data_hilbert_phase;
% %         DV_label = 'theta band phase (4-8 Hz)';
% %         freq_phase = 6;
% %         
% %         
% %     % alpha    
% %     case {'alpha_amp'  'alpha_log_amp'} 
% %         data_MEG = data_hilbert_amp;
% %         
% %         if isLogAmp
% %             DV_label = 'alpha band log amplitude (8-12 Hz)';
% %             
% %             for i_trial = 1:length(data_MEG.trial)
% %                 data_MEG.trial{i_trial} = log( data_MEG.trial{i_trial} );
% %             end
% %             
% %         else
% %             DV_label = 'alpha band amplitude (8-12 Hz)';
% %         end         
% %         
% %         
% %     case 'alpha_phase'
% %         data_MEG = data_hilbert_phase;
% %         DV_label = 'alpha band phase (8-12 Hz)'; 
% %         freq_phase = 10;
% %         
% %         
% %     % beta    
% %     case {'beta_amp'  'beta_log_amp'} 
% %         data_MEG = data_hilbert_amp;
% %         
% %         if isLogAmp
% %             DV_label = 'beta band log amplitude (12-30 Hz)';
% %             
% %             for i_trial = 1:length(data_MEG.trial)
% %                 data_MEG.trial{i_trial} = log( data_MEG.trial{i_trial} );
% %             end
% %             
% %         else
% %             DV_label = 'beta band amplitude (12-30 Hz)';
% %         end         
% %         
% %     case 'beta_phase'
% %         data_MEG = data_hilbert_phase;
% %         DV_label = 'beta band phase (12-30 Hz)';            
% %         freq_phase = 21;
% %         
% %         
% %     % low gamma
% %     case {'low_gamma_amp'  'low_gamma_log_amp'} 
% %         data_MEG = data_hilbert_amp;
% %         
% %         if isLogAmp
% %             DV_label = 'low gamma band log amplitude (30-50 Hz)';
% %             
% %             for i_trial = 1:length(data_MEG.trial)
% %                 data_MEG.trial{i_trial} = log( data_MEG.trial{i_trial} );
% %             end
% %             
% %         else
% %             DV_label = 'low gamma band amplitude (30-50 Hz)';
% %         end   
% %         
% %         
% %     case 'low_gamma_phase'
% %         data_MEG = data_hilbert_phase;
% %         DV_label = 'low gamma band phase (30-50 Hz)';            
% %         freq_phase = 40;
% %         
% %         
% %     % mid gamma
% %     case {'mid_gamma_amp'  'mid_gamma_log_amp'} 
% %         data_MEG = data_hilbert_amp;
% %         
% %         if isLogAmp
% %             DV_label = 'mid gamma band log amplitude (50-80 Hz)';
% %             
% %             for i_trial = 1:length(data_MEG.trial)
% %                 data_MEG.trial{i_trial} = log( data_MEG.trial{i_trial} );
% %             end
% %             
% %         else
% %             DV_label = 'mid gamma band amplitude (50-80 Hz)';
% %         end        
% %         
% %         
% %     case 'mid_gamma_phase'
% %         data_MEG = data_hilbert_phase;
% %         DV_label = 'mid gamma band phase (50-80 Hz)';
% %         freq_phase = 65;
% %         
% %         
% %     % high gamma
% %     case {'high_gamma_amp'  'high_gamma_log_amp'} 
% %         data_MEG = data_hilbert_amp;
% %         
% %         if isLogAmp
% %             DV_label = 'high gamma band log amplitude (80-150 Hz)';
% %             
% %             for i_trial = 1:length(data_MEG.trial)
% %                 data_MEG.trial{i_trial} = log( data_MEG.trial{i_trial} );
% %             end
% %             
% %         else
% %             DV_label = 'high gamma band amplitude (80-150 Hz)';
% %         end
% %         
% %         
% %     case 'high_gamma_phase'
% %         data_MEG = data_hilbert_phase;
% %         DV_label = 'high gamma band phase (80-150 Hz)';           
% %         freq_phase = 115;

end


%% extract trials

data_MEG = extract_trials(data_MEG, sub);


%% filter behavioral data to exclude trials that don't have MEG recording

dat_fields  = fieldnames(data);
for i = 1:length(dat_fields)
    if eval(['length(data.' dat_fields{i} ') == 360'])
        eval(['data.' dat_fields{i} ' = data.' dat_fields{i} '(si.good_trials);']);
    end
end


stim_fields = fieldnames(data.stim);
for i = 1:length(stim_fields)
    if eval(['length(data.stim.' stim_fields{i} ') == 360'])
        eval(['data.stim.' stim_fields{i} ' = data.stim.' stim_fields{i} '(si.good_trials);']);
    end
end


%% append stimulus / behavioral data to MEG data

data_MEG.behav = data;
data_MEG.stim  = data.stim;


%% calculate within-condition correlations for each sensor and experimental condition

nSensors = size( data_MEG.trial{1}, 1 );

%%%% WITHIN STIMULUS CORRELATION %%%%%

toneDur_inSecs  = .3;
fs              = data_MEG.fsample;
nSamplesPerTone = toneDur_inSecs * fs;

t0 = GetSecs;

for i_beta = input_betaID % 1:5

    i_beta
    
    sigIDs = [0 1];
    for i_sig = 1:2
        
        predIDs = [-1 0 1];      
        for i_pred = 1:3
                        
            % extract data for trials at a specific level of beta and
            % predictionID and sigma_e
            filter = data_MEG.stim.betaID == i_beta  ...
                        & data_MEG.stim.predID == predIDs(i_pred) ...
                        & data_MEG.stim.sigma_e_ID == sigIDs(i_sig);

            field_length = length(data_MEG.trial);
                    
            d       = myft_trialfilter(data_MEG, filter, field_length);
            d.behav = myft_trialfilter(d.behav, filter, field_length);
            d.stim  = myft_trialfilter(d.stim, filter, field_length);            
                
                        
            % now that the data for this stimulus set is selected,
            % conduct all possible correlations at each sensor...
            
            nTrials = length(d.trial);            
            nPerm   = nchoosek(nTrials, 2);
            
            corr_within{i_beta, i_sig, i_pred}.r = zeros(nSensors, nPerm);
            corr_within{i_beta, i_sig, i_pred}.p = zeros(nSensors, nPerm);
            corr_within{i_beta, i_sig, i_pred}.trial_pair = zeros(nPerm, 2);
            corr_within{i_beta, i_sig, i_pred}.DV_label = DV_label;
            
            ind = 0;
            for i_trial1 = 1 : nTrials-1
                for i_trial2 = i_trial1+1 : nTrials

                    ind = ind + 1;
                    
                    corr_within{i_beta, i_sig, i_pred}.trial_pair(ind, :) = [i_trial1,  i_trial2];

                    for i_sensor = 1 : nSensors
                        
                        % compute correlation after removing the final tone
                        s1 = d.trial{i_trial1}(i_sensor, 1:end-nSamplesPerTone);
                        s2 = d.trial{i_trial2}(i_sensor, 1:end-nSamplesPerTone);

                        if isPhase
                            % phase_option values
                            % 0 --> do circular correlation on phase
                            % 1 --> unwrap phase, detrend, and do linear correlation
                            % 2 --> unwrap phase, subtract expected phase, and do linear correlation                            
                            switch phase_option
                                case 0
                                    [r p] = circ_corrcc(s1', s2');
                                    
                                case 1
                                    s1 = detrend(unwrap(s1));
                                    s2 = detrend(unwrap(s2));
                                    [r p] = corr(s1', s2', 'type', 'Pearson');
                                    
                                case 2
                                    t  = d.time{1}(1,:);  % time values for series
                                    m  = 2*pi*freq_phase; % slope of expected value for unwrapped phase
                                   
                                    s1 = unwrap(s1);
                                    b1 = s1(1);
                                    s1_exp = m*t + b1;
                                    s1_exp = s1_exp(1:length(s1)); % remove time points corresponding to final tone
                                    s1 = s1 - s1_exp;
                                    
                                    s2 = unwrap(s2);
                                    b2 = s2(1);
                                    s2_exp = m*t + b2;
                                    s2_exp = s2_exp(1:length(s2));                                    
                                    s2 = s2 - s2_exp;
                                    
                                    [r p] = corr(s1', s2', 'type', 'Pearson');
                                    
                                    
                            end
                            
                        else
                            [r p] = corr(s1', s2', 'type', 'Pearson');
                            
                        end
                        
                        corr_within{i_beta, i_sig, i_pred}.r(i_sensor, ind) = r;
                        corr_within{i_beta, i_sig, i_pred}.p(i_sensor, ind) = p;

                    end

                end
            end            
            
        end
    end
end

tf_wtn = GetSecs - t0


%% calculate across-condition correlations for each sensor and experimental condition

%%%% ACROSS STIMULUS CORRELATION %%%%%


t0 = GetSecs;

for i_beta = input_betaID % 1:5

    i_beta
    
    sigIDs = [0 1];
    for i_sig = 1:2
        
        predIDs = [-1 0 1];      
        for i_pred_wtn = 1:2

            for i_pred_acr = i_pred_wtn + 1 : 3
                
                % extract data for trials at a specific level of beta and
                % predictionID and sigma_e
                filter = data_MEG.stim.betaID == i_beta  ...
                            & data_MEG.stim.predID == predIDs(i_pred_wtn) ...
                            & data_MEG.stim.sigma_e_ID == sigIDs(i_sig);

                field_length = length(data_MEG.trial);

                d_wtn       = myft_trialfilter(data_MEG, filter, field_length);
                d_wtn.behav = myft_trialfilter(d_wtn.behav, filter, field_length);
                d_wtn.stim  = myft_trialfilter(d_wtn.stim, filter, field_length);                   
                
                
                            
                % extract data for trials at a specific level of beta and
                % predictionID and sigma_e
                filter = data_MEG.stim.betaID == i_beta  ...
                            & data_MEG.stim.predID == predIDs(i_pred_acr) ...
                            & data_MEG.stim.sigma_e_ID == sigIDs(i_sig);

                field_length = length(data_MEG.trial);

                d_acr       = myft_trialfilter(data_MEG, filter, field_length);
                d_acr.behav = myft_trialfilter(d_acr.behav, filter, field_length);
                d_acr.stim  = myft_trialfilter(d_acr.stim, filter, field_length);                    
                
            
                
                nTrials_wtn = length(d_wtn.trial);
                nTrials_acr = length(d_acr.trial);
                
                nPerm   = nTrials_wtn * nTrials_acr;

                corr_across{i_beta, i_sig, i_pred_wtn, i_pred_acr}.r = zeros(nSensors, nPerm);
                corr_across{i_beta, i_sig, i_pred_wtn, i_pred_acr}.p = zeros(nSensors, nPerm);
                corr_across{i_beta, i_sig, i_pred_wtn, i_pred_acr}.trial_pair = zeros(nPerm, 2);
                corr_across{i_beta, i_sig, i_pred_wtn, i_pred_acr}.DV_label = DV_label;
                
                ind = 0;
                for i_trial_wtn = 1 : nTrials_wtn
                    for i_trial_acr = 1 : nTrials_acr

                        ind = ind + 1;

                        corr_across{i_beta, i_sig, i_pred_wtn, i_pred_acr}.trial_pair(ind, :) = ...
                            [i_trial_wtn,  i_trial_acr];

                        for i_sensor = 1 : nSensors
                            
                            % compute correlation after removing the final tone
                            s1 = d_wtn.trial{i_trial_wtn}(i_sensor, 1:end-nSamplesPerTone);
                            s2 = d_acr.trial{i_trial_acr}(i_sensor, 1:end-nSamplesPerTone);
                        
                            
                            if isPhase
                                % phase_option values
                                % 0 --> do circular correlation on phase
                                % 1 --> unwrap phase, detrend, and do linear correlation
                                % 2 --> unwrap phase, subtract expected phase, and do linear correlation                            
                                switch phase_option
                                    case 0
                                        [r p] = circ_corrcc(s1', s2');

                                    case 1
                                        s1 = detrend(unwrap(s1));
                                        s2 = detrend(unwrap(s2));
                                        [r p] = corr(s1', s2', 'type', 'Pearson');

                                    case 2
                                        t  = d.time{1}(1,:);  % time values for series
                                        m  = 2*pi*freq_phase; % slope of expected value for unwrapped phase

                                        s1 = unwrap(s1);
                                        b1 = s1(1);
                                        s1_exp = m*t + b1;
                                        s1_exp = s1_exp(1:length(s1)); % remove time points corresponding to final tone                                        
                                        s1 = s1 - s1_exp;

                                        s2 = unwrap(s2);
                                        b2 = s2(1);
                                        s2_exp = m*t + b2;
                                        s2_exp = s2_exp(1:length(s2)); % remove time points corresponding to final tone                                        
                                        s2 = s2 - s2_exp;

                                        [r p] = corr(s1', s2', 'type', 'Pearson');


                                end
                                
                            else
                                [r p] = corr(s1', s2', 'type', 'Pearson');
                                
                            end
                            

                            corr_across{i_beta, i_sig, i_pred_wtn, i_pred_acr}.r(i_sensor, ind) = r;
                            corr_across{i_beta, i_sig, i_pred_wtn, i_pred_acr}.p(i_sensor, ind) = p;

                        end

                    end
                end            

            end
        end
        
    end
 end

tf_acr = GetSecs - t0

%% do ANOVA to test correlations for within vs across, including stimulus condition factors

zvals  = [];
Beta   = [];
Sigma  = [];
Pred   = [];
Within = [];


for i_sensor = 1:271

zvals_i  = [];
Beta_i   = [];
Sigma_i  = [];
% Pred_i   = [];
Within_i = [];
    
for i_beta = input_betaID  % 1:5
    for i_sig = 1:2
        for i_pred = 1:3


            % within-stimulus correlations
            rs_within = corr_within{i_beta, i_sig, i_pred}.r(i_sensor, :)';
            zvals_i = [zvals_i; r2z( rs_within )];

            Beta_i   = [Beta_i;   i_beta * ones(size(rs_within))];
            Sigma_i  = [Sigma_i;  i_sig * ones(size(rs_within))];
%             Pred_i   = [Pred_i;   i_pred * ones(size(rs_within))];
            Within_i = [Within_i; ones(size(rs_within))];


            % across-stimulus correlations
            rs_across = [];
            switch i_pred
                case 1
%                     rs_across = [ corr_across{i_beta, i_sig, 1, 2}.r(i_sensor,:) ...
%                                   corr_across{i_beta, i_sig, 1, 3}.r(i_sensor,:) ]';
                    rs_across = [ corr_across{i_beta, i_sig, 1, 2}.r(i_sensor,:)  ]';

                              
                case 2
%                     rs_across = [ corr_across{i_beta, i_sig, 1, 2}.r(i_sensor,:) ...
%                                   corr_across{i_beta, i_sig, 2, 3}.r(i_sensor,:) ]';
                    rs_across = [ corr_across{i_beta, i_sig, 1, 3}.r(i_sensor,:)  ]';


                case 3
%                     rs_across = [ corr_across{i_beta, i_sig, 1, 3}.r(i_sensor,:) ...
%                                   corr_across{i_beta, i_sig, 2, 3}.r(i_sensor,:) ]';
                    rs_across = [ corr_across{i_beta, i_sig, 2, 3}.r(i_sensor,:)  ]';

            end

            zvals_i = [zvals_i; r2z( rs_across )];

            Beta_i   = [Beta_i;   i_beta * ones(size(rs_across))];
            Sigma_i  = [Sigma_i;  i_sig * ones(size(rs_across))];
%             Pred_i   = [Pred_i;   i_pred * ones(size(rs_across))];
            Within_i = [Within_i; zeros(size(rs_across))];
                          

        end
    end
end

zvals(i_sensor, :)  = zvals_i';
Beta(i_sensor, :)   = Beta_i';
Sigma(i_sensor, :)  = Sigma_i';
% Pred(i_sensor, :)   = Pred_i';
Within(i_sensor, :) = Within_i';

% CHANGE FROM VERSION 1
% [p, atab{i_sensor}] = anovan(zvals_i, ...
%                 {Beta_i Sigma_i Pred_i Within_i}, ...
%                 'model', 4, 'sstype', 2, 'display', 'off', ...
%                 'varnames', strvcat('Beta', 'Sigma', 'Pred', 'Within'));


% [p, atab{i_sensor}] = anovan(zvals_i, ...
%                 {Beta_i Sigma_i Within_i}, ...
%                 'model', 3, 'sstype', 2, 'display', 'off', ...
%                 'varnames', strvcat('Beta', 'Sigma', 'Within'));            


[p, atab{i_sensor}] = anovan(zvals_i, ...
                {Sigma_i Within_i}, ...
                'model', 2, 'sstype', 2, 'display', 'off', ...
                'varnames', strvcat('Sigma', 'Within'));       

anova_p(i_sensor, :) = p';


end


nFactors = 3;
for i_factor = 1 : nFactors
    factor_labels{i_factor} = atab{1}{i_factor+1, 1};
end


%% do statistical tests on correlations for within vs across

% for i_beta = 1:5
%     for i_sig = 1:2
%         for i_pred = 1:3
%             
%             for i_sensor = 1:nSensors
%                 
%                 rs_within = corr_within{i_beta, i_sig, i_pred}.r(i_sensor,:);
%                 
%                 rs_across = [];
%                 switch i_pred
%                     case 1
%                         rs_across = [ corr_across{i_beta, i_sig, 1, 2}.r(i_sensor,:) ...
%                                       corr_across{i_beta, i_sig, 1, 3}.r(i_sensor,:) ];
%                                   
%                     case 2
%                         rs_across = [ corr_across{i_beta, i_sig, 1, 2}.r(i_sensor,:) ...
%                                       corr_across{i_beta, i_sig, 2, 3}.r(i_sensor,:) ];
%                                   
%                     case 3
%                         rs_across = [ corr_across{i_beta, i_sig, 1, 3}.r(i_sensor,:) ...
%                                       corr_across{i_beta, i_sig, 2, 3}.r(i_sensor,:) ];
%                 end
%                 
%                 [h, p, ci, stats] = ttest2( r2z( rs_within )', r2z( rs_across )' );
%                 
%                 corr_ttest{i_beta, i_sig, i_pred}.t(i_sensor)  = stats.tstat;
%                 corr_ttest{i_beta, i_sig, i_pred}.df(i_sensor) = stats.df;
%                 corr_ttest{i_beta, i_sig, i_pred}.sd(i_sensor) = stats.sd;
%                 corr_ttest{i_beta, i_sig, i_pred}.p(i_sensor)  = p;
%                 corr_ttest{i_beta, i_sig, i_pred}.DV_label     = DV_label;
%                 
%             end
%             
%         end
%     end
% end


%% save the correlation data for later analysis


disp('saving...')
% save(savefile, 'corr_within', 'corr_across', 'corr_ttest', '-v7.3');
% save(savefile, 'corr_within', 'corr_across', ...
%                'anova_p', 'atab', ...
%                'zvals', 'Beta', 'Sigma', 'Pred', 'Within', '-v7.3');
% save(savefile, 'corr_within', 'corr_across', ...
%                'anova_p', 'atab', 'factor_labels', ...
%                'zvals', 'Beta', 'Sigma', 'Within', '-v7.3');

save(savefile, 'corr_within', 'corr_across', ...
               'anova_p', 'atab', 'factor_labels', ...
               'zvals', 'Sigma', 'Within', '-v7.3');

disp('done.')

end
