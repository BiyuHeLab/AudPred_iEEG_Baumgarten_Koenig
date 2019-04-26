function sfa_expt4_HisTrack_CombOddEvenSets(sub, tonedur_text, baseline_correct)
%TJB: New version of tone_history_CV_plot.m
%Aim: Combine results from both even and odd traning set regressions 
%(i.e., results from K-value computation (via sfa_expt3_HisTrack_CompSaveKval_SLIDE16.m)
%Implementaion:
%1) Compute k-values for oth train and test set (SLIDE16 = last 16 tones, sliding time window)
%2) Save output (matrix k-values)

%% 0.1) Specify vars, paths, and setup fieldtrip
addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
% addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')

sfa_expt4_setVars
paths_sfa_expt4 = sfa_expt4_paths(vars.location);

%Add base dir and own script dir
addpath(genpath(paths_sfa_expt4.BaseDir));
addpath(genpath(paths_sfa_expt4.ScriptsDir));

%% 0.2) Determine subject- and condition- specific parameters (whole-recording)
sfa_expt4_subjectinfo %load subject info file (var: si)
subs_PreProcSettings = sfa_expt4_subjectPreProcSettings; %load in file with individual preproc infos

%Tone duration condition for data load-in
if strcmp(tonedur_text,'0.2')
    tonedur_title = '200';
elseif  strcmp(tonedur_text,'0.4')
    tonedur_title = '400';
end

if baseline_correct == 1
    path_load = [paths_sfa_expt4.ECoGdata_HisTrack sub '/baseline_corrected/Exp/'];
    path_save = [paths_sfa_expt4.ECoGdata_HisTrack sub '/baseline_corrected/Exp/'];
    BC_text = 'BC';

elseif baseline_correct == 0
    path_load = [paths_sfa_expt4.ECoGdata_HisTrack sub '/not_baseline_corrected/Exp/'];
    path_save = [paths_sfa_expt4.ECoGdata_HisTrack sub '/not_baseline_corrected/Exp/'];
    BC_text = 'NBC';
end


%% 1) EVEN set analysis
%1.1 load regstats data
 load([path_load sub '_HisTrackSLIDE16_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur_trainEven.mat']);

%1.2 Define analysis parameters
fsample = 512;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

win_size    = 25;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / fsample) * 1000;

nSensors = size(stats_linear,3); %all sensors, not only selected ones
nWindows = length(windows);
n_k      = size(stats_linear,2) %number model order;

k_list = 1:n_k;

%1.3 Restructure regstats data and read out relevant stats for model selection from even set
for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        
        for i_k = 1:n_k
            
%             linear_AIC(i_k) = stats_linear{i_win, i_k, i_sensor}.AIC_train; %n*log(sig2_train) + 2*K
%             %AIC = Akaike information criterion - esitmator of relative quality of statistical models for given set of data; means for model selection
%             linear_BIC(i_k) = stats_linear{i_win, i_k, i_sensor}.BIC_train; %n*log(sig2_train) + K*log(n)
%             %BIC = Bayesian information criterion - criterion for model selection
            sig2s(i_k) = stats_linear{i_win, i_k, i_sensor}.sig2_test;  %sum of squared residuals divided by length of model order
            R2s(i_k) = stats_linear{i_win, i_k, i_sensor}.R2_test; %r-square for test set based on estimated data     
            
        end
        
%1.4 Compute criteria for model selection from even set        
        % linear AIC
%         linear_AIC_w_even{i_win}(i_sensor,:) = akaikeWeights(linear_AIC);
%         [max_w, max_w_ind] = max( linear_AIC_w_even{i_win}(i_sensor,:) );
%         % linear_AIC_k{i_win}(i_sensor) = sum(k_list .* linear_AIC_w{i_win}(i_sensor,:));
%         linear_AIC_k_even{i_win}(i_sensor) = k_list(max_w_ind);
%         linear_AIC_k_w_even{i_win}(i_sensor) = max_w;
%         
%         % linear BIC
%         linear_BIC_w_even{i_win}(i_sensor,:) = akaikeWeights(linear_BIC);
%         [max_w, max_w_ind] = max( linear_BIC_w_even{i_win}(i_sensor,:) );
%         linear_BIC_k_even{i_win}(i_sensor) = k_list(max_w_ind);
%         linear_BIC_k_w_even{i_win}(i_sensor) = max_w;

        % cross validation k and VAR_res
        sig2_test_even{i_win}(i_sensor, :) = sig2s; %sum of squared residuals divided by length of model order
        [min_s, min_s_ind] = min( sig2s );%read out minimum squared residuals divided by length of model order and the respective position in array
        sig2_test_k_even{i_win}(i_sensor) = k_list(min_s_ind); %k-prime-value (i.e., preferred model order/Number of tones back best explaining MEG data)
        sig2_test_k_sig_even{i_win}(i_sensor) = min_s; %minimum squared residuals divided by length of model order of K-value (i.e., preferred model order/Number of tones back)
                
        % cross validation R^2
        R2_test_even{i_win}(i_sensor, :) = R2s;%r-square for test set based on estimated data   
        R2_test_k_even{i_win}(i_sensor) = R2s(min_s_ind);
                
    end
end
%1.5 Clean-up
clear stats_linear linear_AIC linear_BIC sig2s R2s

%% 2) ODD set analysis
%2.1 load regstats data
load([path_load sub '_HisTrackSLIDE16_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur_trainOdd.mat']);

%1.2 Define analysis parameters
win_size    = 25;  %25 data points with fs 512 = ~49ms as mentioned in manuscript
s_per_win = (1/fsample)*win_size;
win_overlap = 0; %as mentioned in manuscript

fsample = 512;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

windows = [1 win_size];
while windows(end,end) < nSamplesPerTone
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end

if windows(end,end) > nSamplesPerTone
    windows(end,:) = [];
end

windows_inms = (windows / fsample) * 1000;

nSensors = size(stats_linear,3); %all sensors, not only selected ones
nWindows = length(windows);
n_k      = size(stats_linear,2) %number model order;

k_list = 1:n_k;

%2.3 Restructure regstats data and read out relevant stats for model selection from odd set
for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        
        for i_k = 1:n_k
            
%             linear_AIC(i_k) = stats_linear{i_win, i_k, i_sensor}.AIC_train;
%             linear_BIC(i_k) = stats_linear{i_win, i_k, i_sensor}.BIC_train;

            sig2s(i_k) = stats_linear{i_win, i_k, i_sensor}.sig2_test;
            R2s(i_k) = stats_linear{i_win, i_k, i_sensor}.R2_test;
            
        end
        
%2.4 Compute criteria for model selection from odd set        
%         % linear AIC
%         linear_AIC_w_odd{i_win}(i_sensor,:) = akaikeWeights(linear_AIC);
%         [max_w, max_w_ind] = max( linear_AIC_w_odd{i_win}(i_sensor,:) );
%         % linear_AIC_k{i_win}(i_sensor) = sum(k_list .* linear_AIC_w{i_win}(i_sensor,:));
%         linear_AIC_k_odd{i_win}(i_sensor) = k_list(max_w_ind);
%         linear_AIC_k_w_odd{i_win}(i_sensor) = max_w;
%         
%         % linear BIC
%         linear_BIC_w_odd{i_win}(i_sensor,:) = akaikeWeights(linear_BIC);
%         [max_w, max_w_ind] = max( linear_BIC_w_odd{i_win}(i_sensor,:) );
%         linear_BIC_k_odd{i_win}(i_sensor) = k_list(max_w_ind);
%         linear_BIC_k_w_odd{i_win}(i_sensor) = max_w;

        % cross validation k
        sig2_test_odd{i_win}(i_sensor, :) = sig2s;
        [min_s, min_s_ind] = min( sig2s );
        sig2_test_k_odd{i_win}(i_sensor) = k_list(min_s_ind);
        sig2_test_k_sig_odd{i_win}(i_sensor) = min_s;
                
        % cross validation R^2
        R2_test_odd{i_win}(i_sensor, :) = R2s;
        R2_test_k_odd{i_win}(i_sensor) = R2s(min_s_ind);
        
    end
end
%2.5 Clean-up
clear stats_linear linear_AIC linear_BIC sig2s R2s


%% 3) Combine EVEN & ODD sets
%i.e., average model selection criteria by averaging across sets

for i_win = 1:nWindows
    for i_sensor = 1:nSensors
        
%         % linear AIC
%         linear_AIC_w{i_win}(i_sensor,:) = mean([linear_AIC_w_even{i_win}(i_sensor,:);  linear_AIC_w_odd{i_win}(i_sensor,:)]);
%         linear_AIC_k{i_win}(i_sensor) = mean([linear_AIC_k_even{i_win}(i_sensor);  linear_AIC_k_odd{i_win}(i_sensor)]);
%         linear_AIC_k_w{i_win}(i_sensor) = mean([linear_AIC_k_w_even{i_win}(i_sensor);  linear_AIC_k_w_odd{i_win}(i_sensor)]);
% 
%         
%         % linear BIC
%         linear_BIC_w{i_win}(i_sensor,:) = mean([linear_BIC_w_even{i_win}(i_sensor,:);  linear_BIC_w_odd{i_win}(i_sensor,:)]);
%         linear_BIC_k{i_win}(i_sensor) = mean([linear_BIC_k_even{i_win}(i_sensor);  linear_BIC_k_odd{i_win}(i_sensor)]);
%         linear_BIC_k_w{i_win}(i_sensor) = mean([linear_BIC_k_w_even{i_win}(i_sensor);  linear_BIC_k_w_odd{i_win}(i_sensor)]);
        
        
        % cross validation k
        sig2_test{i_win}(i_sensor, :) = mean([sig2_test_even{i_win}(i_sensor,:);  sig2_test_odd{i_win}(i_sensor,:)]); 
        sig2_test_k{i_win}(i_sensor) = mean([sig2_test_k_even{i_win}(i_sensor);  sig2_test_k_odd{i_win}(i_sensor)]); 
        %k-prime values (i.e., number of model with minimum sum of squared residuals divided by length of model order
        %CAVE: Notes Brian: Matlab data k-prime is offset by +1 cmpared to manuscript k-prime (due to intercept term)
        sig2_test_k_sig{i_win}(i_sensor) = mean([sig2_test_k_sig_even{i_win}(i_sensor);  sig2_test_k_sig_odd{i_win}(i_sensor)]);
       
        
        % cross validation R^2
        R2_test{i_win}(i_sensor, :) = mean([R2_test_even{i_win}(i_sensor,:);  R2_test_odd{i_win}(i_sensor,:)]);
        R2_test_k{i_win}(i_sensor) = mean([R2_test_k_even{i_win}(i_sensor); R2_test_k_odd{i_win}(i_sensor)]);        
        
    end
end



%% 4) Save output variables
savefile = [path_save sub '_HisTrackSLIDE16_regstatsEXP_' BC_text '_' tonedur_title 'msToneDur_kModelComp_CombSets.mat'];

save(savefile, 'sig2_test', 'sig2_test_k', 'sig2_test_k_sig', ...
               'R2_test', 'R2_test_k');

% save(savefile, 'linear_AIC_w', 'linear_BIC_w',  ... %w = weights
%                'linear_AIC_k', 'linear_BIC_k',  ... %k = k-values
%                'linear_AIC_k_w', 'linear_BIC_k_w', ... %k_w = w of selected k
%                'sig2_test', 'sig2_test_k', 'sig2_test_k_sig', ...
%                'R2_test', 'R2_test_k');
        
end