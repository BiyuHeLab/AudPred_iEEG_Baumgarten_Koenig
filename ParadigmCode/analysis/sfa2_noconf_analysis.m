clear 

%% define script paths

addpath('/data/gogodisk1/brian/scripts/');
addpath(genpath('/data/gogodisk1/brian/scripts/mopttb/'))


%% load behavioral data

path = '/data/gogodisk1/amy/data/sfa_expt2_noconf/';

subs = {'S2' 'S3' 'S4' 'S5' 'S6' 'S7' 'S8' 'S9' 'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18'};

bad = {'S5' ...   % rho = .05
       'S12' ...  % rho = .07
       };

subs = setdiff(subs,bad);

sub = 'S2';

switch sub
    case 'S2',  load([path 'S2/behavioral/sfa_expt2 2014-09-10 14-27-52.mat']);
    case 'S3',  load([path 'S3/behavioral/sfa_expt2 2014-09-16 10-50-03.mat']);
    case 'S4',  load([path 'S4/behavioral/sfa_expt2 2014-09-19 11-59-16.mat']);
    case 'S5',  load([path 'S5/behavioral/sfa_expt2 2014-09-19 14-53-12.mat']);
    case 'S6',  load([path 'S6/behavioral/sfa_expt2 2014-09-24 13-46-35.mat']);
    case 'S7',  load([path 'S7/behavioral/sfa_expt2 2014-09-25 10-43-38.mat']);
    case 'S8',  load([path 'S8/behavioral/sfa_expt2 2014-09-25 13-49-54.mat']);
    case 'S9',  load([path 'S9/behavioral/sfa_expt2 2014-09-26 10-22-35.mat']);
    case 'S10', load([path 'S10/behavioral/sfa_expt2 2014-10-03 14-51-59.mat']);
    case 'S11', load([path 'S11/behavioral/sfa_expt2 2014-10-08 14-36-36.mat']);
    case 'S12', load([path 'S12/behavioral/sfa_expt2 2014-10-24 11-22-42.mat']);
    case 'S13', load([path 'S13/behavioral/sfa_expt2 2014-10-28 10-41-21.mat']);
    case 'S14', load([path 'S14/behavioral/sfa_expt2 2014-11-06 11-17-11.mat']);
    case 'S15', load([path 'S15/behavioral/sfa_expt2 2014-11-10 09-58-27.mat']);
    case 'S16', load([path 'S16/behavioral/sfa_expt2 2014-11-12 10-56-04.mat']);
    case 'S17', load([path 'S17/behavioral/sfa_expt2 2014-11-22 15-39-30.mat']);
    case 'S18', load([path 'S18/behavioral/sfa_expt2 2014-11-24 11-24-09.mat']);

end



%% define filters

% general filter - ensure all behavioral responses were entered properly
filter_resp = data.resp_beta >= 0 & data.resp_prob >= 0 & data.trialNum > 0;
filter_rt   = data.r1_RT > 0 & data.r2_RT > 0;


% manually defined filters
filter_manual = ones(1,360);

switch sub
%     case 'LS'
%         % manual filter for LS
%         % - omit trials in block 2 where r1_RT < .1 (button got stuck)
%         % - omit trials in block 3 (button was stuck for most of block)
%         manualFilter = manualFilter & data.r1_RT > .1;
%         manualFilter = manualFilter & (data.trialNum < 61 | data.trialNum > 90);
end

filter = filter_resp & filter_rt & filter_manual;



%% apply filters
% filtered data will be applied to newly defined structs "dataf" and "stimf"

% definte dataf
dfn = fieldnames(data);
for j = 1:length(dfn)
    if eval(['length(data.' dfn{j} ') == 360'])
        eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);
    end
end


% define stimf
sfn = fieldnames(stim);
for j = 1:length(sfn)
    if eval(['length(stim.' sfn{j} ') == 360'])
        eval(['stimf.' sfn{j} ' = stim.' sfn{j} '(filter);']);
    end
end


%% model selection

del_logf_pen  = abs( stimf.logf_final - stimf.logf_pen );
del_logf_pred = abs( stimf.logf_final - stimf.logf_pred );
del_SD_pen    = abs( stimf.logf_final - stimf.logf_pen ) ./ stimf.sigma_e;
del_SD_pred   = abs( stimf.logf_final - stimf.logf_pred ) ./ stimf.sigma_e;

dv = 5-dataf.resp_prob;

% dv(stimf.finalID < 0) = -1*dv(stimf.finalID < 0);

iv_list = {del_logf_pen,  del_logf_pred,  del_SD_pen,  del_SD_pred};
iv_labels2 = {'logf pen', 'logf pred', 'SD pen', 'SD pred'};
iv_labels = {'|log f_{t+1} - log f_t|', '|log f_{t+1} - E[log f_{t+1}]|', ...
             '|log f_{t+1} - log f_t| / \sigma_\epsilon', '|log f_{t+1} - E[log f_{t+1}]| / \sigma_\epsilon'};
% iv_labels2 = {'|log f_{t+1} - log f_t|', '|log f_{t+1} - E[log f_{t+1}|', ...
%              '|log f_{t+1} - log f_t| / s_e', '|log f_{t+1} - E[log f_{t+1}]| / s_e'};

         
h=figure;
fs = 15;
set(h,'DefaultAxesFontSize',fs);
for i_iv = 1:4
    iv = iv_list{i_iv};

    subplot(2,2,i_iv);
    plot(iv, dv, 'b.');
    ylim([.5 5.5])
    xlabel(iv_labels{i_iv});
    ylabel('prob rating')

end


for i_iv = 1:4
    iv = iv_list{i_iv};
    
    m = [dv' iv'];
    in.DATA = m;
    out = MATLAB_Ordered_Probit_Estimate(in);

    logL(i_iv) = out.Likelihood.LLV;

    % recover model probabilities
    x = iv*out.Beta;

    probit_x{i_iv} = x;
    probit_beta{i_iv} = out.Beta;
    probit_cutpoints{i_iv} = out.Cut_Points;

    for i_trial = 1:length(x)
        prob_model{i_iv}(i_trial) = 1 + sum( x(i_trial) > out.Cut_Points );
    end

% % %     [beta mu l] = fit_ordered_probit_MLE(iv', dv', 4);
% % %     logL(i_iv) = l;
% % %     x = iv*beta;
% % % 
% % %     probit_x{i_iv} = x;
% % %     probit_beta{i_iv} = beta;
% % %     probit_cutpoints{i_iv} = mu;
% % %     
% % %     for i_trial = 1:length(x)
% % %         prob_model{i_iv}(i_trial) = 1 + sum( x(i_trial) > mu );
% % %     end

   
% % %     [mu sd l] = fit_ordered_probit_MLE2(iv', dv', 4);
% % %     logL(i_iv) = l;
% % %     x = iv;
% % % 
% % %     probit_x{i_iv} = x;
% % %     probit_sd{i_iv} = sd;
% % %     probit_cutpoints{i_iv} = mu;
% % %     
% % %     for i_trial = 1:length(x)
% % %         prob_model{i_iv}(i_trial) = 1 + sum( x(i_trial) > mu );
% % %     end
    
    
    for k = 1:5
        p_b(k) = sum(dv==k)/length(dv);
        p_m{i_iv}(k) = sum(prob_model{i_iv}==k)/length(prob_model{i_iv});
    end
end

figure;
subplot(211);
b = [p_b' p_m{1}' p_m{2}' p_m{3}' p_m{4}'];
bar(b);

subplot(212);
bdiff = [p_b'-p_m{1}' p_b'-p_m{2}' p_b'-p_m{3}' p_b'-p_m{4}'];
bar(bdiff);

K = 5;
N = length(dv);

BIC = -2*logL + K.*log(N);
AIC = -2*logL + 2.*K.*N./(N-K-1);

BICw = exp( -.5*(BIC-min(BIC)) ) / sum( exp( -.5*(BIC-min(BIC)) ));
AICw = exp( -.5*(AIC-min(AIC)) ) / sum( exp( -.5*(AIC-min(AIC)) ))

figure;
fs = 15;

bar(AICw);
set(gca, 'XTickLabel', iv_labels2);
ylabel('Akaike weight', 'FontSize', fs)
set(gca, 'FontSize', fs);



