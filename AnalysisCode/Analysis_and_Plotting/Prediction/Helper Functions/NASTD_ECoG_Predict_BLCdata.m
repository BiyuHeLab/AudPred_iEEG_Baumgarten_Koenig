function BLC_data = NASTD_ECoG_Predict_BLCdata(BLwindow, input_data)

%0. Find samples coresponding to BLwindow markers
Dist = abs(input_data.time{1} - BLwindow(1));
minDist = min(Dist);
BLstartsample = find(Dist == minDist);
clear Dist minDist

Dist = abs(input_data.time{1} - BLwindow(2));
minDist = min(Dist);
BLstopsample = find(Dist == minDist);
clear Dist minDist 

%1. Demean and baseline-correct (500 ms prestim) for low freqs
disp(['Baseline Correction for :' num2str(BLwindow(1)) ' to ' num2str(BLwindow(2)) ' sec'])

cfg = [];
cfg.demean = 'yes';
cfg.method = 'trial';
cfg.baselinewindow = [BLwindow(1) BLwindow(2)];
BLC_data = ft_preprocessing(cfg,input_data);

%2. Copy special subfields
BLC_data.cfg.info_elec      = input_data.cfg.info_elec;
BLC_data.cfg.info_ref       = input_data.cfg.info_ref;
BLC_data.cfg.info_trigger   = input_data.cfg.info_trigger;
BLC_data.behav              = input_data.behav;
BLC_data.stim               = input_data.stim;


%3. Also baseline-correct for Gamma-Amplitude
for i_trial = 1:length(input_data.trial)
    for i_chan = 1:length(input_data.label)
        
        %3.1 Demean
        proxy_avgGammaAmp_perTrial = mean(input_data.GammaAmp_Norm{i_trial}(i_chan,:));    
        BLC_data.GammaAmp{i_trial}(i_chan,:) = input_data.GammaAmp_Norm{i_trial}(i_chan,:) - proxy_avgGammaAmp_perTrial; 
      
        proxy_avgLogGammaAmp_perTrial = mean(input_data.LogGammaAmp_Norm{i_trial}(i_chan,:));    
        BLC_data.LogGammaAmp{i_trial}(i_chan,:) = input_data.LogGammaAmp_Norm{i_trial}(i_chan,:) - proxy_avgLogGammaAmp_perTrial;
        
%         figure;
%         plot(preprocData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));
%         hold on;
%         plot(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));
        
        %3.2 BLC
        proxy_avgGammaAmp_perBLW = mean(BLC_data.GammaAmp{i_trial}(i_chan,BLstartsample:BLstopsample));    
        BLC_data.GammaAmp{i_trial}(i_chan,:) = BLC_data.GammaAmp{i_trial}(i_chan,:) - proxy_avgGammaAmp_perBLW;
        
%         hold on;
%         plot(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));
%         legend('original','demeaned','prestim BLC')
%         lineBL_pos = zeros(1,length(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:)));
%         lineBL_pos(BLstartsample:BLstopsample) = max(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));         
% %         lineBL_avgAmp = zeros(1,length(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:)));
% %         plot(lineBL_pos,'k','LineWidth',3)
%         lineBL_avgAmp(BLstartsample:BLstopsample) = proxy_avgGammaAmp_perBLW;
%          plot(lineBL_avgAmp,'g','LineWidth',3)
       
        proxy_avgLogGammaAmp_perBLW = mean(BLC_data.LogGammaAmp{i_trial}(i_chan,BLstartsample:BLstopsample));    
        BLC_data.LogGammaAmp{i_trial}(i_chan,:) = BLC_data.LogGammaAmp{i_trial}(i_chan,:) - proxy_avgLogGammaAmp_perBLW;       

    end
end
 
end