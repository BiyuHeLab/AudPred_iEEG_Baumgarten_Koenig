function BLCData_perTD = NASTD_ECoG_TLA_BLCdataperTD(BLwin, LPfreq, preprocData_perTD)

%0. Find samples coresponding to BLwindow markers
Dist = abs(preprocData_perTD.time{1} - BLwin(1));
minDist = min(Dist);
BLstartsample = find(Dist == minDist);
clear Dist minDist

Dist = abs(preprocData_perTD.time{1} - BLwin(2));
minDist = min(Dist);
BLstopsample = find(Dist == minDist);
clear Dist minDist 

%1. Demean and baseline-correct (500 ms prestim) for low freqs
cfg = [];
cfg.demean = 'yes';

cfg.baselinewindow = [BLwin(1) BLwin(2)];
disp(['Baseline Correction for :' num2str(BLwin(1)) ' to ' num2str(BLwin(2)) ' sec'])

% cfg.hpfilter = 'yes';
% cfg.hpfreq = 0.1;

cfg.lpfilter = 'yes';
cfg.lpfreq = LPfreq;

BLCData_perTD = ft_preprocessing(cfg,preprocData_perTD);

%2. Copy special subfields
BLCData_perTD.cfg.info_elec = preprocData_perTD.cfg.info_elec;
BLCData_perTD.cfg.info_ref = preprocData_perTD.cfg.info_ref;
BLCData_perTD.cfg.info_trigger = preprocData_perTD.cfg.info_trigger;
BLCData_perTD.behav = preprocData_perTD.behav;

%3. Also baseline-correct for Gamma-Amplitude
for i_trial = 1:length(preprocData_perTD.trial)
    for i_chan = 1:length(preprocData_perTD.label)
        
        %3.1 Demean
        proxy_avgGammaAmp_perTrial = mean(preprocData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));    
        BLCData_perTD.GammaAmp{i_trial}(i_chan,:) = preprocData_perTD.GammaAmp_Norm{i_trial}(i_chan,:) - proxy_avgGammaAmp_perTrial; 
      
        proxy_avgLogGammaAmp_perTrial = mean(preprocData_perTD.LogGammaAmp_Norm{i_trial}(i_chan,:));    
        BLCData_perTD.LogGammaAmp{i_trial}(i_chan,:) = preprocData_perTD.LogGammaAmp_Norm{i_trial}(i_chan,:) - proxy_avgLogGammaAmp_perTrial;
        
%         figure;
%         plot(preprocData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));
%         hold on;
%         plot(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));
        
        %3.2 BLC
        proxy_avgGammaAmp_perBLW = mean(BLCData_perTD.GammaAmp{i_trial}(i_chan,BLstartsample:BLstopsample));    
        BLCData_perTD.GammaAmp{i_trial}(i_chan,:) = BLCData_perTD.GammaAmp{i_trial}(i_chan,:) - proxy_avgGammaAmp_perBLW;
        
%         hold on;
%         plot(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));
%         legend('original','demeaned','prestim BLC')
%         lineBL_pos = zeros(1,length(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:)));
%         lineBL_pos(BLstartsample:BLstopsample) = max(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:));         
% %         lineBL_avgAmp = zeros(1,length(BLCData_perTD.GammaAmp_Norm{i_trial}(i_chan,:)));
% %         plot(lineBL_pos,'k','LineWidth',3)
%         lineBL_avgAmp(BLstartsample:BLstopsample) = proxy_avgGammaAmp_perBLW;
%          plot(lineBL_avgAmp,'g','LineWidth',3)
       
        proxy_avgLogGammaAmp_perBLW = mean(BLCData_perTD.LogGammaAmp{i_trial}(i_chan,BLstartsample:BLstopsample));    
        BLCData_perTD.LogGammaAmp{i_trial}(i_chan,:) = BLCData_perTD.LogGammaAmp{i_trial}(i_chan,:) - proxy_avgLogGammaAmp_perBLW;
        


    end
end
 
end