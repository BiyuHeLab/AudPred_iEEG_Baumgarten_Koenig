function stim = sfa_expt4_makeTrainingStim(param, pahandle)

% addpath(genpath('C:\Users\localbiyu\Desktop\Brian\sfa\sfa_expt4\'))
% addpath(genpath('C:/Users/localbiyu/Desktop/Thomas/Experiments/Auditory History Tracking/sfa_expt4_Adeen_July18/')); %Path experimental laptop
% addpath(genpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\AuditoryPrediction\iEEG\sfa_expt4_Adeen_Current\')); %Path desktop/gogo
% addpath(genpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\NaturalisticAuditorySequences_ToneDuration(NAS_TD)\ECoG\Experiment_Versions\sfa_expt4_Adeen_Current\stim_creation\TestThomas_CreateTrainingStim\'));
addpath('/isilon/LFMI/VMdrive/Thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/ECoG/Experiment_Versions/sfa_expt4_Adeen_Current/stim_creation/TestThomas_CreateTrainingStim/')
load('seriesselection4training_9x2seqs_3pertrend_beta1.5_TJB.mat'); %file with training sequences

%% 1)  Gather together the stimulus variables and create stimulus var

y = 1:18; %individual sequences
seqTable = reshape(y, [3 3 2]); %Create matrix determining following unique condition combinations:
%3 math. expected final tone pitch (p*34: low, middle, high)
%2 unique sequences
%1 tone duration conditions

i_b = 3; %presumably index for beta selection

%Create array for trial-wise p34 determination
%1 = 5.3936; 2 = 5.7402; 3 = 6.4333; 4 = 6.7799 logHz
  tf = 1:18; 
  tf(1) = 2; %trial 1: low p*34(logHz) = 5.9135; p34(logHz) = 5.7402; closest - highly likely
  tf(2) = 4; %trial 2: low p*34(logHz) = 5.9713; p34(logHz) = 6.7799; furthest - highly UNlikely
  tf(3) = 2; %trial 3: low p*34(logHz) = 5.9135; p34(logHz) = 5.7402; closest - highly likely
  tf(4) = 2; %trial 4: medium p*34(logHz) = 6.0868; p34(logHz) = 5.7402; closest - highly likely
  tf(5) = 2; %trial 5: medium p*34(logHz) = 6.0868; p34(logHz) = 5.7402; closest - highly likely
  tf(6) = 4; %trial 6: medium p*34(logHz) = 6.0868; p34(logHz) = 6.7799; furthest - highly UNlikely
  tf(7) = 3; %trial 7: high p*34(logHz) = 6.2601; p34(logHz) = 6.4333; closest - highly likely
  tf(8) = 3; %trial 8: high p*34(logHz) = 6.2023; p34(logHz) = 6.4333; closest - highly likely
  tf(9) = 1; %trial 9: high p*34(logHz) = 6.2601; p34(logHz) = 5.3936; furthest - highly UNlikely

  tf(10) = 2; %trial 10: low p*34(logHz) = 5.9135; p34(logHz) = 5.7402; closest - highly likely
  tf(11) = 2; %trial 11: low p*34(logHz) = 5.9713; p34(logHz) = 5.7402; closest - highly likely
  tf(12) = 4; %trial 12: low p*34(logHz) = 5.9135; p34(logHz) = 6.7799; furthest - highly UNlikely
  tf(13) = 4; %trial 13: medium p*34(logHz) = 6.0868; p34(logHz) = 6.7799; furthest - highly UNlikely
  tf(14) = 2; %trial 14: medium p*34(logHz) = 6.0868; p34(logHz) = 5.7402; closest - highly likely
  tf(15) = 2; %trial 15: medium p*34(logHz) = 6.0868; p34(logHz) = 5.7402; closest - highly likely
  tf(16) = 1; %trial 16: high p*34(logHz) = 6.2601; p34(logHz) = 5.3936; furthest - highly UNlikely
  tf(17) = 3; %trial 17: high p*34(logHz) = 6.2023; p34(logHz) = 6.4333; closest - highly likely
  tf(18) = 3; %trial 18: high p*34(logHz) = 6.2601; p34(logHz) = 6.4333; closest - highly likely

  match = 1:18; %quantifies if p*34 and p34 match; 0 = no match (max p*34-p34 distance), 1 = match (min p*34-p34 distance)
  match(1) = 1;
  match(2) = 0;
  match(3) = 1;
  match(4) = 1;
  match(5) = 1;
  match(6) = 0;
  match(7) = 1;
  match(8) = 1;
  match(9) = 0;
  match(10) = 1;
  match(11) = 1;
  match(12) = 0;
  match(13) = 0;
  match(14) = 1;
  match(15) = 1;
  match(16) = 0;
  match(17) = 1;
  match(18) = 1;
  
ind = 0;
for i_e = 1:3 %trend (high, med, low)
    for i_w = 1:6 %2 unique sequences (per trend)
                ind = ind + 1;

                i_tf = tf(ind); %select trial-wise p34 determinant
                
                seqID(ind) = seqTable(i_e, i_w); %select unique p*34/sequence/tonedur combination
                trialID(ind) = ind; %trial number serving as trial ID

                stim.beta(ind)    = 1.5; %fixed beta 
                stim.betaID(ind)  = 3;
                
                stim.uSeqID(ind) = i_e*6-6 + i_w; %determine unique sequence ID (i.e., sequence independent on p34)

                stim.toneDur(ind) = 0.4; %determine tone duration in s
                stim.toneDurID(ind) = 2;

                stim.logf_pen(ind)   = series_logHz{i_b}(i_e*6-6 + i_w,end); %read out logf of penultimate tone from tone data matrix (always 440Hz)
                
                stim.match(ind) = match(ind);
                
                f_semitone = series_Hz_tf_Hz{3}(1,:);
                logf_semitone = log(f_semitone); 

                stim.logf_final(ind) = logf_semitone(i_tf); %read out logf of presented final tone (p34) (220-880Hz)

                    switch i_tf
                        case 1, stim.finalID(ind) = -2; % for downward trend, choose lowest value for final tone
                        case 2, stim.finalID(ind) = -1; %for medium trend, choosemiddle value for final tone
                        case 3, stim.finalID(ind) =  1; %for medium trend, choosemiddle value for final tone
                        case 4, stim.finalID(ind) =  2; % for upward trend, choose highest value for final tone
                    end
                    
                stim.series_f{ind}  = [series_Hz{i_b}(i_e*6-6 + i_w,:) f_semitone(i_tf)]; 
                %determine frequencies for tone sequence by combining
                %frequency values for first 33 tones from loaded matrix and
                %newly detremined frequecnies for final tone pitch/p34


                %Determine predicted freq of final tone (i.e., freq of p*34)
                stim.logf_pred(ind)  = series_pred_logHz{i_b}(i_e*6-6 + i_w);
                switch i_e
                    case 1, stim.predID(ind) = -1; % min predicted frequency
                    case 2, stim.predID(ind) = 0; % max predicted frequency
                    case 3, stim.predID(ind) =  +1; % median predicted frequency
                end


%                   stim.soundwave{ind} = series2soundwave(stim.series_f{ind}, param.toneDur_inSecs(2), param.f_sample_inHz);
        
        end
    end
                stim.nTrials    = length(stim.betaID);
                stim.nBlocks    = stim.nTrials / param.nTrialsPerBlock;
                stim.nBlocks    = 2;
                stim.pahandle   = pahandle;
                stim.ind_order  = [1:stim.nTrials];


%Plot overview for all training series
k = 33;
figure;
col_down(1:length(seqTable)) = {'c'};
col_med(1:length(seqTable)) = {'m'};
col_up(1:length(seqTable)) = {'y'};
col = [col_down, col_med, col_up col_down, col_med, col_up];
% col = {'r' 'r' 'r' 'r' 'g' 'g' 'g' 'g' 'b' 'b' 'b' 'b' 'b'};
suptitle(['Training Sequences in logHz (all beta = 1.5)'])       

for j = 1:stim.nTrials
subplot((stim.nTrials/3),3,j); hold on;
    hold on
    %plot tone 1:33 of indiv sequence
    plot(1:k, series_logHz{i_b}(j,:), 'k','Marker', '.','MarkerSize',30,'LineWidth',2);
    plot(1:k, series_logHz{i_b}(j,:), [col{j} '-'],'LineWidth',2);
    %plot all possible p34
    plot(k+1, log(series_Hz_tf_Hz{i_b}(1,1)), ['k' 's'], 'MarkerSize',8,'MarkerFaceColor','k');
    plot(k+1, log(series_Hz_tf_Hz{i_b}(1,2)), ['k' 's'], 'MarkerSize',8,'MarkerFaceColor','k');
    plot(k+1, log(series_Hz_tf_Hz{i_b}(1,3)), ['k' 's'], 'MarkerSize',8,'MarkerFaceColor','k');
    plot(k+1, log(series_Hz_tf_Hz{i_b}(1,4)), ['k' 's'], 'MarkerSize',8,'MarkerFaceColor','k');   
    %plot p*34 for this seq
    plot(k+1, series_pred_logHz{i_b}(j), ['b' 'd'],'MarkerSize',10,'MarkerFaceColor','b');
    %plot p34 for this trial
    plot(k+1, log(stim.series_f{j}(end)), ['r' '<'], 'MarkerSize',10, 'MarkerFaceColor','r');    
    %plot middle line
    plot([0 34], [log(440) log(440)], '--', 'color', [0.4 0.4 0.4])
    
    if match(j) == 1
        label_match = 'MIN p*34-p34';
    else
        label_match = 'MAX p*34-p34';
    end
    
    if series_pred_Hz{i_b}(j) < 440
    title({['trial: ' num2str(stim.uSeqID(j)) ', down'],['p*34 = ' num2str(series_pred_logHz{i_b}(j)) 'logHz, p34 =' num2str(log(stim.series_f{j}(end))) 'logHz ' label_match]})
    elseif series_pred_Hz{i_b}(j) == 440
    title({['trial: ' num2str(stim.uSeqID(j)) ', medium'],['p*34 = ' num2str(series_pred_logHz{i_b}(j)) ', p34 =' num2str(log(stim.series_f{j}(end))) 'logHz ' label_match]})
    elseif series_pred_Hz{i_b}(j) > 440
    title({['trial: ' num2str(stim.uSeqID(j)) ', up'],['p*34 = ' num2str(series_pred_logHz{i_b}(j)) ', p34 =' num2str(log(stim.series_f{j}(end))) 'logHz ' label_match]})
    end
end

end