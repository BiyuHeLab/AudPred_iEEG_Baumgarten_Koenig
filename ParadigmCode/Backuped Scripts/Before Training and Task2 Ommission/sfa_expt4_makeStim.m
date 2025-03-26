function stim = sfa_expt4_makeStim(param, pahandle)

% addpath(genpath('C:\Users\localbiyu\Desktop\Brian\sfa\sfa_expt4\'))
% addpath(genpath('C:/Users/localbiyu/Desktop/Thomas/Experiments/Auditory History Tracking/sfa_expt4_Adeen_July18/')); %Path experimental laptop
% addpath(genpath('\\gogo.sb.nyumc.org\data\gogodisk4\thomas\AuditoryPrediction\iEEG\sfa_expt4_Adeen_Current\')); %Path desktop/gogo

load series_selection_scale2Hz_jenn
%loads series_selection_scale2Hz_jenn.mat from subfolder stim_creation

%% 1)  Gather together the stimulus variables and create stimulus var

y = 1:30;
seqTable = reshape(y, [3 5 2]); %Create matrix determining following unique condition combinations:
%3 math. expected final tone pitch (p*34: low, middle, high)
%5 unique sequences
%2 tone duration conditions

i_b = 3; %presumably index for beta selection

%Create array for trial-wise p34 determination
% tf = repmat([1:6]',20,1); %120 trials total, varying between 1:6 (p34/ presented final tone pitch)
%  tf = repmat([1:4]',30,1); %120 trials total, varying between 1:4 (p34/ presented final tone pitch)
                %critical: every unique trial repeats every 30 repetitions;
                %the output of logf_semitone (1:4) thus need to be
                %allocated to trial (e.g.,) 1, 31, 61, 91, we thus create
                %tf as a 30*1/30*2/30*3/30*4 array
  tf = 1:120; 
  tf(1:30) = 1;
  tf(31:60) = 2;
  tf(61:90) = 3;
  tf(91:120) =4;
  tf = tf';

ind = 0;
for i_block = 1:4 %4 blocks
for i_toneDur = 1:2 %2tone durations
    for i_e = 1:3 %trend (high, med, low)
        for i_w = 1:5 %5 unique sequences
                ind = ind + 1;

                i_tf = tf(ind); %select trial-wise p34 determinant
                
                seqID(ind) = seqTable(i_e, i_w, i_toneDur); %select unique p*34/sequence/tonedur combination
                trialID(ind) = ind; %trial number serving as trial ID

                stim.beta(ind)    = 1.5; %fixed beta 
                stim.betaID(ind)  = 3;
                
                stim.uSeqID(ind) = i_e*5-5 + i_w; %determine unique sequence ID (i.e., sequence independent on p34)

                stim.toneDur(ind) = param.toneDur_inSecs(i_toneDur); %determine tone duration in s
                stim.toneDurID(ind) = i_toneDur;

                stim.logf_pen(ind)   = series_logHz{i_b}(i_e*5-5 + i_w,end); %read out logf of penultimate tone from tone data matrix (always 440Hz)
                %i_e*5-5 + i_w determines row (sequence), column (tone) always #33   
                %i_e*5-5 + i_w = 1/6/11 for i_e =1/2/3, i_w =1; 2/7/12 for i_e =1/2/3 for
                %i_w =2...
                
%                 %Original determination of final tone pitches for 6 possible final tone pitches
%                 %4,8,or 12 semi tones above/below 440 Hz
%                     stim.logf_final(ind) = series_logHz_tf_Hz{i_b}(i_e*5-5 + i_w, i_tf); %read out logf of presented final tone (p34) (220-880Hz)
%                     %i_e*5-5 + i_w = row doesn't matter, all same; column
%                     %determines distance from 440Hz - i.e., 4,8,r 12 semi tones above/below 440 Hz
%                     %here 220 Hz, 277 Hz, 349 Hz, 554 Hz, 698 Hz, 880 Hz
%                     switch i_tf
%                         case 1, stim.finalID(ind) = -3; % lowest value for final tone
%                         case 2, stim.finalID(ind) = -2;
%                         case 3, stim.finalID(ind) = -1;
%                         case 4, stim.finalID(ind) =  1;
%                         case 5, stim.finalID(ind) =  2;
%                         case 6, stim.finalID(ind) =  3; % highest value for final tone
%                     end
%                 
%                 stim.series_f{ind}  = [series_Hz{i_b}(i_e*5-5 + i_w,:) series_Hz_tf_Hz{i_b}(i_e*5-5 + i_w, i_tf)];
    
                %New determination of final tone pitches for 4 possible final tone pitches
                %6 or 12 semi tones above/below 440 Hz
                %TJB: calculation of semitones:
                base_freq = 440; %Hz
                semitones = [-12 -6 6 12];
                for i_semitone = 1:length(semitones)
                f_semitone(i_semitone) = base_freq *2^(semitones(i_semitone)/12);
                logf_semitone(i_semitone) = log(f_semitone(i_semitone));                
                end
                

                stim.logf_final(ind) = logf_semitone(i_tf); %read out logf of presented final tone (p34) (220-880Hz)

                    switch i_tf
                        case 1, stim.finalID(ind) = -2; % lowest value for final tone
                        case 2, stim.finalID(ind) = -1;
                        case 3, stim.finalID(ind) =  1;
                        case 4, stim.finalID(ind) =  2; % highest value for final tone
                    end
                    
                stim.series_f{ind}  = [series_Hz{i_b}(i_e*5-5 + i_w,:) f_semitone(i_tf)]; 
                %determine frequencies for tone sequence by combining
                %frequency values for first 33 tones from loaded matrix and
                %newly detremined frequecnies for final tone pitch/p34


                %Determine predicted freq of final tone (i.e., freq of p*34)
                stim.logf_pred(ind)  = series_pred_logHz{i_b}(i_e*5-5 + i_w);
                switch i_e
                    case 1, stim.predID(ind) = -1; % min predicted frequency
                    case 2, stim.predID(ind) = +1; % max predicted frequency
                    case 3, stim.predID(ind) =  0; % median predicted frequency
                end

                stim.soundwave{ind} = series2soundwave(stim.series_f{ind}, param.toneDur_inSecs(i_toneDur), param.f_sample_inHz);

            end
        end
    end
end

%% 2) Randomize trial/stimulus order
% if ~param.restarting
%     rand_ind = randperm(ind);
% else
%     rand_ind = param.stim_ind_order;
% end

% randomize the stimuli, with the constraint that each block holds one
% instance of beta (3) x p*34 (3) x toneDur (3)

if isfield(param, 'stim_ind_order')
    rand_ind = param.stim_ind_order;
    
else
    seqID2   = seqID;
    trialID2 = trialID;
    rand_ind = [];
    rand_seq_ind = [];
    for i_block = 1:4

        block_trial_ind = [];
        for i_trial = 1:30

            found = 0;
            while ~found %continue to shuffle until found condition below is met
                [seqID_r, ind] = Shuffle(seqID2); %shuffle seqID and keep their original index
                trialID_r = trialID2(ind); %reorder trial number according to shuffled index

                if seqID_r(1) == i_trial %if seqID matches with i_trial (i.e., at some random point within shuffling)
                    found = 1;
                    block_trial_ind(i_trial) = trialID_r(1); %copy first trial ID to block-specific trial ID
                    block_seq_ind(i_trial) = seqID_r(1);%label trials within block in ascending order (1:30)

                    trialID2 = trialID_r(2:end); %put remaining trials back into pool to reshuffle
                    seqID2 = seqID_r(2:end);
                end
            end

        end

        [block_trial_ind, ind] = Shuffle(block_trial_ind); %shuffle trial again within block
        block_seq_ind = block_seq_ind(ind); %reorder trials within block according to shuffle

        rand_ind = [rand_ind  block_trial_ind]; %fuse indices over block again
        rand_seq_ind = [rand_seq_ind  block_seq_ind];
    end

end

stim.uSeqID     = stim.uSeqID(rand_ind);
stim.beta       = stim.beta(rand_ind);
stim.betaID     = stim.betaID(rand_ind);
stim.toneDur    = stim.toneDur(rand_ind);
stim.toneDurID  = stim.toneDurID(rand_ind);
stim.logf_pen   = stim.logf_pen(rand_ind);
stim.logf_final = stim.logf_final(rand_ind);
stim.finalID    = stim.finalID(rand_ind);
stim.logf_pred  = stim.logf_pred(rand_ind);
stim.predID     = stim.predID(rand_ind);
stim.series_f   = stim.series_f(rand_ind);
stim.soundwave  = stim.soundwave(rand_ind);

stim.nTrials    = length(stim.betaID);
stim.nBlocks    = stim.nTrials / param.nTrialsPerBlock;
stim.pahandle   = pahandle;
stim.ind_order  = rand_ind;

%% Check trial numbers
%I.e., check if for each uSeq/p*34/tonedur combi (30) 4 trials are present
%and these 4 trials have 4 different final tone pitches (p34)
   
for i_uSeq = 1:15 %15 unique sequences (3 p*34 * 5 uSeq)
    for i_tonedur = 1:2 %2 tone durations   
            f_uSeq = stim.uSeqID == i_uSeq; 
            f_tonedur = stim.toneDurID == i_tonedur; 
            f_comb = f_uSeq & f_tonedur;
            
            check_balance_length = length(stim.finalID(f_comb)) == 4;
            if ~check_balance_length
                disp(['Error in stim generation for i_uSeq ', num2str(i_uSeq), ' and tonedur ' num2str(i_tonedur)])
                pause
            end
            check_balance_sum = sum(stim.finalID(f_comb)) == 0;
            if ~check_balance_sum
                disp(['Error in stim generation for i_uSeq ', num2str(i_uSeq), ' and tonedur ' num2str(i_tonedur)])
                pause
            end
    end
end

end