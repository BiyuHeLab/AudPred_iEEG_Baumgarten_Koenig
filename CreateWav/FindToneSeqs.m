%% Read out used tone sequence information %%

%% 1. Set up paths and load data
addpath('/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/'); %main project path from server/gogo

addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/'); %main project path from desktop
addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/'); %main project path from desktop

path_rawdata = '/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/Data/Behavior_rawdata/';
path_rawdata = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/Data/Behavior_rawdata/';

sub = 'S1'; %Doesn't matter since identical sequences for all subjects

load([path_rawdata 'subjects/sfa_expt3_' sub '.mat']); %load single-subject behavioral data

%% 2. Read out individual sequences

index_p34 = [-3 -2 -1 1 2 3]; %p34 options 
index_pred34 = [-1 0 1]; %p*34 options

for i_beta = 1:3 %beta level
    for i_dur = 1:3 %tone dur
        for i_pred = 1:3 %p*34
            for i_final = 1:6 %p34              
                    
                    %set up filter
                    f_beta = stim.betaID == i_beta;
                    f_dur = stim.toneDurID == i_dur;
                    f_pred = stim.predID == index_pred34(i_pred);
                    f_final = stim.finalID == index_p34(i_final);
                    f_all = f_beta&f_dur&f_pred&f_final;
                    
                    %find sequence corresponding to filter
                    inds = find(f_all); 

                    %read out log freq of corresponding sequence
                    log_seq{i_beta,i_dur,i_pred,i_final} = log(stim.series_f{inds(1)});
                        %All individual sequences, ordered by beta level,
                        %p*34 and p34.
                     
                    %create a playable soundwave for each seq    
                    sound_seq{i_beta, i_dur, i_pred,i_final} = series2soundwave(stim.series_f{inds(1)}, stim.toneDur(inds(1)), 44100);
                      
                    %read out log freq of p*34
                    pred34{i_beta,i_dur,i_pred,i_final} = stim.logf_pred(inds(1));
                    
            end
        end
    end
end

%% 3. Play audio of selected sequence
%3.1 Select parameters defining sequence
i_beta = 1;
i_dur = 1;
i_pred = 1; 
i_final = 1;

%3.2 Play slected sequence
play_sec = audioplayer(sound_seq{i_beta, i_dur, i_pred, i_final}, 44100);
play(play_sec)

%3.3. Save selected sequence as .wav
path_save = '//gogo.sb.nyumc.org/data/gogodisk4/thomas/NaturalisticAuditorySequences_ToneDuration(NAS_TD)/MEG/NASTD_MEG_Matlab/';
filename = ['ToneSeq_Beta' num2str(i_beta) '_ToneDur' num2str(i_dur) '_predP34' num2str(i_pred) '_realP34' num2str(i_final) '.wav'];
audiowrite([path_save filename], sound_seq{i_beta, i_dur, i_pred, i_final}, 44100)

%% 4. Plot individual sequences as dot-line combintaion
%Note: Plot agrees with above to-be-played sequence - i_beta defined row,
%i_pred defines column

vec_finaltones = unique(stim.logf_final);
x_fin = [34 34 34 34 34 34];

y_pen = log(440);
x_pen = 33;
y_fin = [5.9135 6.0868 6.2601];

count = 0;
h = figd(26,3,40);
set(gcf,'units','normalized','outerposition',[0 0 1 1])

for i_beta = 1:3
        for i_pred = 1:3

            i_dur = 1; %doesn't mater, since same sequences for all tone dur
            i_final = 1;%doesn't mater, since tone 1-33 same across all p34
            count = count + 1;

            subplot(3,3,count)
            
            plot([0 34], [log(440) log(440)], '--', 'color', [0.4 0.4 0.4])
            hold on
            plot(log_seq{i_beta,i_dur,i_pred,i_final}(1:33), 'Marker', '.', 'Color', 'k')
            hold on
            scatter(x_fin,vec_finaltones, 100,'o', 'filled', 'b');
            hold on
            scatter(x_fin(1),y_fin(i_pred),100, '+', 'r', 'LineWidth', 4)


            xlim([1 34])
            ylim([log(220) log(880)])
            
            ax = gca;
            
            set(gca,'Xtick',[])
            set(gca,'Ytick',[])
            set(gca,'YTickLabel', {'220','440','880'});
            
            switch count
                case 1
                    title('p*_3_4 = low')
                    ylabel('\beta = 0.5')
                    set(gca,'YTick',[log(220), log(440), log(880)]);
                case 2
                    title('p*_3_4 = med')
                case 3
                    title('p*_3_4 = high')
                case 4
                    ylabel('\beta = 1.0')
                    set(gca,'YTick',[log(220), log(440), log(880)]);
                case 7
                    ylabel('\beta = 1.5')
                    set(gca,'YTick',[log(220), log(440), log(880)]);
                    set(gca,'XTick',[10, 20, 30]);
                case 8
                    set(gca,'XTick',[10, 20, 30]);
                case 9
                    set(gca,'XTick',[10, 20, 30]);
            end

    end
end

