function NASTD_ECoG_Predict_ComparePredEffects_TWvsSample...
    (sub, FuncInput_DataType, FuncInput_ToneDur_text,  ...
    param,...
    FuncInput_InputData, ...
    save_poststepFigs, paths_NASTD_ECoG)

%Aims:
%Plot neural signal during p33 as a function of predp34 and highlight both
%sign. time windows determined via TW-approach and sign. samples determined
%via sample-based approach.

%1. Load in TW and sample base data
%% 0.1) Specify vars, paths, and setup fieldtrip

%Add base dir and own script dir
addpath(genpath(paths_NASTD_ECoG.BaseDir));
addpath(genpath(paths_NASTD_ECoG.ScriptsDir));

path_inputdata_TW = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/' param.Label_TW 'sTW/'];
path_inputdata_Sample = [paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Data/Samplewise/'];

path_fig = ([paths_NASTD_ECoG.ECoGdata_Prediction '/PredEffects/' ...
    sub '/Figs/ERFs/PredEffect/' FuncInput_DataType '/']);
if (~exist(path_fig, 'dir')); mkdir(path_fig); end

SampleFreq  = FuncInput_InputData.fsample;
nSensors    = size(FuncInput_InputData.trial{1},1);

%% 1. Load in data

%1.1 TW data
load([path_inputdata_TW sub '_PredEffects_' FuncInput_DataType '_' FuncInput_ToneDur_text  'sTD.mat'], ...
    'PredEffect');
PredEffect_TW = PredEffect;
clear PredEffect

%1.2 Sample-wise data
load([path_inputdata_Sample sub '_PredEffectsCluster_' FuncInput_DataType '_' FuncInput_ToneDur_text 'sTD.mat'], ...
    'PredEffect', 'Data_p33');
PredEffect_Sample = PredEffect;
clear PredEffect

%% 2. Separate trials into predp34-conditions
label_Predp34 = [-1 0 1];% low, medium, high;
for i_Predp34 = 1:length(label_Predp34)
    TrialFilt_Predp34 = FuncInput_InputData.behav.stim.predID == label_Predp34(i_Predp34);
    Index_Predp34{i_Predp34} = find(TrialFilt_Predp34 == 1);
end
clear TrialFilt_Predp34

%% 3. Define TWs
win_size    = param.SamplesTW; 
win_overlap = 0;
ToneDur_Sec  = str2num(FuncInput_ToneDur_text);
nSamples_perTone = ToneDur_Sec * SampleFreq;

windows = [1 win_size];
while windows(end,end) < min(nSamples_perTone)
    windows = [windows; windows(end,:) + (win_size - win_overlap)];
end
if windows(end,end) > min(nSamples_perTone) 
    windows(end,:) = [];
end

%% 4. Plot neural signal during p33
%Compute AVG + STD across trials per Predp34
for i_Predp34 = 1:length(label_Predp34)
    Avgp33{i_Predp34} = mean(Data_p33(:,:,Index_Predp34{i_Predp34}),3);
    STDp33{i_Predp34} = std(Data_p33(:,:,Index_Predp34{i_Predp34}),1,3);
end

%Plot ERF for all selected elecs
NumFig = ceil(nSensors/30); %20 elecs per plot
plot_param.color = ...
    {[0, 0.4470, 0.7410, 0.7], ...
    [0, 0.75, 0.75, 0.7],[0.8500, ...
    0.3250, 0.0980, 0.7]};

%ERF for p33
i_elec = 0;
for i_fig = 1:NumFig
    h = figure('visible','on'); %ensures figure doesn't pup during plotting
%     h = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    for i_subplot = 1:30
        if i_elec < nSensors
            subplot(5,6,i_subplot)
            i_elec = i_elec+1;
            
            %Plot signal traces
            plot(1:length(Avgp33{1}(i_elec,:)),Avgp33{1}(i_elec,:),...
                'color',plot_param.color{1},'LineWidth', 2)
            hold on;
            plot(1:length(Avgp33{2}(i_elec,:)),Avgp33{2}(i_elec,:),...
                'color',plot_param.color{2},'LineWidth', 2)
            hold on;
            plot(1:length(Avgp33{3}(i_elec,:)),Avgp33{3}(i_elec,:),...
                'color',plot_param.color{3},'LineWidth', 2)
            axis tight
            a = axis;
            xlabel('samples (p33)')
            
            %Highlight sign. samples/clusters
            maxY = max([max(Avgp33{1}(i_elec,:)), max(Avgp33{2}(i_elec,:)), max(Avgp33{3}(i_elec,:))]);
            minY = min([min(Avgp33{1}(i_elec,:)), min(Avgp33{2}(i_elec,:)), min(Avgp33{3}(i_elec,:))]);
            if ~isnan(PredEffect_Sample.clusterstat{i_elec}.cluster_pval)
                for i_cluster = 1:length(PredEffect_Sample.clusterstat{i_elec}.cluster_pval)
                    if PredEffect_Sample.clusterstat{i_elec}.cluster_pval(i_cluster) < param.alpha
                        sign_samples = ...
                            find(PredEffect_Sample.clusterstat{i_elec}.cluster_timecourse ...
                            == i_cluster);
                        plot(sign_samples, ...
                            ones(length(sign_samples)) * max([minY maxY])+abs(max([minY maxY])*0.05),...
                            'linewidth',7,'Color',[175 175 175]./255)
                        text(sign_samples(1),max([minY maxY])+abs(max([minY maxY])*0.05), ...
                            [num2str(round(...
                            PredEffect_Sample.clusterstat{i_elec}.cluster_pval(i_cluster),3))], ...
                            'FontSize',6);
                    end
                end
            end
            
            %Highlight sign. TWs
            signTW = [];
            for i_win = 1:size(windows,1)
                %Add demarcation lines for TW (50ms)
                ind_TWstart = windows(i_win, 1);
                ind_TWtop = windows(i_win, 2);
                hold on;
                plot([ind_TWstart ind_TWstart],...
                     [min([minY maxY])-+abs(max([minY maxY])*0.35) ...
                    max([minY maxY])+abs(max([minY maxY])*0.35)], ...
                    'k--','LineWidth', 0.2)
                
                %Add shading
                if PredEffect_TW.pval{1}{i_win}(i_elec) < param.alpha
                    signTW = [signTW 1];
                    hold on;
                    plot(ind_TWstart:ind_TWtop, ...
                        ones(length(ind_TWstart:ind_TWtop)) * max([minY maxY])+abs(max([minY maxY])*0.35),...
                        'linewidth',7,'Color',[150 150 150]./255)
                    text(ind_TWstart,max([minY maxY])+abs(max([minY maxY])*0.35), ...
                        [num2str(round(PredEffect_TW.pval{1}{i_win}(i_elec),3))], ...
                        'FontSize',6);
                else
                    signTW = [signTW 0];                  
                end
            end
            
            %Scale y-axis to show sign. text if present
            if any(PredEffect_Sample.clusterstat{i_elec}.cluster_pval < param.alpha) || sum(signTW) > 0
                ylim([min([minY maxY])-+abs(max([minY maxY])*0.4) ...
                    max([minY maxY])+abs(max([minY maxY])*0.4)]); %With space for sign. info 
            else
                ylim([min([minY maxY]) max([minY maxY])]); %Without space for sign. info
            end
            
            %Title
            title(['Elec: ' FuncInput_InputData.label{i_elec}])            
        end
    end
    
    %Add subtitle
    sgtitle({[sub ' - ' FuncInput_DataType ' - ' FuncInput_ToneDur_text ...
        'msTD'], ['Pred Effect (p33 ERF as a function of p*34, p < ' num2str(param.alpha) '), [' ...
        num2str(i_fig) '/' num2str(NumFig) ']']},'Interpreter','none')
    
    if save_poststepFigs == 1
        filename     = ['p33ERF_AllElecs_TWvsClusterStat_' sub '_' FuncInput_DataType '_' ...
            FuncInput_ToneDur_text 'sTD_' num2str(i_fig) '.png'];
        figfile      = [path_fig filename];
        saveas(gcf, [figfile], 'png'); %save png version
        close
    end
    
end

end