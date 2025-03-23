function data_ECoGfiltref_blocks = NASTD_ECoG_Preproc_CommonAvgRef...
    (sub, numstd_thresh, data_ECoGfilt_blocks, subs_PreProcSettings, ...
    paths_NASTD_ECoG, plot_poststepFigs, save_poststepFigs)

%Aim: Re-reference grid & strip electrodes using a common-average reference
%computed across all clean grid & strip electrodes

%1. Determine number of grid+strip ECoG-recording channels and compute
%filter that reads out selected electrodes from all electrodes
num_selECoGelecs    = length(data_ECoGfilt_blocks.cfg.info_elec.selected.Label);
filter_selECoGelecs = zeros(1,length(data_ECoGfilt_blocks.label));
filter_selECoGelecs(data_ECoGfilt_blocks.cfg.info_elec.selected.index4EDF') = 1;
filter_selECoGelecs = logical(filter_selECoGelecs);
index_selECoGelecs  = data_ECoGfilt_blocks.cfg.info_elec.selected.index4EDF';
label_selECoGelecs  = data_ECoGfilt_blocks.cfg.info_elec.selected.Label;

%2. Compute std across samples for each selected channel, then average
%across std across channels
for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)
    std_perchan(i_block,:) = std(data_ECoGfilt_blocks.trial{i_block}...
        (filter_selECoGelecs,:),[],2)';
    avgstd_perchan(i_block,:) = mean(std_perchan(i_block,:));
end

%3. Detect channels with STD above thresh
%3.1 Read out number, index, and amount of supra-thresh chan per block
index_suprathreshchan_allblocks = [];
label_suprathreshchan_allblocks = [];
for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)
    index_suprathreshchan{i_block} = ...
        find(std_perchan(i_block,:) > (avgstd_perchan(i_block,:)*numstd_thresh));
    label_suprathreshchan{i_block} = ...
        label_selECoGelecs(index_suprathreshchan{i_block});
    
%     index_suprathreshchan_allblocks = ...
%         [index_suprathreshchan_allblocks; index_suprathreshchan{i_block}'];
    label_suprathreshchan_allblocks = ...
        [label_suprathreshchan_allblocks; strcat(label_suprathreshchan{i_block})];
end

%3.2 Determine unique above-thresh channels across blocks
%and also remove electrodes that are rejected at the end of preprocessing
label_suprathreshchan_allblocks = ...
    unique([subs_PreProcSettings.(sub).rejectedChan_label'; ...
    unique(label_suprathreshchan_allblocks)]);

index_suprathreshchan_allblocks = [];
for i_elec = 1:length(label_suprathreshchan_allblocks)
    index_suprathreshchan_allblocks = ...
        [index_suprathreshchan_allblocks; ...
        find(strcmp(label_selECoGelecs,label_suprathreshchan_allblocks{i_elec}))];
end
index_suprathreshchan_allblocks = sort(index_suprathreshchan_allblocks);

% index_suprathreshchan_allblocks = unique(index_suprathreshchan_allblocks); %only superthresh elecs
% label_suprathreshchan_allblocks = unique(label_suprathreshchan_allblocks);


if plot_poststepFigs == 1
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
    line_vec = ((1:num_selECoGelecs)*0)+1;
    for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)
        subplot(length(data_ECoGfilt_blocks.sampleinfo)/4,...
            length(data_ECoGfilt_blocks.sampleinfo)/3,i_block);
        
        plot(1:num_selECoGelecs,std_perchan(i_block,:),'b-')
        hold on
        plot(1:num_selECoGelecs,(avgstd_perchan(i_block,:)*numstd_thresh)*line_vec,'k--')
        fig_title = {['Block:' num2str(i_block)], ...
            ['Thresh = ' num2str(numstd_thresh) 'STD - '...
            num2str(length(index_suprathreshchan{i_block})) ' Chan > Thresh']};
        title(fig_title)
        
    end
    Figtitle = [sub ' STD per ECoG chan - ' ...
        num2str(length(index_suprathreshchan_allblocks)) 'selectChan > Thresh'];
    sgtitle(Figtitle)
    
    if save_poststepFigs == 1
        path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/Preproc/CommonAvgRef/' ];
        mkdir(path_fig);
        
        filename     = strcat([sub '_SDTperChanperBlock_STDThresh' num2str(numstd_thresh) '.png']);
        figfile      = [path_fig filename];
        
        saveas(gcf, [figfile], 'png'); %save png version
    end
    close
end

%4. remove selected channels from dataset (to use cleaned dataset for referencing)
%4.1 Create channel_index with only subthresh channels
channel_selection_label = label_selECoGelecs;
channel_selection_label(index_suprathreshchan_allblocks) = [];

%4.2 Create temporary datset with only subthresh channels to compute the
%Common Average Reference
cfg = [];
cfg.channel = channel_selection_label;
ref_data1 = ft_preprocessing(cfg,data_ECoGfilt_blocks);

data_ECoGfiltref_blocks = data_ECoGfilt_blocks; %copy struct for new output file

%4.3 For all strip & grid EEG electrodes (not only STDthresh = selected, but not for DC
%and Trigger), subtract CAR determined by selected channels for each sample per block
for i_block = 1:length(ref_data1.sampleinfo)
    
    Chanavg_data = mean(ref_data1.trial{i_block}); 
    %common average across selected electrodes for each sample per block
    
    data_ECoGfiltref_blocks.trial{i_block}(filter_selECoGelecs,:) = ...
        data_ECoGfilt_blocks.trial{i_block}(filter_selECoGelecs,:) ...
        - repmat(Chanavg_data,[num_selECoGelecs 1]);
end

%4.4. Compute cross-correlation betwen channels for un-referenced and
%referenced (by mean across subthresh channels) datasets and plot to check
%dataset-quality
if plot_poststepFigs == 1
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
    for i_block = 1:length(ref_data1.sampleinfo)
        %For each block, plot unreferenced and referenced cross-correlation
        %between channels
        subplot(2,length(data_ECoGfilt_blocks.sampleinfo),i_block);
        imagesc(corr(data_ECoGfilt_blocks.trial{i_block}(filter_selECoGelecs,:)'),[-1 1]) 
        %correlation across channels for unreferenced data
        if i_block == length(ref_data1.sampleinfo)
            colorbar
        end
        grid on
        set(gca,'xtick',length(data_ECoGfilt_blocks.label(filter_selECoGelecs,:)));
        set(gca,'ytick',length(data_ECoGfilt_blocks.label(filter_selECoGelecs,:)));
        fig_title = {['Block:' num2str(i_block)], ['No CAR']};
        title(fig_title)
        
        subplot(2,length(data_ECoGfiltref_blocks.sampleinfo),i_block+length(data_ECoGfiltref_blocks.sampleinfo))
        imagesc(corr(data_ECoGfiltref_blocks.trial{i_block}(filter_selECoGelecs,:)'),[-1 1])
        %correlation across channels for referenced data
        set(gca,'xtick',length(data_ECoGfiltref_blocks.label(filter_selECoGelecs,:)));
        set(gca,'ytick',length(data_ECoGfiltref_blocks.label(filter_selECoGelecs,:)));
        if i_block == length(ref_data1.sampleinfo)
            colorbar
        end        
        grid on
        fig_title = {['Block:' num2str(i_block)], ['SelChan CAR']};
        title(fig_title)
    end
    Figtitle = [sub ' CrossCorrMatrix for non-ref and CARref data (CAR based on Chan <' num2str(numstd_thresh) '*STDThresh)'];
    sgtitle(Figtitle)
    
    if save_poststepFigs == 1
        path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub '/Figs/Preproc/CommonAvgRef/' ];
        mkdir(path_fig);
        
        filename     = strcat([sub '_CrossCorrMatrix_NonrefvsCARref_CARchansel' num2str(numstd_thresh) 'STDThresh.png']);
        figfile      = [path_fig filename];
        
        saveas(gcf, [figfile], 'png'); %save png version
    end
    close
end

%5. Add thresh and selected-channel-info to output file
data_ECoGfiltref_blocks.cfg.info_ref = struct;
data_ECoGfiltref_blocks.cfg.info_ref.numstd_thresh = numstd_thresh;

for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)
    data_ECoGfiltref_blocks.cfg.info_ref.index_suprathreshchan{i_block} ...
        = index_selECoGelecs(index_suprathreshchan{i_block});
    %Note: Here we transfer indeces based on only Strip&Grid electrodes
    %back to all electrodes, as saved in preproc data
    data_ECoGfiltref_blocks.cfg.info_ref.label_suprathreshchan{i_block} ...
        = label_suprathreshchan{i_block};
end

data_ECoGfiltref_blocks.cfg.info_ref.index_suprathreshchan_allblocks...
    = index_selECoGelecs(index_suprathreshchan_allblocks);
data_ECoGfiltref_blocks.cfg.info_ref.label_suprathreshchan_allblocks...
    = label_suprathreshchan_allblocks;

end