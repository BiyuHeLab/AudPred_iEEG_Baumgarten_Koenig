function data_ECoGfiltref_blocks = sfa_expt4_ECoGrawdata_CommonAvgReference...
    (sub, numstd_thresh, data_ECoGfilt_blocks, subs_PreProcSettings, paths_sfa_expt4, plot_poststepFigs, save_poststepFigs)

%6.1 Determine number of ECoG-recording channels
    num_ECoGelecs = subs_PreProcSettings.(sub).number_ECoGchan; %Determine number of ECoG channels
    
%6.2 Compute variance/std of signal for all channels
%compute std across channels, then select all channels that are > X std
for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)  
    std_perchan(i_block,:) = std(data_ECoGfilt_blocks.trial{i_block}(1:num_ECoGelecs,:),[],2)';
    avgstd_perchan(i_block,:) = mean(std_perchan(i_block,:));
end

%6.3 Detect channels with STD above thresh
%6.3.1 Read out number, index, and amount of supra-thresh chan per block
    index_suprathreshchan_allblocks = [];
    label_suprathreshchan_allblocks = [];
for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)  
    index_suprathreshchan{i_block} = find(std_perchan(i_block,:) > (avgstd_perchan(i_block,:)*numstd_thresh));
    label_suprathreshchan{i_block} = data_ECoGfilt_blocks.label(index_suprathreshchan{i_block});
    
    index_suprathreshchan_allblocks = [index_suprathreshchan_allblocks; index_suprathreshchan{i_block}'];    
    label_suprathreshchan_allblocks = [label_suprathreshchan_allblocks; strcat(label_suprathreshchan{i_block})];
end
%6.3.2 Determine unique above-thresh channels across blocks
    index_suprathreshchan_allblocks = unique(index_suprathreshchan_allblocks);
    label_suprathreshchan_allblocks = unique(label_suprathreshchan_allblocks);
%6.3.3 Determine which of the unique above-thresh channels belong to our
%selected channels
    index_suprathreshchan4selectedchan_allblocks = data_ECoGfilt_blocks.cfg.info_elec.selected.index4EDF...
        (find(ismember(data_ECoGfilt_blocks.cfg.info_elec.selected.index4EDF, index_suprathreshchan_allblocks)));
    label_suprathreshchan4selectedchan_allblocks = data_ECoGfilt_blocks.cfg.info_elec.selected.Label...
        (find(ismember(data_ECoGfilt_blocks.cfg.info_elec.selected.index4EDF, index_suprathreshchan_allblocks)));

    if plot_poststepFigs == 1
        figure
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
        line_vec = ((1:num_ECoGelecs)*0)+1;
        for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo)  
            subplot(length(data_ECoGfilt_blocks.sampleinfo)/4,length(data_ECoGfilt_blocks.sampleinfo)/3,i_block);    
            plot(1:num_ECoGelecs,std_perchan(i_block,:),'b-')
            hold on
            plot(1:num_ECoGelecs,(avgstd_perchan(i_block,:)*numstd_thresh)*line_vec,'k--')
            fig_title = {['Block:' num2str(i_block)], ...
                ['Thresh = ' num2str(numstd_thresh) 'STD - ' num2str(length(index_suprathreshchan{i_block})) ' Chan > Thresh']};
            title(fig_title)

        end
        Figtitle = [sub ' STD per ECoG chan - ' num2str(length(index_suprathreshchan4selectedchan_allblocks)) 'selectChan > Thresh'];
        suptitle(Figtitle)

        if save_poststepFigs == 1            
            path_fig = [paths_sfa_expt4.Fig_ECoGdataraw_SingleSubs sub '/' 'Preproc/CommonAvgRef/' ];
            mkdir(path_fig);

            filename     = strcat([sub '_SDTperChanperBlock_STDThresh' num2str(numstd_thresh) '.png']);            
            figfile      = [path_fig filename];                

            saveas(gcf, [figfile], 'png'); %save png version  
        end
    close
    end

%6.4 remove selected channels from dataset (to use cleaned datset for referencing)
%6.4.1 Create channel_index with only subthresh channels
    channel_selection = 1:num_ECoGelecs;
    channel_selection(index_suprathreshchan_allblocks) = [];
%6.4.2 Create temporary datset with only subthresh channels to compute the
%Common Average Reference
    cfg = [];
    cfg.channel = channel_selection;
    ref_data1 = ft_preprocessing(cfg,data_ECoGfilt_blocks);
    
    num_selectECoGelecs = length(channel_selection); %determine number of selected electrodes

    data_ECoGfiltref_blocks = data_ECoGfilt_blocks; %copy struct for new output file

%6.4.3 For all EEG electrodes (not only STDthresh=selected, but not for DC 
%and Trigger), subtract CAR determined by selected channels for each sample per block
    for i_block = 1:length(ref_data1.sampleinfo)
        Chanavg_data = mean(ref_data1.trial{i_block}); %common average across selected electrodes for each sample per block

        data_ECoGfiltref_blocks.trial{i_block}(1:num_ECoGelecs,:) = ...
            data_ECoGfilt_blocks.trial{i_block}(1:num_ECoGelecs,:) - repmat(Chanavg_data,[num_ECoGelecs 1]);
    end
    
%6.4.4. Compute cross-correlation betwen channels for un-referenced and
%referenced (by mean across subthresh channels) datasets and plot to check 
%dataset-quality
if plot_poststepFigs == 1
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen%
    for i_block = 1:length(ref_data1.sampleinfo)
        %For each block, plot unreferenced and referenced cross-correlation
        %between channels
        subplot(2,length(data_ECoGfilt_blocks.sampleinfo),i_block);
        imagesc(corr(data_ECoGfilt_blocks.trial{i_block}(1:num_selectECoGelecs,:)'),[-1 1]) %correlation across channels for unreferenced data
        grid on
        set(gca,'xtick',length(data_ECoGfilt_blocks.label(1:num_selectECoGelecs,:)));
        set(gca,'ytick',length(data_ECoGfilt_blocks.label(1:num_selectECoGelecs,:)));
        fig_title = {['Block:' num2str(i_block)], ['No CAR']};
        title(fig_title)   

        subplot(2,length(data_ECoGfiltref_blocks.sampleinfo),i_block+length(data_ECoGfiltref_blocks.sampleinfo))
        imagesc(corr(data_ECoGfiltref_blocks.trial{i_block}(1:num_selectECoGelecs,:)'),[-1 1])%correlation across channels for referenced data
        set(gca,'xtick',length(data_ECoGfiltref_blocks.label(1:num_selectECoGelecs,:)));
        set(gca,'ytick',length(data_ECoGfiltref_blocks.label(1:num_selectECoGelecs,:)));
        grid on
        fig_title = {['Block:' num2str(i_block)], ['SelChan CAR']};
        title(fig_title)
    end
        Figtitle = [sub ' CrossCorrMatrix for non-ref and CARref data (CAR based on Chan <' num2str(numstd_thresh) '*STDThresh)'];
        suptitle(Figtitle)

    if save_poststepFigs == 1            
        path_fig = [paths_sfa_expt4.Fig_ECoGdataraw_SingleSubs sub '/' 'Preproc/CommonAvgRef/' ];
        mkdir(path_fig);

        filename     = strcat([sub '_CrossCorrMatrix_NonrefvsCARref_CARchansel' num2str(numstd_thresh) 'STDThresh.png']);            
        figfile      = [path_fig filename];                

        saveas(gcf, [figfile], 'png'); %save png version  
    end
    close
end

%6.5. Add thresh and selected-channel-info to output file
    data_ECoGfiltref_blocks.cfg.info_ref = struct;
    data_ECoGfiltref_blocks.cfg.info_ref.numstd_thresh = numstd_thresh;
    for i_block = 1:length(data_ECoGfilt_blocks.sampleinfo) 
        data_ECoGfiltref_blocks.cfg.info_ref.index_suprathreshchan{i_block} = index_suprathreshchan{i_block};
        data_ECoGfiltref_blocks.cfg.info_ref.label_suprathreshchan{i_block} = label_suprathreshchan{i_block};
    end
    data_ECoGfiltref_blocks.cfg.info_ref.index_suprathreshchan_allblocks = index_suprathreshchan_allblocks;
    data_ECoGfiltref_blocks.cfg.info_ref.label_suprathreshchan_allblocks = label_suprathreshchan_allblocks;
    data_ECoGfiltref_blocks.cfg.info_ref.index_suprathreshchan4selectedchan_allblocks = index_suprathreshchan4selectedchan_allblocks;
    data_ECoGfiltref_blocks.cfg.info_ref.label_suprathreshchan4selectedchan_allblocks = label_suprathreshchan4selectedchan_allblocks;
 
end
