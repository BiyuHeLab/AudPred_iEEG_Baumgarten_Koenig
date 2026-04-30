function [trialdata_filtinterp, info_filtinterp] = ...
    NASTD_ECoG_FiltNaNinterp_HP05toLP30Hz...
    (sub, TD_label, ...
    input_data, ...
    plot_poststepFigs, save_poststepFigs, paths_NASTD_ECoG)

%Aim: Filter (0.5-30Hz) data containing interspersed NaNs and interpolate
%NaN values based on filtered neighbours. Circumvents problems with filter
%due to NaN entries.
%Since we can't use NaN entries for filter operations and
%interpolation prior to filtering introduces artifacts, we
%1) cut the signal into uninterrupted segements prior to filtering
%2) compute filters on these segments
%3) inpolate the remaining NaNs using filtered data
    
%1) Set filter parameters
FilterOrder     = 3;
HighPassFreq    = 0.5;
LowPassFreq     = 30;
[lpb_LP,lpa_LP] = butter(FilterOrder,(LowPassFreq*2)/input_data.fsample,'low');
[lpb_HP,lpa_HP] = butter(FilterOrder,(HighPassFreq*2)/input_data.fsample,'high');

min_epochlength = 50; %in samples

%2) Cut data into uninterrupted segements prior to filtering and then
%filter said segments
datastruct = input_data;

for i_trial = 1:length(input_data.trial)    

    clear temp*
    for i_elec = 1:length(input_data.label)
        
        %determine noNaN and NaN samples
        ind_noNaNsamples{i_trial}{i_elec} = ...
            find(~isnan(input_data.trial{i_trial}(i_elec,:)));
        ind_NaNsamples{i_trial}{i_elec} = ...
            find(isnan(input_data.trial{i_trial}(i_elec,:)));
        
        if isempty(ind_NaNsamples{i_trial}{i_elec})
            %if no NaNs are present, filter using the whole trial
            
            datastruct.trial{i_trial}(i_elec,:) = ...
                filtfilt(lpb_LP,lpa_LP,input_data.trial{i_trial}(i_elec,:)')';
            datastruct.trial{i_trial}(i_elec,:) = ...
                filtfilt(lpb_HP,lpa_HP,datastruct.trial{i_trial}(i_elec,:)')';            
            
        else%if NaNs are present,
            iteration_counter = 1;
            while iteration_counter < 3 
                %ensure additional iteration if any too short segments are changed to NaN
                %Set up empty structs
                temp_num_noNaNchanges = []; temp_num_noNaNsamples_epoch = [];
                temp_ind_noNaNchanges = []; temp_ind_noNaNsamples_epoch = [];
                %Redo determination of noNaN and NaN samples
                ind_noNaNsamples{i_trial}{i_elec} = ...
                    find(~isnan(input_data.trial{i_trial}(i_elec,:)));
                ind_NaNsamples{i_trial}{i_elec} = ...
                    find(isnan(input_data.trial{i_trial}(i_elec,:)));
                %determine number and location of each change in noNaN samples
                %to determine where separate noNaNepochs lie
                temp_num_noNaNchanges = ...
                    length(find(ischange(ind_noNaNsamples{i_trial}{i_elec},...
                    'linear','Threshold',2)));
                temp_ind_noNaNchanges = ...
                    find(ischange(ind_noNaNsamples{i_trial}{i_elec},...
                    'linear','Threshold',2));
                
                if temp_num_noNaNchanges == 0 %if there are no changes in noNaN samples = 1 continuous noNaNepoch
                    temp_num_noNaNsamples_epoch = ...
                        length(ind_noNaNsamples{i_trial}{i_elec}(1:end));
                    temp_ind_noNaNsamples_epoch{1} = ...
                        ind_noNaNsamples{i_trial}{i_elec}(1:end);
                else %if there are changes = multiple noNaNepochs
                    %determine length and indices of each non-NaN-epoch to ensure sufficient data for filtering
                    for i_changes = 1:(temp_num_noNaNchanges+1) %+1 because X changes result in X+1 noNaNepochs
                        if i_changes == 1 %first change
                            temp_ind_noNaNsamples_epoch{i_changes} = ...%get samples from start to before first change
                                ind_noNaNsamples{i_trial}{i_elec}...
                                (1:(temp_ind_noNaNchanges(i_changes)-1));
                        elseif i_changes == (temp_num_noNaNchanges+1) %last change
                            temp_ind_noNaNsamples_epoch{i_changes} = ...%get samples from after last change to end
                                ind_noNaNsamples{i_trial}{i_elec}...
                                (temp_ind_noNaNchanges(i_changes-1):end);
                        else %every other change
                            temp_ind_noNaNsamples_epoch{i_changes} = ...%get samples from previous change to before current change
                                ind_noNaNsamples{i_trial}{i_elec}...
                                (temp_ind_noNaNchanges(i_changes-1):temp_ind_noNaNchanges(i_changes)-1);
                        end
                        temp_num_noNaNsamples_epoch(i_changes) = ...
                            length(temp_ind_noNaNsamples_epoch{i_changes});
                    end
                end
                
                %for each noNaNepoch, check if number of samples is too short
                if any(temp_num_noNaNsamples_epoch < min_epochlength)
                    for i_noNaNepoch = 1:length(temp_ind_noNaNsamples_epoch)
                        if length(temp_ind_noNaNsamples_epoch{i_noNaNepoch}) ...
                                < min_epochlength
                            %if it is too short, NaN those samples and restart process
                            input_data.trial{i_trial}...
                                (i_elec,temp_ind_noNaNsamples_epoch{i_noNaNepoch}) ...
                                = nan;
                        end
                    end
                    iteration_counter = iteration_counter + 1; %set final iteration
                else
                    iteration_counter = iteration_counter + 2; %end iteration
                end
            end
            
            %Filter segments
            for i_noNaNepoch = 1:length(temp_ind_noNaNsamples_epoch)        
                temp_filtereddata{i_elec}{i_noNaNepoch}(1,:) = ...
                    filtfilt(lpb_LP,lpa_LP,...
                    input_data.trial{i_trial}...
                    (i_elec,temp_ind_noNaNsamples_epoch{i_noNaNepoch})')';
                temp_filtereddata{i_elec}{i_noNaNepoch}(1,:) = ...
                    filtfilt(lpb_HP,lpa_HP,...
                    temp_filtereddata{i_elec}{i_noNaNepoch}(1,:)')';
                temp_filtereddata{i_elec}{i_noNaNepoch}(2,:) = ...
                    temp_ind_noNaNsamples_epoch{i_noNaNepoch};
            end
            
            %Place filtered data from separate noNaNepochs in correct place in common trial struct
            datastruct.trial{i_trial}(i_elec,:) = ...
                nan(1,length(input_data.trial{i_trial}));
            for i_noNaNepoch = 1:length(temp_filtereddata{i_elec})
                datastruct.trial{i_trial}...
                    (i_elec,temp_ind_noNaNsamples_epoch{i_noNaNepoch}) = ...
                    temp_filtereddata{i_elec}{i_noNaNepoch}(1,:);
            end
        end
    end
end

%3) Compute summary info data about Nan entries
num_NaNelecs_pertrial = zeros(1,length(ind_NaNsamples));
avgNaNsamplenum_pertrial = zeros(1,length(ind_NaNsamples));
perc_NaNsample_pertrial = zeros(1,length(ind_NaNsamples));

for i_trial = 1:length(ind_NaNsamples)
    temp_numNaNsamples = [];
    temp_NaNpresent = 0;
    for i_elec = 1:length(ind_NaNsamples{1})
        temp_numNaNsamples = [temp_numNaNsamples, ...
            length(ind_NaNsamples{i_trial}{i_elec})];
        if length(ind_NaNsamples{i_trial}{i_elec}) > 0
            temp_NaNpresent = temp_NaNpresent + 1;
        end
    end
    num_NaNelecs_pertrial(i_trial) = temp_NaNpresent;
    avgNaNsamplenum_pertrial(i_trial) = nanmean(temp_numNaNsamples);
    perc_NaNsample_pertrial(i_trial) = nanmean(temp_numNaNsamples) ...
        ./ length(input_data.trial{i_trial});
end

% %4A) Interpolate by connecting both edge values
% cfg = [];
% cfg.method = 'linear'; %'nearest','linear','spline','pchip','cubic','v5cubic'
% cfg.prewindow = 0.1; %in s
% cfg.postwindow = 0.1;
% datastruct_interpolated = ...
%     ft_interpolatenan(cfg, datastruct);

%4B) Interpolate by substituting the average of neighbouring noNaN samples
datastruct_interpolated = datastruct;
interpwindow_samples = 50;
for i_trial = 1:length(datastruct.trial)
    for i_elec = 1:size(datastruct.trial{i_trial},1)
        for i_NaNsample = find(isnan(datastruct.trial{i_trial}(i_elec,:)))  
            
            %prew window
            interpwindow_pre_values = [];
            iteration_counter = 0;
            counter = 1;
            while iteration_counter < 1
                if i_NaNsample - counter < 1
                    iteration_counter = 1;
                else
                    if ~isnan(datastruct.trial{i_trial}(i_elec,i_NaNsample - counter))
                        interpwindow_pre_values = ...
                            [interpwindow_pre_values ...
                            datastruct.trial{i_trial}(i_elec,i_NaNsample - counter)];
                    end
                    if length(interpwindow_pre_values) == interpwindow_samples
                        iteration_counter = 1;
                    end
                    counter = counter + 1;        
                end
            end
            
            %post window
            interpwindow_post_values = [];
            iteration_counter = 0;
            counter = 1;
            while iteration_counter < 1
                if i_NaNsample + counter > length(datastruct.trial{i_trial})
                    iteration_counter = 1;
                else
                    if ~isnan(datastruct.trial{i_trial}(i_elec,i_NaNsample + counter))
                        interpwindow_post_values = ...
                            [interpwindow_post_values ...
                            datastruct.trial{i_trial}(i_elec,i_NaNsample + counter)];
                    end
                    if length(interpwindow_post_values) == interpwindow_samples
                        iteration_counter = 1;
                    end
                    counter = counter + 1;        
                end
            end
            
            datastruct_interpolated.trial{i_trial}(i_elec,i_NaNsample) = ...
                mean([interpwindow_pre_values interpwindow_post_values]);
        end             
    end
end

%5) Visually check if NaNepoch placement and filtering are correct
% ft_databrowser([],datastruct_interpolated)
if plot_poststepFigs == 1
    for i_trial = 1:length(datastruct.trial)
        index_nanelects = [];
        for i_elec = 1:length(datastruct.label)
            if ~isempty(ind_NaNsamples{i_trial}{i_elec})
                index_nanelects = [index_nanelects i_elec];
            end
        end
        
        if ~isempty(index_nanelects)
            subplot_counter = 0;
            figure('visible','off');
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %fullscreen
            for i_elec = index_nanelects
                subplot_counter = subplot_counter + 1;
                subplot(round(sqrt(length(index_nanelects))),...
                    ceil(sqrt(length(index_nanelects))),subplot_counter);

                hold on;
                %Plot NaNepoch
                plot_nan = zeros(1,length(datastruct.trial{i_trial}));
                plot_nan(ind_NaNsamples{i_trial}{i_elec}) = ...
                    max(datastruct_interpolated.trial{i_trial}(i_elec,:));
                a1  = area(plot_nan);
                a1.FaceAlpha = 0.2; a1.FaceColor = [1 0 0];
                a1.EdgeAlpha = 0.2; a1.EdgeColor = [1 0 0];
                plot_nan = zeros(1,length(datastruct.trial{i_trial}));
                plot_nan(ind_NaNsamples{i_trial}{i_elec}) = ...
                    min(datastruct_interpolated.trial{i_trial}(i_elec,:));                
                a2  = area(plot_nan);
                a2.FaceAlpha = 0.2; a2.FaceColor = [1 0 0];
                a2.EdgeAlpha = 0.2; a2.EdgeColor = [1 0 0];
                %plot interploated data 
                plot(datastruct_interpolated.trial{i_trial}(i_elec,:),'r')
                %plot not interpolated data
                plot(datastruct.trial{i_trial}(i_elec,:),'k')
                
                %plot vertical line indicating trial beginning and end
                plot([datastruct_interpolated.fsample ...
                    datastruct_interpolated.fsample],...
                    [min(datastruct_interpolated.trial{i_trial}(i_elec,:)) ...
                    max(datastruct_interpolated.trial{i_trial}(i_elec,:))], ...
                    'k--')
                plot([length(datastruct_interpolated.trial{i_trial}) - datastruct_interpolated.fsample ...
                    length(datastruct_interpolated.trial{i_trial}) - datastruct_interpolated.fsample],...
                    [min(datastruct_interpolated.trial{i_trial}(i_elec,:)) ...
                    max(datastruct_interpolated.trial{i_trial}(i_elec,:))], ...
                    'k--')
                
                %plot title and suptitle
                title(['Elec ' datastruct.label(i_elec)])
            end
            sgtitle([sub ' - TD: ' TD_label ' - trial: ' num2str(i_trial)])

            if save_poststepFigs == 1
                path_fig = [paths_NASTD_ECoG.Preproc_ECoGdata sub ...
                    '/Figs/Preproc/FiltInterp/LP2to30Hz/' ];
                if (~exist(path_fig, 'dir')); mkdir(path_fig); end

                filename     = strcat([sub '_FiltInterp_LP2to30Hz_TD' TD_label ...
                    's_trial' num2str(i_trial) '.png']);            
                figfile      = [path_fig filename];
                saveas(gcf, [figfile], 'png'); %save png version  
                close
            end
        end
    end    
end

%6) Place filtered data and info in final output struct
trialdata_filtinterp                        = datastruct_interpolated.trial;
info_filtinterp.NaNind                      = ind_NaNsamples;
info_filtinterp.num_NaNelecs_pertrial       = num_NaNelecs_pertrial;
info_filtinterp.avgNaNsamplenum_pertrial    = avgNaNsamplenum_pertrial;
info_filtinterp.perc_NaNsample_pertrial     = perc_NaNsample_pertrial;

end