cd 'D:\Work\Kongresse&Workshops\2019\OCNC Okinawa Summer School\Project\Data\Analysis\Preproc_ECoGdata\NY688\'
load('NY688_reducedData_200ms.mat')
load('NY688_reducedData_200ms_matlab.mat')


% Apply trial labels to SPM data
for i_trial = 1:length(D.trials)
    
    if data_ECoG_reduced.trialinfo(i_trial,3) == 1
        D.trials(i_trial).label = 'likely';
    elseif data_ECoG_reduced.trialinfo(i_trial,3) == -1
        D.trials(i_trial).label = 'unlikely';
    end
end
path_save = 'D:\Work\Kongresse&Workshops\2019\OCNC Okinawa Summer School\Project\Data\Analysis\Preproc_ECoGdata\NY688\';
savefile_spm = [path_save 'NY688_reducedData_200ms.mat'];
save(savefile_spm, 'D','-v7.3');
