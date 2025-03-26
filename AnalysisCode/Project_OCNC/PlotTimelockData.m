%% 1) Specify SPM-input data in FT

%0) Set up paths and input vars%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0.1) Specify paths and setup fieldtrip
restoredefaultpath
addpath('D:\Programs\MATLAB\fieldtrip-20190314\') %path FT
ft_defaults

addpath(genpath('D:\Work\Kongresse&Workshops\2019\OCNC Okinawa Summer School\Project\Data\')) %path data
path_save = 'D:\Work\Kongresse&Workshops\2019\OCNC Okinawa Summer School\Project\Data\Analysis\Preproc_ECoGdata\NY688\';

%0.2) Specify vars
% toneIndex = 33; %focused tone
tonedur_text = '0.2';

%0.3) Specify subject (exemplary sub1)
sub = 'NY688';
subs_PreProcSettings = NASTD_SubjectPreProcSettings; %load in file with individual preproc infos

%0.4) Load data corresponding to SOM-input data
cd D:\Work\Kongresse&Workshops\2019\'OCNC Okinawa Summer School'\Project\Data\Analysis\Preproc_ECoGdata\NY688\
tic
disp('loading...')
    load('NY688_data_200ms_matlab.mat')                  %load preprossed ECoG data
    load('sfa_expt4 2018-10-19 13-43 data TJB.mat') %load behavioral data
disp(['done loading in ' num2str(toc) ' sec'])


%% Trial selection: Compute ERP and plot output (all and focussed on last 2 tones)
%Determine trial groups (i.e., large vs. small p*34-p34 distance)
for i_trial = 1:length(data_ECoG.trial)
    logdiff_predvspresent(i_trial) = data_ECoG.stim.logf_pred(i_trial) - data_ECoG.stim.logf_final(i_trial);
end
% figure; hold on
%     bar(1:length(unique(logdiff_predvspresent)), sort(abs(unique(logdiff_predvspresent))),'FaceColor',[0.7 0.7 0.7])
%     plot(1:length(unique(logdiff_predvspresent)),ones(1,9)*mean(abs(logdiff_predvspresent)),'r')
%     plot(1:length(unique(logdiff_predvspresent)),ones(1,9)*median(abs(logdiff_predvspresent)),'k')

f_likely = find(abs(logdiff_predvspresent) < median(abs(logdiff_predvspresent))); %lower than median p*34-p34 distance
f_unlikely = find(abs(logdiff_predvspresent) > median(abs(logdiff_predvspresent)));
f_unsure = find(abs(logdiff_predvspresent) == median(abs(logdiff_predvspresent))); %higher than median p*34-p34 distance

f_likely_resp = find(data_ECoG.behav.resp_prob > 3);
f_unsure_resp = find(data_ECoG.behav.resp_prob == 3);
f_unlikely_resp = find(data_ECoG.behav.resp_prob < 3);


%% Temporal parameters
nTrials  = length(data_ECoG.trial);
% nSensors = size( data_input.trial{1}, 1 ); %all channels
nSensors = length(data_ECoG.cfg.info_elec.selected.index4EDF); %selected channels

fsample         = data_ECoG.fsample;
toneDur_inSecs  = str2num(tonedur_text);
nSamplesPerTone = toneDur_inSecs * fsample;

nSamplesPerSeries = length(data_ECoG.trial{1}(1,:));
nSecsPerSeries = nSamplesPerSeries/fsample;
nSecsPerSeriesinTheory = toneDur_inSecs*34;


%% ERP/Timelock analysis
%Compute ERP
cfg = [];
cfg.latency = [data_ECoG.time{1}(end-floor(nSamplesPerTone)) data_ECoG.time{1}(end)]; %only final tone (p34)
cfg.keeptrials = 'yes';
cfg.removemean = 'no';

cfg.trials = f_likely;
ERP_likely = ft_timelockanalysis(cfg, data_ECoG);
cfg.trials = f_unlikely;
ERP_unlikely = ft_timelockanalysis(cfg, data_ECoG);
cfg.trials = f_likely_resp;
ERP_likely_resp = ft_timelockanalysis(cfg, data_ECoG);
cfg.trials = f_unlikely_resp;
ERP_unlikely_resp = ft_timelockanalysis(cfg, data_ECoG);

%determine 1) A1 channels, 2) STG channels, 3) IFG/Frontal channels
chan = struct;
chan.A1 = {'G63'};
%T1 coordinates:
% G39 -64.159805 -50.935322 19.088045 100.00% cSTG 
% G46 -62.884319 -39.909706 27.642714 100.00% cSTG 
% G47 -65.450333 -40.965710 19.031191 100.00% cSTG 
% G54 -64.122597 -30.033758 27.670658 100.00% cSTG 
% G55 -66.271805 -31.305193 18.915991 100.00% cSTG 
% G62 -63.625118 -20.027098 27.078493 100.00% mSTG 
% G63 -65.397896 -21.648180 18.049524 100.00% mSTG 
% T06 -64.231972 -7.734340 11.056163 100.00% mSTG 
% T07 -65.498558 -15.902032 15.346996 100.00% mSTG 
% T08 -65.058823 -24.614094 21.837561 100.00% cSTG 

% A1 for MNI coordinates: -52 -19 7
%MNI coordinates:
% G39 -69.500000 -46.500000 15.500000 G 
% G46 -70.000000 -34.000000 21.000000 G 
% G47 -72.000000 -37.333333 12.666667 G 
% G54 -72.000000 -26.000000 17.000000 G 
% G55 -73.200000 -29.200000 8.400000 G 
% G62 -72.666667 -17.000000 13.333333 G 
% G63 -74.000000 -19.333333 5.333333 G 
% T06 -73.142857 -7.142857 -5.142857 S 
% T07 -74.000000 -14.666667 0.666667 S 
% T08 -73.333333 -22.000000 8.666667 S 

chan.STG = {'G46'};
%T1 coordinates:
% G39 -64.159805 -50.935322 19.088045 100.00% cSTG 
% G46 -62.884319 -39.909706 27.642714 100.00% cSTG 
% G47 -65.450333 -40.965710 19.031191 100.00% cSTG 
% G54 -64.122597 -30.033758 27.670658 100.00% cSTG 
% G55 -66.271805 -31.305193 18.915991 100.00% cSTG 
% G62 -63.625118 -20.027098 27.078493 100.00% mSTG 
% G63 -65.397896 -21.648180 18.049524 100.00% mSTG 
% T06 -64.231972 -7.734340 11.056163 100.00% mSTG 
% T07 -65.498558 -15.902032 15.346996 100.00% mSTG 
% T08 -65.058823 -24.614094 21.837561 100.00% cSTG 

%MNI coordinates:
% G39 -69.500000 -46.500000 15.500000 G 
% G46 -70.000000 -34.000000 21.000000 G 
% G47 -72.000000 -37.333333 12.666667 G 
% G54 -72.000000 -26.000000 17.000000 G 
% G55 -73.200000 -29.200000 8.400000 G 
% G62 -72.666667 -17.000000 13.333333 G 
% G63 -74.000000 -19.333333 5.333333 G 
% T06 -73.142857 -7.142857 -5.142857 S 
% T07 -74.000000 -14.666667 0.666667 S 
% T08 -73.333333 -22.000000 8.666667 S 

chan.IFG = {'IF1'};
%T1
% IF01 -42.187981 34.914268 34.117516 100.00% parstriangularis 
% IF02 -46.726345 26.285583 36.303383 100.00% parstriangularis 
% IF03 -51.203339 16.831656 37.808311 100.00% parsopercularis 
% IF04 -52.773922 8.375415 41.319572 100.00% parsopercularis 
% IF05 -54.814228 -1.799271 42.940117 100.00% precentral 
% IF06 -56.709232 -11.409325 44.624222 100.00% precentral 

%MNI
% IF01 -57.500000 41.500000 11.000000 S 
% IF02 -60.666667 34.000000 15.333333 S 
% IF03 -64.666667 23.000000 18.666667 S 
% IF04 -64.500000 15.500000 26.500000 S 
% IF05 -66.000000 6.666667 31.333333 S 
% IF06 -66.000000 -2.666667 34.666667 S 

%Plot ERP
cfg = [];
cfg.ylim = [-25 25];

cfg.channel  = chan.A1;
figure;
    ft_singleplotER(cfg, ERP_likely, ERP_unlikely);
    title('ERP A1')
    legend('Likely','Unlikely')
cfg.channel  = chan.STG;
figure;
    ft_singleplotER(cfg, ERP_likely, ERP_unlikely);
    title('ERP STG')
    legend('Likely','Unlikely')
cfg.channel  = chan.IFG;
figure;
    ft_singleplotER(cfg, ERP_likely, ERP_unlikely);
    title('ERP IFG')
    legend('Likely','Unlikely')
