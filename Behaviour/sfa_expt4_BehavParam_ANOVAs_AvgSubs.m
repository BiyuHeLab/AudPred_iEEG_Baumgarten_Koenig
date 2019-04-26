function [stats_anova2] = sfa_expt4_BehavParam_ANOVAs_AvgSubs(subs, tonedur_text, location)

%% Group-level ANOVAs for FTPLestimation
%Script computes 1) across subject average ANOVAs for FTPLestimates and 2)
%uSeq identification against chance

% subs = {'TestThomas', 'TestLark','TestElla','TestJonathan'}; %healthy test subjects normal task2
% tonedur_text = {'0.2' '0.4' 'all'};
% location = 'server'; %gago/gogo


%% 0) Specify paths
if strcmp(location,'server'); %gago/gogo
    addpath('/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
elseif strcmp(location,'desktop'); %local on office desktop
    addpath('//gogo.sb.nyumc.org/data/gogodisk4/thomas/AuditoryPrediction/iEEG/')
end

paths_sfa_expt4 = sfa_expt4_paths(location);

     global ft_defaults;
     ft_default.trackusage = 'no';

 %% 1) load behavioral data, specify paths
for i_sub = 1:length(subs)
     
    sub = subs{i_sub};

    sfa_expt4_subjectinfo %load subject info file (var: si)

load([si.path_behavioraldata_sub]);%Load raw behavioral data 
% load([si.path_stimorder_sub]);%Load stimulus presentation order 

for i_tonedur = 1:length(tonedur_text)

tonedur = tonedur_text{i_tonedur};

%% 2) check counterbalancing (i.e., if all stim combinations present per block)

trials_by_block = reshape(1:120, [10 12]);
for i_block = 1:12
    
    trials = trials_by_block(:, i_block);
    
    count = zeros(2,3);
    
    for i_trial = trials
        count(data.stim.toneDurID(i_trial),... %tone duration conditions
            data.stim.predID(i_trial)+2) = ... %label for math. expected FTP (-1 = low, 0 = medium, 1 = high)  
            count(data.stim.toneDurID(i_trial),...
            data.stim.predID(i_trial)+2) + 1;
    end
    
    if all(count(:) == 1)
        block_check(i_block) = 1;
    else
        block_check(i_block) = 0;
        disp('bad counterbalancing')
    end
end
    

%% 3.1) define filters for trial selection

    % 3.1.1 general filter - ensure all behavioral responses were entered properly
    filter_resp = data.uSeq.respABC >= 0 & data.resp_prob >= 0 & data.trialNum > 0;
    %This takes care of the -1 responses in the final tone likelihood estimates
    filter_rt   = data.r1_RT > 0 & data.r2_RT > 0;

    % 3.1.2 Tone duration filter 
    switch tonedur
        case '0.2'
            filter_toneDur = data.stim.toneDur == 0.2;

        case '0.4'
            filter_toneDur = data.stim.toneDur == 0.4;
 
        otherwise
            filter_toneDur = ones(1,120);
    end

    % 3.1.3 manually defined filters
    filter_manual = ones(1,120);
    
    % 3.1.4 combined filter
    filter = filter_resp & filter_rt & filter_toneDur & filter_manual; %filter var for trial selection

    %% 3.2) apply filters
    % filtered data will be applied to newly defined structs "dataf" and
    % "stimf", since these are the only relevant subfields in behavioral data
   
    %define new struct 'dataf' containing all data subfields with behavioral responses but only for above-filtered trials
    %from data subfields
    dfn = fieldnames(data);
    for j = 1:length(dfn)
        if eval(['length(data.' dfn{j} ') == 120'])
           eval(['dataf.' dfn{j} ' = data.' dfn{j} '(filter);']);            
        end      
    end
    
    %from data.uSeq (i.e., task2) subfields
    dfn = fieldnames(data.uSeq);
    for j = 1:length(dfn)       
        if eval(['length(data.uSeq.' dfn{j} ') == 120'])
            if strcmp(dfn{j},'disp3') %special because 3 columns instead of one
            eval(['dataf.' dfn{j} ' = data.uSeq.' dfn{j} '(filter,:);']);            
            else   
            eval(['dataf.' dfn{j} ' = data.uSeq.' dfn{j} '(filter);']);
            end
        end  
    end
    dataf = rmfield(dataf,{'resp_beta','conf_beta', 'correct_beta', 'diff_beta', 'r3_RT'});

    %define new struct 'stimf' containing all data subfields with all stimulus-wise parameters but only for above-filtered trials
    sfn = fieldnames(data.stim);
    for j = 1:length(sfn)
        if eval(['length(data.stim.' sfn{j} ') == 120'])
            eval(['stimf.' sfn{j} ' = data.stim.' sfn{j} '(filter);']);
        end
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 4. Create data struct for ANOVA
    %% 4.1 data struct for 1way ANOVA (dep var: FTPLrating, factors/IV: presented FTP (p34))
    %empty placeholder vectors, later filled single subject data
    vecPredRating = []; %subject response/final tone likelihood rating (how well does the
    %final tone pitch fit in given the previously presented sequence 
    %(i.e., distance between presented and math. expected final tone;
    %1 - very unlikely, 5 - very likely))
    vecActual = []; %actually presented final tone frequency 
    %([-3:3],negative values = < 440 Hz, positive values = > 440 Hz)
    
    if strcmp(subs(i_sub),'NY688')
        fprintf('Stop - SubNY688 has different experimental parameters - not comparable with other subjects')
        pause
    else
        
        actual_FTP = {'-2' '-1' '1' '2'}; %actually presented final tone pitches
            if strcmp(tonedur,'all')
            %rowcount determined by [1 - (number of trials per condition))] 
            rowCount = -29; % all tone dur = 30 trials/condition (4 conditions)
            else 
            rowCount = -14; %for 0.2/0.4s tone dur = 15 trials/condition
            end
        
    end
               
            for L = 1:length(actual_FTP) %p34 loop
                actual = actual_FTP{L};   
        
                switch tonedur %filter tone duration
                    case '0.2'
                        filter_toneDur = data.stim.toneDur == 0.2;
                    case '0.4'
                        filter_toneDur = data.stim.toneDur == 0.4;
                    otherwise
                        filter_toneDur = ones(1,120);
                end               

                switch actual %filter actually presented final tone
                    %Note: Only NY688 has 6 options, all other have ony 4
                    case '-3'
                        filter_actual = data.stim.finalID == -3;  
                    case '-2'
                        filter_actual = data.stim.finalID == -2;
                    case '-1'
                        filter_actual = data.stim.finalID == -1;
                    case '1'
                        filter_actual = data.stim.finalID == 1;
                    case '2'
                        filter_actual = data.stim.finalID == 2;
                    case '3'
                        filter_actual = data.stim.finalID == 3;
                end 

            filter = filter_toneDur & filter_actual;

            %Place trial selection in new structure used for MATLAB ANOVA
            vecPredRating = [vecPredRating; data.resp_prob(filter ==1)']; 
            %vector reading out the subject's final tone likelihood rating for each trial in
            %order determined by above loops (i.e., hierarchical order)

            nTrialsPerCond = length(data.resp_prob(filter ==1));

            %Update row
            rowCount = rowCount + nTrialsPerCond;
            lastRow = rowCount + nTrialsPerCond - 1;

            %Copy order of loop vars into vectors (i.e., hierarchical sorting for
            %1) tone duration, 2) math. expected final tone, 3) actually presented
            %final tone
            [vecToneDur(rowCount: lastRow)] = deal(i_tonedur);
            [vecActual(rowCount: lastRow)] = deal(L);
            end       
        
            %Output: single-subject cell containing vectors with all trials,
            %hiearchically sorted by 1) tone dur (1:3), 2) Expected tone (1:3), 3) Actual tone (1:6)
            FTPLrating_anova1{i_sub,i_tonedur}.PredRating = vecPredRating;
            FTPLrating_anova1{i_sub,i_tonedur}.ToneDur = vecToneDur';
            FTPLrating_anova1{i_sub,i_tonedur}.Actual = vecActual';

            %Question: vecPredRating should only contain [1:5], but in some enties
            %there are -1. This is mentioned in th Group-level ANOVA script,
            %where -1 entries are then replaced with 1. Presumably, -1 entries
            %resemble missed responses. Therefore, I now replace them with NaNs.
            ii = find(FTPLrating_anova1{i_sub,i_tonedur}.PredRating == -1);
            FTPLrating_anova1{i_sub,i_tonedur}.PredRating(ii) = NaN;            

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 4.2 Data Struct for 2way ANOVA (dep var: FTPLrating, factors/IV: presented FTP (p34),
    %math. expected FTP (p*34))
    vecPredRating = []; %subject response/final tone likelihood rating (how well does the
        %final tone pitch fit in given the previously presented sequence 
        %(i.e., distance between presented and math. expected final tone;
        %1 - very unlikely, 5 - very likely))
    vecExpected = []; %math. expected final tone (-1 = low, 0 = medium, 1 = high)
    vecActual = []; %actually presented final tone frequency 
    %([-3:3],negative values = < 440 Hz, positive values = > 440 Hz)

    %empty placeholder vectors, later filled single subject data
    expected_FTP = {'low' 'med' 'high'}; %math. expected final tone pitches
    if strcmp(subs(i_sub),'NY688')
        fprintf('Stop - SubNY688 has different experimental parameters - not comparable with other subjects')
        pause
          
    else
        actual_FTP = {'-2' '-1' '1' '2'}; %actually presented final tone pitches
        
        if strcmp(tonedur,'all')
        rowCount = -9; %for all tone dur = 10 trials/condition
        else 
        rowCount = -4; %for 0.2/0.4s tone dur = 5 trials/condition
        end    
    
    end
      
        for k = 1:length(expected_FTP) %p*34 loop
            expected = expected_FTP{k};
            
            for L = 1:length(actual_FTP) %p34 loop
                actual = actual_FTP{L};
        
                switch tonedur %filter tone duration
                    case '0.2'
                        filter_toneDur = data.stim.toneDur == 0.2;
                    case '0.4'
                        filter_toneDur = data.stim.toneDur == 0.4;
                    otherwise
                        filter_toneDur = ones(1,120);
                end
                
                switch expected %filter math. expected final tone
                    case 'low'
                        filter_expected = data.stim.predID == -1;
                    case 'med'
                        filter_expected = data.stim.predID == 0;
                    case 'high'
                        filter_expected = data.stim.predID == 1;
                end

                switch actual %filter actually presented final tone
                    %Note: Only NY688 has 6 options, all other have ony 4
                    case '-3'
                        filter_actual = data.stim.finalID == -3;  
                    case '-2'
                        filter_actual = data.stim.finalID == -2;
                    case '-1'
                        filter_actual = data.stim.finalID == -1;
                    case '1'
                        filter_actual = data.stim.finalID == 1;
                    case '2'
                        filter_actual = data.stim.finalID == 2;
                    case '3'
                        filter_actual = data.stim.finalID == 3;
                end 

            filter = filter_toneDur & filter_expected & filter_actual;

            %Place trial selection in new structure used for MATLAB ANOVA
            vecPredRating = [vecPredRating; data.resp_prob(filter ==1)']; 
            %vector reading out the subject's final tone likelihood rating for each trial in
            %order determined by above loops (i.e., hierarchical order)

            nTrialsPerCond = length(data.resp_prob(filter ==1));

            %Update row
            rowCount = rowCount + nTrialsPerCond;
            lastRow = rowCount + nTrialsPerCond - 1;

            %Copy order of loop vars into vectors (i.e., hierarchical sorting for
            %1) tone duration, 2) math. expected final tone, 3) actually presented
            %final tone
            [vecToneDur(rowCount: lastRow)] = deal(i_tonedur);
            [vecExpected(rowCount: lastRow)] = deal(k);
            [vecActual(rowCount: lastRow)] = deal(L);
            end

        end
        
            %Output: single-subject cell containing vectors with all trials,
            %hiearchically sorted by 1) tone dur (1:3), 2) Expected tone (1:3), 3) Actual tone (1:6)
            FTPLrating_anova2{i_sub,i_tonedur}.PredRating = vecPredRating;
            FTPLrating_anova2{i_sub,i_tonedur}.ToneDur = vecToneDur';
            FTPLrating_anova2{i_sub,i_tonedur}.Expected = vecExpected';
            FTPLrating_anova2{i_sub,i_tonedur}.Actual = vecActual';

            %Question: vecPredRating should only contain [1:5], but in some enties
            %there are -1. This is mentioned in th Group-level ANOVA script,
            %where -1 entries are then replaced with 1. Presumably, -1 entries
            %resemble missed responses. Therefore, I now replace them with NaNs.
            ii = find(FTPLrating_anova2{i_sub,i_tonedur}.PredRating == -1);
            FTPLrating_anova2{i_sub,i_tonedur}.PredRating(ii) = NaN;            

            clear vec* subDurTbl rowCount2 lastRow2

end
end

%% 5 Concatenate trials across subjects per tone duration condition and compute respective ANOVA
addpath(paths_sfa_expt4.ScriptsDir) %path fo anova funcitons

% 5.1 1way ANOVA (dep var: FTPLrating, factors/IV: presented FTP (p34))
for i_tone = 1:length(tonedur_text) %Separate RM ANOVA for each tone duration
    
    %Create empty proxies
    tonePredRating{i_tone} = [];
    toneActual{i_tone} = [];
    subVec{i_tone} = [];
    
    for i_sub = 1:length(subs) %Concatenate all trials for all subs per tone condition and fill proxies        
        tonePredRating{i_tone}  = [tonePredRating{i_tone};  FTPLrating_anova1{i_sub,i_tone}.PredRating];
        toneActual{i_tone}      = [toneActual{i_tone};      FTPLrating_anova1{i_sub,i_tone}.Actual];
        ntrials                 = length(FTPLrating_anova2{i_sub,i_tone}.PredRating);
        subVec{i_tone}          = [subVec{i_tone};          repelem(i_sub,ntrials)'];
    end
    %Repeated Measures Two-way Analysis of Variance Test
    %dependent var: final tone likelihood rating;
    %factors: math. expected FTP, presented FTP
     RMAOV1([tonePredRating{i_tone}, toneActual{i_tone}, subVec{i_tone}]);
    %IV1: actual FTP  
  
end

% 5.2 2way ANOVA (dep var: FTPLrating, factors/IV: expected FTP (p*34), presented FTP (p34),
for i_tone = 1:length(tonedur_text) %Separate RM ANOVA for each tone duration
    
    %Create empty proxies
    tonePredRating{i_tone} = [];
    toneExpected{i_tone} = [];
    toneActual{i_tone} = [];
    subVec{i_tone} = [];
    
    for i_sub = 1:length(subs) %Concatenate all trials for all subs per tone condition and fill proxies        
        tonePredRating{i_tone}  = [tonePredRating{i_tone};  FTPLrating_anova2{i_sub,i_tone}.PredRating];
        toneExpected{i_tone} 	= [toneExpected{i_tone};    FTPLrating_anova2{i_sub,i_tone}.Expected];
        toneActual{i_tone}      = [toneActual{i_tone};      FTPLrating_anova2{i_sub,i_tone}.Actual];
        ntrials                 = length(FTPLrating_anova2{i_sub,i_tone}.PredRating);
        subVec{i_tone}          = [subVec{i_tone};          repelem(i_sub,ntrials)'];
    end
    %Repeated Measures Two-way Analysis of Variance Test
    %dependent var: final tone likelihood rating;
    %factors: math. expected FTP, presented FTP
    RMAOV2([tonePredRating{i_tone}, toneExpected{i_tone}, toneActual{i_tone}, subVec{i_tone}]);
    stats_anova2{i_tone} = rm_anova2(tonePredRating{i_tone}, subVec{i_tone}, toneExpected{i_tone}, toneActual{i_tone}, {'math.expected tone','presented tone'});    
    %IV1: math. expected FTP
    %IV2: actual FTP      
   
end

end
