function MP_GRABRL_RPE_MLR(dataIndex)

% add the cutpoint later

nFiles = size(dataIndex,1);
%% reference: Sul et al.
% s(t) = a0 + a1C(t) + a2R(t) + a3X(t) + a4deltaQ(t) + a5Qc(t) + A(t)
% no A(t) first
%% fit the model per session
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    fn_latent = fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_FQRPECKlatentV.mat']);
    load(fn_latent);
    % load pupil files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    
    if length(fn_fluo) == 1
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo/LinearRegr');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
        
        
        load(fullfile(fn_fluo.folder,fn_fluo.name));
        
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        %saveMLRmatpath = fullfile(savematpath,[fn_beh.name(1:end-7),'regRL_CK.mat']);
        saveRegName = fullfile(savematpath,'regRL_RPE_norm.mat');  % regression for fluo changehange
        % saveMLRmatpath_outcome = fullfile(dataIndex.BehPath{ii},[fn_beh.name(1:end-7),'regRL_lag0_outcome_cut_fitall.mat']);
        
        if exist(saveRegName)
            load(saveRegName)
        end
        
        tic
        posIndex = stats.r > 0;
        negIndex = stats.r == 0;
        
        params=[];
        
        choice = NaN(size(trials.left));
        if strcmp(dataIndex.RecordingSite{ii},'left')
            choice(trialData.response == 2) = 0;
            choice(trialData.response == 3) = 1;
        else
            choice(trialData.response == 2) = 1;
            choice(trialData.response == 3) = 0;
        end
        % dummycode left: 0, right 1
        %C(n)
        params.trigEvent = choice;
        
        
        params.trigEvent2= [NaN; params.trigEvent(1:end-1,1)]; % c(n-1)
        
        %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];
        reward = NaN(size(trials.go));
        reward(trials.reward) = 1;
        reward(trials.noreward) = 0;
        params.trigEvent3= [NaN;reward(1:end-1)];   % r(n-1)
        
        % dQ
        if strcmp(dataIndex.RecordingSite{ii},'left')
            params.trigEvent4 =stats_new.qr-stats_new.ql;
        else
            params.trigEvent4 =stats_new.ql-stats_new.qr;
        end
        
        %RPE
        params.trigEvent5= stats_new.rpe;
        % dK
        if strcmp(dataIndex.RecordingSite{ii},'left')
            params.trigEvent6 =stats_new.ckr-stats_new.ckl;
        else
            params.trigEvent6 =stats_new.ckl-stats_new.ckr;
        end
        % CKE
        % 1-choice choice
        params.trigEvent7 = NaN(size(trials.go));
        params.trigEvent7(stats_new.c==-1) = 1-stats_new.ckl(stats_new.c==-1); % choosing left
        params.trigEvent7(stats_new.c==1) = 1-stats_new.ckr(stats_new.c==1);
        
        params.trigEvent8 = NaN(size(trials.go));
        for kk = 1:length(trials.left)
            if kk <= 20
                params.trigEvent8(kk) = sum(trials.reward(1:kk))/kk;
            else
                params.trigEvent8(kk) = sum(trials.reward(kk-19:kk))/20;
            end
        end
        
        params.trigEvent9=NaN(size(trials.left));
        for kk = 1:length(trials.left)
            params.trigEvent9(kk) = sum(trials.reward(1:kk));
        end
        params.trigEvent9 = params.trigEvent9/sum(trials.reward);
        
        
        % make matrix
        RL_event = concat_event(params);
        
        % when align pupil signal to cue
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.nback = 0;       %how many trials back to regress against
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        params.interaction = false;
        params.ifplot = 0;
        %only perform analysis on this subset of trials
        fieldname={'go'};
        trialMask = getMask(trials,fieldname);
        tlabel={'C(n)','C(n-1)','R(n-1)','deltaQ','RPE', 'deltaK','CKE','Reward Rate', 'Cumulavtive reward'};
        
        %reg_cr=linear_regr( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params );
        %         reg_cr_ctrl = linear_regr_ctrl( pupil.dia, pupil.t, RL_event, params.trigTime, trialMask, params, tlabel);
        %         save(saveMLRmatpath, 'reg_cr','reg_cr_ctrl');
        
        params.xtitle = {'Time from cue (s)'};
        LR_t = cells.t;
        if isfield(cells,'normdFF')
            LR_dFF = cells.normdFF;
        else
            LR_dFF = cells.dFF;
        end
        %if ~exist('reg_cr','var')
            parfor j=1:numel(LR_dFF )
                if length(LR_t) > length(LR_dFF {1})
                    reg_cr{j}=linear_regr( LR_dFF {j}, LR_t(1:length(LR_dFF {1})), RL_event, params.trigTime, trialMask, params );
                else
                    reg_cr{j}=linear_regr( LR_dFF {j}(1:length(LR_t)), LR_t, RL_event, params.trigTime, trialMask, params );
                end
            end
           %end 
            MP_plot_regr(reg_cr,[],params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','MLR-valueupdating_norm');    %png format
            saveas(gcf, 'MLR-valueupdating_norm', 'fig');
        
        %% positive RPE trials
       % if ~exist('reg_cr_pos','var')
            paramsPos.trigEvent = stats.c(posIndex,1);
            % dummycode left: 0, right 1
            paramsPos.trigEvent(paramsPos.trigEvent == -1) = 0;
            
            paramsPos.trigEvent2 = stats_new.ql(posIndex)-stats_new.qr(posIndex); % delta q
            
            preChoice = [NaN; params.trigEvent(1:end-1,1)];
            paramsPos.trigEvent3= preChoice(posIndex); % c(n-1)
            
            %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];
            preReward = [NaN;stats.r(1:end-1)];
            paramsPos.trigEvent4= preReward(posIndex);   % r(n-1)
            %second predictor is action value right
            %params.trigEvent5=[stats_sim.qr(2:end) + stats_sim.ql(2:end);NaN];
            paramsPos.trigEvent5= stats_new.rpe(posIndex);
            
            % make matrix
            RL_event_pos = concat_event(paramsPos);
            
            params.trigTime = trialData.cueTimes(posIndex);
            % when align pupil signal to cue
            
            
            %only perform analysis on this subset of trials
            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
            
            tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'posRPE'};
            params.xtitle = {'Time from cue (s)'};
            LR_t = cells.t;
            if isfield(cells,'normdFF')
                LR_dFF = cells.normdFF;
            else
                LR_dFF = cells.dFF;
            end
            
            parfor j=1:numel(LR_dFF )
                if length(LR_t) > length(LR_dFF {1})
                    reg_cr_pos{j}=linear_regr( LR_dFF {j}, LR_t(1:length(LR_dFF {1})), RL_event_pos, params.trigTime, trialMask, params );
                else
                    reg_cr_pos{j}=linear_regr( LR_dFF {j}(1:length(LR_t)), LR_t, RL_event_pos, params.trigTime, trialMask, params );
                end
            end
            %end
            params.xtitle = {'Time from cue (s)'};
            MP_plot_regr(reg_cr_pos,[],params.pvalThresh,tlabel,params.xtitle);
            
            
            print(gcf,'-dpng','MLR-RL_posRPE_fitall');    %png format
            saveas(gcf, 'MLR-RL_posRPE_fitall', 'fig');
        
        %% negative RPE trials
       % if ~exist('reg_cr_neg','var')
            paramsNeg.trigEvent = stats.c(negIndex,1);
            % dummycode left: 0, right 1
            paramsNeg.trigEvent(paramsNeg.trigEvent == -1) = 0;
            
            paramsNeg.trigEvent2 = stats_new.ql(negIndex)-stats_new.qr(negIndex); % delta q
            
            
            paramsNeg.trigEvent3= preChoice(negIndex); % c(n-1)
            
            %params.trigEvent4=[stats_sim.ql(2:end)-stats_sim.qr(2:end);NaN];
            
            paramsNeg.trigEvent4= preReward(negIndex);   % r(n-1)
            %second predictor is action value right
            %params.trigEvent5=[stats_sim.qr(2:end) + stats_sim.ql(2:end);NaN];
            paramsNeg.trigEvent5= stats_new.rpe(negIndex);
            
            % make matrix
            RL_event_neg = concat_event(paramsNeg);
            params.trigTime = trialData.cueTimes(negIndex);
            
            % when align pupil signal to cue
            
            %only perform analysis on this subset of trials
            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
            
            tlabel={'C(n)','QL-QR','C(n-1)','R(n-1)', 'negRPE'};
            params.xtitle = {'Time from cue (s)'};
            LR_t = cells.t;
            if isfield(cells,'normdFF')
                LR_dFF = cells.normdFF;
            else
                LR_dFF = cells.dFF;
            end
            
            parfor j=1:numel(LR_dFF )
                if length(LR_t) > length(LR_dFF {1})
                    reg_cr_neg{j}=linear_regr( LR_dFF {j}, LR_t(1:length(LR_dFF {1})), RL_event_neg, params.trigTime, trialMask, params );
                else
                    reg_cr_neg{j}=linear_regr( LR_dFF {j}(1:length(LR_t)), LR_t, RL_event_neg, params.trigTime, trialMask, params );
                end
            end
         %end   
            params.xtitle = {'Time from cue (s)'};
            MP_plot_regr(reg_cr_neg,[],params.pvalThresh,tlabel,params.xtitle);
            
            print(gcf,'-dpng','MLR-RL_negRPE_fitall');    %png format
            saveas(gcf, 'MLR-RL_RPE_negfitall', 'fig');
        
       
        save(saveRegName, 'reg_cr','reg_cr_pos','reg_cr_neg');
        close all;
        clear reg_cr reg_cr_pos reg_cr_neg
        toc
    end
end
end