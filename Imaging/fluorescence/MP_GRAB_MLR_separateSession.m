function MP_GRAB_MLR_separateSession(dataIndex)

%% reference to Sul et al.2011
% running multilinear regression
% z-score = b0 + b1C(n) + b2C(n-1) + b3C(n-1) + b4R(n) + b5R(n-1) +
% b6R(n-2) + b7X(n) + b8X(n-1) + b9X(n-2);
% C: choice, R: reward, X:interaction


% separate a session into first/second half, run linear regression
% separately
nFiles = size(dataIndex,1);


%% load the behavior files first to get the session length
sesLength = [];

for ii = 1:nFiles
     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    % load fluo files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1
        
        if trialData.cutPoint == 0
            sesLength = [sesLength,length(trialData.cue)];
        else
            sesLength = [sesLength,trialData.cutPoint];
        end

       
    end
end


% check the 33, 66 percentile point of every session


%% go through every session, running MLR separately
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
        
        
        load(fullfile(fn_fluo.folder,fn_fluo.name));
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        
        
        %saveRegName = fullfile(savematpath,[fn_beh.name(1:end-7),'regCR.mat']);
        saveRegName = fullfile(savematpath,'regCR_sep.mat');  % regression for fluo change
        % saveRegName_ITI =  fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);
        
        %% linear regression with C(n+1), observable variables
        
        %  C(n-2) - C(n+1)
        %  R(n-2) - R(n+1)
        %   interaction term
        %if ~exist(saveRegName)
            
            % convert left/right choices into ipsi/contra choices
            % dummy code: ipsi: 0; contra: 1;
            
            choice = NaN(size(trials.left));
            if strcmp(dataIndex.RecordingSite{ii},'left')
                choice(trialData.response == 2) = 0;
                choice(trialData.response == 3) = 1;
            else
                choice(trialData.response == 2) = 1;
                choice(trialData.response == 3) = 0;
            end
            
            params_future.trigEvent1 = [choice(2:end);NaN];  % C(n+1)
            params_future.trigEvent2 = choice; %C(n)
            params_future.trigEvent3 = [NaN;choice(1:end-1)];   % C(n-1)
            params_future.trigEvent4 = [NaN;NaN;choice(1:end-2);];  % C(n-2)
            
            %second predictor is outcome; dummy-code: reward=1, error=0, miss=NaN
            %params_future.trigEvent4 = [NaN;NaN;params_future.trigEvent(1:end-2)]; % C(n-2)
            
            reward = NaN(size(trials.left));
            reward(trialData.response ~= 0 & trials.reward == 0) = 0;
            reward(trialData.response ~= 0 & trials.reward == 1) = 1;
            
            params_future.trigEvent5 = [reward(2:end);NaN];  % R(n+1)
            params_future.trigEvent6 = reward; % R(n)
            params_future.trigEvent7 = [NaN;reward(1:end-1)]; % R(n-1)
            params_future.trigEvent8 = [NaN;NaN;reward(1:end-2)]; % R(n-2)
            
            % interaction
            params_future.trigEvent9 = params_future.trigEvent1 .* params_future.trigEvent5;
            params_future.trigEvent10 = params_future.trigEvent2 .* params_future.trigEvent6;
            params_future.trigEvent11 = params_future.trigEvent3 .* params_future.trigEvent7;
            params_future.trigEvent12 = params_future.trigEvent4 .* params_future.trigEvent8;
            
            % average reward rate on 20 trials window
            params_future.trigEvent13 = NaN(size(trials.go));
            for kk = 1:length(trials.left)
                if kk <= 20
                    params_future.trigEvent13(kk) = sum(trials.reward(1:kk))/kk;
                else
                    params_future.trigEvent13(kk) = sum(trials.reward(kk-19:kk))/20;
                end
            end
            
            % cumulative reward
            params_future.trigEvent14=NaN(size(trials.left));
            for kk = 1:length(trials.left)
                params_future.trigEvent14(kk) = sum(trials.reward(1:kk));
            end
            params_future.trigEvent14 = (params_future.trigEvent14)/sum(trials.reward);
            
            
            
            
            future_event = concat_event(params_future);
            % params.trigEvent2(trials.doublereward) = 2;
            
            fieldname={'go'};
            trialMask = getMask(trials,fieldname);
            
            
            params.xtitle = 'Time from cue (s)';
            params.window = [-3:0.1:5];
            params.nback = 0;       %how many trials back to regress against
            params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
            params.interaction = false;
            params.ifplot = 0;
            params.trigTime = trialData.cueTimes;
            %only perform analysis on this subset of trials
            
            tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n+1)','R(n)', 'R(n-1)','R(n-2)',...
                'C(n+1)*R(n+1)','C(n)*R(n)','C(n-1)*R(n-1)','C(n-2)*R(n-2)','Reward Rate','Cumulative Reward'};
            
            trialMask1  = trialMask;
            trialMask1(round(length(trialMask)/2)+1:end) = 0;
            trialMask2 = trialMask;
            trialMask2(1:round(length(trialMask)/2)) = 0;
            
            % reg_cr_future=linear_regr( fluo.dia, fluo.t, future_event, params.trigTime, trialMask, params );
            reg_t = cells.t;
            if isfield(cells,'normdFF')
                reg_dFF = cells.normdFF;
            else
                reg_dFF = cells.dFF;
            end
            parfor j=1:numel(reg_dFF)
                if length(reg_t) > length(reg_dFF{1})
                    %reg_cr{j}=linear_regr( cells.dFF{j}, cells.t(1:length(cells.dFF{1})), future_event, params.trigTime, trialMask, params );
                    reg_cr1{j}=linear_regr( reg_dFF{j}, reg_t(1:length(cells.normdFF{1})), future_event, params.trigTime, trialMask1, params );
                else
                    %reg_cr{j}=linear_regr( cells.dFF{j}, cells.t, future_event, params.trigTime, trialMask, params );
                    reg_cr1{j}=linear_regr(reg_dFF{j}, reg_t, future_event, params.trigTime, trialMask1, params );
                end
            end
            parfor j=1:numel(reg_dFF)
                if length(reg_t) > length(reg_dFF{1})
                    %reg_cr{j}=linear_regr( cells.dFF{j}, cells.t(1:length(cells.dFF{1})), future_event, params.trigTime, trialMask, params );
                    reg_cr2{j}=linear_regr( reg_dFF{j}, reg_t(1:length(cells.normdFF{1})), future_event, params.trigTime, trialMask2, params );
                else
                    %reg_cr{j}=linear_regr( cells.dFF{j}, cells.t, future_event, params.trigTime, trialMask, params );
                    reg_cr2{j}=linear_regr(reg_dFF{j}, reg_t, future_event, params.trigTime, trialMask2, params );
                end
            end
            MP_plot_regr(reg_cr1,[],params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','MLR-norm-choiceoutcome_1half');    %png format
            saveas(gcf, 'MLR-norm-choiceoutcome_1half', 'fig');
            
            MP_plot_regr(reg_cr2,[],params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','MLR-norm-choiceoutcome_2half');    %png format
            saveas(gcf, 'MLR-norm-choiceoutcome_2half', 'fig');
            %%plot the field of view, with significance of choice/reward in every grid
            
            
            % plot_regr_fluo(reg_cr,params.pvalThresh,tlabel,params.xtitle);
            % print(gcf,'-dpng','MLR-choiceoutcome_sig_choice0_1');    %png format
            % saveas(gcf, 'MLR-choiceoutcome_sigchoice0_1', 'fig');
            % MP_plot_regrcoef_fluo(reg_cr_future_change,params.pvalThresh,tlabel,params.xtitle);
            
            %         print(gcf,'-dpng','MLR-choiceoutcome_cut_c(n+1)');    %png format
            %         saveas(gcf, 'MLR-choiceoutcome_cut_c(n+1)', 'fig');
            %
            %% separate ITI between C(n) and C(n-1)
            %
            %             iti_time_1 = NaN(1, length(trialData.cueTimes));
            %
            %             for tt=2:length(trialData.cueTimes)
            %                 iti_time_1(tt) = trialData.cueTimes(tt) - trialData.outcomeTimes(tt-1);
            %             end
            %             trialIndex_1 = 1:length(trialData.cueTimes);
            %
            %             prcIndex1_1 = trialIndex_1(iti_time_1'<prc1);
            %             prcIndex2_1 = trialIndex_1(iti_time_1'<prc2&iti_time_1'>prc1);
            %             prcIndex3_1 = trialIndex_1(iti_time_1'>prc2);
            %
            %             % the previous trials need to be manually set
            %             params1_1 = MP_get_previous(prcIndex1_1, trials, trialData);
            %
            %             params2_1 = MP_get_previous(prcIndex2_1, trials, trialData);
            %
            %             params3_1 = MP_get_previous(prcIndex3_1, trials, trialData);
            %                 % params.trigEvent2(trials.omitreward) = 0;
            %                 % params.trigEvent2(trials.doublereward) = 2;
            %             ITI_event_1 = concat_event(params1_1);
            %             ITI_event_2 = concat_event(params2_1);
            %             ITI_event_3 = concat_event(params3_1);
            %
            %             fieldname={'go'};
            %             trialMask = getMask(trials,fieldname);
            %
            %             % reg_cr1_change=linear_regr( fluo.resp, fluo.respT, [params1.trigEvent_1' params1.trigEvent_2' params1.trigEvent_3' params1.trigEvent2_1' params1.trigEvent2_2' params1.trigEvent2_3' params1.trigEvent3_1' params1.trigEvent3_2' params1.trigEvent3_3'], params1.trigTime1, trialMask, params1 );
            %             %reg_cr1_change_1=linear_regr( fluo.resp, fluo.respT, [params1_1.trigEvent' params1_1.trigEvent_1' params1_1.trigEvent_2' params1_1.trigEvent_3' params1_1.trigEvent2' params1_1.trigEvent2_1' params1_1.trigEvent2_2' params1_1.trigEvent2_3' params1_1.trigEvent3' params1_1.trigEvent3_1' params1_1.trigEvent3_2' params1_1.trigEvent3_3'], params1_1.trigTime1, trialMask, params1_1 );
            %             reg_cr1_change_1=linear_regr_PR( fluo.resp, fluo.respT,ITI_event_1, params1_1.trigTime, trialMask, params1_1, trialData.cueTimes );
            %             reg_cr2_change_1=linear_regr_PR( fluo.resp, fluo.respT, ITI_event_2, params2_1.trigTime, trialMask, params2_1, trialData.cueTimes );
            %             reg_cr3_change_1=linear_regr_PR( fluo.resp, fluo.respT, ITI_event_3,  params3_1.trigTime,trialMask, params3_1, trialData.cueTimes );
            %
            %
            %% running control multilinear regression
            % shuffle every factor one by one, keeping other factors intact
            
            % check
            
            % construct every regression factor
            %             params_ctrl.trigEvent1 = params.trigEvent;  % C(n)
            %             params_ctrl.trigEvent2 = [NaN;params_ctrl.trigEvent1(1:end-1)];  %C(n-1)
            %             params_ctrl.trigEvent3 = [NaN; NaN; params_ctrl.trigEvent1(1:end-2)];  %c(n-2);
            %             params_ctrl.trigEvent4 = params.trigEvent2;  % R(n);
            %             params_ctrl.trigEvent5 = [NaN; params_ctrl.trigEvent4(1:end-1)]; % R(n-1)
            %             params_ctrl.trigEvent6 = [NaN; NaN; params_ctrl.trigEvent4(1:end-2)];  %R(n-2)
            %             params_ctrl.trigEvent7 = params_ctrl.trigEvent1 .* params_ctrl.trigEvent4;  %X(n)
            %             params_ctrl.trigEvent8 = params_ctrl.trigEvent2 .* params_ctrl.trigEvent5;  %X(n-1);
            %             params_ctrl.trigEvent9 = params_ctrl.trigEvent3 .* params_ctrl.trigEvent6;  %X(n-1);
            %
            %             % concatenate it into a matrix
            %             params_ctrlMat = [];
            %             fields = fieldnames(params_ctrl);
            %             for jj = 1:length(fields)
            %                 params_ctrlMat = [params_ctrlMat params_ctrl.(fields{jj})];
            %             end
            %
            %             % shul
            %             fieldname={'go'};
            %             trialMask = getMask(trials,fieldname);
            %
            %
            %             params.xtitle = 'Time from cue (s)';
            %             params.window = [-3:0.1:5];
            %             params.nback = 0;       %how many trials back to regress against
            %             params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
            %             params.interaction = false;
            %             params.trigTime = trialData.cueTimes;
            %             params.ifplot = 0;
            %             %only perform analysis on this subset of trials
            %
            %             % iterate through all 9 factors, shuffle name one by one to get
            %             % the control regression for every factor
            %             tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};
            %             reg_cr_ctrl = linear_regr_ctrl(fluo.dia, fluo.t, params_ctrlMat, params.trigTime, trialMask, params, tlabel);
            %             reg_cr_change_ctrl = linear_regr_ctrl(fluo.resp, fluo.respT, params_ctrlMat, params.trigTime, trialMask, params, tlabel);
            %
            %             % MP_plot_regrcoef_fluo(reg_cr_shuffleR,params.pvalThresh,tlabel,params.xtitle);
            %
            %             % control for future choice
            %             params_ctrlMat_future = [];
            %             fields = fieldnames(params_future);
            %             for jj = 1:length(fields)
            %                 params_ctrlMat_future = [params_ctrlMat_future params_future.(fields{jj})];
            %             end
            %
            %             tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n+1)','R(n)', 'R(n-1)','R(n-2)',...
            %                     'C(n+1)*R(n+1)','C(n)*R(n)','C(n-1)*R(n-1)','C(n-2)*R(n-2)','Reward Rate','Cumulative Reward'};
            %
            %             reg_cr_future_ctrl = linear_regr_ctrl(fluo.dia, fluo.t, params_ctrlMat_future, params.trigTime, trialMask, params, tlabel);
            %             reg_cr_future_change_ctrl = linear_regr_ctrl(fluo.resp, fluo.respT, params_ctrlMat_future,  params.trigTime, trialMask, params, tlabel );
            %
            
            %% save all regression into one mat file
            % save all these in a structure
            %save(saveRegName, 'reg_cr', 'reg_cr1','reg_cr2','reg_cr3','reg_cr_future','reg_cr_future_ctrl','reg_cr_ctrl');
            save(saveRegName,'reg_cr1','reg_cr2');
            %save(saveRegName_ITI,'reg_cr1_change_1','reg_cr2_change_1','reg_cr3_change_1');
            close all;
%         else
%             display('Regression already done');
%         end
    end
end
end
