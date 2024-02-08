function MP_GRABRL_sumQ_MLR(dataIndex)

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
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
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
        saveRegName = fullfile(savematpath,'regRL_sumQ_norm.mat');  % regression for fluo changehange
        % saveMLRmatpath_outcome = fullfile(dataIndex.BehPath{ii},[fn_beh.name(1:end-7),'regRL_lag0_outcome_cut_fitall.mat']);
        
        %if ~exist(saveRegName)
            params=[];
            
             choice = NaN(size(trials.left));
            if strcmp(dataIndex.RecordingSite{ii},'left')
                choice(trialData.response == 2) = 0;
                choice(trialData.response == 3) = 1;
            else
                choice(trialData.response == 2) = 1;
                choice(trialData.response == 3) = 0;
            end
            % dummycode ipsi: 0, contra 1
            %C(n)
            params.trigEvent = choice;
            %R(n)
            reward = NaN(size(trials.go));
            reward(trials.reward) = 1;
            reward(trials.noreward) = 0;
            params.trigEvent2 = reward;
            %C(n)*R(n)
            params.trigEvent3 = params.trigEvent .* params.trigEvent2; % interaction term
            % C(n-1)
            params.trigEvent4= [NaN; choice(1:end-1,1)];
            % R(n-1)
            params.trigEvent5 = [NaN;reward(1:end-1,1)];
            % C(n-1)*R(n-1)
            params.trigEvent6 = params.trigEvent4 .* params.trigEvent5;
            
            % delta Q
           params.trigEvent7 = NaN(size(trials.left));
            if strcmp(dataIndex.RecordingSite{ii},'left')
                params.trigEvent7 =stats_new.qr-stats_new.ql;
            else
                params.trigEvent7 =stats_new.ql-stats_new.qr;
            end
          
          
            % sum Q
            params.trigEvent8 = stats_new.ql+stats_new.qr;
%             % delta choice kernel
%             params.trigEvent9 = stats_new.ckl - stats_new.ckr;
%             
%             % chosen choice kernel
%             params.trigEvent10 = NaN(length(stats_new.ql),1);
%             params.trigEvent10(stats.c(:,1)==-1) = stats_new.ckl(stats.c(:,1)==-1);
%             params.trigEvent10(stats.c(:,1) == 1) = stats_new.ckr(stats.c(:,1) == 1);
%             
            % average reward rate on 20 trials window
            params.trigEvent9 = NaN(size(trials.go));
            for kk = 1:length(trials.left)
                if kk <= 20
                    params.trigEvent9(kk) = sum(trials.reward(1:kk))/kk;
                else
                    params.trigEvent9(kk) = sum(trials.reward(kk-19:kk))/20;
                end
            end
            
            % cumulative reward
            params.trigEvent10=NaN(size(trials.left));
            for kk = 1:length(trials.left)
                params.trigEvent10(kk) = sum(trials.reward(1:kk));
            end
            params.trigEvent10 = params.trigEvent10/sum(trials.reward);
            
            
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
            tlabel={'c(n)','r(n)','c(n)xr(n)','c(n-1)','r(n-1)', 'c(n-1)xr(n-1)','dQ', 'sumQ','Reward Rate', 'Cumulavtive reward'};
            
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
            parfor j=1:numel(LR_dFF )
                if length(LR_t) > length(LR_dFF {1})
                    reg_cr{j}=linear_regr( LR_dFF {j}, LR_t(1:length(LR_dFF {1})), RL_event, params.trigTime, trialMask, params );
                else
                    reg_cr{j}=linear_regr( LR_dFF {j}, LR_t, RL_event, params.trigTime, trialMask, params );
                end
            end
            MP_plot_regr(reg_cr,[],params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','MLR-sumQ-norm');    %png format
            saveas(gcf, 'MLR-sumQ-norm', 'fig');
            
            %% running control multilinear regression
            % shuffle every factor one by one, keeping other factors intact
            
            % check
            
            %             % construct every regression factor
            %             params_ctrl = params;
            %
            %             % concatenate it into a matrix
            %             params_ctrlMat = [];
            %             fields = fieldnames(params_ctrl);
            %             for jj = 1:length(fields)
            %                 if contains(fields{jj},'trigEvent')
            %                     params_ctrlMat = [params_ctrlMat params_ctrl.(fields{jj})];
            %                 end
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
            %             %params.trigTime = trialData.cueTimes;
            %             params.ifplot = 0;
            %             %only perform analysis on this subset of trials
            %
            %             % iterate through all 9 factors, shuffle name one by one to get
            %             % the control regression for every factor
            %             tlabel={'C(n)','R(n)','C(n)xR(n)','C(n-1)','R(n-1)', 'C(n-1)xR(n-1)','QL-QR', 'ChosenQ', 'Reward Rate', 'Cumulavtive reward'};
            %
            %
            %             reg_cr_ctrl = linear_regr_ctrl(pupil.dia, pupil.t, params_ctrlMat, params.trigTime, trialMask, params, tlabel);
            %             reg_cr_change_ctrl = linear_regr_ctrl(pupil.resp, pupil.respT, params_ctrlMat,  params.trigTime, trialMask, params, tlabel );
            %
            %             save(saveMLRmatpath, 'reg_cr','reg_cr_ctrl');
            save(saveRegName, 'reg_cr');
            close all;
        %else
        %    display('Regression already done');
        %end
        
    end
end
end