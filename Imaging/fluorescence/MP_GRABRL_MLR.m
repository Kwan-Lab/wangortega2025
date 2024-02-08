function MP_GRABRL_MLR(dataIndex)

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
        saveRegName = fullfile(savematpath,'regRL_norm.mat');  % regression for fluo changehange
        % saveMLRmatpath_outcome = fullfile(dataIndex.BehPath{ii},[fn_beh.name(1:end-7),'regRL_lag0_outcome_cut_fitall.mat']);

       % if ~exist(saveRegName)
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

            % Q and K need to be ipsi/contra
            % delta Q, Q(contra) - Q(ipsi)

            params.trigEvent7 = NaN(size(trials.left));
            if strcmp(dataIndex.RecordingSite{ii},'left')
                params.trigEvent7 =stats_new.qr-stats_new.ql;
            else
                params.trigEvent7 =stats_new.ql-stats_new.qr;
            end


            % chosen Q
            params.trigEvent8 = NaN(length(stats_new.ql),1);
            params.trigEvent8(stats.c(:,1)==-1) = stats_new.ql(stats.c(:,1)==-1);
            params.trigEvent8(stats.c(:,1) == 1) = stats_new.qr(stats.c(:,1) == 1);

            % delta choice kernel
            params.trigEvent9 = NaN(size(trials.left));
            if strcmp(dataIndex.RecordingSite{ii},'left')
                params.trigEvent9 =stats_new.ckr-stats_new.ckl;
            else
                params.trigEvent9 =stats_new.ckl-stats_new.ckr;
            end

            % chosen choice kernel
            params.trigEvent10 = NaN(length(stats_new.ql),1);
            params.trigEvent10(stats.c(:,1)==-1) = stats_new.ckl(stats.c(:,1)==-1);
            params.trigEvent10(stats.c(:,1) == 1) = stats_new.ckr(stats.c(:,1) == 1);

            % average reward rate on 20 trials window
            params.trigEvent11 = NaN(size(trials.go));
            for kk = 1:length(trials.left)
                if kk <= 20
                    params.trigEvent11(kk) = sum(trials.reward(1:kk))/kk;
                else
                    params.trigEvent11(kk) = sum(trials.reward(kk-19:kk))/20;
                end
            end

            % cumulative reward
            params.trigEvent12=NaN(size(trials.left));
            for kk = 1:length(trials.left)
                params.trigEvent12(kk) = sum(trials.reward(1:kk));
            end
            params.trigEvent12 = params.trigEvent12/sum(trials.reward);


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
            tlabel={'c(n)','r(n)','c(n)xr(n)','c(n-1)','r(n-1)', 'c(n-1)xr(n-1)','dQ', 'ChosenQ', 'dK', 'ChosenK','Reward Rate', 'Cumulavtive reward'};

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
                    reg_cr{j}=linear_regr( LR_dFF {j}(1:length(LR_t)), LR_t, RL_event, params.trigTime, trialMask, params );
                end
            end

            % get bootstrap
            all_coeff = [];
            for rr = 1:length(reg_cr)
                all_coeff = cat(3,all_coeff, reg_cr{rr}.coeff);
            end
            reg_cr_all.coeff= all_coeff;

            % use bootstrp to get coefficient
            reg_cr_all = getBootstrp(reg_cr_all, 0, 0.05);

            reg_cr_all.regr_time = reg_cr{1}.regr_time;
            reg_cr_all.numPredictor = reg_cr{1}.numPredictor;
            reg_cr_all.nback = reg_cr{1}.nback;
            reg_cr_all.interaction = reg_cr{1}.interaction;
            reg_cr_all.pvalThresh= 0.01;

            MP_plot_regrcoef_fluo(reg_cr_all,params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','MLR-norm_choiceselection');    %png format
            saveas(gcf, 'MLR-norm_choiceselection', 'fig');
            saveas(gcf, 'MLR-norm_choiceselection','svg');

            MP_plot_regr(reg_cr,[],params.pvalThresh,tlabel,params.xtitle);
            print(gcf,'-dpng','MLR-choiceselection-norm');    %png format
            saveas(gcf, 'MLR-choiceselection-norm', 'fig');

            save(saveRegName, 'reg_cr');
            close all;
%         else
%             display('Regression already done');
%             load(saveRegName);
%             params.xtitle = {'Time from cue (s)'};
%             tlabel={'c(n)','r(n)','c(n)xr(n)','c(n-1)','r(n-1)', 'c(n-1)xr(n-1)','dQ', 'ChosenQ', 'dK', 'ChosenK','Reward Rate', 'Cumulavtive reward'};
%             params.pvalThresh = 0.01;
% 
%             % get bootstrap
%             all_coeff = [];
%             for rr = 1:length(reg_cr)
%                 all_coeff = cat(3,all_coeff, reg_cr{rr}.coeff);
%             end
%             reg_cr_all.coeff= all_coeff;
% 
%             % use bootstrp to get coefficient
%             reg_cr_all = getBootstrp(reg_cr_all, 0, 0.05);
% 
%             reg_cr_all.regr_time = reg_cr{1}.regr_time;
%             reg_cr_all.numPredictor = reg_cr{1}.numPredictor;
%             reg_cr_all.nback = reg_cr{1}.nback;
%             reg_cr_all.interaction = reg_cr{1}.interaction;
%             reg_cr_all.pvalThresh= 0.01;
% 
%             MP_plot_regrcoef_fluo(reg_cr_all,params.pvalThresh,tlabel,params.xtitle);
%             print(gcf,'-dpng','MLR-norm_choiceselection');    %png format
%             saveas(gcf, 'MLR-norm_choiceselection', 'fig');
%             saveas(gcf, 'MLR-norm_choiceselection','svg');
% 
%             MP_plot_regr(reg_cr,[],params.pvalThresh,tlabel,params.xtitle);
%             print(gcf,'-dpng','MLR-choiceselection-norm');    %png format
%             saveas(gcf, 'MLR-choiceselection-norm', 'fig');
% 
%             close all;
%         end

    end
end
end