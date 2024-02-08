function MP_GRAB_RF_RL(dataIndex)

%% build a random forest regressor to estimate the variable importance

% 1. OOB-predictor importance -- how to compare them. mean/std. make sense
% to compare them
% 2. add/delete predictors: based on the pred importance.
% 3. multiple random forest to see how the predictor importance varied: done.
% The estimation is stable
% 4. add more variable (c(n+1), r(n+1),cumulative reward)
nFiles = size(dataIndex,1);




%% go through every session, running RF separately
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    fn_latent = fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_FQRPECKlatentV.mat']);
    load(fn_latent);
    
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
        saveRFName = fullfile(savematpath,'RF_AS_100.mat');  % random forest, action-selection
        % saveRegName_ITI =  fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);
        %
        
        if ~exist(saveRFName)
            fn_ROI = dir(fullfile(dataIndex.LogFilePath{ii},'cell*.mat'));
            meanProjPath = dir(fullfile(dataIndex.BehPath{ii},'*Mean*.tif'));
            meanProj = loadtiffseq([],fullfile(meanProjPath.folder,meanProjPath.name));
            
            %get the mean intensity in ROI
            meanIntensity = zeros(1, size(fn_ROI,1));
            varIntensity = zeros(1,size(fn_ROI,1));
            for uu = 1:size(fn_ROI,1)
                roi = load(fullfile(fn_ROI(uu).folder,fn_ROI(uu).name));
                meanIntensity(uu) = mean(meanProj(roi.bw==1));
                varIntensity(uu) = nanstd(cells.cellf{uu});
            end
            
            %% random forest
            % use c(n), c(n-1),r(n),r(n-1),x(n),x(n-1), average r for now
            
            %  C(n-2) - C(n+1)
            %  R(n-2) - R(n+1)
            %   interaction term
            %         if ~exist(saveRFName)
            tic
            nPred = 10; nTrials = length(trials.left);
            event = zeros(nTrials,nPred);
            choice = NaN(size(trials.left));
            choice(trialData.response == 2) = 0;
            choice(trialData.response == 3) = 1;
            
            %params_future.trigEvent1 = [choice(2:end);NaN];  % C(n+1)
            %params_future.trigEvent2 = choice; %C(n)
            %params_future.trigEvent3 = [NaN;choice(1:end-1)];   % C(n-1)
            %params_future.trigEvent4 = [NaN;NaN;choice(1:end-2);];  % C(n-2)
            event(:,1) = categorical(choice); %C(n)
            event(:,2) = categorical([NaN;choice(1:end-1)]);   % C(n-1)
            %second predictor is outcome; dummy-code: reward=1, error=0, miss=NaN
            %params_future.trigEvent4 = [NaN;NaN;params_future.trigEvent(1:end-2)]; % C(n-2)
            
            reward = NaN(size(trials.left));
            reward(trialData.response ~= 0 & trials.reward == 0) = 0;
            reward(trialData.response ~= 0 & trials.reward == 1) = 1;
            
            event(:,3) =  categorical(reward); % R(n)
            event(:,4) =  categorical([NaN;reward(1:end-1)]); % R(n-1)
            %             params_future.trigEvent5 = [reward(2:end);NaN];  % R(n+1)
            %             params_future.trigEvent6 = reward; % R(n)
            %             params_future.trigEvent7 = [NaN;reward(1:end-1)]; % R(n-1)
            %             params_future.trigEvent8 = [NaN;NaN;reward(1:end-2)]; % R(n-2)
            
            % interaction
            %event(:,5) = event(:,1) .* event(:,3);
            %event(:,6) = event(:,2) .* event(:,4);
            %             params_future.trigEvent9 = params_future.trigEvent1 .* params_future.trigEvent5;
            %             params_future.trigEvent10 = params_future.trigEvent2 .* params_future.trigEvent6;
            %             params_future.trigEvent11 = params_future.trigEvent3 .* params_future.trigEvent7;
            %             params_future.trigEvent12 = params_future.trigEvent4 .* params_future.trigEvent8;
            
           
            
            % c(n+1)
             % delta Q
            
            event(:,5) = stats_new.ql-stats_new.qr;  % dQ
            % chosenQ
             event(:,6) = NaN(length(stats_new.ql),1);
            event(stats.c(:,1)==-1,6) = stats_new.ql(stats.c(:,1)==-1);
            event(stats.c(:,1) == 1,6) = stats_new.qr(stats.c(:,1) == 1);
              % dK
              event(:,7) = stats_new.ckl - stats_new.ckr;
              % chosenK
              event(:,8) = NaN(length(stats_new.ql),1);
            event(stats.c(:,1)==-1,8) = stats_new.ckl(stats.c(:,1)==-1);
            event(stats.c(:,1) == 1,8) = stats_new.ckr(stats.c(:,1) == 1);
             % average reward rate on 20 trials window, results shows not
            % important
           
             event(:,9) = NaN(size(trials.go));
            for kk = 1:length(trials.left)
                if kk <= 20
                    event(kk,9) = sum(trials.reward(1:kk))/kk;
                else
                    event(kk,9) = sum(trials.reward(kk-19:kk))/20;
                end
            end
            
            %cumulative reward
            event(:,10)=NaN(size(trials.left));
            for kk = 1:length(trials.left)
                event(kk,10) = sum(trials.reward(1:kk));
            end
            event(:,10) = event(:,10)/sum(trials.reward);
            %
            X = array2table(event);
            % first 4 predictors are categorical predictors
            for uu = [1:4]
                X.(uu) = categorical(X.(uu));
            end
            
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
            params.bin = [];
            %only perform analysis on this subset of trials
            %
            %             tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n+1)','R(n)', 'R(n-1)','R(n-2)',...
            %                 'C(n+1)*R(n+1)','C(n)*R(n)','C(n-1)*R(n-1)','C(n-2)*R(n-2)','Reward Rate','Cumulative Reward'};
            %
            % reg_cr_future=linear_regr( fluo.dia, fluo.t, future_event, params.trigTime, trialMask, params );
            tic
            RF_t = cells.t;
            RF_dFF = cells.dFF;

            warning('off','all')

            parfor j=1:numel(RF_dFF)
                %                 for mm = 1:20
                if length(RF_t) > length(RF_dFF{1})
                    rf_cr{j}=random_forest( RF_dFF{j}, RF_t(1:length(RF_dFF{1})), X, params.trigTime, trialMask, params );
                else
                    rf_cr{j}=random_forest( RF_dFF{j}, RF_t, X, params.trigTime, trialMask, params );
                end
                %                 end
            end
            toc
            warning('on','all')
            % check variation, the estimation is stable
            %             rsquare = zeros(80,20);
            %             predImp = zeros(80,7,20);
            %             for mm = 1:20
            %                 rsquare(:,mm) = rf_cr{mm}.rsquare;
            %                 predImp(:,:,mm) = rf_cr{mm}.predImp;
            %             end
            %
            %             tic
            %             parfor j=1:numel(10)
            %                 if length(cells.t) > length(cells.dFF{1})
            %                     rf_cr{j}=random_forest( cells.dFF{j}, cells.t(1:length(cells.dFF{1})), X, params.trigTime, trialMask, params );
            %                 else
            %                     rf_cr{j}=random_forest( cells.dFF{j}, cells.t, X, params.trigTime, trialMask, params );
            %                 end
            %             end
            %             toc
            %
            meanrSquare = zeros(1, length(rf_cr));
            for nn = 1:length(rf_cr)
                meanrSquare(nn) = nanmean(rf_cr{nn}.rsquare);
            end
            
            figure;plot(meanrSquare); % negative correlation with meanIntensity
            
            yyaxis right; hold on; plot(meanIntensity);
            hold on; plot(varIntensity);
            legend('R2','meanInt','varInt');
            print(gcf,'-dpng',fullfile(savefluofigpath,'RF-actionselection-r2meanvar'));
            saveas(gcf, fullfile(savefluofigpath,'actionselection-r2meanvar'), 'fig');
            
            % get mean predictor importance from all 196 ROIs
            meanPredImp = zeros(length(rf_cr{1}.regr_time),10);
            
            for npred = 1:10
                for tim = 1:length(rf_cr{1}.regr_time)
                    sumPredImp = zeros(1,length(rf_cr));
                    for nroi = 1:length(rf_cr)
                        sumPredImp(nroi) = rf_cr{nroi}.predImp(tim,npred);
                    end
                    meanPredImp(tim,npred) = nanmean(sumPredImp);% average across ROIs
                end
            end
            %
%                         figure;
%                         for mm = 1:10
%                             subplot(2,5,mm)
%                             plot(rf_cr{1}.regr_time, meanPredImp(:,mm))
%                         end
            xtitle = {'c(n)','c(n-1)','r(n)','r(n-1)','dQ','chosenQ','dK','chosenK','Reward rate','Cum. reward'};
            maxValue = max(max(meanPredImp));
            minValue = min(min(meanPredImp));
            figure;
            for mm = 1:10
                subplot(2,5,mm)
                plot(rf_cr{1}.regr_time, meanPredImp(:,mm))
                title(xtitle{mm});
                set(gca,'box','off');
                ylim([minValue-0.05,maxValue+0.05]);
            end
            print(gcf,'-dpng',fullfile(savefluofigpath,'RF-actionselection-avepredImp'));

            saveas(gcf, fullfile(savefluofigpath,'RF-actionselectionavepredImp'), 'fig');
            toc

            saveas(gcf, fullfile(savefluofigpath,'RF-actionselection-avepredImp'), 'fig');
            

            
%             figure;
%             for mm = 1:10
%                 subplot(2,5,mm)
%                 plot(rf_cr{1}.regr_time, rf_cr{1}.predImp(:,mm))
%                 title(xtitle{mm});
%                 set(gca,'box','off');
%                 ylim([minValue-0.05,maxValue+0.05]);
%             end
            %% save all regression into one mat file
            % save all these in a structure
            %save(saveRegName, 'reg_cr', 'reg_cr1','reg_cr2','reg_cr3','reg_cr_future','reg_cr_future_ctrl','reg_cr_ctrl');
            save(saveRFName,'rf_cr');
            %save(saveRegName_ITI,'reg_cr1_change_1','reg_cr2_change_1','reg_cr3_change_1');
            close all;
            %         else
            %             display('Random forest already done');
            %         end
        else
            display('Random forest already computed');
        end
    end
end
end
