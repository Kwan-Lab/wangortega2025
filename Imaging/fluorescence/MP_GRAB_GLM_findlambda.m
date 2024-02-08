function MP_GRAB_GLM_findlambda(dataIndex)

%% manuelly find lambda by looking at the average correlation plot

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
        saveRFName = fullfile(savematpath,'GLM_CR.mat');  % random forest
        % saveRegName_ITI =  fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);
        %
        
        %if ~exist(saveRFName)
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
           choice = NaN(size(trials.left));
            choice(trialData.response == 2) = 0;
            choice(trialData.response == 3) = 1;
            
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
            %
            %             tlabel={'C(n+1)','C(n)','C(n-1)','C(n-2)','R(n+1)','R(n)', 'R(n-1)','R(n-2)',...
            %                 'C(n+1)*R(n+1)','C(n)*R(n)','C(n-1)*R(n-1)','C(n-2)*R(n-2)','Reward Rate','Cumulative Reward'};
            %
            % reg_cr_future=linear_regr( fluo.dia, fluo.t, future_event, params.trigTime, trialMask, params );
            
            glm_t = cells.t;
            glm_dFF = cells.dFF;
            %parfor j=1:numel(glm_dFF)
            % randomly pick 20 ROIs for determining lambda
            runningInd = randsample(numel(glm_dFF),40);
            tic
            parfor j=1:numel(runningInd)
                %                 for mm = 1:20
                if length(glm_t) > length(glm_dFF{1})
                    greg_cr{j}=generalized_linear_lambda( glm_dFF{runningInd(j)}, glm_t(1:length(glm_dFF{1})), future_event, params.trigTime, trialMask, params );
                else
                    greg_cr{j}=generalized_linear_lambda( glm_dFF{runningInd(j)}, glm_t, future_event, params.trigTime, trialMask, params );
                end
                %                 end
            end
            toc
            testLambda=[0.05:0.05:0.6];
            meanCorr = zeros(numel(testLambda),1);
            for j = 1:numel(runningInd)
%                 
%                 figure;
%                 plot(testLambda,mean(greg_cr{j},2))
                meanCorr = (meanCorr + mean(greg_cr{j},2)/numel(runningInd));
            end
            
            figure;plot(testLambda,meanCorr);
xlabel('Lambda');
ylabel('Average correlation of y and yhat');
print(gcf,'-dpng','ridge-lambda');    %png format

            
        %end
    end
end
end
