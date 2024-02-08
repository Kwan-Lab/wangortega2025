function MP_GRAB_tonic(dataIndex)

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
        saveRFName = fullfile(savematpath,'tonic.mat');  % random forest
        % saveRegName_ITI =  fullfile(savematpath,[fn_beh.name(1:end-7),'regCR_ITI.mat']);
        %
        
        if ~exist(saveRFName)
            
            
            % examine the relationship of tonic NE fluorescence and
            % reward rate
            binSize = 10:10:100;
            
             % average reward rate on 20 trials window, results shows not
            % important
           
            averageReward = NaN(length(trials.go),length(binSize));
            for rr = 1:length(binSize)
                for kk = 1:length(trials.left)
                    if kk <= binSize(rr)
                        averageReward(kk,rr) = sum(trials.reward(1:kk))/kk;
                    else
                       averageReward(kk,rr) = sum(trials.reward(kk-binSize(rr)+1:kk))/binSize(rr);
                    end
                end
            end
            
            dffBaseline = zeros(length(cells.dFF),length(trialData.cue));
            for tt = 1:length(trialData.cue)
                for gg = 1:length(cells.dFF)
                    dffBaseline(gg,tt) = mean(cells.dFF{gg}(cells.t<trialData.cueTimes(tt) & cells.t>trialData.cueTimes(tt)-2));
                end
            end
            
            RF_t = cells.t;
            RF_dFF = cells.dFF;
            for j=1:numel(RF_dFF)
                %                 for mm = 1:20
                % calculate the correlation 
                %                 end
                figure;
                plot(trialData.cueTimes,averageReward(:,4));
                yyaxis right
                hold on;
                plot(trialData.cueTimes,smooth(dffBaseline(1,:),10))
                figure;
                scatter(averageReward(:,4),mean(dffBaseline,1));
                [R,P] = corrcoef(averageReward(:,4),mean(dffBaseline,1),'rows','pairwise')
                plot(cells.t,lowpass(cells.cellf{1},0.5,15))
                xlim([0 1000])
            end
            
       
            
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
