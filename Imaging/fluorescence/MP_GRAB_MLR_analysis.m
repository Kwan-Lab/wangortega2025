function MP_GRAB_MLR_analysis(dataIndex)

%% reference to Sul et al.2011
% analyze the MLR results
% 1. plot field of view, mark ROI - significant variance 
% 2. amount of variance explained
% 3. Patriachi et. al. Science 2018. plot df/f 



nFiles = size(dataIndex,1);

%% go through every session, run analysis
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    load(fullfile(fn_fluo.folder,fn_fluo.name));
    
     savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
    fn_MLR = dir(fullfile(savematpath,'regCR.mat'));  % regression for fluo change
    
     params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = NaN;  %confidence interval
        params.minNumTrial = 50;
        
    if length(fn_MLR) == 1
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
     
        load(fullfile(fn_MLR.folder,fn_MLR.name));
      
        % save figure path
        %% plot the ROI grid with significance (consider only c(n), r(n), x(n))
        % get the ROI data
        fn_ROI = dir(fullfile(dataIndex.LogFilePath{ii},'cell*.mat'));
        
        % load the field of view
        fn_FOV = dir(fullfile(dataIndex.BehPath{ii},'reg_stackMean.tif'));
        Img_FOV = imread(fullfile(fn_FOV.folder, fn_FOV.name));
        
        %% plot the variance explained by the full model, against the mean/var of the ROI intensity
        % R^2: variance explained by the full model
        % increment R^2 or anova to get variance explained by a single
        % model
        nBins = size(reg_cr{1}.coeff,1);
        sumSq = zeros(nBins, length(fn_ROI));
        sumSq_total = zeros(nBins, length(fn_ROI));
        mseList = zeros(1, length(fn_ROI));
        for tt = 1:length(fn_ROI)
            sumSq(:,tt) = sum(reg_cr{tt}.SumSq(:,2:end),2);
            sumSq_total(:,tt) = sum(reg_cr{tt}.SumSq(:,1:end),2);
            mseList(tt) = reg_cr{tt}.mse;
        end
        
        % get mean/var of projection
        meanProjPath = dir(fullfile(dataIndex.BehPath{ii},'*Mean*.tif'));
        meanProj = loadtiffseq([],fullfile(meanProjPath.folder,meanProjPath.name));
    
        % get the mean intensity in ROI
        meanIntensity = zeros(1, size(fn_ROI,1));
        varIntensity = zeros(1,size(fn_ROI,1));
        for uu = 1:size(fn_ROI,1)
            roi = load(fullfile(fn_ROI(uu).folder,fn_ROI(uu).name));
            meanIntensity(uu) = mean(meanProj(roi.bw==1));
            varIntensity(uu) = nanstd(cells.cellf{uu});
        end
        
        %% get mean coefficient
        all_coeff = [];
        for cc = 1:length(reg_cr)
            all_coeff = cat(3,all_coeff, reg_cr{cc}.coeff);
        end
        reg_cr_c.coeff= all_coeff;
        reg_cr_c.regr_time = reg_cr{1}.regr_time;
        reg_cr_c.numPredictor = 14;
        reg_cr_c.nback = 0;
        reg_cr_c.interaction = 0;
    % use bootstrp to get coefficient
        reg_cr_c = getBootstrp(reg_cr_c, 0, 0.05);
        xtitle='Time from cue (s)';
        tlabel={'c(n+1)','c(n)','c(n-1)','c(n-2)','r(n+1)','r(n)', 'r(n-1)','r(n-2)',...
        'c(n+1)*r(n+1)','c(n)*r(n)','c(n-1)*r(n-1)','c(n-2)*r(n-2)','Reward Rate','Cumulative Reward'};
        MP_plot_regrcoef_pupil(reg_cr_c, 0.05, tlabel,xtitle);
        
        sigAlpha = 0.01;
        
        sigCn = zeros(1,length(fn_ROI));
        sigRn = zeros(1,length(fn_ROI));
        sigXn = zeros(1,length(fn_ROI));
        for jj = 1:length(fn_ROI)
            roi = load(fullfile(fn_ROI(jj).folder,fn_ROI(jj).name));
            
            %check Cn
            sigCn(jj) = sum(sum(reg_cr{jj}.pval(:,3)<sigAlpha));
            sigRn(jj) = sum(sum(reg_cr{jj}.pval(:,7)<sigAlpha));
            sigXn(jj) = sum(sum(reg_cr{jj}.pval(:,11)<sigAlpha));
            % check if c/r/x is significant (choice/reward/interaction)
            % get mean df/f for different trial types, run wilcoxon ranksum
            % test to determine whether it is modulated by certain variable
            %check pval
            %save(saveRegName_ITI,'reg_cr1_change_1','reg_cr2_change_1','reg_cr3_change_1');
%             fieldname{1}={'left'};fieldname{2}= {'right'};
%             % try wilcoxon ranksum test?
%             for kk=1:numel(fieldname)
%                     trialMask = getMask(trials,fieldname{kk});
%                     psth_panel.sig{kk} = get_psth( cells.dFF{jj}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
%             end
%             [p_choice,~] = ranksum(psth_panel.sig{1}.signal, psth_panel.sig{2}.signal);
%             
%              fieldname{1}={'reward'};fieldname{2}= {'noreward'};
%             % try wilcoxon ranksum test?
%             for kk=1:numel(fieldname)
%                     trialMask = getMask(trials,fieldname{kk});
%                     psth_panel.sig{kk} = get_psth( cells.dFF{jj}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
%             end
%             [p_reward,~] = ranksum(psth_panel.sig{1}.signal, psth_panel.sig{2}.signal);
%             
        end
        
        % use 5 time points as threshold for now (if there are more than 5
        % time points (0.5s) with significant correlation coefficient, the 
        
        % black: modulate only by choice; red: modulate only by reward;
        % yellow: modulate by both\
        threshold = 5;
        figure; imagesc(meanProj);
        axis equal
        set(gca,'box','off', 'Xcolor', 'none', 'Ycolor', 'none');
        
        % also ignore the ROIs with mean intensity below 45 percentile
        for tt = 1:length(fn_ROI)
            if meanIntensity(tt) > prctile(meanIntensity,45)
             roi = load(fullfile(fn_ROI(tt).folder,fn_ROI(tt).name));
             [B,L]=bwboundaries(roi.bw,'noholes');
             if sigCn(tt) >= threshold && sigRn(tt) < threshold
                 hold on; plot(B{1}(:,2),B{1}(:,1),'k','LineWidth',1);
             elseif sigCn(tt) < threshold && sigRn(tt) >= threshold
                 hold on; plot(B{1}(:,2),B{1}(:,1),'r','LineWidth',1);
             elseif sigCn(tt) >= threshold && sigRn(tt) >= threshold
                 hold on; plot(B{1}(:,2),B{1}(:,1),'w','LineWidth',1);
             end
            end
        end
        print(gcf,'-dpng',fullfile(savefluofigpath,'FOV-ROICR'));
            saveas(gcf, fullfile(savefluofigpath,'FOV-ROICR'), 'fig');
        
            figure; imagesc(meanProj);
        axis equal
        set(gca,'box','off', 'Xcolor', 'none', 'Ycolor', 'none');
        
        % also ignore the ROIs with mean intensity below 33 percentile
        for tt = 1:length(fn_ROI)
            if meanIntensity(tt) > prctile(meanIntensity,45)
             roi = load(fullfile(fn_ROI(tt).folder,fn_ROI(tt).name));
             [B,L]=bwboundaries(roi.bw,'noholes');
            
                 hold on; plot(B{1}(:,2),B{1}(:,1),'k','LineWidth',1);
             
            end
        end
        print(gcf,'-dpng',fullfile(savefluofigpath,'FOV-ROIincluded'));
            saveas(gcf, fullfile(savefluofigpath,'FOV-ROIincluded'), 'fig');
        
            close all;

        end
    end
end

