function MP_GRAB_clustering(dataIndex);

% try: hierachical clustering; spectral clustering


nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    %load(fullfile(fn_beh.folder,fn_beh.name));
    
    
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
        
        warning('off','all')
        
        clusterMat = zeros(length(cells.normdFF),length(cells.normdFF{1}));
        for cc = 1:length(cells.normdFF)
            clusterMat(cc,:) = cells.normdFF{cc};
        end
        
        maxclust = 5;
        T = clusterdata(clusterMat,'Linkage','ward','SaveMemory','on','Maxclust',maxclust,'distance','correlation');
        
        % calculate correlation matrix
        corrMat = zeros(length(cells.normdFF));
         %regroup the cells with cluster results
         cellInd = 1:length(cells.normdFF);
         clustInd = [];
         clusterNum = zeros(1,maxclust);
         for gg=1:maxclust
             clustInd = [clustInd,cellInd(T==gg)];
             clusterNum(gg) = sum(T==gg);
         end
         
        for gg=1:length(clustInd)
            for yy = gg:length(clustInd)
                coeff = corrcoef(cells.normdFF{clustInd(gg)},cells.normdFF{clustInd(yy)});
                corrMat(gg,yy) = coeff(1,2);
            end
        end
        
        for yy=1:length(clustInd)
            for gg = yy+1:length(clustInd)
                
                corrMat(gg,yy) = corrMat(yy,gg);
            end
        end
        
         colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-0.2 0.2];
        
        z  =linkage(T,'ward');
        figure;image(corrMat,'CDataMapping','scaled')
        %hold on;dendrogram(z)
         axis square;
    colormap(colors);
    caxis([colorRange(1) colorRange(2)]);
    % plot the cluster
    hold on;
    for cc = 1:length(clusterNum)
        startInd = 1+sum(clusterNum(1:cc-1));
        endInd = sum(clusterNum(1:cc));
        plot([startInd endInd endInd startInd startInd],[startInd startInd endInd endInd startInd],'k:');
    end
     
% load regression
 saveCRName = fullfile(savematpath,'regCR_norm.mat');
         if exist(saveCRName)
            reg1 = load(saveCRName);
         end
        saveRLName = fullfile(savematpath,'regRL_norm.mat');
          if exist(saveRLName)
            reg2 = load(saveRLName);
         end
        saveRPEName = fullfile(savematpath,'regRL_RPE_norm.mat');
          if exist(saveRPEName)
            reg3 = load(saveRPEName);
          end
 
          for TT = 1:maxclust
              cluInd = cellInd(T==TT);
              for rr = 1:length(cluInd)
                 
                reg_cr{rr}.numPredictor = 15;
              reg_cr{rr}.nback = reg1.reg_cr{1}.nback;
               reg_cr{rr}.interaction = reg1.reg_cr{1}.interaction;
                reg_cr{rr}.regr_time = reg1.reg_cr{1}.regr_time;
                reg_cr{rr}.coeff = zeros(length(reg1.reg_cr{1}.regr_time),1);reg_cr{rr}.pval = zeros(length(reg1.reg_cr{1}.regr_time),1);
                % regression1
                % cn+1,cn,cn-1,rn,rn-1,xn,xn-1,ave_r,cum_r
                reg_cr{rr}.coeff = cat(2,reg_cr{rr}.coeff,reg1.reg_cr{cluInd(rr)}.coeff(:,[2,3,4,7,8,11,12,14,15]));
                reg_cr{rr}.pval = cat(2,reg_cr{rr}.pval,reg1.reg_cr{cluInd(rr)}.pval(:,[2,3,4,7,8,11,12,14,15]));
             % regression 2
             % dQ,chosenQ,dK,chosenK
              reg_cr{rr}.coeff = cat(2,reg_cr{rr}.coeff,reg2.reg_cr{cluInd(rr)}.coeff(:,8:11));
                reg_cr{rr}.pval = cat(2,reg_cr{rr}.pval,reg2.reg_cr{cluInd(rr)}.pval(:,8:11));
                % regression 3
                %RPE CKE
                 reg_cr{rr}.coeff = cat(2,reg_cr{rr}.coeff,reg3.reg_cr{cluInd(rr)}.coeff(:,[6,8]));
                reg_cr{rr}.pval = cat(2,reg_cr{rr}.pval,reg3.reg_cr{cluInd(rr)}.pval(:,[6,8]));
              end
              pvalThresh = 0.01;xtitle = 'Time from cue(s)';
                tlabel={'C(n+1)','C(n)','C(n-1)','R(n)', 'R(n-1)',...
                'C(n)*R(n)','C(n-1)*R(n-1)','Reward Rate','Cumulative Reward',...
                'dQ','chosenQ','dK','chosenK',...
                'RPE','CKE'};
            
              MP_plot_regr(reg_cr,[],pvalThresh,tlabel,xtitle);
          end
        %end
        
        % outcome regression mask
       % if ~exist('outcomeRegData','var') % if choice regression mask not computed
            label = 'outcome';
            outcomeInd = 7;
            outcome = selectRg2D(reg1.reg_cr,outcomeInd,[]);
            sigOutcome=outcome.coeff;
            sigOutcome(outcome.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,2,:) = sigOutcome;
        
        % x(n) regression mask
        %if ~exist('xnRegData','var') % if choice regression mask not computed
            label = 'interaction';
            xnInd = 11;
             interaction = selectRg2D(reg1.reg_cr,xnInd,[]);
            siginteraction=interaction.coeff;
            siginteraction(interaction.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,3,:) = siginteraction;
        %end
        
        % average reward regression mask
        %if ~exist('ave_rRegData','var') % if choice regression mask not computed
            label = 'average reward';
            ave_rInd = 14;
            ave_r = selectRg2D(reg1.reg_cr,ave_rInd,[]);
            sigave_r=ave_r.coeff;
            sigave_r(ave_r.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,4,:) = sigave_r;
        %end
        
        % cumulative reward regression mask
        %if ~exist('cum_rRegData','var') % if choice regression mask not computed
            label = 'cumulative reward';
            cum_rInd = 15;
            cum_r = selectRg2D(reg1.reg_cr,cum_rInd,[]);
            sigcum_r=cum_r.coeff;
            sigcum_r(cum_r.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,5,:) = sigcum_r;
        %end
     MP_plot_regr(reg_cr,[],params.pvalThresh,tlabel,params.xtitle);   
        %% regression 2-----------------------------------------------------------------------
        
        %if ~exist('dQRegData','var') % if choice regression mask not computed
            label = 'delta Q';
            dQInd = 8;
            dQ = selectRg2D(reg2.reg_cr,dQInd,[]);
            sigdQ=dQ.coeff;
            sigdQ(dQ.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,6,:) = sigdQ;
    %end
        %if ~exist('chosenQRegData','var') % if choice regression mask not computed

        %if ~exist('dKRegData','var') % if choice regression mask not computed
            label = 'delta K';
            dKInd = 10;
            dK = selectRg2D(reg2.reg_cr,dKInd,[]);
            sigdK=dK.coeff;
            sigdK(dK.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,7,:) = sigdK;
        %end
        
        
        %% regression 3-----------------------------------------------------------------------
        
        %if ~exist('RPERegData','var') % if choice regression mask not computed
            label = 'RPE';
            RPEInd = 6;
            RPE = selectRg2D(reg3.reg_cr,RPEInd,[]);
            sigRPE=RPE.coeff;
            sigRPE(RPE.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,8,:) = sigRPE;
        %end
        
        %if ~exist('CKERegData','var') % if choice regression mask not computed
            label = 'CKE';
            CKEInd = 8;
            CKE = selectRg2D(reg3.reg_cr,CKEInd,[]);
            sigCKE=CKE.coeff;
            sigRPE(CKE.pval>=sigThresh.alpha)= 0;
            clusterData(:,:,9,:) = sigCKE;
        %end
        
           %% set threshold
        Thresh.pos = 0.2; Thresh.neg = -0.2; % threshold to find area with pos/neg modulation
        sigThresh.value = 0.05; sigThresh.alpha = 0.01;  % threshold to find area with significant modulation
        
        tic
        
        %% check other variables using regression results (should focus on significance rather than coefficient)
        
        % variable list:
        % regression 1: cn,rn,cn+1,cn-1,rn-1,cn*rn,average_r,cumulative_r
        % regression 2: dQ,dK,chosenQ,chosenK
        % regression 3: RPE, CKE,
        saveregclusterpath = fullfile(savematpath,'cluster.mat');  % file path to save the results
        
        if exist(saveregclusterpath)
            load(saveregclusterpath);
        end
        
        % need to get the boundaries of possible coefficient value of all
        % regressions
       
         
          % data should be save in 28*28*n*80
          % or 748*n*80?
        clusterData = zeros(28,28,9,80); 
        % using choice,reward,interaction,aveR,cumR,dQ,dK,RPE,CKE for now
        %% regression 1 -------------------------------------------------------------------------
             % get the choice from linear regression
        %if ~exist('choiceRegData','var') % if choice regression mask not computed
           
        HCluster = reshape(clusterData,784,720);
        T = clusterdata(HCluster,'Linkage','ward','SaveMemory','on','Maxclust',20);
        
        % plot average coefficient for different clusters
        ClusterMat = zeros(28,28);
        for tt = 1:length(T)
            ind1 = floor((tt-1)/28)+1;ind2 =mod(tt,28);
            if ind2==0
                ind2 = 28;
            end
            ClusterMat(ind1,ind2) = T(tt);
        end
        figure;imagesc(ClusterMat)
        
        tree = linkage(HCluster,'average');

        figure;dendrogram(tree);
        
        rng('default');  % For reproducibility
X = rand(20000,3);
 T = clusterdata(X,'Linkage','ward','SaveMemory','on','Maxclust',4,'distance','correlation');
        
        %% save the results
        save(saveregmaskpath,'outcomeRegData','choiceRegData','cn_1RegData','cn__1RegData',...
                             'rn__1RegData', 'xnRegData', 'ave_rRegData', 'cum_rRegData',...
                             'dQRegData', 'chosenQRegData','dKRegData','chosenKRegData',...
                             'RPERegData', 'CKERegData')
    end
    
        close all
        clearvars -except analysispath dataIndex logfilepath root_path save_path_fluo
    toc
end
end