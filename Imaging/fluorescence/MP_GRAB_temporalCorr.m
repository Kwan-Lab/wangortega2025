function MP_GRAB_temporalCorr(dataIndex)

% calculate temporal correlation of different varible coefficients within
% grids

% check if coefficient sign changes over time

nFiles = size(dataIndex,1);

for ii = 1:nFiles

    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));

    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1

        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end

        cd(savefluofigpath);
   
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end

        %% set threshold
        sigThresh.value = 0.05; sigThresh.alpha = 0.01;  % threshold to find area with significant modulation

        tic
    
        %% check other variables using regression results (should focus on significance rather than coefficient)

        % variable list:
        % regression 1: cn,rn,cn+1,cn-1,rn-1,cn*rn,average_r,cumulative_r
        % regression 2: dQ,dK,chosenQ,chosenK
        % regression 3: RPE, CKE,
        saveregmaskpath = fullfile(savematpath,'regressionTempCorr.mat');  % file path to save the results

        if exist(saveregmaskpath)
            load(saveregmaskpath);
        end

        % need to get the boundaries of possible coefficient value of all
        % regressions
        saveCRName = fullfile(savematpath,'regCR_norm.mat');
         if exist(saveCRName)
            reg_cr{1} = load(saveCRName);
         end
        saveRLName = fullfile(savematpath,'regRL_norm.mat');
          if exist(saveRLName)
            reg_cr{2} = load(saveRLName);
         end
        saveRPEName = fullfile(savematpath,'regRL_RPE_norm.mat');
          if exist(saveRPEName)
            reg_cr{3} = load(saveRPEName);
          end  
         % get corresponding coefficient
         numGrids = numel(reg_cr{1}.reg_cr);
          varName = {'cn+1','cn','cn-1','rn','rn-1','xn','xn-1','ave_r','Cum_r','dQ','chosenQ','dK','chosenK','RPE','CKE'};
        varInd = [2,3,4,7,8,11,12,14,15,8,9,10,11,6,8];
        varReg = [1,1,1,1,1,1,1,1,1,2,2,2,2,3,3];
        zerolag = length(reg_cr{1}.reg_cr{1}.regr_time);
       
        tempCorrMat = cell(length(varInd),length(varInd));
        corrCoeff0 = zeros(length(varInd),length(varInd),numGrids);
        corrMax = zeros(length(varInd),length(varInd),numGrids);
        corrLag = zeros(length(varInd),length(varInd),numGrids);
        % calculate crosscorrelation
        %if ~exist(saveregVarcorrpath)
        for uu = 1:length(varName)
            for vv=uu:length(varName)

                tempCorrMat{uu,vv} = getRegTempCorr(reg_cr{varReg(uu)},reg_cr{varReg(vv)},varName{uu},varName{vv},varInd(uu),varInd(vv), sigThresh, savefluofigpath);
                corrCoeff0(uu,vv,:) = squeeze(tempCorrMat{uu,vv}(zerolag,1,:));
                [corrMax(uu,vv,:),maxInd] = max(squeeze(tempCorrMat{uu,vv}(:,1,:)),[],1);
                lagVec = squeeze(tempCorrMat{uu,vv}(:,2,:));
                corrLag(uu,vv,:) = lagVec(maxInd);
            end
        %end
        end
        
        % save the data
        save(saveregmaskpath,'tempCorrMat','corrCoeff0','corrMax','corrLag');

%% with in variables, calculate temporal correlation with a same reference grid
    saveregpath = fullfile(savematpath,'regressionVarTempCorr.mat');  % file path to save the results

        %% regression 1 -------------------------------------------------------------------------
             % get the choice from linear regression
        %if ~exist('choicetempData','var') % if choice regression mask not computed
            label = 'choice';
            choiceInd = 3;

            choicetempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,choiceInd,[0 3],sigThresh,savefluofigpath);
         % end

        % outcome regression mask
        
        if ~exist('outcometempData','var') % if choice regression mask not computed
            label = 'outcome';
            outcomeInd = 7;
            outcometempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,outcomeInd,[0,3],sigThresh,savefluofigpath);
          tlabel1='Outcome coefficient';
             xtitle='Time from cue(s)';colorRange = [-0.05 0.05];
            sigGrid =~isnan(outcometempData.tempCorrCoeff);
            outcomeSel.t = reg_cr{1}.reg_cr{1}.regr_time;
            %outcomeSel.coeff= outcometempData.coeff(sigGrid(:),:);
            %outcomeSel.lag = outcometempData.tempCorrLag(sigGrid);
            outcomeSel.coeff= outcometempData.coeff;
            %outcomeSel.lag = outcometempData.tempCorrLag;
       %plot_coeff_sort(outcomeSel,[0 3],tlabel1,xtitle,colorRange,savefluofigpath);

        %cluster grid based on the correlation within 0.5-4s
        allInd = 1:size(outcometempData.coeff,1);
        sigInd = allInd(sigGrid(:));
        outcomeSel.sigInd = sigInd;
        outcomeSel.sigMat = sigGrid;
        [clustId, oriInd] = regCoef_cluster(outcomeSel, [-3, 5], 2, tlabel1,xtitle,savefluofigpath);
         
%           outcome.t = reg_cr{1}.reg_cr{1}.regr_time;
%             outcome.coeff= outcometempData.coeff;
%            clustId = regCoef_cluster(outcomeSel, [-3, 5], 3, tlabel1,xtitle,savefluofigpath);

        saveDataOutcome.clustID = clustId;
        saveDataOutcome.oriInd = oriInd;
        saveDataOutcome.T = zeros(1,length(oriInd));
        clusterNum = [sum(saveDataOutcome.clustID==1), sum(saveDataOutcome.clustID==2), sum(saveDataOutcome.clustID==3)];
         saveDataOutcome.T(1:clusterNum(1)) = 1;
         saveDataOutcome.T(clusterNum(1)+1:sum(clusterNum(1:2))) = 2;
         saveDataOutcome.T(sum(clusterNum(1:2))+1:end) = 3;
        saveDataOutcome.t = outcomeSel.t;
        saveDataOutcome.coeff = outcomeSel.coeff;
        saveDataOutcome.sigInd = outcomeSel.sigInd;
%    
            %end
        % save the result

         % c(n+1) regression mask
%         if ~exist('cn_1tempData','var') % if choice regression mask not computed
            label = 'upcoming choice';
            cn_1Ind = 2;
            cn_1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,cn_1Ind,[-3 0],sigThresh,savefluofigpath);
%         end
% 
%         % r(n+1) regression mask
%         if ~exist('rn_1tempData','var') % if choice regression mask not computed
            label = 'upcoming reward';
            rn_1Ind = 6;
            rn_1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,rn_1Ind,[-3 0],sigThresh,savefluofigpath);

% 
%         % c(n-1) regression mask
%         if ~exist('cn__1tempData','var') % if choice regression mask not computed
            label = 'previous choice';
            cn__1Ind = 4;
            cn__1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,cn__1Ind,[-3 0],sigThresh,savefluofigpath);
 %       end
  label = 'previous+1 choice';
            cn__2Ind = 5;
            cn__2tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,cn__1Ind,[-3 0],sigThresh,savefluofigpath);

% 
%         % r(n-1) regression mask
%         if ~exist('rn__1tempData','var') % if choice regression mask not computed
            label = 'previous reward';
            rn__1Ind = 8;
            rn__1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,rn__1Ind,[-3 0],sigThresh,savefluofigpath);
%         end
label = 'previous+1 reward';
            rn__2Ind = 9;
            rn__2tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,rn__1Ind,[-3 0],sigThresh,savefluofigpath);

% 
%         % x(n) regression mask
 %       if ~exist('xntempData','var') % if choice regression mask not computed
            label = 'interaction';
            
            xnInd = 11;
            xntempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,xnInd,[0 3],sigThresh,savefluofigpath);
% 
%          if ~exist('xn__1tempData','var')
 label = 'upcoming interaction';
            xn_1Ind = 10;
            xn_1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,xn_1Ind,[-3 0],sigThresh,savefluofigpath);


         label = 'previous interaction';
            xn__1Ind = 12;
            xn__1tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,xn__1Ind,[-3 0],sigThresh,savefluofigpath);
%          end
 label = 'previous+1 interaction';
            xn__2Ind = 13;
            xn__2tempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,xn__2Ind,[-3 0],sigThresh,savefluofigpath);

% 
%         % average reward regression mask
%         if ~exist('ave_rtempData','var') % if choice regression mask not computed
            label = 'average reward';
            ave_rInd = 14;
            ave_rtempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,ave_rInd,[-3 0],sigThresh,savefluofigpath) ;
%         end
% 
%         % cumulative reward regression mask
%         if ~exist('cum_rtempData','var') % if choice regression mask not computed
            label = 'cumulative reward';
            cum_rInd = 15;
            cum_rtempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,cum_rInd,[-3 0],sigThresh,savefluofigpath) ;
%         end
% 
%         %% regression 2-----------------------------------------------------------------------
%         if ~exist('dQtempData','var') % if choice regression mask not computed
            label = 'delta Q';
            dQInd = 8;
            dQtempData = getRegautoCorrData(reg_cr{2}.reg_cr,label,dQInd,[0,3],sigThresh,savefluofigpath) ;
%         end
% 
%         if ~exist('chosenQtempData','var') % if choice regression mask not computed
%             label = 'chosen Q';
%             chosenQInd = 9;
%             chosenQtempData = getRegautoCorrData(reg_cr{2}.reg_cr,label,chosenQInd, sigThresh,savefluofigpath) ;
%         end
% 
%         if ~exist('dKtempData','var') % if choice regression mask not computed
            label = 'delta K';
            dKInd = 10;
            dKtempData = getRegautoCorrData(reg_cr{2}.reg_cr,label,dKInd,[0,3],sigThresh,savefluofigpath) ;
%         end
% 
%         if ~exist('chosenKtempData','var') % if choice regression mask not computed
%             label = 'chosen K';
%             chosenKInd = 11;
%             chosenKtempData = getRegautoCorrData(reg_cr{2}.reg_cr,label,chosenKInd, sigThresh,savefluofigpath) ;
%         end
% 
%         %% regression 3-----------------------------------------------------------------------

% 
%         if ~exist('CKEtempData','var') % if choice regression mask not computed
            label = 'CKE';
            CKEInd = 8;
            CKEtempData = getRegautoCorrData(reg_cr{3}.reg_cr,label,CKEInd,[0,3],sigThresh,savefluofigpath) ;
%         end
% 
%      pos/neg RPE use all grids with coefficient
        label = 'posRPE';
        RPEInd = 6;
        posRPEtempData = getRegautoCorrData(reg_cr{3}.reg_cr_pos,label,RPEInd,[0 3],sigThresh,savefluofigpath) ;
          label = 'negRPE';
        RPEInd = 6;
        negRPEtempData = getRegautoCorrData(reg_cr{3}.reg_cr_neg,label,RPEInd,[0 3],sigThresh,savefluofigpath) ;

%         %% save the results
        save(saveregpath,'cn_1tempData','cn__1tempData','cn__2tempData',...
                             'rn_1tempData','rn__1tempData', 'rn__2tempData',...
                             'xn_1tempData','xn__1tempData', 'xn__2tempData',...
                             'ave_rtempData', 'cum_rtempData');

%% based on the outcome result, plot pos/neg RPE coefficient, choice/interaction of the corresponding results
t =  reg_cr{1}.reg_cr{1}.regr_time;
nCells = length(saveDataOutcome.oriInd);

figure;
subplot(2,3,1)
tlabel = 'Outcome';
image(t,1:nCells,outcometempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel);

subplot(2,3,4)

line1 = outcometempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = outcometempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  outcometempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');

% plot pos RPE
subplot(2,3,2)
tlabel = 'posRPE';
image(t,1:nCells, posRPEtempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel);

subplot(2,3,5)

line1 = posRPEtempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = posRPEtempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  posRPEtempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');

% plot neg RPE
subplot(2,3,3)
tlabel = 'negRPE';
image(t,1:nCells, negRPEtempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel)
subplot(2,3,6)

line1 = negRPEtempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = negRPEtempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  negRPEtempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');

savepath = fullfile(savefluofigpath,'clustering');
if ~exist(savepath)
    mkdir(savepath)
end

print(gcf,'-dpng',fullfile(savepath,[tlabel,' outcome-RPE-cluster']));
saveas(gcf, fullfile(savepath,[tlabel,' outcome-RPE-cluster']), 'fig');
saveas(gcf, fullfile(savepath,[tlabel,' outcome-RPE-cluster']), 'svg');

%% outcome, choice, interaction
figure;
subplot(2,3,1)
tlabel = 'Outcome';
image(t,1:nCells,outcometempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel);

subplot(2,3,4)

line1 = outcometempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = outcometempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  outcometempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');
% plot choice
subplot(2,3,2)
tlabel = 'Choice';
image(t,1:nCells, choicetempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel);
subplot(2,3,5)

line1 = choicetempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = choicetempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  choicetempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');


% plot interaction
subplot(2,3,3)
tlabel = 'Interaction';
image(t,1:nCells,xntempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel);
subplot(2,3,6)

line1 = xntempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = xntempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  xntempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');

print(gcf,'-dpng',fullfile(savepath,[tlabel,' outcome-choice-cluster']));
saveas(gcf, fullfile(savepath,[tlabel,' outcome-choice-cluster']), 'fig');
saveas(gcf, fullfile(savepath,[tlabel,' outcome-choice-cluster']), 'svg');
%         

%% outcome, dQ, dK
%% outcome, choice, interaction
figure;
subplot(2,3,1)
tlabel = 'Outcome';
image(t,1:nCells,outcometempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel);

subplot(2,3,4)

line1 = outcometempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = outcometempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  outcometempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');
% plot dQ
subplot(2,3,2)
tlabel = 'dQ';
image(t,1:nCells, dQtempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel);
subplot(2,3,5)

line1 = dQtempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = dQtempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  dQtempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');


% plot dK
subplot(2,3,3)
tlabel = 'dK';
image(t,1:nCells,dKtempData.coeff(saveDataOutcome.oriInd,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.1 0.1]);      %normalize dF/F heatmap to max of all conditions
title(tlabel);
subplot(2,3,6)

line1 = dKtempData.coeff(saveDataOutcome.oriInd(1:clusterNum(1)),:);
ste = nanstd(line1,0,1)/sqrt(size(line1,1));
plot(t,nanmean(line1,1),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,1)-ste,nanmean(line1,1)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = dKtempData.coeff(saveDataOutcome.oriInd(clusterNum(1)+1:sum(clusterNum(1:2))),:);
ste = nanstd(line2,0,1)/sqrt(size(line2,1));
plot(t,nanmean(line2,1),'k');
errorshade(t,nanmean(line2,1)-ste,nanmean(line2,1)+ste,[0 0 0]/255,0.5);
hold on;
line3 =  dKtempData.coeff(saveDataOutcome.oriInd(sum(clusterNum(1:2))+1:end),:);
ste = nanstd(line3,0,1)/sqrt(size(line3,1));
plot(t,nanmean(line3,1), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,1)-ste,nanmean(line3,1)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');

print(gcf,'-dpng',fullfile(savepath,[tlabel,' outcome-dQ-cluster']));
saveas(gcf, fullfile(savepath,[tlabel,' outcome-dQ-cluster']), 'fig');
saveas(gcf, fullfile(savepath,[tlabel,' outcome-dQ-cluster']), 'svg');
%         
close all;
     save(fullfile(savematpath,'cluster.mat'), 'saveDataOutcome','choicetempData','xntempData','posRPEtempData','negRPEtempData','dQtempData','dKtempData','CKEtempData');
    end

        close all
        clearvars -except analysispath dataIndex logfilepath root_path save_path_fluo
    toc
end
end


