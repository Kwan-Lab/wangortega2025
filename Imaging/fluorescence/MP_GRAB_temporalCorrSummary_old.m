function MP_GRAB_temporalCorrSummary(dataIndex, savesumfigpath,savematsumpath)

% calculate temporal correlation of different varible coefficients within
% grids

% check if coefficient sign changes over time
% for eachgrid, calculate the number of + time point and - time point
% then average

nFiles = size(dataIndex,1);

%% get clustering results, average similar clusters, also rerun clustering on all grids together
rCoeff = cell(0); rClust = cell(0);
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


        %load(fullfile(fn_fluo.folder,fn_fluo.name));

        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end

        % outcome
        tlabel1 = 'Outcome coefficient';
        saveregpath =  fullfile(savematpath,[tlabel1,' cluster.mat']);  % file path to save the results
        if exist(saveregpath)
            load(saveregpath)
        end

        rCoeff{ii} = saveData.coeff;
        rClust{ii} = saveData.clustID;
        rt = saveData.t;

    end
end

% trim the clustering.
% 1. delete cluster with too few grids . <10 grids?
% 2. combine goup with high correlation >0.9?
g1=cell(0); g2 = cell(0); g3 = cell(0);
clustCoeff = []; clustMean = [];
rClust_new = cell(0);
rCoeff_new = cell(0);
newInd = 1;
rHPeakT = [];
rLPeakT = [];
rDecayTau = [];
oriGroup = [];
% save all grids that cannot be grouped into the 3 groups



t0 = 0; t1 = 4;  % time range for correlation calculation

for gg = 1:nFiles
    g1Num = sum(rClust{gg}==1); g2Num = sum(rClust{gg}==2); g3Num = sum(rClust{gg}==3);
    % calculate the average coefficient with in 1-2.5s for each group
    g1Mean =nanmean(rCoeff{gg}(rClust{gg}==1,:),1);  g2Mean = nanmean(rCoeff{gg}(rClust{gg}==2,:),1); g3Mean =nanmean(rCoeff{gg}(rClust{gg}==3,:),1);

    %clustMean = [clustMean;g1Mean;g2Mean;g3Mean];
    %% smooth the mean and find peak and decay tau
    % find 1): highest peak and time;
    % 2): lowest peak and time;
    % 3) exponential decay tau after high peak



    % (islocalmin & islocalmax)
    tInd = 1:length(rt);
    searchInd = tInd(rt>0);
    % group 1
    smoothed1 = smooth(g1Mean);
    [minValue1,minInd] = nanmin(smoothed1(rt>0));
    minInd1 = searchInd(minInd);
    [maxValue1,maxInd] = nanmax(smoothed1(rt>0));
    maxInd1 = searchInd(maxInd);
    rHPeakT = [rHPeakT, rt(maxInd1)];
    rLPeakT = [rLPeakT, rt(minInd1)];
    oriGroup = [oriGroup,1];
    % fit for tau
    %     if maxInd1<length(rt)
    %         f1 = fit(rt(maxInd1:end),smoothed1(maxInd1:end),'exp1');
    %         rDecayTau = [rDecayTau, f1.b];
    %     else
    %         f1 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==1,:);
    newInd = newInd+1;

    % group 2
    smoothed2 = smooth(g2Mean);
    [minValue2,minInd] = nanmin(smoothed2(rt>0));
    minInd2 = searchInd(minInd);
    [maxValue2,maxInd] = nanmax(smoothed2(rt>0));
    maxInd2 = searchInd(maxInd);
    rHPeakT = [rHPeakT, rt(maxInd2)];
    rLPeakT = [rLPeakT, rt(minInd2)];
    oriGroup = [oriGroup,2];
    % fit for tau
    %     if maxInd2<length(rt)
    %         f2 = fit(rt(maxInd2:end),smoothed2(maxInd2:end),'exp1');
    %         rDecayTau = [rDecayTau, f2.b];
    %     else
    %         f2 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==2,:);
    newInd = newInd+1;

    % group 3
    smoothed3 = smooth(g3Mean);
    [minValue3,minInd] = nanmin(smoothed3(rt>0));
    minInd3 = searchInd(minInd);
    [maxValue3,maxInd] = nanmax(smoothed3(rt>0));
    maxInd3 = searchInd(maxInd);
    rHPeakT = [rHPeakT, rt(maxInd3)];
    rLPeakT = [rLPeakT, rt(minInd3)];
    oriGroup = [oriGroup,3];
    % fit for tau
    %     if maxInd3<length(rt)
    %         f3 = fit(rt(maxInd3:end),smoothed3(maxInd3:end),'exp1');
    %         rDecayTau = [rDecayTau, f3.b];
    %     else
    %         f3 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==3,:);
    newInd = newInd+1;

    % plot the fit results
    savefigpath =  fullfile(dataIndex.BehPath{gg},'figs-fluo','clustProp');
    if ~exist(savefigpath)
        mkdir(savefigpath)
    end

    % colors = [[241, 84, 18]/255;[0, 0, 0];[52, 179, 241]/255 ];
    % figure;plot(rt,smooth(g1Mean),'Color',colors(1,:));hold on;plot(rt,g1Mean,'--','Color',colors(1,:));
    % hold on;scatter(rt(minInd1),minValue1,160, 'filled');scatter(rt(maxInd1),maxValue1,160, 'filled');
    % if ~isempty(f1)
    %     plot(f1,rt(maxInd1:end),smoothed1(maxInd1:end)); legend('hide');
    % end
    %
    % hold on;plot(rt,smooth(g2Mean),'Color',colors(2,:));hold on;plot(rt,g2Mean,'--','Color',colors(2,:))
    % hold on;scatter(rt(minInd2),minValue2,160, 'filled');scatter(rt(maxInd2),maxValue2,160, 'filled');
    % if ~isempty(f2)
    % plot(f2,rt(maxInd2:end),smoothed2(maxInd2:end)); legend('hide');
    % end
    %
    %  plot(rt,smooth(g3Mean),'Color',colors(3,:));hold on;plot(rt,g3Mean,'--','Color',colors(3,:))
    %  hold on;scatter(rt(minInd3),minValue3,160, 'filled');scatter(rt(maxInd3),maxValue3,160, 'filled');
    %  if ~isempty(f3)
    % plot(f3,rt(maxInd3:end),smoothed3(maxInd3:end)); legend('hide');
    %  end
    % print(gcf,'-dpng',fullfile(savefigpath,tlabel1));
    % close;


    %% not working, group groups from individual session with a reference
    % ACh, first session as ref group
    %     if gg  == 1
    %         ref1 = g1Mean;
    %         ref2 = g2Mean;
    %         ref3 = g3Mean;
    %         rClust_new{gg} = rClust{gg};
    %         rCoeff_new{gg} = rCoeff{gg};
    %     end
    %     figure;plot(rt,ref1);hold on;
    %     plot(rt,ref2);plot(rt,ref3);
    %     plot(rt,g1Mean)
    %
    %      figure;plot(rt,ref1);hold on;
    %     plot(rt,ref2);plot(rt,ref3);
    %     plot(rt,g3Mean)
    %
    %      figure;plot(rt,ref1);hold on;
    %     plot(rt,ref2);plot(rt,ref3);
    %     plot(rt,g2Mean)
    %     if gg > 1
    %         corr11 = corrcoef(g1Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr12 = corrcoef(g1Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr13 = corrcoef(g1Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
    %          corr21 = corrcoef(g2Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr22 = corrcoef(g2Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr23 = corrcoef(g2Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
    %         corr31 = corrcoef(g3Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr32 = corrcoef(g3Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr33 = corrcoef(g3Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
    %
    %         corr1 = [corr11(1,2),corr12(1,2),corr13(1,2)];
    %         corr2 = [corr21(1,2),corr22(1,2),corr23(1,2)];
    %         corr3 = [corr31(1,2),corr32(1,2),corr33(1,2)];
    %
    %          [maxCorr,ind] = max(corr1);
    %          if maxCorr>0.75 & g1Num > 25
    %              rClust_new{gg} = [rClust_new{gg}; ones(sum(rClust{gg}==1),1)*ind];
    %              rCoeff_new{gg} = [rCoeff_new{gg}; rCoeff{gg}(rClust{gg}==1,:)];
    %          else
    %              rClust_new{gg} = [rClust_new{gg}; NaN(sum(rClust{gg}==1),1)];
    %              outlierGroups{outInd} = rCoeff{gg}(rClust{gg}==1,:);
    %              outInd = outInd+1;
    %          end
    %          [maxCorr,ind] = max(corr2);
    %          if maxCorr>0.75 & g2Num > 25
    %               rClust_new{gg} = [rClust_new{gg}; ones(sum(rClust{gg}==2),1)*ind];
    %              rCoeff_new{gg} = [rCoeff_new{gg}; rCoeff{gg}(rClust{gg}==2,:)];
    %          else
    %                rClust_new{gg} = [rClust_new{gg}; NaN(sum(rClust{gg}==2),1)];
    %              outlierGroups{outInd} = rCoeff{gg}(rClust{gg}==2,:);
    %              outInd = outInd+1;
    %          end
    %          [maxCorr,ind] = max(corr3);
    %          if maxCorr>0.75 & g3Num > 25
    %                rClust_new{gg} = [rClust_new{gg}; ones(sum(rClust{gg}==3),1)*ind];
    %              rCoeff_new{gg} = [rCoeff_new{gg}; rCoeff{gg}(rClust{gg}==3,:)];
    %          else
    %                rClust_new{gg} = [rClust_new{gg}; NaN(sum(rClust{gg}==3),1)];
    %              outlierGroups{outInd} = rCoeff{gg}(rClust{gg}==3,:);
    %              outInd = outInd+1;
    %          end
    %     end

    clustCoeff = [clustCoeff;rCoeff{gg}];
    % how to determine the group?
    % in ACh use the first session as the reference group

end

% average
group1 = []; group2 = []; group3 = []; %group4 = [];

% if the surface area between (0,2) less than 0, group1
% if the surface area between (2,4) less than 0, group4
for ggg = 1:length(rCoeff_new)
    if oriGroup(ggg) == 1
        line11 = nanmean(rCoeff_new{ggg},1);
        areaUnder1 = trapz(rt(rt>0&rt<2),line11(rt>0 & rt<2));
       %areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
        if areaUnder1 < 0 %& areaUnder2 > areaUnder1
            group1 = [group1;rCoeff_new{ggg}];
%         elseif areaUnder2 < 0
%             group4 = [group4;rCoeff_new{ggg}];
        else
            group2 = [group2;rCoeff_new{ggg}];
        end
    elseif oriGroup(ggg) == 2
        line11 = nanmean(rCoeff_new{ggg},1);
        areaUnder1 = trapz(rt(rt>0&rt<2),line11(rt>0 & rt<2));
        %areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
        if areaUnder1 < 0 %& areaUnder2 > areaUnder1
            group1 = [group1;rCoeff_new{ggg}];
%         elseif areaUnder2 < 0
%             group4 = [group4;rCoeff_new{ggg}];
        else
            group2 = [group2;rCoeff_new{ggg}];
        end
    else
%         areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
%         if areaUnder2 < 0
%             group4 = [group4;rCoeff_new{ggg}];
%         else
            group3 = [group3;rCoeff_new{ggg}];
        %end
    end
end

% sort group2 and group3
g2.coeff =group2; g2.t = rt;
g2sortOrd = coeff_sort(g2,[0,3]);
g3.coeff = group3; g3.t = rt;
g3sortOrd = coeff_sort(g3,[0,3]);

temp2 = group2(g2sortOrd,:); group2 = temp2;
temp3 = group3(g3sortOrd,:); group3 = temp3;

% calculate correlation
nCells = size(group1,1)+size(group2,1) + size(group3,1);
 corrMat = zeros(nCells);
         %regroup the cells with cluster results
       
 total = [group1;group2;group3];
         
        for gg=1:nCells
            for yy = gg:nCells
                coeff = corrcoef(total(gg,:),total(yy,:));
                corrMat(gg,yy) = coeff(1,2);
            end
        end
        
        for yy=1:nCells
            for gg = yy+1:nCells
                
                corrMat(gg,yy) = corrMat(yy,gg);
            end
        end

figure;
subplot(1,3,1);

image(rt,1:nCells,[group1;group2;group3],'CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([rt(1) rt(end)],[size(group1,1) size(group1,1)],'Color',[241, 84, 18]/255)
hold on; plot([rt(1) rt(end)],[size(group1,1)+size(group2,1) size(group1,1)+size(group2,1)],'k')
%hold on; plot([rt(1) rt(end)],[size(group1,1)+size(group2,1)+size(group3,1) size(group1,1)+size(group2,1)+size(group3,1)],'k')
colors=cbrewer('div','RdBu',256);
colors=flipud(colors);
colorRange = [-1 1];
colormap(colors);
caxis([-0.05 0.05]);      %normalize dF/F heatmap to max of all conditions
ylabel('Cells');



% plot the cross correlation
clusterNum = [size(group1,1),size(group2,1),size(group3,1)];
 subplot(2,2,2)
     image(corrMat,'CDataMapping','scaled')
        %hold on;dendrogram(z)
         axis square;
            colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
    colormap(colors);
    caxis([-1 1]);
    hold on;
    for cc = 1:length(clusterNum)
        startInd = 1+sum(clusterNum(1:cc-1));
        endInd = sum(clusterNum(1:cc));
        if cc==1
            color = [241, 84, 18]/255;
        elseif cc == 2
            color = [0 0 0];
        elseif cc == 3
            color = [52, 179, 241]/255;
        end
        plot([startInd endInd endInd startInd startInd],[startInd startInd endInd endInd startInd],':','Color',color);
    end
    xticklabels({})
yticklabels({})

subplot(2,2,4)
line1 = group1';
ste = nanstd(line1,0,2)/sqrt(size(line1,2));
plot(rt,nanmean(line1,2),'Color',[241, 84, 18]/255);
hold on;
errorshade(rt,nanmean(line1,2)-ste,nanmean(line1,2)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = group2';
ste = nanstd(line2,0,2)/sqrt(size(line2,2));
plot(rt,nanmean(line2,2),'k');
errorshade(rt,nanmean(line2,2)-ste,nanmean(line2,2)+ste,[0 0 0]/255,0.5);
hold on;
line3 = group3';
ste = nanstd(line3,0,2)/sqrt(size(line3,2));
plot(rt,nanmean(line3,2), 'Color',[52, 179, 241]/255);
errorshade(rt,nanmean(line3,2)-ste,nanmean(line3,2)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');

tlabel = 'Outcome';
print(gcf,'-dpng',fullfile(savesumfigpath,[tlabel,'-cluster-summary']));
saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'fig');
saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'svg');

%% using the average to regroup?
% newGroup = oriGroup;
%
% ref1 = nanmean(group1,1);
% ref2 = nanmean(group2,1);
% ref3 = nanmean(group3,1);
% outInd  =1;
% outlierGroups = cell(0);
%
% for nn = 1:length(rCoeff_new)
%     % ACh, first session as ref group
%
%
%     g1Num = sum(rClust{nn}==1); g2Num = sum(rClust{nn}==2); g3Num = sum(rClust{nn}==3);
%     % calculate the average coefficient with in 1-2.5s for each group
%     g1Mean =nanmean(rCoeff{nn}(rClust{nn}==1,:),1);  g2Mean = nanmean(rCoeff{nn}(rClust{nn}==2,:),1); g3Mean =nanmean(rCoeff{nn}(rClust{nn}==3,:),1);
%
%
%     figure;plot(rt,ref1);hold on;
%     plot(rt,ref2);plot(rt,ref3);
%     plot(rt,g1Mean)
%
%     figure;plot(rt,ref1);hold on;
%     plot(rt,ref2);plot(rt,ref3);
%     plot(rt,g2Mean)
%
%     figure;plot(rt,ref1);hold on;
%     plot(rt,ref2);plot(rt,ref3);
%     plot(rt,g3Mean)
%     corr11 = corrcoef(g1Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr12 = corrcoef(g1Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr13 = corrcoef(g1Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
%     corr21 = corrcoef(g2Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr22 = corrcoef(g2Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr23 = corrcoef(g2Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
%     corr31 = corrcoef(g3Mean(rt>t0 & rt<t1),ref1(rt>t0 & rt<t1)); corr32 = corrcoef(g3Mean(rt>t0 & rt<t1),ref2(rt>t0 & rt<t1)); corr33 = corrcoef(g3Mean(rt>t0 & rt<t1),ref3(rt>t0 & rt<t1));
%
%     corr1 = [corr11(1,2),corr12(1,2),corr13(1,2)];
%     corr2 = [corr21(1,2),corr22(1,2),corr23(1,2)];
%     corr3 = [corr31(1,2),corr32(1,2),corr33(1,2)];
%
%     [maxCorr,ind] = max(corr1);
%     if maxCorr>0.5 & g1Num > 25
%         newGroup(nn) = ind;
%     else
%         newGroup(nn) =4; % 4 means outlier
%         outlierGroups{outInd} = rCoeff{nn}(rClust{nn}==1,:);
%         outInd = outInd+1;
%     end
%     [maxCorr,ind] = max(corr2);
%     if maxCorr>0.5 & g2Num > 25
%         newGroup(nn) = ind
%     else
%         newGroup(nn) =4; % 4 means outlier
%         outlierGroups{outInd} = rCoeff{nn}(rClust{nn}==2,:);
%         outInd = outInd+1;
%     end
%     [maxCorr,ind] = max(corr3);
%     if maxCorr>0.5 & g3Num > 25
%         newGroup(nn) = ind;
%
%     else
%         newGroup(nn) =4; % 4 means outlier
%         outlierGroups{outInd} = rCoeff{nn}(rClust{nn}==3,:);
%         outInd = outInd+1;
%     end
% end



% clust.t = rt;
% clust.coeff = clustCoeff;
% tRange = [-3,5]; numClust=6;
% regCoef_cluster(clust,tRange,numClust,tlabel1,'Time from cue(s)',savefluofigpath);


%% get lag
% cnLagCorr = zeros(1,nFiles); cnLagP =  zeros(1,nFiles); cnDis = []; cnCorr = []; cnLag = [];
% rnLagCorr = zeros(1,nFiles); rnLagP = zeros(1,nFiles); rnDis = []; rnCorr = []; rnLag = [];
% xnLagCorr = zeros(1,nFiles); xnLagP = zeros(1,nFiles); xnDis = []; xnCorr = []; xnLag = [];
%
%
% for ii = 1:nFiles
%
%     % load behavior files
%     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
%     %load(fullfile(fn_beh.folder,fn_beh.name));
%
%
%     % load fluorescent files
%     fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
%     if length(fn_fluo) == 1
%
%         savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
%         if ~exist(savefluofigpath,'dir')
%             mkdir(savefluofigpath);
%         end
%
%         cd(savefluofigpath);
%
%
%         %load(fullfile(fn_fluo.folder,fn_fluo.name));
%
%         % make folders to save analysis and plots
%         savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
%         if ~exist(savematpath,'dir')
%             mkdir(savematpath);
%         end
%
%         saveregpath = fullfile(savematpath,'regressionVarTempCorr.mat');  % file path to save the results
%         if exist(saveregpath)
%             load(saveregpath)
%         end
%
%         % get the data
%         % scatter plot of distance from reference grid and lag
%
%         % get distance
%         cnDistance = []; cnlag = []; cnCoeff = [];
%         for tt = 1:size(choicetempData.tempCorrLag,1)
%            for uu = 1:size(choicetempData.tempCorrLag,2)
%                if ~isnan(choicetempData.tempCorrLag(tt,uu))
%                    cnDistance = [cnDistance,sqrt((tt-choicetempData.maxX)^2+(uu-choicetempData.maxY)^2)];
%                    cnlag = [cnlag, choicetempData.tempCorrLag(tt,uu)];
%                    cnCoeff = [cnCoeff, choicetempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%
%         [R,P] = corrcoef(cnDistance,abs(cnlag));
%         if ~isnan(R)
%             cnLagCorr(ii) = R(1,2); cnLagP(ii) = P(1,2); cnDis = [cnDis,cnDistance]; cnCorr = [cnCorr, cnCoeff]; cnLag = [cnLag, cnlag];
%         end
%             %figure;scatter(cnDistance,cnLag);
%
%         rnDistance = []; rnlag = []; rnCoeff = [];
%         for tt = 1:size(outcometempData.tempCorrLag,1)
%            for uu = 1:size(outcometempData.tempCorrLag,2)
%                if ~isnan(outcometempData.tempCorrLag(tt,uu))
%                    rnDistance = [rnDistance,sqrt((tt-outcometempData.maxX)^2+(uu-outcometempData.maxY)^2)];
%                    rnlag = [rnlag, outcometempData.tempCorrLag(tt,uu)];
%                    rnCoeff = [rnCoeff, outcometempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%            [R,P] = corrcoef(rnDistance,abs(rnlag));
%             if ~isnan(R)
%         rnLagCorr(ii) = R(1,2); rnLagP(ii) = P(1,2); rnDis = [rnDis,rnDistance]; rnCorr = [rnCorr, rnCoeff]; rnLag = [rnLag, rnlag];
%             end
%
%         xnDistance = []; xnlag = []; xnCoeff = [];
%         for tt = 1:size(xntempData.tempCorrLag,1)
%            for uu = 1:size(xntempData.tempCorrLag,2)
%                if ~isnan(xntempData.tempCorrLag(tt,uu))
%                    xnDistance = [xnDistance,sqrt((tt-xntempData.maxX)^2+(uu-xntempData.maxY)^2)];
%                    xnlag = [xnlag, xntempData.tempCorrLag(tt,uu)];
%                    xnCoeff = [xnCoeff, xntempData.tempCorrCoeff(tt,uu)];
%                end
%            end
%         end
%            [R,P] = corrcoef(xnDistance,abs(xnlag));
%
%             if ~isnan(R)
%         xnLagCorr(ii) = R(1,2); xnLagP(ii) = P(1,2); xnDis = [xnDis,xnDistance]; xnCorr = [xnCorr, xnCoeff]; xnLag = [xnLag, xnlag];
%             end
%     end
% end
%
% % save the results
% savematdir = fullfile(savematsumpath, 'tempCorr');
% if ~exist(savematdir)
%     mkdir(savematdir)
% end
% cnTemp.dis = cnDis; cnTemp.corr = cnCorr; cnTemp.lag = cnLag; cnTemp.corrp = cnLagP; cnTemp.corrCoef = cnLagCorr;
% rnTemp.dis = rnDis; rnTemp.corr = rnCorr; rnTemp.lag = rnLag; rnTemp.corrp = rnLagP; rnTemp.corrCoef = rnLagCorr;
% xnTemp.dis = xnDis; xnTemp.corr = xnCorr; xnTemp.lag = xnLag; xnTemp.corrp = xnLagP; xnTemp.corrCoef = xnLagCorr;
% save(fullfile(savematdir,'tempCorrLag.mat'),'cnTemp','rnTemp','xnTemp');
%
% cnLagCorrSig = cnLagCorr(cnLagP<0.05);
% rnLagCorrSig = rnLagCorr(rnLagP<0.05);
% xnLagCorrSig = xnLagCorr(xnLagP<0.05);
%
%  [p,h,stats] = signrank(cnLagCorrSig,0,'tail','right')
%   [p,h,stats] = signrank(rnLagCorrSig,0,'tail','right')
%    [p,h,stats] = signrank(xnLagCorrSig,0,'tail','right')
%
% figure;scatter(cnDis, abs(cnLag),'.');
% [R,P] = corrcoef(cnDis,abs(cnLag))
%
% figure;scatter(rnDis, abs(rnLag),'.');
% [R,P] = corrcoef(rnDis,abs(rnLag))
%
% figure;scatter(xnDis, abs(xnLag),'.');
% [R,P] = corrcoef(xnDis,abs(xnLag))