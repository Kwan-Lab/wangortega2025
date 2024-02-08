function plot_groupSummary(group1,group2, group3, rt, tlabel,savesumfigpath)


%% plot groups separately

%% group1
h1=figure;
sgtitle([tlabel,'-group1'])
subplot(1,3,1);

nCells = size(group1,1);

image(rt,1:nCells,group1,'CDataMapping','scaled');

hold on; plot([0 0],[0 nCells+1],'w');
colors=cbrewer('div','RdBu',256);
colorRange = [-1 1];
colormap(colors);
colors=flipud(colors);
colorRange(1)=-0.1;
colorRange(2)=0.1;
caxis([colorRange(1) colorRange(2)]);
ylabel('Cells');

%subplot
subplot(3,20,48);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);

% calculating correlations
corrMat = zeros(size(group1,1));

for xx = 1:size(group1,1)
    for yy = xx:size(group1,1)
        tempCorr = corrcoef(group1(xx,:),group1(yy,:));
        corrMat(xx,yy) = tempCorr(1,2);
    end
end
for xx = 1:size(group1,1)
    for yy = 1:xx-1
        corrMat(xx,yy) = corrMat(yy,xx);
    end
end
subplot(2,2,2)
image(corrMat,'CDataMapping','scaled')
colormap(colors);
xlim([1, size(group1,1)])
caxis([-1 1]);
set(gca,'dataAspectRatio',[1 1 1])

subplot(2,2,4)
line1 = group1';
ste = nanstd(line1,0,2)/sqrt(size(line1,2));
plot(rt,nanmean(line1,2),'Color',[241, 84, 18]/255);
hold on;
errorshade(rt,nanmean(line1,2)-ste,nanmean(line1,2)+ste,[241, 84, 18]/255,0.5);
% hold on;
% line2 = group2';
% ste = nanstd(line2,0,2)/sqrt(size(line2,2));
% plot(rt,nanmean(line2,2),'k');
% errorshade(rt,nanmean(line2,2)-ste,nanmean(line2,2)+ste,[0 0 0]/255,0.5);
% hold on;
% line3 = group3';
% ste = nanstd(line3,0,2)/sqrt(size(line3,2));
% plot(rt,nanmean(line3,2), 'Color',[52, 179, 241]/255);
% errorshade(rt,nanmean(line3,2)-ste,nanmean(line3,2)+ste,[52, 179, 241]/255,0.5);
% if ~isempty(group4)
%     line4 = group4';
%     ste = nanstd(line4,0,2)/sqrt(size(line4,2));
%     plot(rt,nanmean(line4,2), 'Color',[238, 238, 238]/255);
%     errorshade(rt,nanmean(line4,2)-ste,nanmean(line4,2)+ste,[100,100,100]/255,0.5);
% 
% end
ylim([-0.1 0.1])
yticks([-0.1:0.1:0.1])
set(gca,'box','off');
print(gcf,'-dpng',fullfile(savesumfigpath,[tlabel,'-cluster-summary-group1']));
savefig(h1,fullfile(savesumfigpath,[tlabel,'-cluster-summary-group1.fig']));
%saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'fig','-v7.3');
saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary-group1']), 'svg');

%% group2
h1=figure;
sgtitle([tlabel,'-group2'])
subplot(1,3,1);

nCells = size(group2,1);

image(rt,1:nCells,group2,'CDataMapping','scaled');

hold on; plot([0 0],[0 nCells+1],'w');
colors=cbrewer('div','RdBu',256);
colorRange = [-1 1];
colormap(colors);
colors=flipud(colors);
colorRange(1)=-0.1;
colorRange(2)=0.1;
caxis([colorRange(1) colorRange(2)]);
ylabel('Cells');

%subplot
subplot(3,20,48);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);

corrMat = zeros(size(group2,1));

for xx = 1:size(group2,1)
    for yy = xx:size(group2,1)
        tempCorr = corrcoef(group2(xx,:),group2(yy,:));
        corrMat(xx,yy) = tempCorr(1,2);
    end
end
for xx = 1:size(group2,1)
    for yy = 1:xx-1
        corrMat(xx,yy) = corrMat(yy,xx);
    end
end
subplot(2,2,2)
image(corrMat,'CDataMapping','scaled')
colormap(colors);
xlim([1, size(group2,1)])
caxis([-1 1]);
set(gca,'dataAspectRatio',[1 1 1])

subplot(2,2,4)
%% NE,interaction
% sepPoint = 600; %NE,interaction
% line1 = group2(1:sepPoint,:)';
% ste = nanstd(line1,0,2)/sqrt(size(line1,2));
% plot(rt,nanmean(line1,2));
% hold on;
% errorshade(rt,nanmean(line1,2)-ste,nanmean(line1,2)+ste,0.5);
% plot(rt,nanmean(line1,2));
% hold on;
% line2 = group2(sepPoint+1:end,:)';
% ste = nanstd(line2,0,2)/sqrt(size(line2,2));
% plot(rt,nanmean(line2,2),'k');
% errorshade(rt,nanmean(line2,2)-ste,nanmean(line2,2)+ste,0.5);
% plot(rt,nanmean(line2,2));



%% normal
line2 = group2';
ste = nanstd(line2,0,2)/sqrt(size(line2,2));
plot(rt,nanmean(line2,2),'k');
hold on; errorshade(rt,nanmean(line2,2)-ste,nanmean(line2,2)+ste,[0 0 0]/255,0.5);
ylim([-0.1 0.1])
yticks([-0.1:0.1:0.1])
set(gca,'box','off');
print(gcf,'-dpng',fullfile(savesumfigpath,[tlabel,'-cluster-summary-group2']));
savefig(h1,fullfile(savesumfigpath,[tlabel,'-cluster-summary-group2.fig']));
%saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'fig','-v7.3');
saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary-group2']), 'svg');

%% group3
% h1=figure;
% sgtitle([tlabel,'-group3'])
% subplot(1,3,1);
% 
% nCells = size(group3,1);
% 
% image(rt,1:nCells,group3,'CDataMapping','scaled');
% 
% hold on; plot([0 0],[0 nCells+1],'w');
% colors=cbrewer('div','RdBu',256);
% colorRange = [-1 1];
% colormap(colors);
% colors=flipud(colors);
% colorRange(1)=-0.1;
% colorRange(2)=0.1;
% caxis([colorRange(1) colorRange(2)]);
% ylabel('Cells');
% 
% %subplot
% subplot(3,20,48);
% image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);
% 
% subplot(2,2,4)
% line3 = group3';
% ste = nanstd(line3,0,2)/sqrt(size(line3,2));
% plot(rt,nanmean(line3,2), 'Color',[52, 179, 241]/255);
% hold on; errorshade(rt,nanmean(line3,2)-ste,nanmean(line3,2)+ste,[52, 179, 241]/255,0.5);
% ylim([-0.1 0.1])
% yticks([-0.1:0.1:0.1])
% set(gca,'box','off');
% print(gcf,'-dpng',fullfile(savesumfigpath,[tlabel,'-cluster-summary-group3']));
% savefig(h1,fullfile(savesumfigpath,[tlabel,'-cluster-summary-group3.fig']));
% %saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'fig','-v7.3');
% saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary-group3']), 'svg');
% 
%% group4
   % sort group4
  
if ~isempty(group3)

      g3.coeff=group3;g3.t = rt;
%     g3sortOrd = coeff_sort(g3,[0,3]);
%     group3= group3(g3sortOrd,:);
h1=figure;
sgtitle([tlabel,'-group3'])
subplot(1,3,1);

nCells = size(group3,1);

image(rt,1:nCells,group3,'CDataMapping','scaled');

hold on; plot([0 0],[0 nCells+1],'w');
colors=cbrewer('div','RdBu',256);
colorRange = [-1 1];
colormap(colors);
colors=flipud(colors);
colorRange(1)=-0.1;
colorRange(2)=0.1;
caxis([colorRange(1) colorRange(2)]);
ylabel('Cells');

%subplot
subplot(3,20,48);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);

corrMat = zeros(size(group3,1));

for xx = 1:size(group3,1)
    for yy = xx:size(group3,1)
        tempCorr = corrcoef(group3(xx,:),group3(yy,:));
        corrMat(xx,yy) = tempCorr(1,2);
    end
end
for xx = 1:size(group3,1)
    for yy = 1:xx-1
        corrMat(xx,yy) = corrMat(yy,xx);
    end
end
subplot(2,2,2)
image(corrMat,'CDataMapping','scaled')
colormap(colors);
xlim([1, size(group1,1)])
caxis([-1 1]);
set(gca,'dataAspectRatio',[1 1 1])

subplot(2,2,4)

 % sepPoint = 657; % ACh, choice
% sepPoint = 131; %ACh, xn
%sepPoint = size(group4,1); %NE, choice
sepPoint = 910; %NE,interaction
line1 = group3(1:sepPoint,:)';
ste = nanstd(line1,0,2)/sqrt(size(line1,2));
plot(rt,nanmean(line1,2));
hold on;
errorshade(rt,nanmean(line1,2)-ste,nanmean(line1,2)+ste,0.5);
plot(rt,nanmean(line1,2));
hold on;
line2 = group3(sepPoint+1:end,:)';
ste = nanstd(line2,0,2)/sqrt(size(line2,2));
plot(rt,nanmean(line2,2),'k');
errorshade(rt,nanmean(line2,2)-ste,nanmean(line2,2)+ste,0.5);
plot(rt,nanmean(line2,2));
set(gca,'box','off')
ylim([-0.1 0.1])
yticks([-0.1:0.1:0.1])
set(gca,'box','off');
print(gcf,'-dpng',fullfile(savesumfigpath,[tlabel,'-cluster-summary-group3']));
savefig(h1,fullfile(savesumfigpath,[tlabel,'-cluster-summary-group3.fig']));
%saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'fig','-v7.3');
saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary-group3']), 'svg');

end

% calculate correlation
if isempty(group3)
    nCells = size(group1,1)+size(group2,1);
    total = [group1;group2];
else
%     if ismember(tlabel,{'Interaction-Sig','choice-Sig','posRPE-Sig','negRPE-Sig','dQ-sig','dK-sig','CKE-sig'})
%         %sort g1,g2,g3 as well
% %          g1.coeff=group1;g1.t = rt;
% %         g1sortOrd = coeff_sort(g1,[0,3]);
% %         group1= group1(g1sortOrd,:);
% %          g2.coeff=group2;g2.t = rt;
% %     g2sortOrd = coeff_sort(g2,[0,3]);
% %     group2= group2(g2sortOrd,:);
% %      g3.coeff=group3;g3.t = rt;
% %     g3sortOrd = coeff_sort(g3,[0,3]);
% %     group3= group3(g3sortOrd,:);
%     end
 
    nCells = size(group1,1)+size(group2,1) + size(group3,1);
    total = [group1;group2;group3];
end
corrMat = zeros(nCells);
%regroup the cells with cluster results



for gg=1:nCells
    parfor yy = gg:nCells
        coeff = corrcoef(total(gg,:),total(yy,:));
        corrMat(gg,yy) = coeff(1,2);
    end
end

for yy=1:nCells
    for gg = yy+1:nCells

        corrMat(gg,yy) = corrMat(yy,gg);
    end
end

h1=figure;
sgtitle(tlabel)
subplot(1,3,1);

if isempty(group3)
    image(rt,1:nCells,[group1;group2],'CDataMapping','scaled');
else
    image(rt,1:nCells,[group1;group2;group3],'CDataMapping','scaled');
end
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([rt(1) rt(end)],[size(group1,1) size(group1,1)],'Color',[241, 84, 18]/255)
hold on; plot([rt(1) rt(end)],[size(group1,1)+size(group2,1) size(group1,1)+size(group2,1)],'k')
%hold on; plot([rt(1) rt(end)],[size(group1,1)+size(group2,1)+size(group3,1) size(group1,1)+size(group2,1)+size(group3,1)],'Color',[52, 179, 241]/255)
colors=cbrewer('div','RdBu',256);
colorRange = [-1 1];
colormap(colors);
% if contains(tlabel,'RPE')
%     colorRange(1)=-0.7;
%     colorRange(2)=0.7;
% else
colorRange(1)=-0.1;
colorRange(2)=0.1;
%normalize dF/F heatmap to max of all conditions
% end
caxis([colorRange(1) colorRange(2)]);
ylabel('Cells');

%subplot
subplot(3,20,48);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);

% plot the cross correlation
if isempty(group3)
    clusterNum = [size(group1,1),size(group2,1)];
else
    clusterNum = [size(group1,1),size(group2,1),size(group3,1)];
end
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
    elseif cc == 4
        color = [100,100,100]/255;
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
% hold on;
% line3 = group3';
% ste = nanstd(line3,0,2)/sqrt(size(line3,2));
% plot(rt,nanmean(line3,2), 'Color',[52, 179, 241]/255);
% errorshade(rt,nanmean(line3,2)-ste,nanmean(line3,2)+ste,[52, 179, 241]/255,0.5);
if ~isempty(group3)
    line3 = group3';
    ste = nanstd(line3,0,2)/sqrt(size(line3,2));
    plot(rt,nanmean(line3,2), 'Color',[238, 238, 238]/255);
    errorshade(rt,nanmean(line3,2)-ste,nanmean(line3,2)+ste,[100,100,100]/255,0.5);

end
ylim([-0.1 0.1])
yticks([-0.1:0.1:0.1])
set(gca,'box','off');
print(gcf,'-dpng',fullfile(savesumfigpath,[tlabel,'-cluster-summary']));
savefig(h1,fullfile(savesumfigpath,[tlabel,'-cluster-summary.fig']));
%saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'fig','-v7.3');
saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'svg');

%% go through group 1 and 2 for rise time
% rise time: time from 10% to 90% maximum
% max time: time to get to maximum relative to the cue

% Two-Sample Kolmogorov-Smirnov Test to determine if group 2 and 3 are
% different
aveGroup1= mean(line1,2);aveGroup2 = mean(line2,2);
[h,p] = kstest2(aveGroup1,aveGroup2)



% figure;
% subplot(1,3,1)
% image(rt,1:size(group2,1),group2,'CDataMapping','scaled');
% 
% hold on; plot([0 0],[0 size(group2,1)+1],'w');
% hold on; plot([rt(1) rt(end)],[334 334],'b')
% hold on; plot([rt(1) rt(end)],[884 884],'k')
% colors=cbrewer('div','RdBu',256);
% colorRange = [-1 1];
% colormap(colors);
% colors=flipud(colors);
% % if contains(tlabel,'RPE')
% %     colorRange(1)=-0.7;
% %     colorRange(2)=0.7;
% % else
% colorRange(1)=-0.1;
% colorRange(2)=0.1;
% %normalize dF/F heatmap to max of all conditions
% % end
% caxis([colorRange(1) colorRange(2)]);
% ylabel('Cells');
% title('Interaction coefficient group2')
% subplot(3,20,48);
% image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);
% 
% subplot(2,2,4)
% plot(rt, nanmean(group2(1:334,:)),'b')
% hold on;  ste = nanstd(group2(1:334,:),0,1)/sqrt(size(group2(1:334,:),1));
%     errorshade(rt,nanmean(group2(1:334,:))-ste,nanmean(group2(1:334,:))+ste,'b',0.5);
% 
% hold on;plot(rt, nanmean(group2(335:884,:)),'k')
% hold on;  ste = nanstd(group2(335:884,:),0,1)/sqrt(size(group2(335:884,:),1));
%     errorshade(rt,nanmean(group2(335:884,:))-ste,nanmean(group2(335:884,:))+ste,'k',0.5);
% 
% plot(rt,nanmean(group2(885:end,:)),'r')
% hold on;  ste = nanstd(group2(885:end,:),0,1)/sqrt(size(group2(885:end,:),1));
%     errorshade(rt,nanmean(group2(885:end,:))-ste,nanmean(group2(885:end,:))+ste,'r',0.5);
% ylim([-0.1 0.1])
% yticks([-0.1:0.1:0.1])
% set(gca,'box','off')
% print(gcf,'-dpng',fullfile(savesumfigpath,[tlabel,'-group2']));
% savefig(h1,fullfile(savesumfigpath,[tlabel,'-group2.fig']));
% %saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'fig','-v7.3');
% saveas(gcf, fullfile(savesumfigpath,[tlabel,'-group2']), 'svg');

% figure;subplot(1,3,1)
% image(rt,1:size(group3,1),group3,'CDataMapping','scaled');
% 
% hold on; plot([0 0],[0 size(group3,1)+1],'w');
% hold on; plot([rt(1) rt(end)],[1344 1344],'b')
% colors=cbrewer('div','RdBu',256);
% colorRange = [-1 1];
% colormap(colors);
% colors=flipud(colors);
% % if contains(tlabel,'RPE')
% %     colorRange(1)=-0.7;
% %     colorRange(2)=0.7;
% % else
% colorRange(1)=-0.1;
% colorRange(2)=0.1;
% %normalize dF/F heatmap to max of all conditions
% % end
% caxis([colorRange(1) colorRange(2)]);
% ylabel('Cells');
% title('Interaction coefficient group3')
% subplot(3,20,48);
% image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);

% subplot(2,2,4)
% plot(rt, nanmean(group3(1:1344,:)),'b')
% hold on;  ste = nanstd(group3(1:1344,:),0,1)/sqrt(size(group3(1:1344,:),1));
%     errorshade(rt,nanmean(group3(1:1344,:))-ste,nanmean(group3(1:1344,:))+ste,'b',0.5);
% 
% plot(rt,nanmean(group3(1344:end,:)),'r')
% hold on;  ste = nanstd(group3(1344:end,:),0,1)/sqrt(size(group3(1344:end,:),1));
% errorshade(rt,nanmean(group3(1344:end,:))-ste,nanmean(group3(1344:end,:))+ste,'r',0.5);
% ylim([-0.1 0.1])
% yticks([-0.1:0.1:0.1])
% set(gca,'box','off')
% print(gcf,'-dpng',fullfile(savesumfigpath,[tlabel,'-group3']));
% savefig(h1,fullfile(savesumfigpath,[tlabel,'-group3.fig']));
% %saveas(gcf, fullfile(savesumfigpath,[tlabel,'-cluster-summary']), 'fig','-v7.3');
% saveas(gcf, fullfile(savesumfigpath,[tlabel,'-group3']), 'svg');


% 

close;