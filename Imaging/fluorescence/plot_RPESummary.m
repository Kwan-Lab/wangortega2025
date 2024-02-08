function plot_RPESummary(group1,group2, group3, rt, tlabel,savesumfigpath)


%% plot groups separately

group = [group1;group2;group3];


if ~isempty(group)

      g.coeff=group;g.t = rt;
    gsortOrd = coeff_sort(g,[0,3]);
    group= group(gsortOrd,:);
h1=figure;
sgtitle([tlabel,'-group'])
subplot(1,3,1);

nCells = size(group,1);

image(rt,1:nCells,group,'CDataMapping','scaled');

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
    line4 = group3';
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



close;