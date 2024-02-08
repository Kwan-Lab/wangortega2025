function audi_GRAB_pupil_plots(dataIndex)

% load behavior and dff data
% make some simple plots

nFiles = size(dataIndex,1);

for ii = 1:nFiles

    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    % load dFF files
    load(fullfile(fn_beh.folder,'dff.mat'));
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii},['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    %load(fullfile(fn_pup.folder,fn_pup.name));

    savefigpath = fullfile(fn_beh.folder,'figs-fluo');
    if ~exist(savefigpath,'dir')
        mkdir(savefigpath);
    end
    cd(savefigpath);

        %% Plot dF/F of all the cells

        %% get cue-aligned dF/F for each cell
        % Fig. 3d in paper came from 140605 data set, cell 8 10 37 74
        savepsthname = fullfile(fn_beh.folder,'psthMat.mat');
%         savepupname = fullfile(fn_pup.folder,'psthPup.mat');
        params=[];
            params.trigTime = trialData.cueTimes;
            params.xtitle = 'Time from cue (s)';
            params.window = [-1:0.05:3];
            params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
            params.CI = 0.95;  %confidence interval
            params.minNumTrial = 50;
        if ~exist(savepsthname)

            for j=1:numel(cells.dFF)
                psth_output=[];

                psth_panel(j).sig{1} = get_psth( cells.dFF{j}, cells.t, params.trigTime(2:end), 'df/f', params );

                tlabel = ['Cell #',num2str(j)];
                %plot_psth(psth_panel,tlabel,params.xtitle);
                %print(gcf,'-dpng',fullfile(savefigpath,['cell' int2str(j)]));
                %close;
            end

            save(savepsthname,'psth_panel');
        else
            display('PSTH already computed');
            load(savepsthname);
        end
        % save psth results
        %if ~exist(savepupname)
%             params.window = [-1:0.05:3];
%             psthpupil_panel(1).sig{1} = get_psth(pupil.dia, pupil.t, params.trigTime(2:end-1), 'Pupil change', params );
%             save(savepupname,'psthpupil_panel');
%         else
%             display('PSTH_pup already computed');
%         end
        gray=[0.7 0.7 0.7];
        avePsth = zeros(size(psth_panel(1).sig{1}.signal,1),1);
        figure;
        for tt=1:length(cells.dFF)
            hold on;
            t=psth_panel(tt).sig{1}.t;

             plot(t,psth_panel(tt).sig{1}.signal,'Color',[0.8500 0.3250 0.0980 0.2],'LineWidth',0.5);
             avePsth = avePsth + psth_panel(tt).sig{1}.signal/length(cells.dFF);
             %errorshade(t,psth_panel(tt).sig{1}.bootlow,psth_panel(tt).sig{1}.boothigh,gray);
        end
        ax1 = gca;
        ax1.YAxis(1).Visible = 'off'; % remove y-axis
        ax1.XAxis.Visible = 'off';
        hold on; plot(t, avePsth,'-','Color',[0.8500 0.3250 0.0980]);
        ylim([0.05 0.25]);
          hold on;
             plot([2.5 2.7], [0.18 0.18],'k-');
             hold on;
             plot([2.7 2.7],[0.18 0.2],'k-');
        % get average psth
%         yyaxis right
%         hold on;
%
%        errorshade(psthpupil_panel(1).sig{1}.t,psthpupil_panel(1).sig{1}.bootlow,psthpupil_panel(1).sig{1}.boothigh,[0 0.4470 0.7410], 0.2);
%         % 0.95 CI is too large
%        hold on; plot(psthpupil_panel(1).sig{1}.t,psthpupil_panel(1).sig{1}.signal,'-','Color',[0 0.4470 0.7410]);
%        hold on;
%        ylim([-0.3 0.6]);
%        ax1 = gca;
%        ax1.YAxis(2).Visible = 'off'; % remove y-axis
%        ax1.XAxis.Visible = 'off'; % remove y-axis
%         set(gca,'box','off');
%
%          hold on;
%              plot([2.5 2.7], [0.18 0.18],'k-');
%              hold on;
%              plot([2.7 2.7],[0.18 0.2],'k-');
%              hold on;
%              plot([0 0],[-0.5 0.5],'k--');
                print(gcf,'-dpng',['tone-fluo_ave' ]);
             saveas(gcf, 'tone-fluo_ave', 'fig');
             saveas(gcf, 'tone-fluo_ave', 'svg');

             close;

             %% cluster
nCells = numel(psth_panel);
         pref = zeros(length(t),nCells);
         for gg = 1:nCells
             pref(:,gg) = psth_panel(gg).sig{1}.signal-nanmean(psth_panel(gg).sig{1}.signal(t<0));
         end

clusterMat = pref';


        maxclust = 3;
        T = clusterdata(clusterMat,'Linkage','ward','SaveMemory','on','Maxclust',maxclust,'distance','correlation');

        % calculate correlation matrix
        corrMat = zeros(nCells);
         %regroup the cells with cluster results
         cellInd = 1:nCells;
         clustInd = [];
         clusterNum = zeros(1,maxclust);
         for gg=1:maxclust
             clustInd = [clustInd,cellInd(T==gg)];
             clusterNum(gg) = sum(T==gg);
         end

        for gg=1:length(clustInd)
            for yy = gg:length(clustInd)
                coeff = corrcoef(pref(:,clustInd(gg)),pref(:,clustInd(yy)));
                corrMat(gg,yy) = coeff(1,2);
            end
        end

        for yy=1:length(clustInd)
            for gg = yy+1:length(clustInd)

                corrMat(gg,yy) = corrMat(yy,gg);
            end
        end



        z  =linkage(T,'ward');


           figure;
     subplot(2,3,3)
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

% plot average coefficient by cluster
% error bar
subplot(2,2,4)
line1 = pref(:,cellInd(T==1));
ste = nanstd(line1,0,2)/sqrt(size(line1,2));
plot(t,nanmean(pref(:,cellInd(T==1)),2),'Color',[241, 84, 18]/255);
hold on;
errorshade(t,nanmean(line1,2)-ste,nanmean(line1,2)+ste,[241, 84, 18]/255,0.5);
hold on;
line2 = pref(:,cellInd(T==2));
ste = nanstd(line2,0,2)/sqrt(size(line2,2));
plot(t,nanmean(pref(:,cellInd(T == 2)),2),'k');
errorshade(t,nanmean(line2,2)-ste,nanmean(line2,2)+ste,[0 0 0]/255,0.5);
hold on;
line3 = pref(:,cellInd(T==3));
ste = nanstd(line3,0,2)/sqrt(size(line3,2));
plot(t,nanmean(pref(:,cellInd(T == 3)),2), 'Color',[52, 179, 241]/255);
errorshade(t,nanmean(line3,2)-ste,nanmean(line3,2)+ste,[52, 179, 241]/255,0.5);
set(gca,'box','off');

% check the spatial position of different clusters
    % plot coefficient by cluster
%
% clustMat = NaN(size(input.sigMat,1),size(input.sigMat,2));
%     tempind = 1;
% for xx = 1:size(input.sigMat,1)
%     for zz = 1:size(input.sigMat,2)
%
%         if input.sigMat(xx,zz)==1
%             clustMat(xx,zz) = T(tempind);
%             tempind = tempind+1;
%         end
%     end
% end

% subplot(2,3,2)
% b=image(clustMat,'CDataMapping','scaled'); %,
% set(b,'AlphaData',~isnan(clustMat))
% set(gca, 'Color', [0.7, 0.7, 0.7])
% axis square;
% colors = [[241, 84, 18]/255;[0, 0, 0];[52, 179, 241]/255 ];
% colormap(colors);
% caxis([1, 3]);
% xticklabels({})
% yticklabels({})
%title(titleText);

%plot in pseudocolor


subplot(1,3,1);
image(t,1:nCells,pref(:,clustInd)','CDataMapping','scaled');
hold on; plot([0 0],[0 nCells+1],'w');
hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
   colors=cbrewer('div','RdBu',256);
        colors=flipud(colors);
        colorRange = [-1 1];
colormap(colors);
caxis([-0.05 0.05]);      %normalize dF/F heatmap to max of all conditions
ylabel('Cells');

% xlabel(xtitle);
% title(tlabel);


%              % plot the
%
% figure;
%              image(t,1:nCells,pref','CDataMapping','scaled');
% hold on; plot([0 0],[0 nCells+1],'w');
% %hold on; plot([t(1) t(end)],[clusterNum(1) clusterNum(1)],'Color',[241, 84, 18]/255)
% %hold on; plot([t(1) t(end)],[sum(clusterNum(1:2)) sum(clusterNum(1:2))],'k')
%    colors=cbrewer('div','RdBu',256);
%         colors=flipud(colors);
%         colorRange = [-1 1];
% colormap(colors);
% caxis([-0.01 0.01]);      %normalize dF/F heatmap to max of all conditions
% ylabel('Cells');
% xlabel(xtitle);
% title(tlabel);

          print(gcf,'-dpng',['correlation' ]);
             saveas(gcf, 'correlation', 'fig');
             saveas(gcf, 'correlation', 'svg');

             close
end
