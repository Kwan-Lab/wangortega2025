function fluo_pupil_plots(dataIndex)
% plot spontaneous pupil & fluorescent dynamics

nFiles = size(dataIndex,1);
cxyAll.coeff = [];

 colors=cbrewer('div','RdBu',256);
colors=flipud(colors);
for ii = 1:nFiles
    
    % load behavior files
     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii},['*',date(1:6),'*_pup.mat']);
    dff_name = fullfile(fn_beh.folder,'dff.mat');
    
    fn_pup = dir(pup_name);
    fn_dff = dir(dff_name);
   
    
    if length(fn_pup) == 1 & length(fn_dff) == 1
        load(fullfile(fn_pup.folder,fn_pup.name));
        load(fullfile(fn_dff.folder,fn_dff.name));
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
        %% Plot pupil for all trials
        
        % MP_plot_pupil( pupil, trialData );
        
        %% Plot cue-aligned pupil
        
        % plot a single trace, create a gif
        %% ii=1; ACh plot; ii = 11, NE plot
       
        % tCue = trialData.cueTimes(1):trialData.cueTimes(2);
        pupili = interp1(pupil.t,pupil.dia,cells.t);
        
        % find the ROI with maximum correlation
        corrMat = zeros(1, length(cells.dFF));
        for rr = 1:length(cells.dFF)
            corrVal = corrcoef(cells.dFF{rr},pupili,'rows','pairwise');
            corrMat(rr) = corrVal(1,2);
        end
        [~,indFF] = max(corrMat);
       
        
        % paper plot 
%         avedFF = zeros(length(cells.t),1);
%         for jj  =1:length(cells.dFF)
%             avedFF = avedFF+cells.dFF{jj}/length(cells.dFF);
%         end
        % plot the last 10 mins of recordings
        tstart = pupil.t(end) - 11*60;
        tend = pupil.t(end);
        pupilPlotInd = pupil.t<tend & pupil.t>=tstart;
        figure;
        plot(pupil.t(pupilPlotInd),smooth(pupil.dia(pupilPlotInd),120));
       
        ylim([-2.4,2.6]);
        ax1 = gca;   
        ax1.YAxis.Visible = 'off'; % remove y-axis
        yyaxis right
        fluoPlotInd = cells.t<tend & cells.t>=tstart;
        dFFCurve = cells.dFF{indFF};
        hold on;plot(cells.t(fluoPlotInd), smooth(dFFCurve(fluoPlotInd),120));
        ylim([-0.04 0.46]);
         ax1 = gca;   
            ax1.YAxis(2).Visible = 'off'; % remove y-axis
             ax1.XAxis.Visible = 'off'; % remove y-axis
             set(gca,'box','off');
             
             % plot scale bar
             hold on;
             plot([tend-30 tend], [0.36 0.36],'k-');
             text(tend-30,0.33,'30 s')
             hold on;
             plot([tend tend],[0.36 0.46],'k-');
             text(tend+2,0.43, 'df/f 0.1; z 1')
             
             print(gcf,'-dpng',['spon_pupil_ave_ACh']);
             saveas(gcf, 'spon_pupil_ave_ACh', 'fig');
             saveas(gcf, 'spon_pupil_ave_ACh', 'svg');
        % plot the whole session 
        edgelength = sqrt(numel(cells.dFF));
    corrCoeff = zeros(edgelength);
    corrP = zeros(edgelength);
    for gg = 1:numel(cells.dFF)
        if mod(gg,edgelength) == 0
            Ind2 = edgelength;
        else
            Ind2 = mod(gg,(edgelength));
        end
        if mod(gg,edgelength) == 0
            Ind1 = gg/edgelength;
        else
            Ind1 = floor(gg/edgelength)+1;
        end
        
%            figure('Renderer', 'painters', 'Position', [10 10 2000 600]);
%             plot(pupil.t,smooth(pupil.dia,60)); ylabel('Pupil z-score');
%            
%             %ax1.XAxis.Visible = 'off'; % remove x-axis
%             yyaxis right
%             hold on; plot(cells.t, smooth(cells.dFF{jj},100));
%             ylabel('df/f (smoothed)');
%            
%             set(gca,'box','off');
%             print(gcf,'-dpng',['spon_pupil_grid' int2str(jj)]);
            % calculate correlation
            [cor,p] = corrcoef(pupili,cells.dFF{gg},'Rows','pairwise');
            corrCoeff(Ind1,Ind2) = cor(1,2);
            corrP(Ind1,Ind2) = p(1,2);
            close all;
        end
         savematname='corrPupFluo.mat';
         save(fullfile(savematpath,savematname),'corrCoeff','corrP');
         
        figure;
        colorRange = [-0.5 0.5];
        imagesc(corrCoeff);
         colormap(colors);
         caxis([colorRange(1) colorRange(2)]);
        axis square
        colorbar;
        print(gcf,'-dpng',['pupil correlation']);
         saveas(gcf, ['pupil correlation'], 'fig');
         
        % coherence analysis
        pupInt = interp1(pupil.t,pupil.dia, cells.t);
        % get rid of the NaN
        pupI = pupInt(~isnan(pupInt));
        % sample rate
        fs = 1/mean(diff(cells.t));
        [cxy,fc] = mscohere(pupI, cells.dFF{1}(~isnan(pupInt)),[],[],[],fs);
        [Pxy,F] = cpsd(pupI, cells.dFF{1}(~isnan(pupInt)),[],[],[],fs);

        cxyList = zeros(length(cells.dFF),length(cxy));
        pxyList = zeros(length(cells.dFF),length(Pxy));
        for tt = 1:length(cells.dFF)
        %[cxy,fc] = mscohere(pupI, cells.dFF{1}(~isnan(pupInt)));
            [cxyList(tt,:),fc] = mscohere(pupI, cells.dFF{tt}(~isnan(pupInt)),[],[],[],fs);
            [pxyList(tt,:),F] = cpsd(pupI, cells.dFF{tt}(~isnan(pupInt)),[],[],[],fs);

        end
        cxyBoot.coeff = cxyList;
        cxy = getBootstrp(cxyBoot, 0, 0.05);
        aveCxy = mean(cxyList);
        figure;
        plot(fc(2:end),cxy.coeff_bootave(2:end),'k');
        hold on;
        gray = [0.7, 0.7, 0.7];
        %patch([fc fliplr(fc)], [cxy.boothigh  fliplr(cxy.bootlow)], [0.7 0.7 0.7])
        hold on;
        errorshade(fc(2:end),cxy.bootlow(2:end),cxy.boothigh(2:end),gray);
        hold on;plot(fc(2:end),cxy.coeff_bootave(2:end),'k');
        
        xlim([0.01, 0.2]);
        ylim([0 1]);
        set(gca,'box','off')
        set(gca, 'XScale', 'log')
                     print(gcf,'-dpng',['spon_pupil_coh']);
             saveas(gcf, 'spon_pupil_coh', 'fig');
             saveas(gcf, 'spon_pupil_coh', 'svg');
        
        % cross spectrum
         pxyBoot.coeff = angle(pxyList)/pi;
        pxy = getBootstrp(pxyBoot, 0, 0.05);
        avePxy = mean(pxyList);
        figure;
        plot(F,pxy.coeff_bootave,'k');
        hold on;
        gray = [0.7, 0.7, 0.7];
        %patch([fc fliplr(fc)], [cxy.boothigh  fliplr(cxy.bootlow)], [0.7 0.7 0.7])
        hold on;
        errorshade(F,pxy.bootlow,pxy.boothigh,gray);
        hold on;plot(F,pxy.coeff_bootave,'k');
        
        xlim([0.01, 0.2]);
        %ylim([0 1]);
        set(gca,'box','off')
        set(gca, 'XScale', 'log')
                     print(gcf,'-dpng',['spon_pupil_xspec']);
             saveas(gcf, 'spon_pupil_xspec', 'fig');
             saveas(gcf, 'spon_pupil_xspec', 'svg');
        


        savematname='coherence.mat';
        save(fullfile(savematpath,savematname),'cxy','fc','pxy','F');
        
        %.coeff = [cxyAll.coeff; cxy.coeff];
    end
        
end

% plot total of 588 grids
% cxyAll = getBootstrp(cxyAll, 0, 0.05);
%         
%         figure;
%         plot(fc(2:end),cxyAll.coeff_bootave(2:end),'k');
%         hold on;
%         gray = [0.7, 0.7, 0.7];
%         %patch([fc fliplr(fc)], [cxy.boothigh  fliplr(cxy.bootlow)], [0.7 0.7 0.7])
%         hold on;
%         errorshade(fc(2:end),cxyAll.bootlow(2:end),cxyAll.boothigh(2:end),gray);
%         hold on;plot(fc(2:end),cxyAll.coeff_bootave(2:end),'k');
%         
%         xlim([0.01, 0.2]);
%         set(gca,'box','off')
%         set(gca, 'XScale', 'log')
%                      print(gcf,'-dpng',['spon_pupil_All']);
%              saveas(gcf, 'spon_pupil_All', 'fig');
%              saveas(gcf, 'spon_pupil_All', 'svg');
end