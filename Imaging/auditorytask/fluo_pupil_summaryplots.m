function fluo_pupil_summaryplots(dataIndex)
% plot spontaneous pupil & fluorescent dynamics

nFiles = size(dataIndex,1);
cxy_NE.coeff = [];
cxy_Ach.coeff = [];

 colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
for ii = 1:nFiles
    fn_coh = dir(fullfile(savematpath,savematname),'cxy','fc'););
    load(fullfile(fn_coh.folder,fn_coh.name));
    
    %if ii = 1:3
    % load behavior files
     
   
   
    
  
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
       
        
        cxyAll.coeff = [cxyAll.coeff; cxy.coeff];
    
        
end

% plot total of 588 grids
cxyAll = getBootstrp(cxyAll, 0, 0.05);
        
        figure;
        plot(fc(2:end),cxyAll.coeff_bootave(2:end),'k');
        hold on;
        gray = [0.7, 0.7, 0.7];
        %patch([fc fliplr(fc)], [cxy.boothigh  fliplr(cxy.bootlow)], [0.7 0.7 0.7])
        hold on;
        errorshade(fc(2:end),cxyAll.bootlow(2:end),cxyAll.boothigh(2:end),gray);
        hold on;plot(fc(2:end),cxyAll.coeff_bootave(2:end),'k');
        
        xlim([0.01, 0.2]);
        set(gca,'box','off')
        set(gca, 'XScale', 'log')
                     print(gcf,'-dpng',['spon_pupil_All']);
             saveas(gcf, 'spon_pupil_All', 'fig');
             saveas(gcf, 'spon_pupil_All', 'svg');
end