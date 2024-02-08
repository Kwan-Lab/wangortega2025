function audi_GRAB_plots(dataIndex)

% load behavior and dff data
% make some simple plots

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    % load dFF files
    load(fullfile(fn_beh.folder,'dff.mat'));
    

    savefigpath = fullfile(fn_beh.folder,'figs-fluo');
    if ~exist(savefigpath,'dir')
        mkdir(savefigpath);
    end
    
    if ~exist(fullfile(savefigpath,'cell49.png'))
        %% Plot dF/F of all the cells
        MP_plot_allCells( cells, trialData ,savefigpath);
        
        close all;
        %% check correaltion between different cells
        autoCorr = zeros(length(cells.dFF));
        for uu = 1:length(cells.dFF)
            for vv = 1:length(cells.dFF)
                autoCorr(uu,vv) = corr(cells.dFF{uu}, cells.dFF{vv});
            end
        end
        figure;
        heatmap(autoCorr);
        print(gcf,'-dpng',fullfile(savefigpath,'corr-cells'));    %png format
        saveas(gcf, fullfile(savefigpath,'corr-cells'), 'fig');
        
        %% Plot cue-aligned dF/F for each cell
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-1:0.05:3];
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = 0.95;  %confidence interval
        params.minNumTrial = 50;
        for j=1:numel(cells.dFF)
            psth_output=[];
            
            psth_panel(1).sig{1} = get_psth( cells.dFF{j}, cells.t, params.trigTime(2:end), 'df/f', params );
           
            tlabel = ['Cell #',num2str(j)];
            plot_psth(psth_panel,tlabel,params.xtitle);
            print(gcf,'-dpng',fullfile(savefigpath,['cell' int2str(j)]));
            close;
        end
    end
  
    trialsConsidered = 200;
    aveBase = zeros(trialsConsidered, numel(cells.dFF));
    aveTone = zeros(trialsConsidered, numel(cells.dFF));
    for tt = 1:trialsConsidered
        baseInd = cells.t>=(trialData.cueTimes(tt)-0.2) & cells.t<trialData.cueTimes(tt);
        toneInd = cells.t>=(trialData.cueTimes(tt)+0.05) & cells.t<(trialData.cueTimes(tt)+0.15);
        for cc = 1:numel(cells.dFF)  
            aveBase(tt,cc) = nanmean(cells.dFF{cc}(baseInd));
            aveTone(tt,cc) = nanmean(cells.dFF{cc}(toneInd));
        end
    end
    
    pInc = zeros(1,numel(cells.dFF));
     for cc = 1:numel(cells.dFF)  
        [~,pInc(cc)] = ttest(aveBase(:,cc),aveTone(:,cc),'tail','left');
     end
    % check the average df/f from -0.1 - 0 and 0.1-0.2
     
    aveCellBase = nanmean(aveBase,1);
    aveCellTone = nanmean(aveTone,1);
   
end