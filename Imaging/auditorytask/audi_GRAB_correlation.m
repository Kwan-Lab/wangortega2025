function audi_GRAB_correlation(dataIndex)

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
    cd(savefigpath)
   
     autoCorr = zeros(length(cells.dFF));
     autoCorrP = zeros(length(cells.dFF));
        for uu = 1:length(cells.dFF)
            for vv = 1:length(cells.dFF)
                [autoCorr(uu,vv),autoCorrP(uu,vv)] = corr(cells.dFF{uu}, cells.dFF{vv});
            end
        end
    % plot correlation historgram
    corrPair = triu(autoCorr,1);
    p = signtest(corrPair(corrPair~=0),0);
    meanCorr = mean(corrPair(corrPair~=0));
    figure;histogram(corrPair(corrPair~=0),'FaceColor','black');
    set(gca,'box','off');
    xlabel('Correlation coefficient');
    ylabel('# of pairs');
    print(gcf,'-dpng','Correlation histogram');    %png format
    saveas(gcf, 'Correlation histogram', 'fig');
    
%     binedges = [-0.4:0.01:0.9];
%     figure;histogram(corr3,binedges);
%     hold on;histogram(corr9,binedges)
%     set(gca,'box','off');
%     xlabel('Correlation coefficient');
%     ylabel('# of pairs');
%     legend({'spontaneous','auditory'});
%     print(gcf,'-dpng','Correlation histogram-902');    %png format
%     saveas(gcf, 'Correlation histogram-902', 'fig');
    % cross-correlogram
    frequency = 1/mean(diff(cells.t));
    timeWin = 200;
%     for gg = 1:20
%         [c,lags] = xcorr(cells.dFF{98},cells.dFF{gg},'normalized');
%         lagTime = lags/frequency;
%         figure;
%         stem(lagTime(lagTime>-timeWin & lagTime<timeWin),c(lagTime>-timeWin & lagTime<timeWin));
%         hold on;
%         plot([0 0],[0 1],'-k','LineWidth',0.5);
%         set(gca,'box','off');
%         xlabel('Time lag (s)');
%         ylabel('Normalized corss-correlation');
%         print(gcf,'-dpng',['Cross-correlation between grid98 and grid',num2str(gg)]);    %png format
%         saveas(gcf, ['Cross-correlation between grid98 and grid',num2str(gg)], 'fig');
%         close
% 
%     end
    
    % plot correlation diagram based on pvalue
    sigCorr = zeros(length(cells.dFF));
    sigCorr(autoCorr>0 & autoCorrP<0.05) = 1;
    sigCorr(autoCorr<0 & autoCorrP<0.05) = -1;
    figure;imagesc(sigCorr);
    colorbar;
    print(gcf,'-dpng',['Cross-correlation sig']);    %png format
        saveas(gcf, ['Cross-correlation sig'], 'fig');
    % average correlation of one grid with other grids
    colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);

    edgelength = sqrt(numel(cells.dFF));
    meanCorrSG = zeros(edgelength);
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

        meanCorrSG(Ind1,Ind2) = mean(autoCorr(gg,:));
    end
    
%     figure;imagesc(meanCorrSG);
%     colorRange = [-0.5 0.5];
%     colormap(colors);
%     caxis([colorRange(1) colorRange(2)]);
%     axis square
%     colorbar;
%    print(gcf,'-dpng',['Cross-correlation-average correlation coefficient']);    %png format
%         saveas(gcf, ['Cross-correlation-average correlation coefficient'], 'fig');
%     
                multidFF = zeros(length(cells.dFF{1}),length(cells.dFF));
        for gg = 1:length(cells.dFF)
            multidFF(:,gg) = cells.dFF{gg};
        end
        meandFF = mean(multidFF,2);
        % correlation with mean dFF
        corrMean = zeros(edgelength);  
        for gg = 1:length(cells.dFF)
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
            corrMean(Ind1,Ind2) = corr(cells.dFF{gg}, meandFF);
        end
         figure;imagesc(corrMean);
    colormap(colors);
    axis square
    colorRange = [-1 1];
    colormap(colors);
    caxis([colorRange(1) colorRange(2)]);
    colorbar;
    print(gcf,'-dpng',['Cross-correlation with average dFF']);    %png format
        saveas(gcf, ['Cross-correlation with average dFF'], 'fig');
    
%      frequency = 1/mean(diff(cells.t));
%     timeWin = 200;
%     for gg = 1:numel(cells.dFF)
%         [c,lags] = xcorr(cells.dFF{gg},meandFF,'normalized');
%         lagTime = lags/frequency;
%         figure;
%         stem(lagTime(lagTime>-timeWin & lagTime<timeWin),c(lagTime>-timeWin & lagTime<timeWin));
%         hold on;
%         plot([0 0],[0 1],'-k','LineWidth',0.5);
%         set(gca,'box','off');
%         xlabel('Time lag (s)');
%         ylabel('Normalized corss-correlation');
%         print(gcf,'-dpng',['Cross-correlation between average dFF and grid',num2str(gg)]);    %png format
%         saveas(gcf, ['Cross-correlation between average dFF and grid',num2str(gg)], 'fig');
%         close
% 
%     end
        %% power spectrum
        pSpec = [];
        for gg = 1:numel(cells.dFF)
            [p,f] = pspectrum(cells.dFF{gg},frequency);
            pSpec = [pSpec;p'];
        end
        figure;
        plot(f,pow2db(mean(pSpec,1)),'-k')
        xlabel('Frequency (Hz)');
        ylabel('Power Spectrum (dB)');
        %xlim([0 2]);
        set(gca,'box','off');
        print(gcf,'-dpng',['Single unit average power spectrum']);    %png format
        saveas(gcf, ['Single unit average power spectrum'], 'fig');
    
        % multiunit spectrum
%         multidFF = zeros(length(cells.dFF{1}),length(cells.dFF));
%         for gg = 1:length(cells.dFF)
%             multidFF(:,gg) = cells.dFF{gg};
%         end
%         figure
%         pspectrum(sum(multidFF,2),frequency);
    % correlation to meandFF
    close all;
end