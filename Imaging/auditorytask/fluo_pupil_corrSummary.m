function fluo_pupil_corrSummary(dataIndex,savesumpath)

% summarize cross-correlation results of pupil and GRAB signal

nFiles = size(dataIndex,1);
colors=cbrewer('div','RdBu',256);
colors=flipud(colors);

Ind = 1:size(dataIndex,1);
Ach = Ind(strcmp(dataIndex.GRAB,'ACh'));
NE = Ind(strcmp(dataIndex.GRAB,'NE'));


maxLagSum_Ach = [];
maxCorrSum_Ach = [];

maxLagSum_NE = [];
maxCorrSum_NE = [];


%% get Ach data


for ii = 1:length(Ach)
    savematpath = fullfile(dataIndex.BehPath{Ach(ii)},'analysis-pupil');
    savefilename = fullfile(savematpath,'xcorr-pupil-fluo.mat');
    %fn_corr = dir(savefilename)
    if exist(savefilename)
        load(savefilename);
        % plot the cross-correlation with lags
        
        if ~exist('Lag_Ach','var')
            Lag_Ach = zeros(length(Ach),length(c),length(newlag(newlag<=2 & newlag > -2)));
            Corr_Ach = zeros(length(Ach),length(c),length(newlag(newlag<=2 & newlag > -2)));
        end
        if ~exist('Lag_Achauto','var')
            Lag_Achauto= zeros(length(NE),length(c),length(newlagFluo(newlagFluo<=100 & newlagFluo>=0)));
            Corr_Achauto = zeros(length(NE),length(c),length(newlagFluo(newlagFluo<=100 & newlagFluo>=0)));
        end
        date = num2str(dataIndex.DateNumber(Ach(ii)));
        savefigpath = fullfile(dataIndex.BehPath{Ach(ii)},[date(1:6),'_crosscorrelation']);
        if ~exist(savefigpath)
            mkdir(savefigpath);
        end
        
        maxCorr = zeros(1, length(c));
        maxLag = zeros(1, length(c));
        for rr = 1:length(c)
            % get the maximum correlation and lags
            corrInRange = c{rr}(newlag<=2 & newlag > -2);
            lagInRange = newlag(newlag<=2 & newlag > -2);
            Lag_Ach(ii,rr,:) = lagInRange;
            Corr_Ach(ii,rr,:) = corrInRange;
            
            [~,I] = max(abs(corrInRange));
            maxCorr(rr) = corrInRange(I);
            maxLag(rr) = lagInRange(I);
            
            aCorrF = autoCorrFluo{rr}(newlagFluo<=100 & newlagFluo>=0);
            aLagF= newlagFluo(newlagFluo<=100 & newlagFluo>=0);
            Lag_Achauto(ii,rr,:) = aLagF(1:size(Lag_Achauto,3));
            Corr_Achauto(ii,rr,:) = aCorrF(1:size(Lag_Achauto,3));
            
            % plot cross- and auto-correlation
%             if ~exist(fullfile(savefigpath,['CrossCorrelationRoi-',num2str(length(c)),'.fig']))
%                 figure;stem(lagInRange,corrInRange);
%                 print(gcf,'-dpng',fullfile(savefigpath,['CrossCorrelationRoi-',num2str(rr)]));
%                 saveas(gcf, fullfile(savefigpath,['CrossCorrelationRoi-',num2str(rr)]), 'fig');
%             end
%             if ~exist(fullfile(savefigpath,['autoCorrelationRoi-',num2str(length(c)),'.fig']))
%                 figure;stem(aLagF,aCorrF);
%                 print(gcf,'-dpng',fullfile(savefigpath,['autoCorrelationRoi-',num2str(rr)]));
%                 saveas(gcf, fullfile(savefigpath,['autoCorrelationRoi-',num2str(rr)]), 'fig');
%                 close all;
%             end
        end
        aCorrP = autoCorrPup(newlagPup<=100 & newlagPup>=0);
        aLagP = newlagPup(newlagPup<=100 & newlagPup>=0);
%         figure;stem(aLagP,aCorrP);
%         print(gcf,'-dpng',fullfile(savefigpath,['autoCorrelationPup']));
%         saveas(gcf, fullfile(savefigpath,['autoCorrelationPup']), 'fig');
        close all;
        
        
        %% plot the maxcorr and maxlag in 2d
        tlabel = 'Correlation coefficient';
        colorRange = [-0.3 0.3];
        plot_xcorr(maxCorr,tlabel,colorRange,savefigpath);
        % max lag
        tlabel = 'lags';
        colorRange = [-2 2];
        date = num2str(dataIndex.DateNumber(Ach(ii)));
        savefigpath = fullfile(dataIndex.BehPath{Ach(ii)},[date(1:6),'_figs-pupil']);
%         plot_xcorr(maxLag,tlabel,colorRange,savefigpath);
        %% plot average cross correlation, check how consistent across different ROIs in the same session
        
%         figure;
%         plot(squeeze(mean(Lag_Ach(ii,:,:),2)),squeeze(mean(Corr_Ach(ii,:,:),2)));
%         print(gcf,'-dpng',fullfile(savefigpath,['ave-xcorr']));
%         saveas(gcf, fullfile(savefigpath,['ave-xcorr']), 'fig');
%         % plot individual cross correlation
%         figure;
%         for zz=1:196,hold on;plot(squeeze(Lag_Ach(ii,zz,:)),squeeze(Corr_Ach(ii,zz,:)));end
%         print(gcf,'-dpng',fullfile(savefigpath,['individual-xcorr']));
%         saveas(gcf, fullfile(savefigpath,['individual-xcorr']), 'fig');
        
        % plot average auto correlation
%         figure;
%         plot(squeeze(mean(Lag_Achauto(ii,:,:),2)),squeeze(mean(Corr_Achauto(ii,:,:),2)));
%         print(gcf,'-dpng',fullfile(savefigpath,['ave-autocorr']));
%         saveas(gcf, fullfile(savefigpath,['ave-autocorr']), 'fig');
        % plot individual cross correlation
%         figure;
%         for zz=1:196,hold on;plot(squeeze(Lag_Achauto(ii,zz,:)),squeeze(Corr_Achauto(ii,zz,:)));end
%         print(gcf,'-dpng',fullfile(savefigpath,['individual-autocorr']));
%         saveas(gcf, fullfile(savefigpath,['individual-autocorr']), 'fig');
        % plot distributions of correlation and lag within single session
        
        
        close all;
        maxCorrSum_Ach = cat(1,maxCorrSum_Ach, maxCorr);
        maxLagSum_Ach= cat(1, maxLagSum_Ach, maxLag);
        
        
        
        %% summarize the lags and xcorr
        
    end
end

%% get NE data
for ii = 1:length(NE)
    savematpath = fullfile(dataIndex.BehPath{NE(ii)},'analysis-pupil');
    savefilename = fullfile(savematpath,'xcorr-pupil-fluo.mat');
    %fn_corr = dir(savefilename)
    if exist(savefilename)
        load(savefilename);
        % plot the cross-correlation with lags
        if ~exist('Lag_NE','var')
            Lag_NE= zeros(length(NE),length(c),length(newlag(newlag<=2 & newlag > -2)));
            Corr_NE = zeros(length(NE),length(c),length(newlag(newlag<=2 & newlag > -2)));
        end
        if ~exist('Lag_NEauto','var')
            Lag_NEauto= zeros(length(NE),length(c),length(newlagFluo(newlagFluo<=100 & newlagFluo>=0)));
            Corr_NEauto = zeros(length(NE),length(c),length(newlagFluo(newlagFluo<=100 & newlagFluo>=0)));
        end
        date = num2str(dataIndex.DateNumber(NE(ii)));
        savefigpath = fullfile(dataIndex.BehPath{NE(ii)},[date(1:6),'_crosscorrelation']);
        if ~exist(savefigpath)
            mkdir(savefigpath);
        end
        
        maxCorr = zeros(1, length(c));
        maxLag = zeros(1, length(c));
        for rr = 1:length(c)
            % get the maximum correlation and lags
            corrInRange = c{rr}(newlag<=2 & newlag > -2);
            lagInRange = newlag(newlag<=2 & newlag > -2);
            Lag_NE(ii,rr,:) = lagInRange;
            Corr_NE(ii,rr,:) = corrInRange;
            
            [~,I] = max(abs(corrInRange));
            maxCorr(rr) = corrInRange(I);
            maxLag(rr) = lagInRange(I);
            
            aCorrF = autoCorrFluo{rr}(newlagFluo<=100 & newlagFluo>=0);
            aLagF= newlagFluo(newlagFluo<=100 & newlagFluo>=0);
            Lag_NEauto(ii,rr,:) = aLagF(1:size(Lag_NEauto,3));
            Corr_NEauto(ii,rr,:) = aCorrF(1:size(Lag_NEauto,3));
            
            % plot cross- and auto-correlation
%             if ~exist(fullfile(savefigpath,['CrossCorrelationRoi-',num2str(length(c)),'.fig']))
%                 figure;stem(lagInRange,corrInRange);
%                 print(gcf,'-dpng',fullfile(savefigpath,['CrossCorrelationRoi-',num2str(rr)]));
%                 saveas(gcf, fullfile(savefigpath,['CrossCorrelationRoi-',num2str(rr)]), 'fig');
%             end
%             if ~exist(fullfile(savefigpath,['autoCorrelationRoi-',num2str(length(c)),'.fig']))
%                 figure;stem(aLagF,aCorrF);
%                 print(gcf,'-dpng',fullfile(savefigpath,['autoCorrelationRoi-',num2str(rr)]));
%                 saveas(gcf, fullfile(savefigpath,['autoCorrelationRoi-',num2str(rr)]), 'fig');
%                 close all;
%             end
            
        end
        aCorrP = autoCorrPup(newlagPup<=100 & newlagPup>=0);
        aLagP = newlagPup(newlagPup<=100 & newlagPup>=0);
%         figure;stem(aLagP,aCorrP);
%         print(gcf,'-dpng',fullfile(savefigpath,['autoCorrelationPup']));
%         saveas(gcf, fullfile(savefigpath,['autoCorrelationPup']), 'fig');
%         close all;
        
        
        % plot the maxcorr and maxlag in 2d
        tlabel = 'Correlation coefficient';
        colorRange = [-0.3 0.3];
        plot_xcorr(maxCorr,tlabel,colorRange,savefigpath);
        % max lag
        tlabel = 'lags';
        colorRange = [-2 2];
        date = num2str(dataIndex.DateNumber(NE(ii)));
        savefigpath = fullfile(dataIndex.BehPath{NE(ii)},[date(1:6),'_figs-pupil']);
%         plot_xcorr(maxLag,tlabel,colorRange,savefigpath);
        %% plot average cross correlation, check how consistent across different ROIs in the same session
        
%         figure;
%         plot(squeeze(mean(Lag_NE(ii,:,:),2)),squeeze(mean(Corr_NE(ii,:,:),2)));
%         print(gcf,'-dpng',fullfile(savefigpath,['ave-xcorr']));
%         saveas(gcf, fullfile(savefigpath,['ave-xcorr']), 'fig');
%         % plot individual cross correlation
%         figure;
%         for zz=1:196,hold on;plot(squeeze(Lag_NE(ii,zz,:)),squeeze(Corr_NE(ii,zz,:)));end
%         print(gcf,'-dpng',fullfile(savefigpath,['individual-xcorr']));
%         saveas(gcf, fullfile(savefigpath,['individual-xcorr']), 'fig');
%         % plot distributions of correlation and lag within single session
%         % plot average auto correlation
%             figure;
%             plot(squeeze(mean(Lag_NEauto(ii,:,:),2)),squeeze(mean(Corr_NEauto(ii,:,:),2)));
%             print(gcf,'-dpng',fullfile(savefigpath,['ave-autocorr']));
%             saveas(gcf, fullfile(savefigpath,['ave-autocorr']), 'fig');
%             % plot individual cross correlation
%             figure;
%             for zz=1:196,hold on;plot(squeeze(Lag_NEauto(ii,zz,:)),squeeze(Corr_NEauto(ii,zz,:)));end
%             print(gcf,'-dpng',fullfile(savefigpath,['individual-autocorr']));
%             saveas(gcf, fullfile(savefigpath,['individual-autocorr']), 'fig');
%         
%         close all;
        maxCorrSum_NE = cat(1,maxCorrSum_NE, maxCorr);
        maxLagSum_NE= cat(1, maxLagSum_NE, maxLag);
        
        
        
        %% summarize the lags and xcorr
        
    end
end
display('Sum');

%% plot correlation coefficient
figure;histogram(maxCorrSum_NE(:),'Normalization','probability','FaceColor',[63,167,150]/255,'EdgeColor',[63,167,150]/255);
hold on;histogram(maxCorrSum_Ach(:),'Normalization','probability','FaceColor',[255 189 53]/255,'EdgeColor',[255 189 53]/255);
set(gca,'box','off');
leg = legend('NE','Ach');
set(leg,'Box','off')
title('Distribution of correlation coefficient');
xlabel('Correlation coefficient');
ylabel('Probability');

print(gcf,'-dpng',fullfile(savesumpath,'GRAB-pupil-spon-corrcoeff'));    %png format
saveas(gcf, fullfile(savesumpath,'GRAB-pupil-spon-corrcoeff'), 'fig');
saveas(gcf, fullfile(savesumpath,'GRAB-pupil-spon-corrcoeff'),'svg');

%% plot lag
figure;histogram(maxLagSum_NE(:),'Normalization','probability','FaceColor',[63,167,150]/255,'EdgeColor',[63,167,150]/255,'BinWidth',0.1);
hold on;histogram(maxLagSum_Ach(:),'Normalization','probability','FaceColor',[255 189 53]/255,'EdgeColor',[255 189 53]/255,'BinWidth',0.1);
set(gca,'box','off');
leg = legend('NE','Ach');
set(leg,'Box','off')
title('Distribution of lags of crosscorrelation');
xlabel('Max lag (s)');
ylabel('Probability');

print(gcf,'-dpng',fullfile(savesumpath,'GRAB-pupil-spon-lag'));    %png format
saveas(gcf, fullfile(savesumpath,'GRAB-pupil-spon-lag'), 'fig');
saveas(gcf, fullfile(savesumpath,'GRAB-pupil-spon-lag'),'svg');

%% plot corr and lag within (-2,2)s
corrNE = maxCorrSum_NE(:);lagNE = maxLagSum_NE(:);
corrAch = maxCorrSum_Ach(:); lagAch = maxLagSum_Ach(:);
figure;histogram(corrNE(lagNE>-2 & lagNE<2),'Normalization','probability','FaceColor',[255 189 53]/255,'EdgeColor',[255 189 53]/255);
hold on;histogram(maxCorrSum_Ach(lagAch>-2 & lagAch<2),'Normalization','probability','FaceColor',[63,167,150]/255,'EdgeColor',[63,167,150]/255);
set(gca,'box','off');
leg = legend('NE','Ach');
set(leg,'Box','off')
title('Distribution of correlation coefficient');
xlabel('Correlation coefficient');
ylabel('Probability');

print(gcf,'-dpng',fullfile(savesumpath,'GRAB-pupil-spon-corrcoeff-(lag-2-2)'));    %png format
saveas(gcf, fullfile(savefigpath,'GRAB-pupil-spon-corrcoeff-(lag-2-2)'), 'fig');
saveas(gcf, fullfile(savefigpath,'GRAB-pupil-spon-corrcoeff-(lag-2-2)'),'svg');

%% max lag 0-
corrNE = maxCorrSum_NE(:);lagNE = maxLagSum_NE(:);
corrAch = maxCorrSum_Ach(:); lagAch = maxLagSum_Ach(:);
figure;histogram(corrNE(lagNE>0 & lagNE<2),'Normalization','probability','FaceColor',[255 189 53]/255,'EdgeColor',[255 189 53]/255);
hold on;histogram(maxCorrSum_Ach(lagAch>0 & lagAch<2),'Normalization','probability','FaceColor',[63,167,150]/255,'EdgeColor',[63,167,150]/255);
set(gca,'box','off');
leg = legend('NE','Ach');
set(leg,'Box','off')
title('Distribution of correlation coefficient');
xlabel('Correlation coefficient');
ylabel('Probability');

print(gcf,'-dpng',fullfile(savesumpath,'GRAB-pupil-spon-corrcoeff-(lag0-2)'));    %png format
saveas(gcf, fullfile(savefigpath,'GRAB-pupil-spon-corrcoeff-(lag0-2)'), 'fig');
saveas(gcf, fullfile(savefigpath,'GRAB-pupil-spon-corrcoeff-(lag0-2)'),'svg');


%% plot lag and max corr
figure
scatter(lagNE(lagNE>-2 & lagNE<2),corrNE(lagNE>-2 & lagNE<2),200,[255 189 53]/255,'.');
hold on; scatter(lagAch(lagAch>-2 & lagAch<2),corrAch(lagAch>-2 & lagAch<2),200,[63,167,150]/255,'.');
set(gca,'box','off');
leg = legend('NE','Ach');
set(leg,'Box','off')
[r_Ach,p_Ach] = corr(lagAch(lagAch>-2 & lagAch<2),corrAch(lagAch>-2 & lagAch<2))
[r_NE,p_NE] = corr(lagNE(lagNE>-2 & lagNE<2),corrNE(lagNE>-2 & lagNE<2))
xlabel('Max lag (s)');
ylabel('Correlation coefficient');
print(gcf,'-dpng',fullfile(savesumpath,'corr-coeff-lag-(lag-2-2)'));    %png format
saveas(gcf, fullfile(savefigpath,'corr-coeff-lag-(lag-2-2)'), 'fig');
saveas(gcf, fullfile(savefigpath,'corr-coeff-lag-(lag-2-2)'),'svg');