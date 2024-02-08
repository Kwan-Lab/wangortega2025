function fluo_pupil_corrSummary(dataIndex)

% summarize cross-correlation results of pupil and GRAB signal

nFiles = size(dataIndex,1);
colors=cbrewer('div','RdBu',256);
colors=flipud(colors);

maxLagSum = [];
maxCorrSum = [];
for ii = 1:nFiles
     savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
     savefilename = fullfile(savematpath,'xcorr-pupil-fluo.mat');
     %fn_corr = dir(savefilename)
     if exist(savefilename)
         load(savefilename);
         
         maxCorr = zeros(1, length(c));
         maxLag = zeros(1, length(c));
         for rr = 1:length(c)
             % get the maximum correlation and lags
             corrInRange = c{rr}(newlag{rr}<=50 & newlag{rr} > -50);
             lagInRange = newlag{rr}(newlag{rr}<=50 & newlag{rr} > -50);
             [~,I] = max(abs(corrInRange));
             maxCorr(rr) = corrInRange(I);
             maxLag(rr) = lagInRange(I);
         end
         
         % plot the cross-correlation with lags
          savefigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
          
         figure;stem(newlag{rr},c{rr});
         print(gcf,'-dpng',fullfile(savefigpath,'CrossCorrelation'));
         saveas(gcf, fullfile(savefigpath,'CrossCorrelation'), 'fig');
         % plot the maxcorr and maxlag in 2d
         tlabel = 'Correlation coefficient';
         colorRange = [-0.3 0.3];
         date = num2str(dataIndex.DateNumber(ii));
        
         plot_xcorr(maxCorr,tlabel,colorRange,savefigpath);
         % max lag
          tlabel = 'lags';
         colorRange = [-50 50];
         date = num2str(dataIndex.DateNumber(ii));
         savefigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
         plot_xcorr(maxLag,tlabel,colorRange,savefigpath)
         % plot distributions of correlation and lag within single session
         
         close all;
         for rr = 1:length(reg_cr) 
            maxCorrSum = cat(3,maxCorrSum, reg_cr{rr}.coeff);
            maxLagSum = cat(3, maxLagSum, reg_cr{rr}.pval);
        end
     end
     
     %% summarize the lags and xcorr
     display('Sum');
end



