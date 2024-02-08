function tempData = getRegautoCorrData(reg_cr,label,Ind, timeRange, sigThresh, savefluofigpath)

Mat = selectRg2D(reg_cr,Ind,[]);

tempData.coeff=reshape(Mat.coeff,size(Mat.coeff,1)*size(Mat.coeff,2),size(Mat.coeff,3));
tempData.pval=reshape(Mat.pval,size(Mat.coeff,1)*size(Mat.coeff,2),size(Mat.coeff,3));

% plot absolut coefficient with significant grids only
isSig = 2;
% get significant coefficient only
sigCoeff = Mat.coeff;
% determine the significant grids (threshold: at least 3 consecutive significant point
% or 10% of the time points are significant
sigMat = zeros(size(Mat.coeff,1),size(Mat.coeff,2));

% use the grids with most significant point as reference
maxSig = 0;maxX = 0; maxY = 0;
for ss = 1:size(Mat.coeff,1)
    for uu = 1:size(Mat.coeff,2)
        sigVec = squeeze(Mat.pval(ss,uu,:)<sigThresh.alpha);
        start1 = strfind([0,sigVec'==1],[0 1]);
        end1 = strfind([sigVec'==1,0],[1 0]);
        if sum((end1-start1+1)>=3) | sum(sigVec) >= 8
            sigMat(ss,uu) = 1;
            if sum(sigVec) > maxSig
                maxSig = sum(sigVec);
                maxX = ss; maxY = uu;
            end
        end
    end
end
tempData.maxX = maxX;
tempData.maxY = maxY;

for zz = 1:size(Mat.coeff,3)
    temp = squeeze(Mat.coeff(:,:,zz));
    temp(~sigMat) = NaN;
    sigCoeff(:,:,zz) = temp;
end
coeffBound.min = floor(min(min(min(sigCoeff)))*10)/10;
coeffBound.max = ceil(max(max(max(sigCoeff)))*10)/10;
savepath = fullfile(savefluofigpath,'tempCorr');
if ~exist(savepath)
    mkdir(savepath)
end

% temporal correlation of coefficient, for significant grids only
        t = reg_cr{1}.regr_time; tInd = (t>=timeRange(1) & t<=timeRange(2));

tempData.tempCorrCoeff = NaN(size(Mat.coeff,1),size(Mat.coeff,2));
tempData.tempCorrLag = NaN(size(Mat.coeff,1),size(Mat.coeff,2));
tempData.tempCorrMax = NaN(size(Mat.coeff,1),size(Mat.coeff,2));
tempData.corrcoef = NaN(size(Mat.coeff,1),size(Mat.coeff,2),sum(tInd)*2-1);
tempData.lag = NaN(size(Mat.coeff,1),size(Mat.coeff,2),sum(tInd)*2-1);
index = 1:size(tempData.coeff,1);
tempData.isSig = logical(sigMat(:));

if maxX~=0 
refCoeff = squeeze(Mat.coeff(maxX,maxY,:));
if sum(isnan(refCoeff)) > 0
    refCoeff = fillmissing(refCoeff,'linear');
end
for ss = 1:size(Mat.coeff,1)
    for uu = 1:size(Mat.coeff,2)    
        % fill the NaN values
        if sigMat(ss,uu)==1  % if the grid is significant
            sigCoeff1 = sigCoeff(ss,uu,:);
            if sum(isnan(sigCoeff1))>0

                sigCoeff1 = fillmissing(sigCoeff1,'linear');
            end
            % calculate temporal correlation within a time range
    
            [tempData.corrcoef(ss,uu,:),tempData.lag(ss,uu,:)] = xcorr(sigCoeff1(tInd),refCoeff(tInd),'normalized');
            tempData.tempCorrCoeff(ss,uu) = tempData.corrcoef(ss,uu,sum(tInd));
            % find max
            [~,maxLag] = max(abs(squeeze(tempData.corrcoef(ss,uu,:))));
            tempData.tempCorrLag(ss,uu) = tempData.lag(ss,uu,maxLag)*0.1;
            tempData.tempCorrMax(ss,uu) = tempData.corrcoef(ss,uu,maxLag);
        end
    end
end
end
% plot the lags
 colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
     titleText = ['Cross correlation lags ',label];
     
figure;
subplot(1,2,1)
colorRange=[-0.5 0.5];
b=image(tempData.tempCorrLag,'CDataMapping','scaled'); %,
set(b,'AlphaData',~isnan(tempData.tempCorrLag));
set(gca, 'Color', [0.7, 0.7, 0.7])
axis square;
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
title(titleText);
hold on;
% show the reference grid
plot(maxY,maxX,'g+')
%get a function to plot different Masks
% plotbound(Mask1);plotbound(Mask2);%plotbound(Mask3);%plotbound(Mask4);
% title([label,'-selecvitity-mask']);
subplot(3,20,60);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
print(gcf,'-dpng',fullfile(savepath,titleText));
saveas(gcf, fullfile(savepath,titleText), 'fig');
saveas(gcf, fullfile(savepath,titleText), 'svg');

% plot the maximum correlation
 colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
    colorRange = [-1 1];
     titleText = ['Cross correlation coefficient ',label];
     
figure;
subplot(1,2,1)
b=image(tempData.tempCorrMax,'CDataMapping','scaled'); %,
set(b,'AlphaData',~isnan(tempData.tempCorrMax));
set(gca, 'Color', [0.7, 0.7, 0.7])
axis square;
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
title(titleText);
hold on;
% show the reference grid
plot(maxY,maxX,'g+')
%get a function to plot different Masks
% plotbound(Mask1);plotbound(Mask2);%plotbound(Mask3);%plotbound(Mask4);
% title([label,'-selecvitity-mask']);
subplot(3,20,60);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
print(gcf,'-dpng',fullfile(savepath,titleText));
saveas(gcf, fullfile(savepath,titleText), 'fig');
saveas(gcf, fullfile(savepath,titleText), 'svg');

close all
