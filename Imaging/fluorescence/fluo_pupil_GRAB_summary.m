function fluo_pupil_GRAB_summary(dataIndex,savefigpath)

% summarize the results of pupil/fluorescent correlation
% GRAB_Ach, before MP, 5,8,10,12

Ind = 1:size(dataIndex,1);
ACh = Ind(strcmp(dataIndex.GRAB,'ACh'));
NE = Ind(strcmp(dataIndex.GRAB,'NE'));

nFiles = size(dataIndex,1);
corrCoeffMat_ACh = [];
corrPMat_ACh = [];
corrCoeffMat_NE= [];
corrPMat_NE = [];
for ii = 1:length(Ach)
     fn_beh = dir(fullfile(dataIndex.BehPath{Ach(ii)},'*beh.mat'));
    savematname = fullfile(fn_beh.folder,'analysis-pupil','corrPupFluo.mat');
    load(savematname);
    corrCoeffMat_ACh = [corrCoeffMat_ACh,corrCoeff(:)'];
    corrPMat_Ach = [corrPMat_Ach, corrP(:)'];
end

for ii = 1:length(NE)
     fn_beh = dir(fullfile(dataIndex.BehPath{NE(ii)},'*beh.mat'));
    savematname = fullfile(fn_beh.folder,'analysis-pupil','corrPupFluo.mat');
    load(savematname);
    corrCoeffMat_NE = [corrCoeffMat_NE,corrCoeff(:)'];
    corrPMat_NE = [corrPMat_NE, corrP(:)'];
end

% plot the correlation histogram with significance indicated
maxCorr_Ach = max(corrCoeffMat_Ach);
minCorr_Ach = min(corrCoeffMat_Ach);
binSize = 0.01;
minBin = floor(minCorr_Ach/binSize)*binSize; maxBin = ceil(maxCorr_Ach/binSize)*binSize;
binCenter = minBin+binSize/2:binSize:maxBin-binSize/2;
numBin = length(binCenter);
binMat_Ach = zeros(numBin,1); binSig_Ach = zeros(numBin,1);

for cc = 1:length(corrCoeffMat_Ach)
    ind = floor((corrCoeffMat_Ach(cc)-minBin)/binSize)+1;
    binMat_Ach(ind) = binMat_Ach(ind)+1;
    if corrPMat_Ach(cc) < 0.05
        binSig_Ach(ind) = binSig_Ach(ind)+1;
    end
end

figure;bar(binCenter,binMat_Ach,'EdgeColor',[0.7,0.7,0.7],'FaceColor',[0.7,0.7,0.7])
hold on;bar(binCenter,binSig_Ach,'EdgeColor',[255 189 53]/255,'FaceColor',[255 189 53]/255);
set(gca,'box','off');
xlabel('Correlation coefficient');
ylabel('Number of ROIs');
title('Correlation distribution of spontaneous Ach and pupil');
print(gcf,'-dpng',fullfile(savefigpath,'Ach-pupil'));    %png format
saveas(gcf, fullfile(savefigpath,'Ach-pupil'), 'fig');
saveas(gcf, fullfile(savefigpath,'Ach-pupil'),'svg');

maxCorr_NE = max(corrCoeffMat_NE);
minCorr_NE = min(corrCoeffMat_NE);
binSize = 0.01;
minBin = floor(minCorr_NE/binSize)*binSize; maxBin = ceil(maxCorr_NE/binSize)*binSize;
binCenter = minBin+binSize/2:binSize:maxBin-binSize/2;
numBin = length(binCenter);
binMat_NE = zeros(numBin,1); binSig_NE = zeros(numBin,1);

for cc = 1:length(corrCoeffMat_NE)
    ind = floor((corrCoeffMat_NE(cc)-minBin)/binSize)+1;
    binMat_NE(ind) = binMat_NE(ind)+1;
    if corrPMat_NE(cc) < 0.05
        binSig_NE(ind) = binSig_NE(ind)+1;
    end
end

figure
bar(binCenter,binMat_NE,'EdgeColor',[0.7,0.7,0.7],'FaceColor',[0.7,0.7,0.7])
hold on;bar(binCenter,binSig_NE,'EdgeColor',[63,167,150]/255,'FaceColor',[63,167,150]/255);
set(gca,'box','off');
xlabel('Correlation coefficient');ylabel('Number of ROIs');
title('Correlation distribution of spontaneous NE and pupil');
print(gcf,'-dpng',fullfile(savefigpath,'NE-pupil'));    %png format
saveas(gcf, fullfile(savefigpath,'NE-pupil'), 'fig');
saveas(gcf, fullfile(savefigpath,'NE-pupil'),'svg');
