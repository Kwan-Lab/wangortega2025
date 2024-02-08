% GRABNE_patch

% find the spatial and temporal structure of adrenergic firing

filepath = "F:\GRABNE\tone\tone-pupil-1\mat";

cd(filepath);
matfiles = dir('*.mat');

IM = load(fullfile(matfiles(1).folder,matfiles(1).name));
figure;
for i = 1:912
    imagesc(IM.stack(:,:,i))
    axis equal
    pause(.002);
end

% pixel-wise cross correlation?
% binned cross0-correlation?

% use desampled frames first
desampledPath = 'F:\GRABNE\tone\tone-pupil-1';
desampledFile = 'file_0_DS8.tif';
desampledmat = loadtiffseq(desampledPath,desampledFile);

% frequency of the desampled video is ~ 2Hz
% binned cross-correlation?
meanInt = mean(desampledmat,3);
figure;imagesc(meanInt);axis equal;


binSize = 10; % 10 frames, ~ 5s long
nBin = floor(size(desampledmat,3)/binSize);
range = 100; % check -50 to + 50 pixels
% try pixal (380,335)
corrPix = zeros(range+1,range+1,nBin);
cX = 335; cY = 380;
for xx = cX-range/2:cX+range/2
    for yy = cY-range/2:cY+range/2
        for bb = 1:nBin
            corrMat = corrcoef(double(desampledmat(cX,cY,binSize*(bb-1)+1:binSize*bb)),double(desampledmat(xx,yy,binSize*(bb-1)+1:binSize*bb))); 
            corrPix(xx-cX+range/2+1,yy-cY+range/2+1,bb) = corrMat(1,2);
        end
    end
end

% exclude surronding pixels (16 pixel)
pxExclude = 16;
corrThresh = 0.2;
minPx = 300;  % ROIs less than 300 are not considered
% also don't go through every pixel. Instead do 1:10:end
iterStep = 10;
xIter = pxExclude+1:iterStep:size(desampledmat,1)-pxExclude; yIter = pxExclude+1:iterStep:size(desampledmat,1)-pxExclude;
cellMask = [];
tic
for xx = xIter
%     if mod(xx,10) == 0
%         display(xx)
%     end
    for yy = yIter
        
        % check if given pixel already in the BW 
        if (~isempty(cellMask) & sum(cellMask(xx,yy,:)) == 0) | isempty(cellMask)
            corrPixAll = zeros(size(desampledmat,1)-pxExclude*2,  size(desampledmat,2)-pxExclude*2);
            for xX = pxExclude+1:size(desampledmat,1)-pxExclude
                for yY = pxExclude+1:size(desampledmat,1)-pxExclude
                    corrMat = corrcoef(double(desampledmat(xX,yY,:)),double(desampledmat(xx,yy,:))); 
                    corrPixAll(xX,yY) = corrMat(1,2);
                end
            end
            BW = (corrPixAll>corrThresh); %Get logical mask of pixels exceeding threshold
            BW = ~bwareaopen(~BW,10,8);
            roiMask = bwareafilt(BW,1,4);
            if sum(sum(roiMask==1)) >= minPx
                cellMask = cat(3,cellMask, roiMask); %Keep only the largest region; 'conn'=4
            end
        end
    end
end
toc

% save result
savefilepath = "F:\GRABNE\tone\tone-pupil-1\corrResult.mat";
save(savefilepath, 'cellMask');
% loop through cellMasks, exclude overlap ROIs (keep the larger one)
% overlapping area larger than 70% of the smaller ROIs

load(savefilepath);

delList = zeros(1, size(cellMask,3));
for tt = 1:size(cellMask,3)
    display(tt);
    maskInd = tt+1;
    %ifDel = 0
    while delList(tt) == 0 & maskInd <= size(cellMask,3)
        area1 = sum(sum(cellMask(:,:,tt))); area2 = sum(sum(cellMask(:,:,maskInd)));
        if sum(sum(cellMask(:,:,tt)&cellMask(:,:,maskInd))) > min(area1,area2)*0.7
            if area1 > area2
                delList(maskInd) = 1;
            else
                delList(tt) = 1;
            end
        end
        maskInd = maskInd+1;
    end
end
savedMask = cellMask(:,:,~delList);
% plot cellMask
figure;imagesc(meanInt);
axis equal;
for ii = 1:size(savedMask,3)
    hold on;
    [B,L]=bwboundaries(savedMask(:,:,ii),'noholes');
    plot(B{1}(:,2),B{1}(:,1));
end

% get centor of every shape
centerCoor = zeros(2,size(savedMask,3));
for tt = 1:size(savedMask,3)
  	props = regionprops(savedMask(:,:,tt), 'Centroid');
	centerCoor(1,tt) = props.Centroid(1);
    centerCoor(2,tt) = props.Centroid(2);
end

figure;
for tt = 1:size(savedMask,3)
    hold on;
    scatter(centerCoor(1,tt),centerCoor(2,tt),'.')
end
% calculate nearby area correlation for the center points
range = 100; % check -50 to + 50 pixels
% try pixal (380,335)
corrPix = zeros(range+1,range+1,size(savedMask,3));
for cc = 1:size(centerCoor,2)
    cX = round(centerCoor(1,cc)); cY = round(centerCoor(2,cc));
    for xx = cX-range/2:cX+range/2
        for yy = cY-range/2:cY+range/2
            if xx<size(desampledmat,1) & yy<size(desampledmat,2) & xx>0 & yy>0
                corrMat = corrcoef(double(desampledmat(cX,cY,:)),double(desampledmat(xx,yy,:))); 
                corrPix(xx-cX+range/2+1,yy-cY+range/2+1,cc) = corrMat(1,2);
            else
                corrPix(xx-cX+range/2+1,yy-cY+range/2+1,cc) = nan;
            end
        end
    end
end
% set the autocorrelation to nan
corrPix(corrPix>0.99) = nan;
figure;imagesc(nanmean(corrPix,3))

% correlation as a function of distance
