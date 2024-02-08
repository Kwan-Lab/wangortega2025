function GRAB_pattern(dataIndex)

% look for grab_ne/grab_ach spatial/temporal activity patterns
% raw data needed, use tone-pupil-901 for now

% 1. calculate pxwise correlation for every pixel, get the mean correlation --
% result: too many not correlated pixels

% 2. unsupervised classification of pixel activity 

% 3. spatial correlation of frames
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
    %cd(savefigpath)
   
    % 
    filepath = "F:\GRABNE\tone\tone-pupil-1\mat";

    cd(filepath);
    matfiles = dir('*.mat');

    % use desampled files first
    % use desampled frames first
    desampledPath = 'F:\GRABNE\tone\tone-pupil-1';
    desampledFile = 'file_0_DS8.tif';
    desampledmat = loadtiffseq(desampledPath,desampledFile);
    
    stepSize = 10; % calculate the crosscorrelation
    corrDist = 50; % calulate the xcorre as far as 100 px
    sizeX = floor(size(desampledmat,1)/stepSize);
    sizeY = floor(size(desampledmat,2)/stepSize);
    xcorrMat = zeros(sizeX,sizeY,corrDist*2+1, corrDist*2+1);
    
    for xx = 1:sizeX
        for yy = 1:sizeY
            display(yy);
            xAx = (xx-1)*stepSize+1; yAx = (yy-1)*stepSize+1;
            for ii = xAx-corrDist:xAx+corrDist
                
                for jj = yAx-corrDist:yAx+corrDist
                    if (ii<1) | (jj<1) | (ii>size(desampledmat,1)) | (jj>size(desampledmat,2)) % if out of bounds
                        xcorrMat(xx,yy,ii+corrDist-xAx+1,jj+corrDist-yAx+1) = NaN;
                    else
                        xcorrmat = corrcoef(double(desampledmat(xAx,yAx,:)),double(desampledmat(ii,jj,:)));
                        xcorrMat(xx,yy,ii+corrDist-xAx+1,jj+corrDist-yAx+1) = xcorrmat(1,2);
                    end
                end
            end
        end
    end
    figure;
    %% spatial frequency
    % fft2
    % normalize the figure? 
    aveSpec = zeros(size(desampledmat,1),size(desampledmat,2));
    for ii = 1:size(desampledmat,3)
        normfigure = double(desampledmat(:,:,ii))-mean(mean(desampledmat(:,:,ii)));
        Y = fft2(normfigure);
        %figure;imagesc(desampledmat(:,:,500));
        dd = abs(fftshift(Y));
        aveSpec = aveSpec+dd;
    end
    aveSpec = aveSpec/size(desampledmat,3);
    figure;
    
    imagesc(aveSpec(240:280,240:280))
    figure;surf(aveSpec)
    
    %% 
    %% cross correlation
     corrDist = 50;
     crrMat = zeros(2*corrDist+1, 2*corrDist+1,size(desampledmat,3)-1);
    for ii = 1:size(desampledmat)-1
        crr = xcorr2(desampledmat(:,:,ii),desampledmat(:,:,ii+1));
        [ssr,snd] = max(crr(:));
        [ij,ji] = ind2sub(size(crr),snd);
        midpointX = size(desampledmat,1);midpointY = size(desampledmat,2);
    % center point (512,512)
    % get distance and value
       
        crrMat(:,:,ii) = crr(midpointX-corrDist:midpointX+corrDist,midpointY-corrDist:midpointY+corrDist);
    end
    
    for jj = 261:280
        figure;surf(crrMat(:,:,jj))
    end
end