%% get mean projection

rootdir = "F:\GRABNE\100micron\893\893-02012021";

cd(rootdir)
matFile = dir(fullfile(rootdir,'mat','*.mat'));
stackFile = dir(fullfile(rootdir,'stack_info.mat'));

stackInfo = load(fullfile(stackFile.folder,stackFile.name));

totalFrames = sum(stackInfo.nFrames);


for ii = 1:length(matFile)
    tempImg = load(fullfile(matFile(ii).folder,matFile(ii).name));
    if ii==1
        meanProj = zeros(size(tempImg.stack,1),size(tempImg.stack,2));
    end
    meanProj = meanProj + sum(tempImg.stack,3)/totalFrames;
    clear tempImg
end
meanProj = uint16(meanProj);
%% load tiff tag from existing mean projection
existTiff = "F:\GRABNE\100micron\893\893-01312021\reg_stackMean.tif";
info = imfinfo(existTiff);
tagstruct.ImageLength=info(1).Height;
tagstruct.ImageWidth=info(1).Width;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = info(1).Height;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB';

savetiffname = fullfile(rootdir,'reg_stackMean.tif');

 t = Tiff(savetiffname,'w');
  t.setTag(tagstruct);
        t.write(meanProj); %write first frame
        t.close();