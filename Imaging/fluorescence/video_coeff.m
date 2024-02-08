function video_coeff(input,tlabel,colorRange,savefigpath)

%% make a video about selectivity change

% load the data

    % change left-right to contra-ipsi

colors=cbrewer('div','RdBu',256);
colors=flipud(colors);

% make videos
videoPath = fullfile(savefigpath,[tlabel,'regrCoeff.avi']);

if ~exist(videoPath)
    writerObj = VideoWriter(videoPath);
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for u=1:size(input,3)
        % convert the image to a frame
        figure;
        subplot(1,2,1)
        b=image(input(:,:,u),'CDataMapping','scaled');
        set(b,'AlphaData',~isnan(input(:,:,u)))
        set(gca, 'Color', [0.7, 0.7, 0.7])
        title([tlabel,' sig coeff t=',num2str(u/10-3)]);
        axis square;
        colormap(colors);
        caxis([colorRange(1) colorRange(2)]);
        subplot(3,20,55);
        image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
        colormap(colors);
        caxis([colorRange(1) colorRange(2)]);
        set(gca,'YDir','normal');
        set(gca,'XTick',[]);
        
        
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
        close;
    end
    % close the writer object
    close(writerObj);

else
     display('Video already generated');
 end
