function video_selectivity(input,tlabel,colorRange,savefigpath,recordingsite)

%% make a video about selectivity change

% load the data
nCells=numel(input);
edgelength = sqrt(nCells);
selMat = zeros(edgelength,edgelength,length(input{1}.t));
% plot the grids in 49*49 matrix
for cc = 1:nCells
    if mod(cc,edgelength) == 0
        Ind2 = edgelength;
    else
        Ind2 = mod(cc,(edgelength));
    end
    if mod(cc,edgelength) == 0
        Ind1 = cc/edgelength;
    else
        Ind1 = floor(cc/edgelength)+1;
    end
    
    % change left-right to contra-ipsi
    if strcmp(tlabel,'Choice') || strcmp(tlabel,'Prev choice')
        if recordingsite == 'left'
            selMat(Ind1, Ind2,:) = -input{cc}.signal;
        else
            selMat(Ind1, Ind2,:) = input{cc}.signal;
        end
    else
        selMat(Ind1, Ind2,:) = input{cc}.signal;
    end
end

colors=cbrewer('div','RdBu',256);
colors=flipud(colors);

% make videos
videoPath = fullfile(savefigpath,[tlabel,'sel.avi']);
if ~exist(videoPath)
    writerObj = VideoWriter(videoPath);
    writerObj.FrameRate = 10;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for u=1:size(selMat,3)
        % convert the image to a frame
        figure;
        subplot(1,2,1)
        image(selMat(:,:,u),'CDataMapping','scaled');
        title([tlabel,' selectivity t=',num2str(u/10-3)]);
        axis square;
        colormap(colors);
        caxis([colorRange(1) colorRange(2)]);
        subplot(3,20,55);
        image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
        colormap(colors);
        caxis([colorRange(1) colorRange(2)]);
        if strcmp(tlabel,'Choice') || strcmp(tlabel,'Prev choice')
            title(['(contra-ipsi)/(contra+ipsi)']);
        else
            title(['(R-NR)/(R+NR)']);
        end
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