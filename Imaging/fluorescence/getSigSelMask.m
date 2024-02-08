function sigMask = getSigSelMask(input, alpha,label, savefluofigpath)

%% take in 2D selectivity matrix,
%  get signed selectivity mask
colors=cbrewer('seq','Reds',256);
%colors=flipud(colors);
% scale color range later


edgelength = size(input,1);
sigMask = zeros(edgelength,edgelength,size(input,3));
sigMask(input<alpha) = 1;
sigMask(input>=alpha) = 0;

% figure;
% subplot(1,2,1)
% image(mean(signMask,3),'CDataMapping','scaled');
% axis square;
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);
% title('Average positive choice selecvitity mask');
% subplot(3,20,60);
% image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);
% print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-sign mask'));
% saveas(gcf, fullfile(savefluofigpath,'c selecvitity-sign mask'), 'fig');

%% get pos/neg/merged masks
 
            sigMask = zeros(edgelength,edgelength,size(input,3));
            for tt = 1:size(input,3)
                sigBW = (input(:,:,tt)<alpha); %Get logical mask of pixels exceeding threshold
                sigBW = ~bwareaopen(~sigBW,50,8); %Remove small holes from pixel mask
                sigMask(:,:,tt) = bwareafilt(sigBW,6,4); 
            end
            
            colorRange = [0 1];
            
            figure;
            subplot(1,2,1)
            image(mean(sigMask,3),'CDataMapping','scaled');
            
            axis square;
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]); 
            title(['Average significant ',label, ' selecvitity mask']);
            subplot(3,20,60);
            image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]);
            print(gcf,'-dpng',fullfile(savefluofigpath,[label,'-significance mask']));
            saveas(gcf, fullfile(savefluofigpath,[label,'-significance mask']), 'fig');

            close;
  %% calculate the area with bwareaopen()
%             posMask = zeros(edgelength,edgelength,size(input,3));
%             negMask = zeros(edgelength,edgelength,size(input,3));
%             for tt = 1:size(input,3)
%                 posBW = (input(:,:,tt)>0); %Get logical mask of pixels exceeding threshold
%                 posBW = ~bwareaopen(~posBW,50,8); %Remove small holes from pixel mask
%                 posMask(:,:,tt) = bwareafilt(posBW,6,4);
%                 negBW = (input(:,:,tt)<0); %Get logical mask of pixels exceeding threshold
%                 negBW = ~bwareaopen(~negBW,50,8); %Remove small holes from pixel mask
%                 negMask(:,:,tt) = bwareafilt(negBW,6,4);
%             end
%             figure;
%             subplot(1,2,1)
%             image(mean(posMask,3),'CDataMapping','scaled');
%             axis square;
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]); 
%             title('Average positive choice selecvitity mask');
%             subplot(3,20,60);
%             image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%             print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-pos mask'));
%             saveas(gcf, fullfile(savefluofigpath,'c selecvitity-pos mask'), 'fig');
% 
%             figure;
%             subplot(1,2,1)
%             image(mean(negMask,3),'CDataMapping','scaled');
%             axis square;
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%              title('Average negative choice selecvitity mask');
%             subplot(3,20,60);
%             image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%             print(gcf,'-dpng',fullfile(savefluofigpath,'c selecvitity-neg mask'));
%             saveas(gcf, fullfile(savefluofigpath,'c selecvitity-neg mask'), 'fig');

%% directly averaging the sign of every grid
% pos = 1; neg = -1
%% determine the area with watershed

% meanSignMask = mean(signMask,3);
% %tt = bwareaopen(meanSignMask,10,8);
% figure;imagesc(meanSignMask)
% absMask = abs(meanSignMask);
% figure;imagesc(absMask,'CDataMapping','scaled');
% axis square;
% colormap(colors);caxis([colorRange(1) colorRange(2)]);
% 
% 
% %figure;histogram(meanSignMask(:));
% absMasktt = imhmin(absMask,0.2);
% figure;imagesc(absMasktt,'CDataMapping','scaled');
% axis square;
% colormap(colors);caxis([colorRange(1) colorRange(2)]);
% 
% % convert the figure to binary
% level = graythresh(absMasktt);
% BW = im2bw(absMasktt,level);
% 
% BW = ~bwareaopen(~BW,10,8);
% BW = bwareaopen(BW,10,8);
% D = bwdist(~BW);
% %D(~BW) = -Inf;
% L = watershed(D);
% figure;image(label2rgb(L,'jet','w'))
% imshow(D,[])
% title('Distance Transform of Binary Image')
% I = imhmin(meanSignMask,0.9); %20 is the height threshold for suppressing shallow minima
% figure;image(I,'CDataMapping','scaled')
% axis square;
% colormap(colors);caxis([colorRange(1) colorRange(2)]);
% 
% L = watershed(D);
% figure;imagesc(L)