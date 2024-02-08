function mergeMask = getSignedSelMask(input, alpha, label, savefluofigpath)

%% take in 2D selectivity matrix,
%  get signed selectivity mask
colors=cbrewer('div','RdBu',256);
colors=flipud(colors);
colorRange = [-1 1];

edgelength = size(input.coeff,1);
signMask = zeros(edgelength,edgelength,size(input.coeff,3));
signMask(input.coeff>0 & input.pval<alpha) = 1;
signMask(input.coeff<0 & input.pval<alpha) = -1;
signMask(input.pval>=alpha) = 0;

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
   posMask = zeros(edgelength,edgelength,size(input.coeff,3));
            negMask = zeros(edgelength,edgelength,size(input.coeff,3));
            mergeMask = zeros(edgelength,edgelength,size(input.coeff,3));
            for tt = 1:size(input.coeff,3)
                posBW = double((input.coeff(:,:,tt)>0 & input.pval(:,:,tt)<alpha)); 
                posBW(posBW==0)=NaN;%Get logical mask of pixels exceeding threshold
                %posBW = ~bwareaopen(~posBW,50,8); %Remove small holes from pixel mask
                %posMask(:,:,tt) = bwareafilt(posBW,6,4);
                posMask(:,:,tt) = posBW;
                negBW = double((input.coeff(:,:,tt)<0  & input.pval(:,:,tt)<alpha)); %Get logical mask of pixels exceeding threshold
                negBW(negBW==0)=NaN;
                %negBW = ~bwareaopen(~negBW,50,8); %Remove small holes from pixel mask
                %negMask(:,:,tt) = -bwareafilt(negBW,6,4);
                negMask(:,:,tt) = -negBW;
                
                % merge the pos/neg Mask (pos px = 1; neg px = -1);
                tempmerge = zeros(edgelength,edgelength);
                temppos = posMask(:,:,tt); tempneg = negMask(:,:,tt);
                tempmerge(temppos==1) =1;
                tempmerge(tempneg==-1) = -1;
                tempmerge(temppos==1 & tempneg==-1) = 0;
                mergeMask(:,:,tt) = tempmerge;   
            end

%             figure;
%             subplot(1,2,1)
%             image(mean(posMask,3),'CDataMapping','scaled');
%             axis square;
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]); 
%             title(['Average positive ',label, ' selecvitity mask']);
%             subplot(3,20,60);
%             image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%             print(gcf,'-dpng',fullfile(savefluofigpath,[label, '-selecvitity-pos mask']));
%             saveas(gcf, fullfile(savefluofigpath,[label, '-selecvitity-pos mask']), 'fig');

%             figure;
%             subplot(1,2,1)
%             image(mean(negMask,3),'CDataMapping','scaled');
%             axis square;
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%              title(['Average negative ',label, ' selecvitity mask']);
%             subplot(3,20,60);
%             image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%             colormap(colors);
%             caxis([colorRange(1) colorRange(2)]);
%             print(gcf,'-dpng',fullfile(savefluofigpath,[label, '-selecvitity-neg mask']));
%             saveas(gcf, fullfile(savefluofigpath,[label, '-selecvitity-neg mask']), 'fig');
%             
             figure;
            subplot(1,2,1)
            image(nanmean(mergeMask,3),'CDataMapping','scaled');
            axis square;
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]); 
            title(['Average merged ',label, ' significant coefficient mask']);
            subplot(3,20,60);
            image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
            colormap(colors);
            caxis([colorRange(1) colorRange(2)]);
            print(gcf,'-dpng',fullfile(savefluofigpath,[label, '-sigCoeff-merge mask']));
            saveas(gcf, fullfile(savefluofigpath,[label, '-sigCoeff-merge mask']), 'fig');

            close
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