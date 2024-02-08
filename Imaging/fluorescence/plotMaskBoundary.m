function plotMaskBoundary(mergeMask,Mask1,Mask2, label,savefluofigpath,isSig)

%% isSig: 
%       0:coefficient sign (1 for positive,-1 for negative)
%       1:pval
%       2:coefficient value 
if isSig==1  % signficance mask
    colorRange = [0 1];
    colors=cbrewer('seq','Reds',256);
    titleText = ['Average ', label, ' significance mask'];
elseif isSig == 0   
    colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
    colorRange = [-1 1];
     titleText = ['Average ', label, ' coefficient sign mask'];
elseif isSig == 2
    colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
    
    colorRange = [-0.05, 0.05];
     titleText = ['Average ', label, ' sig coefficient value mask'];
     elseif isSig == 3
    colors=cbrewer('div','RdBu',256);
    colors=flipud(colors);
    
    colorRange = [-0.05, 0.05];
     titleText = ['Average ', label, ' coefficient value mask'];
end


figure;
subplot(1,2,1)
b=image(nanmean(mergeMask,3),'CDataMapping','scaled'); %,
set(b,'AlphaData',~isnan(nanmean(mergeMask,3)))
set(gca, 'Color', [0.7, 0.7, 0.7])
axis square;
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
title(titleText);
hold on;
%get a function to plot different Masks
% plotbound(Mask1);plotbound(Mask2);%plotbound(Mask3);%plotbound(Mask4);
% title([label,'-selecvitity-mask']);
subplot(3,20,60);
image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
colormap(colors);
caxis([colorRange(1) colorRange(2)]);
print(gcf,'-dpng',fullfile(savefluofigpath,titleText));
%saveas(gcf,fullfile(savefluofigpath,titleText),'png');
saveas(gcf, fullfile(savefluofigpath,titleText), 'fig');
saveas(gcf, fullfile(savefluofigpath,titleText), 'svg');

close;
