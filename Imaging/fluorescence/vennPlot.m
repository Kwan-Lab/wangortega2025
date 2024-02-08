function vennPlot(sig1,sig2,sig3,label,savesumfigpath)

% plot venn plot
a1 = length(sig1); a2 = length(sig2); a3 = length(sig3);
a12 = length(intersect(sig1,sig2)); a13 = length(intersect(sig1,sig3));a23 = length(intersect(sig2,sig3));
a123 = length(intersect(intersect(sig1,sig2),sig3));
figure;
[H,S]=venn([a1,a2,a3],[a12,a13,a23,a123],'Labels',label,'ErrMinMode','None');
axis equal
set(gca,'XColor', 'none','YColor','none');
titletext = '';
for tt = 1:length(label)
    titletext = [titletext,'-',label{tt}];
end
title(titletext(2:end));

print(gcf,'-dpng',fullfile(savesumfigpath,['Ratio of significant grids', titletext,'-venn']));
saveas(gcf, fullfile(savesumfigpath,['Ratio of significant grids', titletext,'-venn']), 'fig');
saveas(gcf, fullfile(savesumfigpath,['Ratio of significant grids', titletext,'-venn']), 'svg');
