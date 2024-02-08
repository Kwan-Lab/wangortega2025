function plot_RPE(posRPECoeff_new,negRPECoeff_new,pRPESigInd,nRPESigInd,rt,savesumfigpath)

posRPE = [];
negRPE = [];

for ii = 1:length(posRPECoeff_new)
    posRPE = [posRPE;posRPECoeff_new{ii}];
    negRPE = [negRPE; negRPECoeff_new{ii}];
end

sigBoth = intersect(pRPESigInd,nRPESigInd);
sigP = setdiff(pRPESigInd,nRPESigInd);
sigN = setdiff(nRPESigInd,pRPESigInd);

pRPESigB = posRPE(sigBoth,:);
nRPESigB = negRPE(sigBoth,:);
% sort posRPE
 pR.coeff=pRPESigB;pR.t = rt;
    pRsortOrd = coeff_sort(pR,[0,3]);
    pRPESigB= pRPESigB(pRsortOrd,:);
    nRPESigB = nRPESigB(pRsortOrd,:);

h1=figure;
subplot(1,3,1);

image(rt,1:length(sigBoth),pRPESigB,'CDataMapping','scaled');
hold on; plot([0 0],[0 length(sigBoth)+1],'w');

colors=cbrewer('div','RdBu',256);
colormap(colors);

colorRange(1)=-2;
colorRange(2)=2;
%normalize dF/F heatmap to max of all conditions
% end
caxis([colorRange(1) colorRange(2)]);
ylabel('Cells');

% neg
subplot(1,3,2);

image(rt,1:length(sigBoth),nRPESigB,'CDataMapping','scaled');
hold on; plot([0 0],[0 length(sigBoth)+1],'w');

colors=cbrewer('div','RdBu',256);
colormap(colors);

colorRange(1)=-2;
colorRange(2)=2;
%normalize dF/F heatmap to max of all conditions
% end
caxis([colorRange(1) colorRange(2)]);
ylabel('Cells');

%
h1=figure;
subplot(1,3,1);

image(rt,1:length(sigBoth),pRPESigB,'CDataMapping','scaled');
hold on; plot([0 0],[0 length(sigBoth)+1],'w');

colors=cbrewer('div','RdBu',256);
colormap(colors);

colorRange(1)=-2;
colorRange(2)=2;
%normalize dF/F heatmap to max of all conditions
% end
caxis([colorRange(1) colorRange(2)]);
ylabel('Cells');

%% neg only
h1=figure;
subplot(1,3,1);
    
image(rt,1:length(sigN),negRPE(sigN,:),'CDataMapping','scaled');
hold on; plot([0 0],[0 length(sigN)+1],'w');

colors=cbrewer('div','RdBu',256);
colormap(colors);

colorRange(1)=-1;
colorRange(2)=1;
%normalize dF/F heatmap to max of all conditions
% end
caxis([colorRange(1) colorRange(2)]);
ylabel('Cells');
% %subplot
% subplot(3,20,48);
% image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);