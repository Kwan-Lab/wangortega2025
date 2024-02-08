function MP_GRAB_tempComp(save_path_mat_NE, save_path_mat_ACh,save_path_fluo_ACh)


varList = {'choice', 'outcome', 'interaction', 'dQ', 'dK', 'posRPE', 'negRPE', 'CKE'};

% get rid of sessions with too few c/x significant grids?

%% compare overlapping area of ACh and NE
sigACh = load(fullfile(save_path_mat_ACh,'sigList.mat'),'SigList');
sigNE = load(fullfile(save_path_mat_NE,'sigList.mat'),'SigList');
nFileACh = max(sigACh.SigList.session);
nFileNE = max(sigNE.SigList.session);

cList_ACh = [];oList_ACh = [];xList_ACh=[];
pRPEList_ACh = []; nRPEList_ACh = [];

ocOverlap_ACh = zeros(1,nFileACh);
coOverlap_ACh = zeros(1,nFileACh);
cxOverlap_ACh = zeros(1,nFileACh);
xcOverlap_ACh = zeros(1,nFileACh);
oxOverlap_ACh = zeros(1,nFileACh);
xoOverlap_ACh = zeros(1,nFileACh);
pnRPEOverlap_ACh = zeros(1,nFileACh);
npRPEOverlap_ACh = zeros(1,nFileACh);

% number of pRPE and nRPE
pRPE_ACh = zeros(1,nFileACh);
nRPE_ACh = zeros(1,nFileACh);

minGridRPE = 13;
minGrid = 100;
for aa = 1:nFileACh
    c = sigACh.SigList.c(sigACh.SigList.session==aa);
    cList_ACh = [cList_ACh,c];
    o = sigACh.SigList.r(sigACh.SigList.session==aa);
    oList_ACh = [oList_ACh,o];
    x = sigACh.SigList.x(sigACh.SigList.session==aa);
    xList_ACh = [xList_ACh,x];
    pRPE = sigACh.SigList.pRPE(sigACh.SigList.session==aa);
    pRPEList_ACh = [pRPEList_ACh,pRPE];
    nRPE = sigACh.SigList.nRPE(sigACh.SigList.session==aa);
    nRPEList_ACh = [nRPEList_ACh,nRPE];

    pRPE_ACh(aa) = sum(pRPE);
    nRPE_ACh(aa) = sum(nRPE);

    if sum(nRPE)>minGridRPE | sum(pRPE) > minGridRPE
        pnRPEOverlap_ACh(aa) = sum(pRPE&nRPE)/sum(nRPE);
        npRPEOverlap_ACh(aa) = sum(pRPE&nRPE)/sum(pRPE);
    else
        pnRPEOverlap_ACh(aa) = NaN;
        npRPEOverlap_ACh(aa) = NaN;
    end
    if sum(c)>minGrid
        ocOverlap_ACh(aa) = sum(c&o)/sum(c);
        coOverlap_ACh(aa) = sum(c&o)/sum(o);
    else
        ocOverlap_ACh(aa)=NaN;
        coOverlap_ACh(aa)=NaN;
    end
      if sum(c)>minGrid & sum(x) > minGrid
            cxOverlap_ACh(aa) = sum(c&x)/sum(x);
            xcOverlap_ACh(aa) = sum(c&x)/sum(c);
    else
        cxOverlap_ACh(aa) = nan;
        xcOverlap_ACh(aa) = nan;
    end
    if sum(x) > minGrid
        oxOverlap_ACh(aa) = sum(o&x)/sum(x);
        xoOverlap_ACh(aa) = sum(o&x)/sum(o);
    else
        oxOverlap_ACh(aa) = NaN;
        xoOverlap_ACh(aa) = NaN;
    end
end
ocOverlap_NE = zeros(1,nFileNE);
coOverlap_NE = zeros(1,nFileNE);
cxOverlap_NE = zeros(1,nFileNE);
xcOverlap_NE = zeros(1,nFileNE);
oxOverlap_NE = zeros(1, nFileNE);
xoOverlap_NE = zeros(1, nFileNE);
pnRPEOverlap_NE = zeros(1,nFileNE);
npRPEOverlap_NE = zeros(1,nFileNE);

cList_NE = [];oList_NE = [];xList_NE=[];
pRPEList_NE = []; nRPEList_NE = [];

% number of pRPE and nRPE
pRPE_NE = zeros(1,nFileNE);
nRPE_NE = zeros(1,nFileNE);

for nn = 1:nFileNE
    c = sigNE.SigList.c(sigNE.SigList.session==nn);
       cList_NE = [cList_NE,c];
    o = sigNE.SigList.r(sigNE.SigList.session==nn);
       oList_NE = [oList_NE,o];
    x = sigNE.SigList.x(sigNE.SigList.session==nn);
        xList_NE = [xList_NE,x];
       pRPE = sigNE.SigList.pRPE(sigNE.SigList.session==nn);
    pRPEList_NE = [pRPEList_NE,pRPE];
    nRPE = sigNE.SigList.nRPE(sigNE.SigList.session==nn);
    nRPEList_NE = [nRPEList_NE,nRPE];

     pRPE_NE(nn) = sum(pRPE);
    nRPE_NE(nn) = sum(nRPE);

    if sum(nRPE)>minGridRPE | sum(pRPE) > minGridRPE
        pnRPEOverlap_NE(nn) = sum(pRPE&nRPE)/sum(nRPE);
        npRPEOverlap_NE(nn) = sum(pRPE&nRPE)/sum(pRPE);
    else
        pnRPEOverlap_NE(nn) = NaN;
        npRPEOverlap_NE(nn) = NaN;
    end

    if sum(c)>minGrid
        ocOverlap_NE(nn) = sum(c&o)/sum(c);
        coOverlap_NE(nn) = sum(c&o)/sum(o);
    else
        coOverlap_NE(nn) = NaN;
        ocOverlap_NE(nn) = NaN;
    end
    if sum(c)>minGrid & sum(x) > minGrid
        cxOverlap_NE(nn) = sum(c&x)/sum(x);
        xcOverlap_NE(nn) = sum(c&x)/sum(c);
    else
        cxOverlap_NE(nn) = nan;
        xcOverlap_NE(nn) = nan;
    end
     if sum(x)>minGrid
      oxOverlap_NE(nn) = sum(o&x)/sum(x);
      xoOverlap_NE(nn) = sum(o&x)/sum(o);
     else
         oxOverlap_NE(nn) = NaN;
         xoOverlap_NE(nn) = NaN;
     end
end

%% chi-square test for independence
%% NE
% choice-outcome
pco = x2testOverlap(cList_NE,oList_NE,length(cList_NE))
% choice-interaction
pcx = x2testOverlap(cList_NE,xList_NE,length(cList_NE))
% outcome-interaction
pox = x2testOverlap(oList_NE,xList_NE,length(cList_NE))

   %% ACh
% choice-outcome
pco = x2testOverlap(cList_ACh,oList_ACh,length(cList_ACh))
% choice-interaction
pcx = x2testOverlap(cList_ACh,xList_ACh,length(cList_ACh))
% outcome-interaction
pox = x2testOverlap(oList_ACh,xList_ACh,length(cList_ACh))


pOC = mediantest(ocOverlap_ACh(~isnan(ocOverlap_ACh)), ocOverlap_NE(~isnan(ocOverlap_NE)))
pCO = mediantest(coOverlap_ACh(~isnan(coOverlap_ACh)), coOverlap_NE(~isnan(coOverlap_NE)))
pOX = mediantest(oxOverlap_ACh(~isnan(oxOverlap_ACh)), oxOverlap_NE(~isnan(oxOverlap_NE)))
pXO = mediantest(xoOverlap_ACh(~isnan(xoOverlap_ACh)), xoOverlap_NE(~isnan(xoOverlap_NE)))
pXC = mediantest(xcOverlap_ACh(~isnan(xcOverlap_ACh)), xcOverlap_NE(~isnan(xcOverlap_NE)))
pCX = mediantest(cxOverlap_ACh(~isnan(cxOverlap_ACh)), cxOverlap_NE(~isnan(cxOverlap_NE)))

Group = [ones(1,length(coOverlap_NE)),2*ones(1,length(coOverlap_ACh))];
figure;
subplot(1,6,1);
boxplot([ocOverlap_NE,ocOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
ylabel('Fraction of grids');
title('P(r|c)')
ylim([0 1])
subplot(1,6,2);
boxplot([coOverlap_NE,coOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
ylabel('Fraction of grids');
title('P(c|r)')
ylim([0 1])
subplot(1,6,3)
boxplot([oxOverlap_NE,oxOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(r|x)')
ylim([0 1])
subplot(1,6,4)
boxplot([xoOverlap_NE,xoOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(x|r)')
ylim([0 1])

subplot(1,6,5)
boxplot([cxOverlap_NE,cxOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(c|x)')
ylim([0 1])
subplot(1,6,6)
boxplot([xcOverlap_NE,xcOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(x|c)')
ylim([0 1])
print(gcf,'-dpng',fullfile(save_path_fluo_ACh,'Percentage of overlapping grids with outcome'));
saveas(gcf, fullfile(save_path_fluo_ACh,'Percentage of overlapping grids with outcome'), 'fig');
saveas(gcf, fullfile(save_path_fluo_ACh,'Percentage of overlapping grids with outcome'), 'svg');

% pos/neg RPE
figure;
subplot(1,2,1)
boxplot([pnRPEOverlap_NE,pnRPEOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(pRPE|nRPE)')
%ylim([0 1])
subplot(1,2,2)
boxplot([npRPEOverlap_NE,npRPEOverlap_ACh]',Group,'PlotStyle','compact');
set(gca,'box','off')
title('P(nRPE|pRPE)')
%ylim([0 1])

pPN = mediantest(pnRPEOverlap_ACh(~isnan(pnRPEOverlap_ACh)), pnRPEOverlap_NE(~isnan(pnRPEOverlap_NE)))
pNP = mediantest(npRPEOverlap_ACh(~isnan(npRPEOverlap_ACh)), npRPEOverlap_NE(~isnan(npRPEOverlap_NE)))

%% ACh
sigList_ACh = load(fullfile(save_path_mat_ACh,'sigNum.mat'));
sigList_NE = load(fullfile(save_path_mat_NE,'sigNum.mat'));
boxGroup = [ones(1,size(sigList_ACh.numSigVar,2)),2*ones(1,size(sigList_NE.numSigVar,2))];

figure;
for pp = 1:length(varList)
subplot(2,4,pp)
boxplot([sigList_ACh.numSigVar(pp,:)./sigList_ACh.numSigTotal,sigList_NE.numSigVar(pp,:)./sigList_NE.numSigTotal],boxGroup, ...
    'Labels',{'ACh','NE'},'PlotStyle','compact','LabelOrientation','horizontal', ...
    'Colors', [255, 189, 53;63,167,150]/255);
title(varList{pp})
set(gca,'box','off')
if pp==1 | pp ==5
    ylabel('Fraction of significant grids');
end
p = ranksum(sigList_ACh.numSigVar(pp,:)./sigList_ACh.numSigTotal,sigList_NE.numSigVar(pp,:)./sigList_NE.numSigTotal);
text(0,1.1,num2str(p))
ylim([0 1]);

end

savesumfigpath = fullfile(save_path_fluo_ACh,'ACh_NE_compare');
if ~exist(savesumfigpath)
    mkdir(savesumfigpath)
end
print(gcf,'-dpng',fullfile(savesumfigpath,'Percentage of significant grids var-comp'));
saveas(gcf, fullfile(savesumfigpath,'Percentage of significant grids var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'Percentage of significant grids var-comp'), 'svg');

%% compare the temporal dynamics, outcome
tempNE = fullfile( save_path_mat_NE,'outcome_groupstat.mat');
tempACh = fullfile( save_path_mat_ACh,'outcome_groupstat.mat');
NETemp_o = load(tempNE);
AChTemp_o = load(tempACh);


colors = [63,167,150;255, 189, 53]/255;
data1.NE_medmaxT = []; data1.NE_varmaxT = []; data1.NE_medmaxV = [];
data2.NE_medmaxT = []; data2.NE_varmaxT = []; data2.NE_medmaxV = [];
data1.ACh_medmaxT = []; data1.ACh_varmaxT = []; data1.ACh_medmaxV = [];
data2.ACh_medmaxT = []; data2.ACh_varmaxT = []; data2.ACh_medmaxV = [];

% get median and variance max time
for ss = 1:max(NETemp_o.o1SesInd)
    data1.NE_medmaxT = [data1.NE_medmaxT,nanmedian(NETemp_o.outcome_maxT_1(NETemp_o.o1SesInd==ss))];
     data1.NE_varmaxT = [data1.NE_varmaxT,nanvar(NETemp_o.outcome_maxT_1(NETemp_o.o1SesInd==ss))];
     data1.NE_medmaxV = [data1.NE_medmaxV,nanmedian(NETemp_o.outcome_maxV_1(NETemp_o.o1SesInd==ss))];
     data2.NE_medmaxT = [data2.NE_medmaxT,nanmedian(NETemp_o.outcome_maxT_2(NETemp_o.o2SesInd==ss))];
     data2.NE_varmaxT = [data2.NE_varmaxT,nanvar(NETemp_o.outcome_maxT_2(NETemp_o.o2SesInd==ss))];
     data2.NE_medmaxV = [data2.NE_medmaxV,nanmedian(NETemp_o.outcome_maxV_2(NETemp_o.o2SesInd==ss))];
end

for ss = 1:max(AChTemp_o.o1SesInd)
    data1.ACh_medmaxT = [data1.ACh_medmaxT,nanmedian(AChTemp_o.outcome_maxT_1(AChTemp_o.o1SesInd==ss))];
     data1.ACh_varmaxT = [data1.ACh_varmaxT,nanvar(AChTemp_o.outcome_maxT_1(AChTemp_o.o1SesInd==ss))];
     data1.ACh_medmaxV = [data1.ACh_medmaxV,nanmedian(AChTemp_o.outcome_maxV_1(AChTemp_o.o1SesInd==ss))];
     data2.ACh_medmaxT = [data2.ACh_medmaxT,nanmedian(AChTemp_o.outcome_maxT_2(AChTemp_o.o2SesInd==ss))];
     data2.ACh_varmaxT = [data2.ACh_varmaxT,nanvar(AChTemp_o.outcome_maxT_2(AChTemp_o.o2SesInd==ss))];
     data2.ACh_medmaxV = [data2.ACh_medmaxV,nanmedian(AChTemp_o.outcome_maxV_2(AChTemp_o.o2SesInd==ss))];
end

po1medT = ranksum(data1.ACh_medmaxT, data1.NE_medmaxT)
po1varT = ranksum(data1.ACh_varmaxT, data1.NE_varmaxT)
po1medV = ranksum(data1.ACh_medmaxV, data1.NE_medmaxV)
po2medT = ranksum(data2.ACh_medmaxT, data2.NE_medmaxT)
po2varT = ranksum(data2.ACh_varmaxT, data2.NE_varmaxT)
po2medV = ranksum(data2.ACh_medmaxV, data2.NE_medmaxV)

boxGroup = [ones(1,size(data1.NE_medmaxT,2)),2*ones(1,size(data1.ACh_medmaxT,2))];

figure;
subplot(1,3,1)
boxplot([data1.NE_medmaxT,data1.ACh_medmaxT],boxGroup,'PlotStyle','compact');
ylabel('med max time (s)');
set(gca, 'box','off')
ylim([0 5])
title('Group1');
subplot(1,3,2)
boxplot([data1.NE_varmaxT,data1.ACh_varmaxT],boxGroup,'PlotStyle','compact');
title('Group1')
ylim([0 5])
set(gca, 'box','off')
ylabel('var max time (s)');
subplot(1,3,3)
boxplot([-data1.NE_medmaxV,-data1.ACh_medmaxV],boxGroup,'PlotStyle','compact');
ylabel('max value');
set(gca, 'box','off')
ylim([-0.1 0.1])
title('Group1')
print(gcf,'-dpng',fullfile(savesumfigpath,'outcome group 1 var-comp'));
saveas(gcf, fullfile(savesumfigpath,'outcome group 1  var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'outcome group 1 var-comp'), 'svg');

boxGroup = [ones(1,size(data2.NE_medmaxT,2)),2*ones(1,size(data2.ACh_medmaxT,2))];

figure;
subplot(1,3,1)
boxplot([data2.NE_medmaxT,data2.ACh_medmaxT],boxGroup,'PlotStyle','compact');
ylabel('med max time (s)');
set(gca, 'box','off')
ylim([0 5])
title('Group2');
subplot(1,3,2)
boxplot([data2.NE_varmaxT,data2.ACh_varmaxT],boxGroup,'PlotStyle','compact');
title('Group2')
ylim([0 5])
set(gca, 'box','off')
ylabel('var max time (s)');
subplot(1,3,3)
boxplot([data2.NE_medmaxV,data2.ACh_medmaxV],boxGroup,'PlotStyle','compact');
ylabel('max value');
set(gca, 'box','off')
ylim([-0.1 0.1])
title('Group2')
print(gcf,'-dpng',fullfile(savesumfigpath,'outcome group 2 var-comp'));
saveas(gcf, fullfile(savesumfigpath,'outcome group 2  var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'outcome group 2 var-comp'), 'svg');

%% compare the temporal dynamics, choice
tempNE = fullfile( save_path_mat_NE,'choice_groupstat.mat');
tempACh = fullfile( save_path_mat_ACh,'choice_groupstat.mat');
NETemp_c = load(tempNE);
AChTemp_c = load(tempACh);
data1.NE_medmaxT = []; data1.NE_varmaxT = []; data1.NE_medmaxV = [];
data2.NE_medmaxT = []; data2.NE_varmaxT = []; data2.NE_medmaxV = [];
data1.ACh_medmaxT = []; data1.ACh_varmaxT = []; data1.ACh_medmaxV = [];
data2.ACh_medmaxT = []; data2.ACh_varmaxT = []; data2.ACh_medmaxV = [];

% get median and variance max time
for ss = 1:max(NETemp_c.c1SesInd)
    data1.NE_medmaxT = [data1.NE_medmaxT,nanmedian(NETemp_c.choice_maxT_1(NETemp_c.c1SesInd==ss))];
     data1.NE_varmaxT = [data1.NE_varmaxT,nanvar(NETemp_c.choice_maxT_1(NETemp_c.c1SesInd==ss))];
     data1.NE_medmaxV = [data1.NE_medmaxV,nanmedian(NETemp_c.choice_maxV_1(NETemp_c.c1SesInd==ss))];
     data2.NE_medmaxT = [data2.NE_medmaxT,nanmedian(NETemp_c.choice_maxT_2(NETemp_c.c2SesInd==ss))];
     data2.NE_varmaxT = [data2.NE_varmaxT,nanvar(NETemp_c.choice_maxT_2(NETemp_c.c2SesInd==ss))];
     data2.NE_medmaxV = [data2.NE_medmaxV,nanmedian(NETemp_c.choice_maxV_2(NETemp_c.c2SesInd==ss))];
end

for ss = 1:max(AChTemp_c.c1SesInd)
    data1.ACh_medmaxT = [data1.ACh_medmaxT,nanmedian(AChTemp_c.choice_maxT_1(AChTemp_c.c1SesInd==ss))];
     data1.ACh_varmaxT = [data1.ACh_varmaxT,nanvar(AChTemp_c.choice_maxT_1(AChTemp_c.c1SesInd==ss))];
     data1.ACh_medmaxV = [data1.ACh_medmaxV,nanmedian(AChTemp_c.choice_maxV_1(AChTemp_c.c1SesInd==ss))];
     data2.ACh_medmaxT = [data2.ACh_medmaxT,nanmedian(AChTemp_c.choice_maxT_2(AChTemp_c.c2SesInd==ss))];
     data2.ACh_varmaxT = [data2.ACh_varmaxT,nanvar(AChTemp_c.choice_maxT_2(AChTemp_c.c2SesInd==ss))];
     data2.ACh_medmaxV = [data2.ACh_medmaxV,nanmedian(AChTemp_c.choice_maxV_2(AChTemp_c.c2SesInd==ss))];
end

pc1medT = ranksum(data1.ACh_medmaxT, data1.NE_medmaxT)
pc1varT = ranksum(data1.ACh_varmaxT, data1.NE_varmaxT)
pc1medV = ranksum(-data1.ACh_medmaxV, data1.NE_medmaxV)
pc2medT = ranksum(data2.ACh_medmaxT, data2.NE_medmaxT)
pc2varT = ranksum(data2.ACh_varmaxT, data2.NE_varmaxT)
pc2medV = ranksum(data2.ACh_medmaxV, -data2.NE_medmaxV)

boxGroup = [ones(1,size(data1.NE_medmaxT,2)),2*ones(1,size(data1.ACh_medmaxT,2))];

figure;
subplot(1,3,1)
boxplot([data1.NE_medmaxT,data1.ACh_medmaxT],boxGroup,'PlotStyle','compact');
ylabel('med max time (s)');
set(gca, 'box','off')
ylim([0 5])
title('Group1');
subplot(1,3,2)
boxplot([data1.NE_varmaxT,data1.ACh_varmaxT],boxGroup,'PlotStyle','compact');
title('Group1')
ylim([0 5])
set(gca, 'box','off')
ylabel('var max time (s)');
subplot(1,3,3)
boxplot([data1.NE_medmaxV,-data1.ACh_medmaxV],boxGroup,'PlotStyle','compact');
ylabel('max value');
set(gca, 'box','off')
ylim([-0.1 0.1])
title('Group1')
print(gcf,'-dpng',fullfile(savesumfigpath,'choice group 1 var-comp'));
saveas(gcf, fullfile(savesumfigpath,'choice group 1  var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'choice group 1 var-comp'), 'svg');

boxGroup = [ones(1,size(data2.NE_medmaxT,2)),2*ones(1,size(data2.ACh_medmaxT,2))];


% group 2
figure
subplot(1,3,1)
boxplot([data2.NE_medmaxT,data2.ACh_medmaxT],boxGroup,'PlotStyle','compact');
ylabel('med max time (s)');
set(gca, 'box','off')
ylim([0 5])
title('Group2');
subplot(1,3,2)
boxplot([data2.NE_varmaxT,data2.ACh_varmaxT],boxGroup,'PlotStyle','compact');
title('Group2')
ylim([0 5])
set(gca, 'box','off')
ylabel('var max time (s)');
subplot(1,3,3)
boxplot([-data2.NE_medmaxV,data2.ACh_medmaxV],boxGroup,'PlotStyle','compact');
ylabel('max value');
set(gca, 'box','off')
ylim([-0.1 0.1])
title('Group2')
print(gcf,'-dpng',fullfile(savesumfigpath,'choice group 2 var-comp'));
saveas(gcf, fullfile(savesumfigpath,'choice group 2  var-comp'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'choice group 2 var-comp'), 'svg');



end
