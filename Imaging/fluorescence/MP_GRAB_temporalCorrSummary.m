function MP_GRAB_temporalCorrSummary(dataIndex, savesumfigpath,savematsumpath)

% calculate temporal correlation of different varible coefficients within
% grids

% check if coefficient sign changes over time
% for eachgrid, calculate the number of + time point and - time point
% then average

nFiles = size(dataIndex,1);

%% get clustering results, average similar clusters, also rerun clustering on all grids together
rCoeff = cell(0); rClust = cell(0);

% combine isSig list, for venn plot
SigList.c= [];SigList.r=[];SigList.x=[];
SigList.c_1 = []; SigList.c__1 = []; SigList.c__2 = [];
SigList.r_1 = []; SigList.r__1 = []; SigList.r__2 = [];
SigList.x_1 = []; SigList.x__1 = []; SigList.x__2 = [];
SigList.aR = []; SigList.cR = [];
SigList.dQ=[];SigList.dK=[];SigList.pRPE=[];SigList.nRPE=[];
SigList.CKE=[];
SigList.session=[];
% count the number of significant grids for every variables and sessions
numSigVar = zeros(8,nFiles);
numSigTotal = zeros(1,nFiles);
varList = {'choice', 'outcome', 'interaction', 'dQ', 'dK', 'posRPE', 'negRPE', 'CKE'};
for ii = 1:nFiles

    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    %load(fullfile(fn_beh.folder,fn_beh.name));


    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    if length(fn_fluo) == 1

        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end

        cd(savefluofigpath);


        %load(fullfile(fn_fluo.folder,fn_fluo.name));

        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end

        % LOAD SAVED DATA
        saveregpath = fullfile(savematpath,'cluster.mat');
         saveregpath2 = fullfile(savematpath,'regressionVarTempCorr.mat');  % data with c/r/x from n+1 to n-2

        if exist(saveregpath)
            load(saveregpath)
        end

        if exist(saveregpath2)
            load(saveregpath2)
        end
%         % outcome
%         tlabel1 = 'Outcome coefficient';
%         saveregpath =  fullfile(savematpath,[tlabel1,' cluster.mat']);  % file path to save the results
%         if exist(saveregpath)
%             load(saveregpath)
%         end

        rCoeff{ii} = saveDataOutcome.coeff(saveDataOutcome.oriInd,:);
        rClust{ii} = saveDataOutcome.T;
        rOriID{ii} = saveDataOutcome.oriInd;
        numSigVar(2,ii) = length(saveDataOutcome.oriInd);
        rt = saveDataOutcome.t;

        choiceCoeff{ii} = choicetempData.coeff;choiceisSig{ii} = choicetempData.isSig; numSigVar(1,ii) = sum(choicetempData.isSig);
        xnCoeff{ii} = xntempData.coeff;xnisSig{ii} = xntempData.isSig;numSigVar(3,ii) = sum(xntempData.isSig);
        posRPECoeff{ii} = posRPEtempData.coeff; posRPEisSig{ii} = posRPEtempData.isSig;numSigVar(6,ii) = sum(posRPEtempData.isSig);
        negRPECoeff{ii} = negRPEtempData.coeff; negRPEisSig{ii} = negRPEtempData.isSig;numSigVar(7,ii) = sum(posRPEtempData.isSig);
        %negRPECoeff{ii} = negRPEtempData.coeff; negRPEisSig{ii} = negRPEtempData.isSig;
        dQCoeff{ii} = dQtempData.coeff; dQisSig{ii} = dQtempData.isSig;numSigVar(4,ii) = sum(dQtempData.isSig);
        dKCoeff{ii} = dKtempData.coeff; dKisSig{ii} = dKtempData.isSig;numSigVar(5,ii) = sum(dKtempData.isSig);
        CKECoeff{ii} = CKEtempData.coeff; CKEisSig{ii} = CKEtempData.isSig;numSigVar(8,ii) = sum(CKEtempData.isSig);

        % sig list for venn plot
        SigList.session = [SigList.session,ii*ones(1,length(choicetempData.isSig))];
        SigList.c= [SigList.c,choicetempData.isSig'];
        outcomeSig = zeros(1,size(saveDataOutcome.coeff,1));
        outcomeSig(saveDataOutcome.sigInd)=1;
        SigList.r = [SigList.r, outcomeSig];
        SigList.x = [SigList.x, xntempData.isSig'];
        SigList.dQ = [SigList.dQ, dQtempData.isSig'];
        SigList.dK = [SigList.dK, dKtempData.isSig'];
        SigList.pRPE = [SigList.pRPE, posRPEtempData.isSig'];
        SigList.nRPE = [SigList.nRPE, negRPEtempData.isSig'];
        SigList.CKE = [SigList.CKE, CKEtempData.isSig'];

        % previous/next c/r/x
        SigList.c_1 = [SigList.c_1,cn_1tempData.isSig']; SigList.c__1 = [SigList.c__1,cn__1tempData.isSig']; SigList.c__2 = [SigList.c__2,cn__2tempData.isSig'];
        SigList.r_1 = [SigList.r_1,rn_1tempData.isSig']; SigList.r__1 = [SigList.r__1,rn__1tempData.isSig']; SigList.r__2 = [SigList.r__2,rn__2tempData.isSig'];
        SigList.x_1 = [SigList.x_1,xn_1tempData.isSig']; SigList.x__1 = [SigList.x__1,xn__1tempData.isSig']; SigList.x__2 = [SigList.x__2,xn__2tempData.isSig'];
        SigList.aR = [SigList.aR,ave_rtempData.isSig'];
        SigList.cR = [SigList.cR,cum_rtempData.isSig'];
        %numSigTotal(ii) = sum(choicetempData.isSig|outcomeSig'|xntempData.isSig|cn_1);
    end
end

% save the siglist
if ~exist(savematsumpath)
    mkdir(savematsumpath);
end
save(fullfile(savematsumpath,'sigList.mat'),'SigList');
save(fullfile(savematsumpath,'sigNum.mat'),'numSigVar','numSigTotal');

%% plot percentage of significant grids
sigForAny = SigList.c|SigList.r|SigList.x|SigList.c_1|SigList.c__1|SigList.c__2|SigList.r_1|SigList.r__1|SigList.r__2|SigList.x_1|SigList.x__1|SigList.x__2|SigList.aR|SigList.cR;
numSigAny = sum(sigForAny);
figure;
explode = [0 1];
p=pie([numSigAny,length(sigForAny)-numSigAny],explode);
print(gcf,'-dpng',fullfile(savesumfigpath,'Percentage of significant grids All'));
saveas(gcf, fullfile(savesumfigpath,'Percentage of significant grids All'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'Percentage of significant grids All'), 'svg');

% violin plot for significant grids with different variables
figure;violinplot(numSigVar'./size(choiceCoeff{1},1),varList,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0]);
ylabel('Fraction of significant grids');
print(gcf,'-dpng',fullfile(savesumfigpath,'Fraction of sig variables'));
saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables'), 'fig');
saveas(gcf, fullfile(savesumfigpath,'Fraction of sig variables'), 'svg');

%

%% venn diagram showing overlap

allInd = 1:length(SigList.c);
cSigInd = allInd(logical(SigList.c)); rSigInd = allInd(logical(SigList.r)); xSigInd = allInd(logical(SigList.x));
dQSigInd = allInd(logical(SigList.dQ)); dKSigInd = allInd(logical(SigList.dK));
pRPESigInd = allInd(logical(SigList.pRPE)); nRPESigInd = allInd(logical(SigList.nRPE));
CKESigInd = allInd(logical(SigList.CKE));

%% test for overlapping significance

% chi-square test for independence

% choice and outcome
x = [length(intersect(cSigInd,rSigInd)),length(setdiff(cSigInd,rSigInd));
    length(setdiff(rSigInd,cSigInd)),length(allInd)-length(union(cSigInd,rSigInd))];
[h,chi,p]=chi2ind(x,0.05)
y=[length(intersect(cSigInd,xSigInd)),length(setdiff(cSigInd,xSigInd));
    length(setdiff(xSigInd,cSigInd)),length(allInd)-length(union(cSigInd,xSigInd))];
[h,chi,p]=chi2ind(y,0.05)
z=[length(intersect(rSigInd,xSigInd)),length(setdiff(rSigInd,xSigInd));
    length(setdiff(xSigInd,rSigInd)),length(allInd)-length(union(rSigInd,xSigInd))];
[h,chi,p]=chi2ind(z,0.05)
close all;
%% trim the clustering.
% 1. delete cluster with too few grids . <10 grids?
% 2. combine goup with high correlation >0.9?
g1=cell(0); g2 = cell(0); g3 = cell(0);
clustCoeff = []; clustMean = [];
rClust_new = cell(0);
rCoeff_new = cell(0);
rSes = []; % session number,outcome
cCoeff_new = cell(0);
xCoeff_new = cell(0);
posRPECoeff_new = cell(0);
negRPECoeff_new = cell(0);
dQCoeff_new = cell(0);
dKCoeff_new = cell(0);
CKECoeff_new = cell(0);

% only take significant cells
% cCoeff_Sig = cell(0);
% xCoeff_Sig = cell(0);
% posRPECoeff_Sig = cell(0);
% negRPECoeff_Sig = cell(0);

% sig ROIs that not sig for outcome
cCoeff_nR = [];
xCoeff_nR = [];xnR_sub=[];
posRPECoeff_nR = [];
negRPECoeff_nR = [];
dQCoeff_nR = [];
dKCoeff_nR = [];
CKECoeff_nR = [];

rOriID_new = cell(0);
newInd = 1;
animalList = unique(dataIndex.Animal);
% rHPeakT = [];
% rLPeakT = [];
% rDecayTau = [];
oriGroup = [];
% save all grids that cannot be grouped into the 3 groups



t0 = 0; t1 = 4;  % time range for correlation calculation
ROIId = 1:size(choiceCoeff{1},1);
% for sig grids counting for group4 (not sig for outcome)
numSigVarg4 = zeros(8,nFiles);

for gg = 1:nFiles
    g1Num = sum(rClust{gg}==1); g2Num = sum(rClust{gg}==2); g3Num = sum(rClust{gg}==3);
    % calculate the average coefficient with in 1-2.5s for each group
    g1Mean =nanmean(rCoeff{gg}(rClust{gg}==1,:),1);  g2Mean = nanmean(rCoeff{gg}(rClust{gg}==2,:),1); g3Mean =nanmean(rCoeff{gg}(rClust{gg}==3,:),1);

    %clustMean = [clustMean;g1Mean;g2Mean;g3Mean];
    %% smooth the mean and find peak and decay tau
    % find 1): highest peak and time;
    % 2): lowest peak and time;
    % 3) exponential decay tau after high peak



    % (islocalmin & islocalmax)
    tInd = 1:length(rt);
    searchInd = tInd(rt>0);
    % group 1
%     smoothed1 = smooth(g1Mean);
%     [minValue1,minInd] = nanmin(smoothed1(rt>0));
%     minInd1 = searchInd(minInd);
%     [maxValue1,maxInd] = nanmax(smoothed1(rt>0));
%     maxInd1 = searchInd(maxInd);
%     rHPeakT = [rHPeakT, rt(maxInd1)];
%     rLPeakT = [rLPeakT, rt(minInd1)];
    oriGroup = [oriGroup,1];
    % fit for tau
    %     if maxInd1<length(rt)
    %         f1 = fit(rt(maxInd1:end),smoothed1(maxInd1:end),'exp1');
    %         rDecayTau = [rDecayTau, f1.b];
    %     else
    %         f1 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==1,:);
    rSes = [rSes,gg];
    cCoeff_new{newInd} = choiceCoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    xCoeff_new{newInd} = xnCoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    posRPECoeff_new{newInd} = posRPECoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    negRPECoeff_new{newInd} = negRPECoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    dQCoeff_new{newInd} = dQCoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    dKCoeff_new{newInd} = dKCoeff{gg}(rOriID{gg}(rClust{gg}==1),:);
    CKECoeff_new{newInd} = CKECoeff{gg}(rOriID{gg}(rClust{gg}==1),:);

    group1Sig = zeros(size(choiceCoeff{1},1),1);
    group1SigId = rOriID{gg}(rClust{gg}==1);
    choiceSigId = ROIId(choiceisSig{gg});
    xnSigId = ROIId(xnisSig{gg});
    posRPESigId = ROIId(posRPEisSig{gg});
    negRPESigId = ROIId(negRPEisSig{gg});
    dQSigId = ROIId(dQisSig{gg});
    dKSigId = ROIId(dKisSig{gg});
    CKESigId = ROIId(CKEisSig{gg});

    choiceSig = zeros(1,length(group1SigId));
    xnSig = zeros(1,length(group1SigId));
    posRPESig = zeros(1,length(group1SigId));
    negRPESig = zeros(1,length(group1SigId));
    dQSig = zeros(1,length(group1SigId));
    dKSig = zeros(1,length(group1SigId));
    CKESig = zeros(1,length(group1SigId));

    for rr = 1:length(group1SigId)
        if ismember(group1SigId(rr),choiceSigId)
            choiceSig(rr) = 1;
        end
        if ismember(group1SigId(rr),xnSigId)
            xnSig(rr) = 1;
        end
        if ismember(group1SigId(rr),posRPESigId)
            posRPESig(rr) = 1;
        end
        if ismember(group1SigId(rr),negRPESigId)
            negRPESig(rr) = 1;
        end
        if ismember(group1SigId(rr),dQSigId)
            dQSig(rr) = 1;
        end
        if ismember(group1SigId(rr),dKSigId)
            dKSig(rr) = 1;
        end
        if ismember(group1SigId(rr),CKESigId)
            CKESig(rr) = 1;
        end
    end
    cisSig{newInd} = choiceSig;
    xisSig{newInd} = xnSig;
    pRPEisSig{newInd} = posRPESig;
    nRPEisSig{newInd} = negRPESig;
    dqisSig{newInd} = dQSig;
    dkisSig{newInd} = dKSig;
    ckeisSig{newInd} = CKESig;

%     cCoeff_Sig{newInd} = choiceCoeff{gg}(group1Sig & choiceisSig{gg},:);
%     xCoeff_Sig{newInd} = xnCoeff{gg}(group1Sig & xnisSig{gg},:);
%     posRPECoeff_Sig{newInd} = posRPECoeff{gg}(group1Sig & posRPEisSig{gg},:);
%     negRPECoeff_Sig{newInd} = negRPECoeff{gg}(group1Sig & negRPEisSig{gg},:);
%
    RisSig = zeros(1,size(choiceCoeff{1},1));
    RisSig(rOriID{gg}) = 1;

    cCoeff_nR = [cCoeff_nR; choiceCoeff{gg}(~RisSig'&choiceisSig{gg},:)];numSigVarg4(1,gg) = sum(~RisSig'&choiceisSig{gg});
    xCoeff_nR = [xCoeff_nR; xnCoeff{gg}(~RisSig'&xnisSig{gg},:)];numSigVarg4(3,gg) = sum(~RisSig'&xnisSig{gg});
    xnR_sub = [xnR_sub; ones(sum(~RisSig'&xnisSig{gg}),1)*find(strcmp(animalList,dataIndex.Animal{gg}))];
    posRPECoeff_nR = [posRPECoeff_nR; posRPECoeff{gg}(~RisSig'&posRPEisSig{gg},:)];numSigVarg4(6,gg) = sum(~RisSig'&posRPEisSig{gg});
    negRPECoeff_nR = [negRPECoeff_nR; negRPECoeff{gg}(~RisSig'&negRPEisSig{gg},:)];numSigVarg4(7,gg) = sum(~RisSig'&negRPEisSig{gg});
    dQCoeff_nR = [dQCoeff_nR; dQCoeff{gg}(~RisSig'&dQisSig{gg},:)];numSigVarg4(4,gg) = sum(~RisSig'&dQisSig{gg});
    dKCoeff_nR = [dKCoeff_nR; dKCoeff{gg}(~RisSig'&dKisSig{gg},:)];numSigVarg4(5,gg) = sum(~RisSig'&dKisSig{gg});
    CKECoeff_nR = [CKECoeff_nR; CKECoeff{gg}(~RisSig'&CKEisSig{gg},:)];numSigVarg4(8,gg) = sum(~RisSig'&CKEisSig{gg});

    newInd = newInd+1;

    % group 2
%     smoothed2 = smooth(g2Mean);
%     [minValue2,minInd] = nanmin(smoothed2(rt>0));
%     minInd2 = searchInd(minInd);
%     [maxValue2,maxInd] = nanmax(smoothed2(rt>0));
%     maxInd2 = searchInd(maxInd);
%     rHPeakT = [rHPeakT, rt(maxInd2)];
%     rLPeakT = [rLPeakT, rt(minInd2)];
    oriGroup = [oriGroup,2];
    % fit for tau
    %     if maxInd2<length(rt)
    %         f2 = fit(rt(maxInd2:end),smoothed2(maxInd2:end),'exp1');
    %         rDecayTau = [rDecayTau, f2.b];
    %     else
    %         f2 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==2,:);
    rSes = [rSes,gg];
    cCoeff_new{newInd} = choiceCoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    xCoeff_new{newInd} = xnCoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    posRPECoeff_new{newInd} = posRPECoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    negRPECoeff_new{newInd} = negRPECoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    dQCoeff_new{newInd} = dQCoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    dKCoeff_new{newInd} = dKCoeff{gg}(rOriID{gg}(rClust{gg}==2),:);
    CKECoeff_new{newInd} = CKECoeff{gg}(rOriID{gg}(rClust{gg}==2),:);

    group2Sig = zeros(size(choiceCoeff{1},1),1);
    group2Sig(rOriID{gg}(rClust{gg}==2)) = 1;
    group2SigId = rOriID{gg}(rClust{gg}==2);

    choiceSig = zeros(1,length(group2SigId));
    xnSig = zeros(1,length(group2SigId));
    posRPESig = zeros(1,length(group2SigId));
    negRPESig = zeros(1,length(group2SigId));
    dQSig = zeros(1,length(group2SigId));
    dKSig = zeros(1,length(group2SigId));
    CKESig = zeros(1,length(group2SigId));

    for rr = 1:length(group2SigId)
        if ismember(group2SigId(rr),choiceSigId)
            choiceSig(rr) = 1;
        end
        if ismember(group2SigId(rr),xnSigId)
            xnSig(rr) = 1;
        end
        if ismember(group2SigId(rr),posRPESigId)
            posRPESig(rr) = 1;
        end
        if ismember(group2SigId(rr),negRPESigId)
            negRPESig(rr) = 1;
        end
        if ismember(group2SigId(rr),dQSigId)
            dQSig(rr) = 1;
        end
        if ismember(group2SigId(rr),dKSigId)
            dKSig(rr) = 1;
        end
        if ismember(group2SigId(rr),CKESigId)
            CKESig(rr) = 1;
        end
    end
    cisSig{newInd} = choiceSig;
    xisSig{newInd} = xnSig;
    pRPEisSig{newInd} = posRPESig;
    nRPEisSig{newInd} = negRPESig;
    dqisSig{newInd} = dQSig;
    dkisSig{newInd} = dKSig;
    ckeisSig{newInd} = CKESig;
%      cCoeff_Sig{newInd} = choiceCoeff{gg}(group2Sig & choiceisSig{gg},:);
%      xCoeff_Sig{newInd} = xnCoeff{gg}(group2Sig & xnisSig{gg},:);
%     posRPECoeff_Sig{newInd} = posRPECoeff{gg}(group2Sig & posRPEisSig{gg},:);
%     negRPECoeff_Sig{newInd} = negRPECoeff{gg}(group2Sig & negRPEisSig{gg},:);
%
    newInd = newInd+1;

    % group 3
%     smoothed3 = smooth(g3Mean);
%     [minValue3,minInd] = nanmin(smoothed3(rt>0));
%     minInd3 = searchInd(minInd);
%     [maxValue3,maxInd] = nanmax(smoothed3(rt>0));
%     maxInd3 = searchInd(maxInd);
%     rHPeakT = [rHPeakT, rt(maxInd3)];
%     rLPeakT = [rLPeakT, rt(minInd3)];
    oriGroup = [oriGroup,3];
    % fit for tau
    %     if maxInd3<length(rt)
    %         f3 = fit(rt(maxInd3:end),smoothed3(maxInd3:end),'exp1');
    %         rDecayTau = [rDecayTau, f3.b];
    %     else
    %         f3 = [];
    %         rDecayTau = [rDecayTau,NaN];
    %     end
    rCoeff_new{newInd} = rCoeff{gg}(rClust{gg}==3,:);
     rSes = [rSes,gg];
     cCoeff_new{newInd} = choiceCoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    xCoeff_new{newInd} = xnCoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    posRPECoeff_new{newInd} = posRPECoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    negRPECoeff_new{newInd} = negRPECoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    dQCoeff_new{newInd} = dQCoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    dKCoeff_new{newInd} = dKCoeff{gg}(rOriID{gg}(rClust{gg}==3),:);
    CKECoeff_new{newInd} = CKECoeff{gg}(rOriID{gg}(rClust{gg}==3),:);

     group3Sig = zeros(size(choiceCoeff{1},1),1);
    group3Sig(rOriID{gg}(rClust{gg}==3)) = 1;
     group3SigId = rOriID{gg}(rClust{gg}==3);

    choiceSig = zeros(1,length(group3SigId));
    xnSig = zeros(1,length(group3SigId));
    posRPESig = zeros(1,length(group3SigId));
    negRPESig = zeros(1,length(group3SigId));
    dQSig = zeros(1,length(group3SigId));
    dKSig = zeros(1,length(group3SigId));
    CKESig = zeros(1,length(group3SigId));

    for rr = 1:length(group3SigId)
        if ismember(group3SigId(rr),choiceSigId)
            choiceSig(rr) = 1;
        end
        if ismember(group3SigId(rr),xnSigId)
            xnSig(rr) = 1;
        end
        if ismember(group3SigId(rr),posRPESigId)
            posRPESig(rr) = 1;
        end
        if ismember(group3SigId(rr),negRPESigId)
            negRPESig(rr) = 1;
        end
        if ismember(group3SigId(rr),dQSigId)
            dQSig(rr) = 1;
        end
        if ismember(group3SigId(rr),dKSigId)
            dKSig(rr) = 1;
        end
        if ismember(group3SigId(rr),CKESigId)
            CKESig(rr) = 1;
        end
    end
    cisSig{newInd} = choiceSig;
    xisSig{newInd} = xnSig;
    pRPEisSig{newInd} = posRPESig;
    nRPEisSig{newInd} = negRPESig;
    dqisSig{newInd} = dQSig;
    dkisSig{newInd} = dKSig;
    ckeisSig{newInd} = CKESig;

%      cCoeff_Sig{newInd} = choiceCoeff{gg}(group3Sig & choiceisSig{gg},:);
%      xCoeff_Sig{newInd} = xnCoeff{gg}(group3Sig & xnisSig{gg},:);
%      posRPECoeff_Sig{newInd} = posRPECoeff{gg}(group3Sig & posRPEisSig{gg},:);
%      negRPECoeff_Sig{newInd} = negRPECoeff{gg}(group3Sig & negRPEisSig{gg},:);
%
    newInd = newInd+1;

    % plot the fit results
    savefigpath =  fullfile(dataIndex.BehPath{gg},'figs-fluo','clustProp');
    if ~exist(savefigpath)
        mkdir(savefigpath)
    end


end

% regroup with criteria *area under curve
group1 = []; group2 = []; group3 = []; %group4 = [];
o1SesInd = []; o2SesInd = []; o3SesInd = []; % keep track of the session number
choice1 = []; choice2 = []; choice3 = [];
c1SesInd = []; c2SesInd = []; c3SesInd = [];

cSig1 = []; cSig2 = []; cSig3 = [];
%choice1Sig = [];choice2Sig = [];choice3Sig = [];
xn1 = []; xn2 = []; xn3 = [];
xSig1 = []; xSig2 = []; xSig3 = [];
%xn1Sig = []; xn2Sig = []; xn3Sig = [];
posRPE1 = []; posRPE2 = []; posRPE3 = [];
pRPESig1 = []; pRPESig2 = []; pRPESig3 = [];
%posRPE1Sig = []; posRPE2Sig = []; posRPE3Sig = [];
negRPE1 = []; negRPE2 = []; negRPE3 = [];
nRPESig1 = []; nRPESig2 = []; nRPESig3 = [];
%negRPE1Sig = []; negRPE2Sig = []; negRPE3Sig = [];
dQ1 = []; dQ2 = []; dQ3 = [];
dQSig1 = []; dQSig2 = []; dQSig3 = [];

dK1 = []; dK2 = []; dK3 = [];
dKSig1 = []; dKSig2 = []; dKSig3 = [];

CKE1 = []; CKE2 = []; CKE3 = [];
CKESig1 = []; CKESig2 = []; CKESig3 = [];
% if the surface area between (0,2) less than 0, group1
% if the surface area between (2,4) less than 0, group4

%
numSigVarg1 = zeros(8,nFiles);
numSigVarg2 = zeros(8,nFiles);
numSigVarg3 = zeros(8,nFiles);

for ggg = 1:length(rCoeff_new)
    animalInd = floor((ggg-1)/3)+1;
    if oriGroup(ggg) == 1
        line11 = nanmean(rCoeff_new{ggg},1);
        areaUnder1 = trapz(rt(rt>0&rt<2),line11(rt>0 & rt<2));
       %areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
        if areaUnder1 < 0 %& areaUnder2 > areaUnder1
            group1 = [group1;rCoeff_new{ggg}];
            o1SesInd = [o1SesInd, ones(1,size(rCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg1(2,animalInd) = numSigVarg1(2,animalInd) + size(rCoeff_new{ggg},1);
            choice1 = [choice1; cCoeff_new{ggg}]; cSig1 = [cSig1,cisSig{ggg}];%choice1Sig = [choice1Sig; cCoeff_Sig{ggg}];
             c1SesInd = [c1SesInd, ones(1,size(cCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg1(1,animalInd) = numSigVarg1(1,animalInd) + sum(cisSig{ggg});
            xn1 = [xn1; xCoeff_new{ggg}]; xSig1 = [xSig1,xisSig{ggg}];%xn1Sig = [xn1Sig; xCoeff_Sig{ggg}];
            numSigVarg1(3,animalInd) = numSigVarg1(3,animalInd) + sum(xisSig{ggg});
            posRPE1 = [posRPE1; posRPECoeff_new{ggg}]; pRPESig1 = [pRPESig1,pRPEisSig{ggg}]; %posRPE1Sig = [posRPE1Sig; posRPECoeff_Sig{ggg}];
            numSigVarg1(6,animalInd) = numSigVarg1(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE1 = [negRPE1; negRPECoeff_new{ggg}]; nRPESig1 = [nRPESig1,nRPEisSig{ggg}];%negRPE1Sig = [negRPE1Sig; negRPECoeff_Sig{ggg}];
            numSigVarg1(7,animalInd) = numSigVarg1(7,animalInd) + sum(nRPEisSig{ggg});
            dQ1 = [dQ1; dQCoeff_new{ggg}]; dQSig1 = [dQSig1,dqisSig{ggg}];
            numSigVarg1(4,animalInd) = numSigVarg1(4,animalInd) + sum(dqisSig{ggg});
            dK1 = [dK1; dKCoeff_new{ggg}]; dKSig1 = [dKSig1,dkisSig{ggg}];
            numSigVarg1(5,animalInd) = numSigVarg1(5,animalInd) + sum(dkisSig{ggg});
            CKE1 = [CKE1; CKECoeff_new{ggg}]; CKESig1 = [CKESig1,ckeisSig{ggg}];
            numSigVarg1(8,animalInd) = numSigVarg1(8,animalInd) + sum(ckeisSig{ggg});

        else
            group2 = [group2;rCoeff_new{ggg}];
            o2SesInd = [o2SesInd, ones(1,size(rCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg2(2,animalInd) = numSigVarg2(2,animalInd) + size(rCoeff_new{ggg},1);
            choice2 = [choice2; cCoeff_new{ggg}]; cSig2 = [cSig2,cisSig{ggg}];%choice2Sig = [choice2Sig; cCoeff_Sig{ggg}];
             c2SesInd = [c2SesInd, ones(1,size(cCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg2(1,animalInd) = numSigVarg2(1,animalInd) + sum(cisSig{ggg});
            xn2 = [xn2; xCoeff_new{ggg}]; xSig2 = [xSig2,xisSig{ggg}];%xn2Sig = [xn2Sig; xCoeff_Sig{ggg}];
            numSigVarg2(3,animalInd) = numSigVarg2(3,animalInd) + sum(xisSig{ggg});
            posRPE2 = [posRPE2; posRPECoeff_new{ggg}]; pRPESig2 = [pRPESig2,pRPEisSig{ggg}]; %posRPE2Sig = [posRPE2Sig; posRPECoeff_Sig{ggg}];
            numSigVarg2(6,animalInd) = numSigVarg2(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE2 = [negRPE2; negRPECoeff_new{ggg}]; nRPESig2 = [nRPESig2,nRPEisSig{ggg}];%negRPE2Sig = [negRPE2Sig; negRPECoeff_Sig{ggg}];
            numSigVarg2(7,animalInd) = numSigVarg2(7,animalInd) + sum(nRPEisSig{ggg});
            dQ2 = [dQ2; dQCoeff_new{ggg}]; dQSig2 = [dQSig2,dqisSig{ggg}];
            numSigVarg2(4,animalInd) = numSigVarg2(4,animalInd) + sum(dqisSig{ggg});
             dK2 = [dK2; dKCoeff_new{ggg}]; dKSig2 = [dKSig2,dkisSig{ggg}];
             numSigVarg2(5,animalInd) = numSigVarg2(5,animalInd) + sum(dkisSig{ggg});
            CKE2 = [CKE2; CKECoeff_new{ggg}]; CKESig2 = [CKESig2,ckeisSig{ggg}];
           numSigVarg2(8,animalInd) = numSigVarg2(8,animalInd) + sum(ckeisSig{ggg});
        end
    elseif oriGroup(ggg) == 2
        line11 = nanmean(rCoeff_new{ggg},1);
        areaUnder1 = trapz(rt(rt>0&rt<2),line11(rt>0 & rt<2));
        %areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
        if areaUnder1 < 0 %& areaUnder2 > areaUnder1
            group1 = [group1;rCoeff_new{ggg}];
            o1SesInd = [o1SesInd, ones(1,size(rCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg1(2,animalInd) = numSigVarg1(2,animalInd) + size(rCoeff_new{ggg},1);
%         elseif areaUnder2 < 0
%             group4 = [group4;rCoeff_new{ggg}];
            choice1 = [choice1; cCoeff_new{ggg}]; cSig1 = [cSig1,cisSig{ggg}];%choice1Sig = [choice1Sig; cCoeff_Sig{ggg}];
             c1SesInd = [c1SesInd, ones(1,size(cCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg1(1,animalInd) = numSigVarg1(1,animalInd) + sum(cisSig{ggg});
            xn1 = [xn1; xCoeff_new{ggg}]; xSig1 = [xSig1,xisSig{ggg}];%xn1Sig = [xn1Sig; xCoeff_Sig{ggg}];
            numSigVarg1(3,animalInd) = numSigVarg1(3,animalInd) + sum(xisSig{ggg});
            posRPE1 = [posRPE1; posRPECoeff_new{ggg}]; pRPESig1 = [pRPESig1,pRPEisSig{ggg}];%posRPE1Sig = [posRPE1Sig; posRPECoeff_Sig{ggg}];
            numSigVarg1(6,animalInd) = numSigVarg1(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE1 = [negRPE1; negRPECoeff_new{ggg}]; nRPESig1 = [nRPESig1,nRPEisSig{ggg}];%negRPE2Sig = [negRPE2Sig; negRPECoeff_Sig{ggg}];
            numSigVarg1(7,animalInd) = numSigVarg1(7,animalInd) + sum(nRPEisSig{ggg});
            dQ1 = [dQ1; dQCoeff_new{ggg}]; dQSig1 = [dQSig1,dqisSig{ggg}];
             numSigVarg1(4,animalInd) = numSigVarg1(4,animalInd) + sum(dqisSig{ggg});
            dK1 = [dK1; dKCoeff_new{ggg}]; dKSig1 = [dKSig1,dkisSig{ggg}];
            numSigVarg1(5,animalInd) = numSigVarg1(5,animalInd) + sum(dkisSig{ggg});
            CKE1 = [CKE1; CKECoeff_new{ggg}]; CKESig1 = [CKESig1,ckeisSig{ggg}];
           numSigVarg1(8,animalInd) = numSigVarg1(8,animalInd) + sum(ckeisSig{ggg});

        else
            group2 = [group2;rCoeff_new{ggg}];
            o2SesInd = [o2SesInd, ones(1,size(rCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg2(2,animalInd) = numSigVarg2(2,animalInd) + size(rCoeff_new{ggg},1);
            choice2 = [choice2; cCoeff_new{ggg}]; cSig2 = [cSig2,cisSig{ggg}];%choice2Sig = [choice2Sig; cCoeff_Sig{ggg}];
             c2SesInd = [c2SesInd, ones(1,size(cCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg2(1,animalInd) = numSigVarg2(1,animalInd) + sum(cisSig{ggg});
            xn2 = [xn2; xCoeff_new{ggg}]; xSig2 = [xSig2,xisSig{ggg}];%xn2Sig = [xn2Sig; xCoeff_Sig{ggg}];
            numSigVarg2(3,animalInd) = numSigVarg2(3,animalInd) + sum(xisSig{ggg});
            posRPE2 = [posRPE2; posRPECoeff_new{ggg}]; pRPESig2 = [pRPESig2,pRPEisSig{ggg}];%posRPE2Sig = [posRPE2Sig; posRPECoeff_Sig{ggg}];
            numSigVarg2(6,animalInd) = numSigVarg2(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE2 = [negRPE2; negRPECoeff_new{ggg}]; nRPESig2 = [nRPESig2,nRPEisSig{ggg}];%negRPE2Sig = [negRPE2Sig; negRPECoeff_Sig{ggg}];
            numSigVarg2(7,animalInd) = numSigVarg2(7,animalInd) + sum(nRPEisSig{ggg});
            dQ2 = [dQ2; dQCoeff_new{ggg}]; dQSig2 = [dQSig2,dqisSig{ggg}];
            numSigVarg2(4,animalInd) = numSigVarg2(4,animalInd) + sum(dqisSig{ggg});
            dK2 = [dK2; dKCoeff_new{ggg}]; dKSig2 = [dKSig2,dkisSig{ggg}];
             numSigVarg2(5,animalInd) = numSigVarg2(5,animalInd) + sum(dkisSig{ggg});
            CKE2 = [CKE2; CKECoeff_new{ggg}]; CKESig2 = [CKESig2,ckeisSig{ggg}];
              numSigVarg2(8,animalInd) = numSigVarg2(8,animalInd) + sum(ckeisSig{ggg});
        end
    else
%         areaUnder2 = trapz(rt(rt>3&rt<4),line11(rt>3 & rt<4));
%         if areaUnder2 < 0
%             group4 = [group4;rCoeff_new{ggg}];
%         else
            group3 = [group3;rCoeff_new{ggg}];
            o3SesInd = [o3SesInd, ones(1,size(rCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg3(2,animalInd) = numSigVarg3(2,animalInd) + size(rCoeff_new{ggg},1);
            choice3 = [choice3; cCoeff_new{ggg}]; cSig3 = [cSig3,cisSig{ggg}];
             c3SesInd = [c3SesInd, ones(1,size(cCoeff_new{ggg},1))*rSes(ggg)];
            numSigVarg3(1,animalInd) = numSigVarg3(1,animalInd) + sum(cisSig{ggg});%choice3Sig = [choice3Sig; cCoeff_Sig{ggg}];
            xn3 = [xn3; xCoeff_new{ggg}]; xSig3 = [xSig3,xisSig{ggg}];%xn3Sig = [xn3Sig; xCoeff_Sig{ggg}];
            numSigVarg3(3,animalInd) = numSigVarg3(3,animalInd) + sum(xisSig{ggg});
            posRPE3 = [posRPE3; posRPECoeff_new{ggg}]; pRPESig3 = [pRPESig3,pRPEisSig{ggg}];%posRPE3Sig = [posRPE3Sig; posRPECoeff_Sig{ggg}];
            numSigVarg3(6,animalInd) = numSigVarg3(6,animalInd) + sum(pRPEisSig{ggg});
            negRPE3 = [negRPE3; negRPECoeff_new{ggg}]; nRPESig3 = [nRPESig3,nRPEisSig{ggg}];%negRPE3Sig = [negRPE3Sig; negRPECoeff_Sig{ggg}];
            numSigVarg3(7,animalInd) = numSigVarg3(7,animalInd) + sum(nRPEisSig{ggg});
            dQ3 = [dQ3; dQCoeff_new{ggg}]; dQSig3 = [dQSig3,dqisSig{ggg}];
            numSigVarg3(4,animalInd) = numSigVarg3(4,animalInd) + sum(dqisSig{ggg});
            dK3 = [dK3; dKCoeff_new{ggg}]; dKSig3 = [dKSig3,dkisSig{ggg}];
            numSigVarg3(5,animalInd) = numSigVarg3(5,animalInd) + sum(dkisSig{ggg});
            CKE3 = [CKE3; CKECoeff_new{ggg}]; CKESig3 = [CKESig3,ckeisSig{ggg}];
            numSigVarg3(8,animalInd) = numSigVarg3(8,animalInd) + sum(ckeisSig{ggg});

            %end

    end
end



% sort group1 and group2
g1.coeff=group1;g1.t = rt;
g1sortOrd = coeff_sort(g1,[0,3]);
g2.coeff =group2; g2.t = rt;
g2sortOrd = coeff_sort(g2,[0,3]);
% g3.coeff = group3; g3.t = rt;
% g3sortOrd = coeff_sort(g3,[0,3]);

temp1 = group1(g1sortOrd,:); group1 = temp1; o1SesInd=o1SesInd(g1sortOrd);
temp2 = group2(g2sortOrd,:); group2 = temp2; o2SesInd=o2SesInd(g2sortOrd);
% temp3 = group3(g3sortOrd,:); group3 = temp3; o3SesInd=o3SesInd(g3sortOrd);

choice1 = choice1(g1sortOrd,:);choice2 = choice2(g2sortOrd,:);%choice3 = choice3(g3sortOrd,:);
cSig1 = logical(cSig1(g1sortOrd));cSig2 = logical(cSig2(g2sortOrd)); %cSig3 = logical(cSig3(g3sortOrd));
c1SesInd = o1SesInd(cSig1); c2SesInd = o2SesInd(cSig2);
xn1 = xn1(g1sortOrd,:);xn2 = xn2(g2sortOrd,:);%xn3 = xn3(g3sortOrd,:);
xSig1 = logical(xSig1(g1sortOrd));xSig2 = logical(xSig2(g2sortOrd)); %xSig3 = logical(xSig3(g3sortOrd));
% posRPE1 = posRPE1(g1sortOrd,:);%posRPE2 = posRPE2(g2sortOrd,:);posRPE3 = posRPE3(g3sortOrd,:);
% pRPESig1 = logical(pRPESig1(g1sortOrd));%pRPESig2 = logical(pRPESig2(g2sortOrd)); pRPESig3 = logical(pRPESig3(g3sortOrd));
% negRPE1 = negRPE1(g1sortOrd,:);%negRPE2 = negRPE2(g2sortOrd,:);negRPE3 = negRPE3(g3sortOrd,:);
% nRPESig1 = logical(nRPESig1(g1sortOrd)); %nRPESig2 = logical(nRPESig2(g2sortOrd)); nRPESig3 = logical(nRPESig3(g3sortOrd));
% dQ1 = dQ1(g1sortOrd,:);%dQ2 = dQ2(g2sortOrd,:);dQ3 = dQ3(g3sortOrd,:);
% dQSig1 = dQSig1(g1sortOrd);%dQSig2 = dQSig2(g2sortOrd); dQSig3 = dQSig3(g3sortOrd);
% dK1 = dK1(g1sortOrd,:);%dK2 = dK2(g2sortOrd,:);dK3 = dK3(g3sortOrd,:);
% dKSig1 = dKSig1(g1sortOrd);%dKSig2 = dKSig2(g2sortOrd); dKSig3 = dKSig3(g3sortOrd);
% CKE1 = CKE1(g1sortOrd,:);%CKE2 = CKE2(g2sortOrd,:);CKE3 = CKE3(g3sortOrd,:);
% CKESig1 = CKESig1(g1sortOrd);%CKESig2 = CKESig2(g2sortOrd); CKESig3 = CKESig3(g3sortOrd);




% outcome
tlabel = 'Outcome';
plot_groupSummary(group1,group2, [], rt, tlabel,savesumfigpath)

%% stats for outcome groups
outcome_maxT_1 = zeros(1, size(group1,1)); outcome_maxV_1 = zeros(1,size(group1,1));
outcome_maxT_2 = zeros(1, size(group2,1)); outcome_maxV_2 = zeros(1,size(group2,1));

% group1
for rr = 1:size(group1,1)
    sCurve = smooth(-group1(rr,:));
    sIntp = interp1(rt,sCurve,rt(1):0.01:rt(end));
    tIntp = rt(1):0.01:rt(end);
    [maxV,ind] = max(sIntp(tIntp>0));
    outcome_maxV_1(rr) = maxV;
    tempt = tIntp(tIntp>0);
    outcome_maxT_1(rr) = tempt(ind);


end
% group2
for rr = 1:size(group2,1)
    sCurve = smooth(group2(rr,:));
   sIntp = interp1(rt,sCurve,rt(1):0.01:rt(end));
    tIntp = rt(1):0.01:rt(end);
    [maxV,ind] = max(sIntp(tIntp>0));
    outcome_maxV_2(rr) = maxV;
    tempt = tIntp(tIntp>0);
    outcome_maxT_2(rr) = tempt(ind);

end

savematpath = fullfile(savematsumpath,'outcome_groupstat.mat');
save(savematpath, 'outcome_maxT_1','outcome_maxT_2','outcome_maxV_1','outcome_maxV_2','o1SesInd','o2SesInd');


tlabel = 'Choice';
plot_groupSummary(choice1,choice2, [], rt, tlabel,savesumfigpath)
tlabel = 'Interaction';
plot_groupSummary(xn1,xn2, xn3, [], rt, tlabel,savesumfigpath)
tlabel = 'posRPE';
plot_groupSummary(posRPE1,posRPE2, posRPE3, [], rt, tlabel,savesumfigpath)
tlabel = 'negRPE';
plot_groupSummary(negRPE1,negRPE2, negRPE3, [], rt, tlabel,savesumfigpath)
tlabel = 'dQ';
plot_groupSummary(dQ1,dQ2, dQ3, [], rt, tlabel,savesumfigpath)
tlabel = 'dK';
plot_groupSummary(dK1,dK2, dK3, [], rt, tlabel,savesumfigpath)
tlabel = 'CKE';
plot_groupSummary(CKE1,CKE2, CKE3, [], rt, tlabel,savesumfigpath)

tlabel = 'Choice-Sig';
plot_groupSummary(choice1(logical(cSig1),:),choice2(logical(cSig2),:), cCoeff_nR, rt, tlabel,savesumfigpath)

%% choice stat, for NE data, group 1 and 2 should be negative
% for interactions, rerun cluster for each group
choice_maxV_1 = zeros(1, size(choice1(logical(cSig1),:),1)); choice_maxT_1 = zeros(1, size(choice1(logical(cSig1),:),1));
choice_maxV_2 = zeros(1, size(choice2(logical(cSig2),:),1)); choice_maxT_2 = zeros(1, size(choice2(logical(cSig2),:),1));
%choice_maxV_3 = zeros(1, size(choice3(logical(cSig3),:),1)); choice_maxT_3 = zeros(1, size(choice3(logical(cSig3),:),1));

group1 = choice1(logical(cSig1),:);
group2 = choice2(logical(cSig2),:);
%group3 = choice3(logical(cSig3),:);

for rr = 1:size(group1,1)
    sCurve = -smooth(group1(rr,:));
    sIntp = interp1(rt,sCurve,rt(1):0.01:rt(end));
    tIntp = rt(1):0.01:rt(end);
    [maxV,ind] = max(sIntp(tIntp>0));
    choice_maxV_1(rr) = maxV;
    tempt = tIntp(tIntp>0);
    choice_maxT_1(rr) = tempt(ind);

end
% group2, -group2 for NE since the coefficient is negative
for rr = 1:size(group2,1)
    sCurve = smooth(-group2(rr,:));
   sIntp = interp1(rt,sCurve,rt(1):0.01:rt(end));
    tIntp = rt(1):0.01:rt(end);
    [maxV,ind] = max(sIntp(tIntp>0));
    choice_maxV_2(rr) = maxV;
    tempt = tIntp(tIntp>0);
    choice_maxT_2(rr) = tempt(ind);

end
% group 3, -group3 for NE since the coefficient is negative
for rr = 1:size(group3,1)
    sCurve = smooth(group3(rr,:));
    sIntp = interp1(rt,sCurve,rt(1):0.01:rt(end));
    tIntp = rt(1):0.01:rt(end);
    [maxV,ind] = max(sIntp(tIntp>0));
    choice_maxV_3(rr) = maxV;
    tempt = tIntp(tIntp>0);
    choice_maxT_3(rr) = tempt(ind);

end
savematpath = fullfile(savematsumpath,'choice_groupstat.mat');
save(savematpath, 'choice_maxT_1','choice_maxT_2','choice_maxT_3','choice_maxV_1','choice_maxV_2','choice_maxV_3','c1SesInd','c2SesInd');

tlabel = 'Interaction-Sig';
plot_groupSummary(xn1(logical(xSig1),:),xn2(logical(xSig2),:),  xCoeff_nR, rt, tlabel,savesumfigpath)
%% plot overall pos/negRPE only
tlabel = 'posRPE-Sig';
plot_RPESummary(posRPE1(logical(pRPESig1),:),posRPE2(logical(pRPESig2),:),  posRPECoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'negRPE-Sig';
plot_groupSummary(negRPE1(logical(nRPESig1),:),negRPE2(logical(nRPESig2),:),  negRPECoeff_nR, rt, tlabel,savesumfigpath)

tlabel = 'dQ-sig';
plot_groupSummary(dQ1(logical(dQSig1),:),dQ2(logical(dQSig2),:), dQ3(logical(dQSig3),:), dQCoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'dK-sig';
plot_groupSummary(dK1(logical(dKSig1),:),dK2(logical(dKSig2),:), dK3(logical(dKSig3),:), dKCoeff_nR, rt, tlabel,savesumfigpath)
tlabel = 'CKE-sig';
plot_groupSummary(CKE1(logical(CKESig1),:),CKE2(logical(CKESig2),:), CKE3(logical(CKESig3),:), CKECoeff_nR, rt, tlabel,savesumfigpath)

%% check pos/neg RPE for overlapping space
plot_RPE(posRPECoeff,negRPECoeff,pRPESigInd,nRPESigInd,rt,savesumfigpath)
% check only significant grids (add another group that is not significant
% for outcome but significant for the variable
