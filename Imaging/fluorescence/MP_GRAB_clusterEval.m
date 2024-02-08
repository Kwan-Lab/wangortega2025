function MP_GRAB_clusterEval(dataIndex, savesumfigpath,savematsumpath)

% find the optimal number of clusters using average silhouette methods
% for each number of clusters, calculate the average sulhouette of
% observation

nFiles = size(dataIndex,1);
numClust = 1:5;

aveSil = zeros(length(numClust),nFiles);
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
        
        % load the regression results
        saveCRName = fullfile(savematpath,'regCR_norm.mat');
         if exist(saveCRName)
            reg_cr{1} = load(saveCRName);
         end
         
         label = 'outcome';
         sigThresh.value = 0.05; sigThresh.alpha = 0.01;  % threshold to find area with significant modulation

            outcomeInd = 7;
            outcometempData = getRegautoCorrData(reg_cr{1}.reg_cr,label,outcomeInd,[0,3],sigThresh,savefluofigpath);

            outcomeSel.t = reg_cr{1}.reg_cr{1}.regr_time;
           
            outcomeSel.coeff= outcometempData.coeff;

        %cluster grid based on the correlation within -3-5s
        allInd = 1:size(outcometempData.coeff,1);
        sigGrid =~isnan(outcometempData.tempCorrCoeff);
        sigInd = allInd(sigGrid(:));
        outcomeSel.sigInd = sigInd;
        outcomeSel.sigMat = sigGrid;
        
        
        xtitle = 'Time from cue';
        tlabel1='Outcome coefficient';
        
        nCells=size(outcomeSel.coeff,1);
            
            %preference is (signal_event1 - signal_event2)/(signal_event1 + signal_event2)
            pref=[];
            for j=1:nCells
                pref(:,j)=outcomeSel.coeff(j,:);
            end
             clusterMat = pref(:,outcomeSel.sigInd)';
        for cc = 1:length(numClust)
            maxclust = numClust(cc);
            T = clusterdata(clusterMat,'Linkage','complete','SaveMemory','off','Maxclust',maxclust,'distance','correlation');           
            [s,h] = silhouette(clusterMat,T,'correlation');
            aveSil(cc,ii) = mean(s);
        end
    end
end

display('Done');
figure;boxplot(aveSil');
[h,p]=ttest(aveSil(1,:),aveSil(2,:))
anova1(aveSil')