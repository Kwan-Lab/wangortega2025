function MP_GRAB_PCA(dataIndex)

%% pca analysis
% dPCA
%  also try t-SNE?
% get the first 2/3 PCs, compute the average variance explained
% separate the trials? or not?

nFiles = size(dataIndex,1);

%% go through every session, run analysis
for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    % load fluorescent files
    fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));
    load(fullfile(fn_fluo.folder,fn_fluo.name));
    
    fn_ROI = dir(fullfile(dataIndex.LogFilePath{ii},'cell*.mat'));
    
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
    fn_PCA = dir(fullfile(savematpath,'pca.mat'));  % regression for fluo change
    
    if isempty(fn_PCA)
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
        
        %% get meam intensity
%         meanProjPath = dir(fullfile(dataIndex.BehPath{ii},'*Mean*.tif'));
%         meanProj = loadtiffseq([],fullfile(meanProjPath.folder,meanProjPath.name));
%         
%         % get the mean intensity in ROI
%         meanIntensity = zeros(1, size(fn_ROI,1));
%         %         varIntensity = zeros(1,size(fn_ROI,1));
%         for uu = 1:size(fn_ROI,1)
%             roi = load(fullfile(fn_ROI(uu).folder,fn_ROI(uu).name));
%             meanIntensity(uu) = mean(meanProj(roi.bw==1));
%             %             varIntensity(uu) = nanstd(cells.cellf{uu});
%         end
%         
        % save figure path
        %% run PCA on dFF of all 196 ROIs
        % get the ROI data
        %         pca_observe = zeros(length(cells.dFF{1}),length(cells.dFF));
        %         for jj = 1:length(cells.dFF)
        %             pca_observe(:,jj) = cells.dFF{jj};
        %         end
        %         [coeff,score,latent,tsquared,explained,mu] = pca(pca_observe);
        %
        %         figure;
        % %         biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4'});
        % %             close all;
        %         scatter3(score(:,1),score(:,2),score(:,3))
        %         axis equal
        %         xlabel('1st Principal Component')
        %         ylabel('2nd Principal Component')
        %         zlabel('3rd Principal Component')
        
        %% dPCA
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = [];  %confidence interval
        params.minNumTrial = 50;
        
        if isfield(cells,'normdFF')
                sel_t= cells.t; sel_dFF = cells.normdFF;
            else
                 sel_t= cells.t; sel_dFF = cells.dFF;
            end
        fieldname{1}={'left','reward'};  trialMask{1} = getMask(trials,fieldname{1});
        fieldname{2}={'right','reward'}; trialMask{2} = getMask(trials,fieldname{2});
        fieldname{3}={'left','noreward'}; trialMask{3} = getMask(trials,fieldname{3});
        fieldname{4}={'right','noreward'}; trialMask{4} = getMask(trials,fieldname{4});
% %         
%          fieldname{1}={'left'};  trialMask{1} = getMask(trials,fieldname{1});
%          fieldname{2}={'right'}; trialMask{2} = getMask(trials,fieldname{2});
        % output should be a N*C*O*T matrix
        % N: number of ROIs
        % C: number of choices: ipsi/contra
        % O: number of outcomes
        % T: number of time points
       psth_output=zeros(numel(sel_dFF),2,2,numel(params.window)-1);
       % psth_output=zeros(numel(sel_dFF),2,numel(params.window)-1);
        tic
        parfor uu=1:numel(sel_dFF)
                temp_output = zeros(2,2,numel(params.window)-1);
                for kk=1:numel(fieldname)
                   
                    out = get_psth(sel_dFF{uu}, sel_t, params.trigTime(trialMask{kk}), strjoin(fieldname{kk}), params );
                    if kk == 1
                        temp_output(1,1,:) = out.signal;
                    elseif kk == 2
                        temp_output(2,1,:) = out.signal;
                    elseif kk == 3
                        temp_output(1,2,:) = out.signal;
                    elseif kk == 4
                        temp_output(2,2,:) = out.signal;
                    end
                end
                psth_output(uu,:,:,:) = temp_output;
        end
        toc
        % running dPCA
      combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Choice', 'Outcome', 'Condition-independent', 'C/O Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
        %maxstep = 1e5;
        
       
tic
%[W] = dpca_old(psth_output, 20,[],[]);
 %   'combinedParams', combinedParams);

[W,V,whichMarg] = dpca(psth_output, 20, 'combinedParams', combinedParams);
explVar = dpca_explainedVariance(psth_output, W, V, 'combinedParams', combinedParams);
toc

time = -2.95:0.1:4.95;
dpca_plot(psth_output, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);
dPCA.data = psth_output;
dPCA.encoder = V;
dPCA.decoder = W;

print(gcf,'-dpng',fullfile(savefluofigpath,'dPCA-choiceoutcome'));
saveas(gcf, fullfile(savefluofigpath,'dPCA-choiceoutcome'), 'fig');
saveas(gcf, fullfile(savefluofigpath,'dPCA-choiceoutcome'), 'svg');
   %'timeEvents', timeEvents,               ...
save(fullfile(fn_PCA.folder,fn_PCA.name),'dPCA');

%% plot trajectories of first three components
% trajLR = W'*squeeze(psth_output(:,1,1,:));
% trajRR = W'*squeeze(psth_output(:,2,1,:));
% trajLN = W'*squeeze(psth_output(:,1,2,:));
% trajRN = W'*squeeze(psth_output(:,2,2,:));
% 
% figure;plot3(trajLR(1,:),trajLR(2,:),trajLR(3,:));
% hold on;plot3(trajRR(1,:),trajRR(2,:),trajRR(3,:));
% hold on;plot3(trajLN(1,:),trajLN(2,:),trajLN(3,:));
% hold on;plot3(trajRN(1,:),trajRN(2,:),trajRN(3,:));
% hold on;plot3(trajLR(1,30),trajLR(2,30),trajLR(3,30),'o')
% hold on;plot3(trajLN(1,30),trajLN(2,30),trajLN(3,30),'o')
% 
% hold on;plot3(trajLR(1,1),trajLR(2,1),trajLR(3,1),'+')
% hold on;plot3(trajLN(1,1),trajLN(2,1),trajLN(3,1),'+')

    end
end
end


