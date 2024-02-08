function MP_GRAB_selectivityCorr(dataIndex)


% calculate perframe cross correlation with significance data
nFiles = size(dataIndex,1);

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
        
        
        % load(fullfile(fn_fluo.folder,fn_fluo.name));
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        
        warning('off','all')
        
        sigThresh = 0.01;
        
        tic
       
        %% regression
        saveregXcorrpath = fullfile(savematpath,'regressionXcorr.mat');  % file path to save the results
        saveregVarcorrpath = fullfile(savematpath,'regressionVarcorr.mat');
        if exist(saveregXcorrpath)
            load(saveregXcorrpath);
        end
        
        saveCRName = fullfile(savematpath,'regCR_norm.mat');
        reg_cr{1}=load(saveCRName);
        saveRLName = fullfile(savematpath,'regRL_norm.mat');
        if exist(saveRLName)
            reg_cr{2}=load(saveRLName);
        end
        saveRPEName = fullfile(savematpath,'regRL_RPE_norm.mat');
        if exist(saveRPEName)
            reg_cr{3}=load(saveRPEName);
        end
        
        %% regression 1--------------
        
        if exist(saveCRName)
            
            %if ~exist('outcomeRegData','var')
            label = 'outcome';
            outcomeInd = 7;
            outcomeRegData = getRegSelCorr(reg_cr{1},label,outcomeInd,sigThresh,savefluofigpath);
            %end
            
            %if ~exist('choiceRegData','var')
            label = 'choice';
            choiceInd = 3;
            choiceRegData = getRegSelCorr(reg_cr{1},label,choiceInd,sigThresh,savefluofigpath);
            %end
            % c(n+1) regression mask
            %if ~exist('cn_1RegData','var') % if choice regression mask not computed
            label = 'upcoming choice';
            cn_1Ind = 2;
            cn_1RegData = getRegSelCorr(reg_cr{1},label,cn_1Ind,sigThresh,savefluofigpath);
            %end
            
            % c(n-1) regression mask
            %if ~exist('cn__1RegData','var') % if choice regression mask not computed
            label = 'previous choice';
            cn__1Ind = 4;
            cn__1RegData = getRegSelCorr(reg_cr{1},label,cn__1Ind,sigThresh,savefluofigpath);
            %end
            
            % r(n-1) regression mask
            %if ~exist('rn__1RegData','var') % if choice regression mask not computed
            label = 'previous reward';
            rn__1Ind = 8;
            rn__1RegData = getRegSelCorr(reg_cr{1},label,rn__1Ind,sigThresh,savefluofigpath);
            %end
            
            % r(n+1) regression mask
            %if ~exist('rn_1RegData','var') % if choice regression mask not computed
            label = 'upcoming reward';
            rn_1Ind = 6;
            rn_1RegData = getRegSelCorr(reg_cr{1},label,rn_1Ind,sigThresh,savefluofigpath);
            %end
            
            % x(n) regression mask
            %if ~exist('xnRegData','var') % if choice regression mask not computed
            label = 'interaction';
            xnInd = 11;
            xnRegData = getRegSelCorr(reg_cr{1},label,xnInd,sigThresh,savefluofigpath);
            %end
            
            % average reward regression mask
            %if ~exist('ave_rRegData','var') % if choice regression mask not computed
            label = 'average reward';
            ave_rInd = 14;
            ave_rRegData = getRegSelCorr(reg_cr{1},label,ave_rInd,sigThresh,savefluofigpath) ;
            %end
            
            % cumulative reward regression mask
            %if ~exist('cum_rRegData','var') % if choice regression mask not computed
            label = 'cumulative reward';
            cum_rInd = 15;
            cum_rRegData = getRegSelCorr(reg_cr{1},label,cum_rInd,sigThresh,savefluofigpath) ;
            %end
            
            %% regression 2-----------------------------------------------------------------------
            
            
            %if ~exist('dQRegData','var') % if choice regression mask not computed
            label = 'delta Q';
            dQInd = 8;
            dQRegData = getRegSelCorr(reg_cr{2},label,dQInd,sigThresh,savefluofigpath) ;
            %end
            
            %if ~exist('chosenQRegData','var') % if choice regression mask not computed
            label = 'chosen Q';
            chosenQInd = 9;
            chosenQRegData = getRegSelCorr(reg_cr{2},label,chosenQInd,sigThresh,savefluofigpath) ;
            %end
            
            %if ~exist('dKRegData','var') % if choice regression mask not computed
            label = 'delta K';
            dKInd = 10;
            dKRegData = getRegSelCorr(reg_cr{2},label,dKInd,sigThresh,savefluofigpath) ;
            %end
            
            %if ~exist('chosenKRegData','var') % if choice regression mask not computed
            label = 'chosen K';
            chosenKInd = 11;
            chosenKRegData = getRegSelCorr(reg_cr{2},label,chosenKInd,sigThresh,savefluofigpath) ;
            %end
            
            %% regression 3-----------------------------------------------------------------------
            
            %if ~exist('RPERegData','var') % if choice regression mask not computed
            label = 'RPE';
            RPEInd = 6;
            RPERegData = getRegSelCorr(reg_cr{3},label,RPEInd,sigThresh,savefluofigpath) ;
            %end
            
            %if ~exist('CKERegData','var') % if choice regression mask not computed
            label = 'CKE';
            CKEInd = 8;
            CKERegData = getRegSelCorr(reg_cr{3},label,CKEInd,sigThresh,savefluofigpath) ;
            %end
            
            
            %% cross-correlation
            varName = {'cn+1','cn','cn-1','rn','rn-1','xn','xn-1','ave_r','Cum_r','dQ','chosenQ','dK','chosenK','RPE','CKE'};
            varInd = [2,3,4,7,8,11,12,14,15,8,9,10,11,6,8];
            varReg = [1,1,1,1,1,1,1,1,1,2,2,2,2,3,3];
            
            
            % calculate crosscorrelation
            if ~exist(saveregVarcorrpath)
                xCorrCell = cell(length(varInd),length(varInd));
                for ii = 1:length(varName)
                    for jj=ii:length(varName)
                        
                        xCorrCell{ii,jj}=getRegSelxCorr(reg_cr{varReg(ii)},reg_cr{varReg(jj)},varName{ii},varName{jj},varInd(ii),varInd(jj), sigThresh, savefluofigpath);
                    end
                end
                save(saveregVarcorrpath,'xCorrCell');
            end
            %% save the results
            save(saveregXcorrpath,'outcomeRegData','choiceRegData','cn_1RegData','cn__1RegData',...
                'rn__1RegData','rn_1RegData', 'xnRegData', 'ave_rRegData', 'cum_rRegData',...
                'dQRegData', 'chosenQRegData','dKRegData','chosenKRegData',...
                'RPERegData', 'CKERegData')
            
        end
        
        close all
        clearvars -except analysispath dataIndex logfilepath root_path save_path_fluo
        toc
    end
end
