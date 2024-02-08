function MP_GRAB_selectivitySpatial(dataIndex)

% spatial properties of choice/outcome selectivity
%% probabily need more grids - pixel selectivity


nFiles = size(dataIndex,1);

for ii = 1:nFiles

    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));


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



           %% set threshold
        Thresh.pos = 0.2; Thresh.neg = -0.2; % threshold to find area with pos/neg modulation
        sigThresh.value = 0.05; sigThresh.alpha = 0.01;  % threshold to find area with significant modulation

        tic

        %% check other variables using regression results (should focus on significance rather than coefficient)

        % variable list:
        % regression 1: cn,rn,cn+1,cn-1,rn-1,cn*rn,average_r,cumulative_r
        % regression 2: dQ,dK,chosenQ,chosenK
        % regression 3: RPE, CKE,
        saveregmaskpath = fullfile(savematpath,'regressionMask.mat');  % file path to save the results

        if exist(saveregmaskpath)
            load(saveregmaskpath);
        end

        % need to get the boundaries of possible coefficient value of all
        % regressions
        saveCRName = fullfile(savematpath,'regCR_norm.mat');
         if exist(saveCRName)
            reg1 = load(saveCRName);
         end
        saveRLName = fullfile(savematpath,'regRL_norm.mat');
          if exist(saveRLName)
            reg2 = load(saveRLName);
         end
        saveRPEName = fullfile(savematpath,'regRL_RPE_norm.mat');
          if exist(saveRPEName)
            reg3 = load(saveRPEName);
          end

         
        %% regression 1 -------------------------------------------------------------------------
             % get the choice from linear regression
        %if ~exist('choiceRegData','var') % if choice regression mask not computed
            label = 'choice';
            choiceInd = 3;

            choiceRegData = getRegSelData(reg1.reg_cr,label,choiceInd,Thresh, sigThresh,savefluofigpath);
        %end

        % outcome regression mask
        
        %if ~exist('outcomeRegData','var') % if choice regression mask not computed
            label = 'outcome';
            outcomeInd = 7;
            outcomeRegData = getRegSelData(reg1.reg_cr,label,outcomeInd,Thresh, sigThresh,savefluofigpath);
        %end
        % save the result

         % c(n+1) regression mask
        %if ~exist('cn_1RegData','var') % if choice regression mask not computed
            label = 'upcoming choice';
            cn_1Ind = 2;
            cn_1RegData = getRegSelData(reg1.reg_cr,label,cn_1Ind,Thresh, sigThresh,savefluofigpath);
        %end

        % r(n+1) regression mask
        %if ~exist('rn_1RegData','var') % if choice regression mask not computed
            label = 'upcoming reward';
            rn_1Ind = 6;
            rn_1RegData = getRegSelData(reg1.reg_cr,label,rn_1Ind,Thresh, sigThresh,savefluofigpath);
        %end

        % c(n-1) regression mask
        %if ~exist('cn__1RegData','var') % if choice regression mask not computed
            label = 'previous choice';
            cn__1Ind = 4;
            cn__1RegData = getRegSelData(reg1.reg_cr,label,cn__1Ind,Thresh, sigThresh,savefluofigpath);
        %end

        % r(n-1) regression mask
        %if ~exist('rn__1RegData','var') % if choice regression mask not computed
            label = 'previous reward';
            rn__1Ind = 8;
            rn__1RegData = getRegSelData(reg1.reg_cr,label,rn__1Ind,Thresh, sigThresh,savefluofigpath);
        %end

        % x(n) regression mask
        %if ~exist('xnRegData','var') % if choice regression mask not computed
            label = 'interaction';
            xnInd = 11;
            xnRegData = getRegSelData(reg1.reg_cr,label,xnInd,Thresh, sigThresh,savefluofigpath);
        %end

         %if ~exist('xn__1RegData','var')
         label = 'previous interaction';
            xnInd = 12;
            xn__1RegData = getRegSelData(reg1.reg_cr,label,xnInd,Thresh, sigThresh,savefluofigpath);
         %end

        % average reward regression mask
        %if ~exist('ave_rRegData','var') % if choice regression mask not computed
            label = 'average reward';
            ave_rInd = 14;
            ave_rRegData = getRegSelData(reg1.reg_cr,label,ave_rInd,Thresh, sigThresh,savefluofigpath) ;
        %end

        % cumulative reward regression mask
        %if ~exist('cum_rRegData','var') % if choice regression mask not computed
            label = 'cumulative reward';
            cum_rInd = 15;
            cum_rRegData = getRegSelData(reg1.reg_cr,label,cum_rInd,Thresh, sigThresh,savefluofigpath) ;
        %end

        %% regression 2-----------------------------------------------------------------------
        %if ~exist('dQRegData','var') % if choice regression mask not computed
            label = 'delta Q';
            dQInd = 8;
            dQRegData = getRegSelData(reg2.reg_cr,label,dQInd,Thresh, sigThresh,savefluofigpath) ;
        %end

        %if ~exist('chosenQRegData','var') % if choice regression mask not computed
            label = 'chosen Q';
            chosenQInd = 9;
            chosenQRegData = getRegSelData(reg2.reg_cr,label,chosenQInd,Thresh, sigThresh,savefluofigpath) ;
        %end

        %if ~exist('dKRegData','var') % if choice regression mask not computed
            label = 'delta K';
            dKInd = 10;
            dKRegData = getRegSelData(reg2.reg_cr,label,dKInd,Thresh, sigThresh,savefluofigpath) ;
        %end

        %if ~exist('chosenKRegData','var') % if choice regression mask not computed
            label = 'chosen K';
            chosenKInd = 11;
            chosenKRegData = getRegSelData(reg2.reg_cr,label,chosenKInd,Thresh, sigThresh,savefluofigpath) ;
        %end

        %% regression 3-----------------------------------------------------------------------

        %if ~exist('RPERegData','var') % if choice regression mask not computed
            label = 'RPE';
            RPEInd = 6;
            RPERegData = getRegSelData(reg3.reg_cr,label,RPEInd,Thresh, sigThresh,savefluofigpath) ;
        %end
            label = 'posRPE';
            RPEInd = 6;
            posRPERegData = getRegSelData(reg3.reg_cr_pos,label,RPEInd,Thresh, sigThresh,savefluofigpath) ;

            label = 'negRPE';
            RPEInd = 6;
            negRPERegData = getRegSelData(reg3.reg_cr_neg,label,RPEInd,Thresh, sigThresh,savefluofigpath) ;
        %if ~exist('CKERegData','var') % if choice regression mask not computed
            label = 'CKE';
            CKEInd = 8;
            CKERegData = getRegSelData(reg3.reg_cr,label,CKEInd,Thresh, sigThresh,savefluofigpath) ;
        %end

        %% save the results
        save(saveregmaskpath,'outcomeRegData','choiceRegData','cn_1RegData','cn__1RegData',...
                             'rn_1RegData','rn__1RegData', 'xnRegData', 'xn__1RegData', 'ave_rRegData', 'cum_rRegData',...
                             'dQRegData', 'chosenQRegData','dKRegData','chosenKRegData',...
                             'RPERegData','posRPERegData','negRPERegData', 'CKERegData')
    end

        close all
        clearvars -except analysispath dataIndex logfilepath root_path save_path_fluo
    toc
end
end
