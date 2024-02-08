function MP_GRAB_RPE_MLR_acrossAnimals(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
setup_figprop
set(0, 'DefaultFigureRenderer', 'painters');
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);

% load the first file to get some parameters

for aa = 1:numel(animalList)
    sessionInclude = [];
    for tt = 1:nFiles
        if dataIndex.Animal{tt}==animalList{aa}
            sessionInclude = [sessionInclude, tt];
        end
    end
    all_coeff = [];
all_pval= [];
all_coeff_pos = [];
all_coeff_neg = [];
all_pval_pos = [];
all_pval_neg = [];
    
    for ii = 1:numel(sessionInclude)
        ind = sessionInclude(ii);
        savematpath = fullfile(dataIndex.BehPath{ind},'analysis-fluo');
        % load behavior files
        fn_beh = dir(fullfile(dataIndex.BehPath{ind},'beh_cut.mat'));
        
       saveRegName = fullfile(savematpath,'regRL_RPE_norm.mat');  % regression for fluo change
        
        if exist(saveRegName)
            load(saveRegName)
            % get subject mask
            
            % load choice and reward regression
            %          reg_all.regr_time = reg_cr_change.regr_time;
            %         reg_all.numPredictor = reg_cr_change.numPredictor;
            %         reg_all.nback = reg_cr_change.nback;
            %         reg_all.interaction = reg_cr_change.interaction;
            %           all_coeff = cat(3,all_coeff, reg_cr_change.coeff);
            %         all_pval = cat(3,all_pval, reg_cr_change.pval);
            
            
            % load the MLR with C(n+1)
            for rr = 1:length(reg_cr)
                all_coeff = cat(3,all_coeff, reg_cr{rr}.coeff);
            all_pval = cat(3, all_pval, reg_cr{rr}.pval);
             
            all_coeff_pos = cat(3,all_coeff_pos, reg_cr_pos{rr}.coeff);
            all_pval_pos = cat(3,all_pval_pos, reg_cr_pos{rr}.pval);
        
            all_coeff_neg = cat(3,all_coeff_neg, reg_cr_neg{rr}.coeff);
            all_pval_neg = cat(3,all_pval_neg, reg_cr_neg{rr}.pval);
            end
            reg_all.regr_time = reg_cr{1}.regr_time;
            reg_all.numPredictor = reg_cr{1}.numPredictor;
            reg_all.nback = reg_cr{1}.nback;
            reg_all.interaction = reg_cr{1}.interaction;
            % load the ITI regression (n+1 and n)
            
            
            %         all_coeff_iti1 = cat(3,all_coeff_iti1, reg_cr1_change.coeff);
            %         all_pval_iti1 = cat(3,all_pval_iti1, reg_cr1_change.pval);
            %
            %         all_coeff_iti2 = cat(3,all_coeff_iti2, reg_cr2_change.coeff);
            %         all_pval_iti2 = cat(3,all_pval_iti2, reg_cr2_change.pval);
            %
            %         all_coeff_iti3 = cat(3,all_coeff_iti3, reg_cr3_change.coeff);
            %         all_pval_iti3 = cat(3,all_pval_iti3, reg_cr3_change.pval);
            %
            % load the ITI regression (n-1 and n)
            
            
          
        end
    end
    
    
    
    
    %% original linear regression
   
    
    % go to the save path
    if ~exist(savefigpath)
        mkdir(savefigpath)
    end
    cd(savefigpath);
    
  
    %
    %%  linear regression with RPE
    reg_cr_all.coeff= all_coeff;
    
    % use bootstrp to get coefficient
    reg_cr_all = getBootstrp(reg_cr_all, 0, 0.05);
    
    if ~exist(savefigpath)
        mkdir(savefigpath)
    end
    cd(savefigpath);
    reg_cr_all.regr_time = reg_all.regr_time;
    reg_cr_all.numPredictor = reg_all.numPredictor;
    reg_cr_all.nback = reg_all.nback;
    reg_cr_all.interaction = reg_all.interaction;
    reg_cr_all.pvalThresh= 0.01;
    
  
    xtitle='Time from cue (s)';
    tlabel={'C(n)','C(n-1)','R(n-1)','deltaQ','RPE', 'deltaK','CKE','Reward Rate', 'Cumulavtive reward'};

    pvalThresh=0.01;
    MP_plot_regrcoef_fluo(reg_cr_all,pvalThresh,tlabel,xtitle);
    print(gcf,'-dpng',['MLR-RPE_averageSession_subject_',animalList{aa}]);    %png format
    saveas(gcf, ['MLR-RPE_averageSession_subject_',animalList{aa}], 'fig');
    saveas(gcf, ['MLR-RPE_averageSession_subject_',animalList{aa}],'svg');
    
    % plot the figure as number of session that is significant
    reg_sig.coeff = all_coeff;
    reg_sig.pval = all_pval;
    reg_sig.regr_time = reg_all.regr_time;
    reg_sig.numPredictor = reg_all.numPredictor;
    reg_sig.nback = reg_all.nback;
    reg_sig.interaction = reg_all.interaction;
    reg_sig.pvalThresh= 0.01;
    
    % preprocessing the control coefficient and pval: bootstrap the pval to
    % determine a baseline percentage of significant session
    %reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);
    
    if exist('all_pval_future_ctrl')
        reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);
        
        % MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
        MP_plot_regr_fluo(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
    else
        MP_plot_regr_fluo(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
    end
    print(gcf,'-dpng',['MLR-RPE_sigCell_subject_',animalList{aa}]);    %png format
    saveas(gcf, ['MLR-RPE_sigCell_subject_',animalList{aa}], 'fig');
    saveas(gcf, ['MLR-RPE_sigCell_subject_',animalList{aa}],'svg');
    
    %% pos RPE
    reg_cr_all.coeff= all_coeff_pos;
    
    % use bootstrp to get coefficient
    reg_cr_all = getBootstrp(reg_cr_all, 0, 0.05);
    
    if ~exist(savefigpath)
        mkdir(savefigpath)
    end
    cd(savefigpath);
    reg_cr_all.regr_time = reg_all.regr_time;
    reg_cr_all.numPredictor = 5;
    reg_cr_all.nback = reg_all.nback;
    reg_cr_all.interaction = reg_all.interaction;
    reg_cr_all.pvalThresh= 0.01;
    
  
    xtitle='Time from cue (s)';
    tlabel={'C(n)','dQ','C(n-1)','R(n-1)','posRPE'};
    pvalThresh=0.01;
    MP_plot_regrcoef_fluo(reg_cr_all,pvalThresh,tlabel,xtitle);
    print(gcf,'-dpng',['MLR-posRPE_averageSession_subject_',animalList{aa}]);    %png format
    saveas(gcf, ['MLR-posRPE_averageSession_subject_',animalList{aa}], 'fig');
    saveas(gcf, ['MLR-posRPE_averageSession_subject_',animalList{aa}],'svg');
    
    % plot the figure as number of session that is significant
    reg_sig.coeff = all_coeff_pos;
    reg_sig.pval = all_pval_pos;
    reg_sig.regr_time = reg_all.regr_time;
    reg_sig.numPredictor = 5;
    reg_sig.nback = reg_all.nback;
    reg_sig.interaction = reg_all.interaction;
    reg_sig.pvalThresh= 0.01;
    
    % preprocessing the control coefficient and pval: bootstrap the pval to
    % determine a baseline percentage of significant session
    %reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);
    
    if exist('all_pval_future_ctrl')
        reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);
        
        % MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
        MP_plot_regr_fluo(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
    else
        MP_plot_regr_fluo(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
    end
    print(gcf,'-dpng',['MLR-posRPE_sigCell_subject_',animalList{aa}]);    %png format
    saveas(gcf, ['MLR-posRPE_sigCell_subject_',animalList{aa}], 'fig');
    saveas(gcf, ['MLR-posRPE_sigCell_subject_',animalList{aa}],'svg');
    
    %% neg RPE
    reg_cr_all.coeff= all_coeff_neg;
    
    % use bootstrp to get coefficient
    reg_cr_all = getBootstrp(reg_cr_all, 0, 0.05);
    
    if ~exist(savefigpath)
        mkdir(savefigpath)
    end
    cd(savefigpath);
    reg_cr_all.regr_time = reg_all.regr_time;
    reg_cr_all.numPredictor =5;
    reg_cr_all.nback = reg_all.nback;
    reg_cr_all.interaction = reg_all.interaction;
    reg_cr_all.pvalThresh= 0.01;
    
  
    xtitle='Time from cue (s)';
    tlabel={'C(n)','dQ','C(n-1)','R(n-1)','negRPE'};

    pvalThresh=0.01;
    MP_plot_regrcoef_fluo(reg_cr_all,pvalThresh,tlabel,xtitle);
    print(gcf,'-dpng',['MLR-negRPE_averageSession_subject_',animalList{aa}]);    %png format
    saveas(gcf, ['MLR-negRPE_averageSession_subject_',animalList{aa}], 'fig');
    saveas(gcf, ['MLR-negRPE_averageSession_subject_',animalList{aa}],'svg');
    
    % plot the figure as number of session that is significant
    reg_sig.coeff = all_coeff_neg;
    reg_sig.pval = all_pval_neg;
    reg_sig.regr_time = reg_all.regr_time;
    reg_sig.numPredictor = 5;
    reg_sig.nback = reg_all.nback;
    reg_sig.interaction = reg_all.interaction;
    reg_sig.pvalThresh= 0.01;
    
    % preprocessing the control coefficient and pval: bootstrap the pval to
    % determine a baseline percentage of significant session
    %reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);
    
    if exist('all_pval_future_ctrl')
        reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);
        
        % MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
        MP_plot_regr_fluo(reg_sig,reg_pval_future_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
    else
        MP_plot_regr_fluo(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
    end
    print(gcf,'-dpng',['MLR-negRPE_sigCell_subject_',animalList{aa}]);    %png format
    saveas(gcf, ['MLR-negRPE_sigCell_subject_',animalList{aa}], 'fig');
    saveas(gcf, ['MLR-negRPE_sigCell_subject_',animalList{aa}],'svg');
end

%%
close all

end