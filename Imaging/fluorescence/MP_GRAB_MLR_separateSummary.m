function MP_GRAB_MLR_separateSummary(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
setup_figprop
set(0, 'DefaultFigureRenderer', 'painters');
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);

% load the first file to get some parameters
all_coeff1_future = [];
all_pval1_future = [];

all_coeff2_future = [];
all_pval2_future = [];

for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
    % load behavior files
       fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    
     saveRegName = fullfile(savematpath,'regCR_sep.mat');   % regression for fluo change
   
    if exist(saveRegName)
        load(saveRegName)
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
        
    
    % load the MLR with C(n+1)
        for rr = 1:length(reg_cr1) 
            all_coeff1_future = cat(3,all_coeff1_future, reg_cr1{rr}.coeff);
            all_pval1_future = cat(3, all_pval1_future, reg_cr1{rr}.pval);
            
             all_coeff2_future = cat(3,all_coeff2_future, reg_cr2{rr}.coeff);
            all_pval2_future = cat(3, all_pval2_future, reg_cr2{rr}.pval);
        end
        reg_all.regr_time = reg_cr1{1}.regr_time;
        reg_all.numPredictor = reg_cr1{1}.numPredictor;
        reg_all.nback = reg_cr1{1}.nback;
        reg_all.interaction = reg_cr1{1}.interaction;
     % load the ITI regression (n+1 and n)
 
    end
end
        



%% original linear regression
% 1. use bootstrap to get the average and 95% CI for each factor, plot the
% bar plot

%% other things can be done for pupil response:
% correlation: pupil response - latent variable
%                                                     - response time

% go to the save path
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);


%%  linear regression with C(n+1), plot the first half
reg_cr_all1.coeff= all_coeff1_future;

% use bootstrp to get coefficient
reg_cr_all1 = getBootstrp(reg_cr_all1, 0, 0.05);

if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);
reg_cr_all1.regr_time = reg_all.regr_time;
reg_cr_all1.numPredictor = reg_all.numPredictor;
reg_cr_all1.nback = reg_all.nback;
reg_cr_all1.interaction = reg_all.interaction;
reg_cr_all1.pvalThresh= 0.01;

xtitle='Time from cue (s)';
 tlabel={'c(n+1)','c(n)','c(n-1)','c(n-2)','r(n+1)','r(n)', 'r(n-1)','r(n-2)',...
        'c(n+1)*r(n+1)','c(n)*r(n)','c(n-1)*r(n-1)','c(n-2)*r(n-2)','Reward Rate','Cumulative Reward'};
pvalThresh=0.01;
MP_plot_regrcoef_fluo(reg_cr_all1,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-norm-choiceoutcome_future_averageSession1');    %png format
saveas(gcf, 'MLR-norm-choiceoutcome_future_averageSession1', 'fig');
saveas(gcf, 'MLR-norm-choiceoutcome_future_averageSession1','svg');

% plot the figure as number of session that is significant
reg_sig1.coeff = all_coeff1_future;
reg_sig1.pval = all_pval1_future;
reg_sig1.regr_time = reg_all.regr_time;
reg_sig1.numPredictor = reg_all.numPredictor;
reg_sig1.nback = reg_all.nback;
reg_sig1.interaction = reg_all.interaction;
reg_sig1.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);

if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr_fluo(reg_sig1,reg_pval_future_ctrl, reg_sig1.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr_fluo(reg_sig1,[], reg_sig1.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR-norm-choiceoutcome_future_sigSession1');    %png format
saveas(gcf, 'MLR-norm-choiceoutcome_future_sigSession1', 'fig');
saveas(gcf, 'MLR-norm-choiceoutcome_future_sigSession1','svg');

%% plot the second half
reg_cr_all2.coeff= all_coeff2_future;

% use bootstrp to get coefficient
reg_cr_all2 = getBootstrp(reg_cr_all2, 0, 0.05);


reg_cr_all2.regr_time = reg_all.regr_time;
reg_cr_all2.numPredictor = reg_all.numPredictor;
reg_cr_all2.nback = reg_all.nback;
reg_cr_all2.interaction = reg_all.interaction;
reg_cr_all2.pvalThresh= 0.01;

xtitle='Time from cue (s)';
 tlabel={'c(n+1)','c(n)','c(n-1)','c(n-2)','r(n+1)','r(n)', 'r(n-1)','r(n-2)',...
        'c(n+1)*r(n+1)','c(n)*r(n)','c(n-1)*r(n-1)','c(n-2)*r(n-2)','Reward Rate','Cumulative Reward'};
pvalThresh=0.01;
MP_plot_regrcoef_fluo(reg_cr_all2,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-norm-choiceoutcome_future_averageSession2');    %png format
saveas(gcf, 'MLR-norm-choiceoutcome_future_averageSession2', 'fig');
saveas(gcf, 'MLR-norm-choiceoutcome_future_averageSession2','svg');

% plot the figure as number of session that is significant
reg_sig2.coeff = all_coeff2_future;
reg_sig2.pval = all_pval2_future;
reg_sig2.regr_time = reg_all.regr_time;
reg_sig2.numPredictor = reg_all.numPredictor;
reg_sig2.nback = reg_all.nback;
reg_sig2.interaction = reg_all.interaction;
reg_sig2.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);

if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr_fluo(reg_sig2,reg_pval_future_ctrl, reg_sig2.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr_fluo(reg_sig2,[], reg_sig1.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR-norm-choiceoutcome_future_sigSession2');    %png format
saveas(gcf, 'MLR-norm-choiceoutcome_future_sigSession2', 'fig');
saveas(gcf, 'MLR-norm-choiceoutcome_future_sigSession2','svg');



%%
close all

end
