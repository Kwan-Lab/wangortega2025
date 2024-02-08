function MP_GRABRL_RPE_MLR_acrossSessions(dataIndex, savefigpath)

%% reference to Sul et al.2011
% averaged within subject
setup_figprop
set(0, 'DefaultFigureRenderer', 'painters');
nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);

% load the first file to get some parameters
all_coeff = [];
all_pval= [];
all_coeff_pos = [];
all_coeff_neg = [];
all_pval_pos = [];
all_pval_neg = [];


for ii = 1:nFiles
    savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
    % load behavior files
       fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    
     saveRegName = fullfile(savematpath,'regRL_RPE_norm.mat');  % regression for fluo change
   
    if exist(saveRegName)
        load(saveRegName)
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex.Animal{ii}));
        
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
    
    end
end
        
% go to the save path
if ~exist(savefigpath)
    mkdir(savefigpath)
end
cd(savefigpath);

%% combine sessions from same animal together
% for ii = 1:length(animalList)
%     newanimalList{ii} = animalList{ii}(1:3);
% end
% newanimalList = unique(newanimalList);
% 
% newsubject_mask = subject_mask;
% for ii=1:length(animalList)
%     if ii == 4
%         newsubject_mask(subject_mask==ii) = 3;
%     elseif ii == 5 || ii == 6 || ii == 7
%         newsubject_mask(subject_mask == ii) = 4;
%     elseif ii == 8
%         newsubject_mask(subject_mask == ii) = 5;
%     end
% end

% %% plot the combined animals
% posSigCn = zeros(size(reg_sig.coeff,1), length(newanimalList));
% negSigCn = zeros(size(reg_sig.coeff,1), length(newanimalList));
% for tt = 1:length(newanimalList)
%     if sum(newsubject_mask == tt) > 0   % if pupil data for certain animal exists
%         reg_sub = reg_sig;
%         reg_sub.coeff = reg_sig.coeff(:,:,newsubject_mask == tt);
%         reg_sub.pval = reg_sig.pval(:,:, newsubject_mask == tt);
%         posSigCn(:,tt) = sum((reg_sub.pval(:,2,:)<reg_sub.pvalThresh&reg_sub.coeff(:,2,:)>0),3)/sum(newsubject_mask == tt);
%         negSigCn(:,tt) = sum((reg_sub.pval(:,2,:)<reg_sub.pvalThresh&reg_sub.coeff(:,2,:)<0),3)/sum(newsubject_mask == tt);
%     end
% end
% 
% figure;
% for uu = 1:length(newanimalList)
%     subplot(5,1,uu)
%     %plot(reg_sig.regr_time, posSig(:,1),'r');
%     H1=area(reg_sig.regr_time,posSigCn(:,uu));
%     set(H1(1),'FaceColor',[1 0.5 0],'EdgeColor',[1,1,0]);
%     hold on;
%     %plot(reg_sig.regr_time, -negSig(:,1));
%     H2=area(reg_sig.regr_time,-negSigCn(:,uu));
%     set(H2(1),'FaceColor',[0 1 0],'EdgeColor',[0,1,0]);
%     hold on;
%     plot([reg_sig.regr_time(1) reg_sig.regr_time(end)], [0 0],'k','LineWidth',1.5);
%     hold on;
%     plot([0 0],[-1 1],'k','LineWidth', 1.5);
%     set(gca,'box','off');
%     set(gca,'XColor','none','YColor','none')
%     pos = get(gca, 'Position');
%     pos(4) = 0.18;
%     set(gca, 'Position', pos)
% end
% print(gcf,'-dpng','MLR_change-posneg5_sigSession');    %png format
% saveas(gcf, 'MLR_change-posneg5_sigSession', 'fig');
% saveas(gcf, 'MLR-change_posneg5_sigSession','svg');
% 
%%  linear regression with C(n+1)
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
print(gcf,'-dpng','MLR-norm_valueupdating_averageSession');    %png format
saveas(gcf, 'MLR-norm_valueupdating_averageSession', 'fig');
saveas(gcf, 'MLR-norm_valueupdating_averageSession','svg');

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
print(gcf,'-dpng','MLR_norm-valueupdating_sigSession');    %png format
saveas(gcf, 'MLR_norm-valueupdating_sigSession', 'fig');
saveas(gcf, 'MLR-norm_valueupdating_sigSession','svg');

%% plot pos/neg RPE

%% plot Positive RPE
reg_sig_pos.coeff = all_coeff_pos;
reg_sig_pos.pval = all_pval_pos;
reg_sig_pos.regr_time = reg_all.regr_time;
reg_sig_pos.numPredictor = 5;
reg_sig_pos.nback = reg_all.nback;
reg_sig_pos.interaction = reg_all.interaction;
reg_sig_pos.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);
 tlabel={'C(n)','deltaQ','C(n-1)','R(n-1)','posRPE'};
if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr_fluo(reg_sig_pos,reg_pval_future_ctrl, reg_sig_pos.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr_fluo(reg_sig_pos,[], reg_sig_pos.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR_norm-valueupdating_posRPEsigSession');    %png format
saveas(gcf, 'MLR_norm-valueupdating_posRPEsigSession', 'fig');
saveas(gcf, 'MLR-norm_valueupdating_posRPEsigSession','svg');

%% plot negative RPE
reg_sig_neg.coeff = all_coeff_neg;
reg_sig_neg.pval = all_pval_neg;
reg_sig_neg.regr_time = reg_all.regr_time;
reg_sig_neg.numPredictor = 5;
reg_sig_neg.nback = reg_all.nback;
reg_sig_neg.interaction = reg_all.interaction;
reg_sig_neg.pvalThresh= 0.01;

% preprocessing the control coefficient and pval: bootstrap the pval to
% determine a baseline percentage of significant session
%reg_pval_futurectrl = getBootstrp(all_pval_futurectrl, 0.01);
 tlabel={'C(n)','deltaQ','C(n-1)','R(n-1)','negRPE'};
if exist('all_pval_future_ctrl')
    reg_pval_future_ctrl = getBootstrp(all_pval_future_ctrl, 0.01, 0.05);

% MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    MP_plot_regr_fluo(reg_sig_neg,reg_pval_future_ctrl, reg_sig_neg.pvalThresh,tlabel,xtitle);
else
    MP_plot_regr_fluo(reg_sig_neg,[], reg_sig_neg.pvalThresh,tlabel,xtitle);
end
print(gcf,'-dpng','MLR_norm-valueupdating_negRPEsigSession');    %png format
saveas(gcf, 'MLR_norm-valueupdating_negRPEsigSession', 'fig');
saveas(gcf, 'MLR-norm_valueupdating_negRPEsigSession','svg');

% plot pos/neg RPE

reg_all_pos = reg_all;  % get other parameters from reg_all
reg_all_pos.coeff= all_coeff_pos;
reg_all_pos = getBootstrp(reg_all_pos,0,0.05);

reg_all_neg = reg_all;
reg_all_neg.coeff= all_coeff_neg;
reg_all_neg = getBootstrp(reg_all_neg, 0, 0.05);

reg_all_sign = reg_all_pos;
reg_all_sign.numPredictor = 2;
% the first one doesnt matter
reg_all_sign.coeff = [reg_all_pos.coeff(:,1,:),reg_all_pos.coeff(:,end,:),reg_all_neg.coeff(:,end,:)];
reg_all_sign.bootSig = [reg_all_pos.bootSig(:,1,:),reg_all_pos.bootSig(:,end,:),reg_all_neg.bootSig(:,end,:)];
reg_all_sign.coeff_bootave = [reg_all_pos.coeff_bootave(:,1),reg_all_pos.coeff_bootave(:,end),reg_all_neg.coeff_bootave(:,end)];
reg_all_sign.bootlow = [reg_all_pos.bootlow(:,1),reg_all_pos.bootlow(:,end), reg_all_neg.bootlow(:,end)];
reg_all_sign.boothigh = [reg_all_pos.boothigh(:,1),reg_all_pos.boothigh(:,end), reg_all_neg.boothigh(:,end)];
tlabel = {'posRPE', 'negRPE'};

MP_plot_regrcoef_GRAB(reg_all_sign,pvalThresh,tlabel,xtitle);
print(gcf,'-dpng','MLR-RL_allRPE_lag0_averageSession_cut');    %png format
saveas(gcf, 'MLR-RL_allRPE_lag0_averageSession_cut', 'fig');

savematpath = fullfile(savefigpath,'Result');
if ~exist(savematpath)
    mkdir(savematpath);
end
save(fullfile(savematpath,'regression3Sum'),'reg_cr_all','reg_sig','reg_sig_pos','reg_sig_neg','reg_all_sign');
%%
close all

end