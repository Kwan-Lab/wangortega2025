function MP_GRAB_MLR_acrossAnimals(dataIndex, savefigpath)

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
    all_coeff_future = [];
    all_pval_future = [];
    
    % all_coeff_iti_1 = [];
    % all_pval_iti_1 = [];
    % all_coeff_iti_2 = [];
    % all_pval_iti_2 = [];
    % all_coeff_iti_3 = [];
    % all_pval_iti_3 = [];
    
    for ii = 1:numel(sessionInclude)
        ind = sessionInclude(ii);
        savematpath = fullfile(dataIndex.BehPath{ind},'analysis-fluo');
        % load behavior files
        fn_beh = dir(fullfile(dataIndex.BehPath{ind},'beh_cut.mat'));
        
        saveRegName = fullfile(savematpath,'regCR_norm.mat');  % regression for fluo change
        
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
                all_coeff_future = cat(3,all_coeff_future, reg_cr{rr}.coeff);
                all_pval_future = cat(3, all_pval_future, reg_cr{rr}.pval);
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
            
            
            %%  load control regression for choice and reward regression
            %no control for now
            %
            %         if ~exist('all_coeff_ctrl') && exist('reg_cr_change_ctrl')
            %             % initialize the control matrix
            %             fieldsName = fieldnames(reg_cr_change_ctrl);
            %             for tt = 1:length(fieldsName)
            %                 all_coeff_ctrl.(fieldsName{tt}) = [];
            %                 all_pval_ctrl.(fieldsName{tt}) = [];
            %             end
            %
            %         end
            %         if exist('reg_cr_change_ctrl')
            %         % load the regression results
            %             fieldsName = fieldnames(reg_cr_change_ctrl);
            %             for uu = 1:length(fieldsName)
            %                 all_coeff_ctrl.(fieldsName{uu}) = cat(3, all_coeff_ctrl.(fieldsName{uu}), reg_cr_change_ctrl.(fieldsName{uu}).coeff);
            %                 all_pval_ctrl.(fieldsName{uu}) = cat(3, all_pval_ctrl.(fieldsName{uu}), reg_cr_change_ctrl.(fieldsName{uu}).pval);
            %             end
            %
            %         end
            %
            %         % load future control
            %         if ~exist('all_coeff_future_ctrl')  && exist('reg_cr_future_change_ctrl')
            %             % initialize the control matrix
            %             fieldsName_future = fieldnames(reg_cr_future_change_ctrl);
            %             for tt = 1:length(fieldsName_future)
            %                 all_coeff_future_ctrl.(fieldsName_future{tt}) = [];
            %                 all_pval_future_ctrl.(fieldsName_future{tt}) = [];
            %             end
            %         end
            %         % load the regression results
            %
            %         if exist('reg_cr_future_change_ctrl')
            %             fieldsName_future = fieldnames(reg_cr_future_change_ctrl);
            %             for uu = 1:length(fieldsName_future)
            %                 all_coeff_future_ctrl.(fieldsName_future{uu}) = cat(3, all_coeff_future_ctrl.(fieldsName_future{uu}), reg_cr_future_change_ctrl.(fieldsName_future{uu}).coeff);
            %                 all_pval_future_ctrl.(fieldsName_future{uu}) = cat(3, all_pval_future_ctrl.(fieldsName_future{uu}), reg_cr_future_change_ctrl.(fieldsName_future{uu}).pval);
            %             end
            %         end
        end
    end
    
    
    
    
    %% original linear regression
    % 1. use bootstrap to get the average and 95% CI for each factor, plot the
    % bar plot
    
    %% other things can be done for pupil response:
    % correlation: pupil response - latent variable
    %                                                     - response time
    
    % reg_all.coeff= all_coeff;
    %
    % % use bootstrp to get coefficient
    % reg_all = getBootstrp(reg_all, 0, 0.05);
    %
    
    % go to the save path
    if ~exist(savefigpath)
        mkdir(savefigpath)
    end
    cd(savefigpath);
    
    % %figure 1
    % figure;
    % xAxis = 1:10;
    % bar(xAxis(2:10),reg_all.coeff_bootave(2:10),'FaceColor',[0.7,0.7,0.7])
    %
    % hold on
    % erneg = reg_all.coeff_bootave(2:10)-reg_all.bootlow(2:10);
    % erpos = reg_all.boothigh(2:10)-reg_all.coeff_bootave(2:10);
    % er = errorbar(xAxis(2:10),reg_all.coeff_bootave(2:10),erneg,erpos);
    %
    % er.Color = [0 0 0];
    % er.LineStyle = 'none';
    % set(gca,'xticklabel',{'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'})
    % hold off
    % pvalThresh = NaN;
    % xtickangle(45)
    % ylabel('Coefficients (a.u.)');
    % title('Coefficient for pupil change - choice and reward');
    % pvalThresh = NaN;
    % xtitle = 'Time from cue (s)';
    % tlabel={'c(n)','c(n-1)','c(n-2)','r(n)','r(n-1)','r(n-2)','c(n)xr(n)','c(n-1)xr(n-1)','c(n-2)xr(n-2)'};
    % MP_plot_regrcoef_pupil(reg_all,pvalThresh,tlabel,xtitle);
    %
    % print(gcf,'-dpng','MLR-change_choiceoutcome_averageSession');    %png format
    % saveas(gcf, 'MLR-change_choiceoutcome_averageSession', 'fig');
    % saveas(gcf, 'MLR-change_choiceoutcome_averageSession','svg');
    %
    % % plot the figure as number of session that is significant
    % reg_sig.coeff = all_coeff;
    % reg_sig.pval = all_pval;
    % reg_sig.regr_time = reg_all.regr_time;
    % reg_sig.numPredictor = reg_all.numPredictor;
    % reg_sig.nback = reg_all.nback;
    % reg_sig.interaction = reg_all.interaction;
    % reg_sig.pvalThresh= 0.01;
    %
    % if exist('all_coeff_ctrl')
    %     reg_pval_ctrl = getBootstrp(all_pval_ctrl, 0.01, 0.01);
    %
    % % MP_plot_regr(reg_sig,reg_sig.pvalThresh,tlabel,xtitle,[all_pval_sC; all_pval_sR]);
    %     MP_plot_regr(reg_sig,reg_pval_ctrl, reg_sig.pvalThresh,tlabel,xtitle);
    % else
    %     MP_plot_regr(reg_sig,[], reg_sig.pvalThresh,tlabel,xtitle);
    % end
    %
    %
    % print(gcf,'-dpng','MLR_change-choiceoutcome_sigSession');    %png format
    % saveas(gcf, 'MLR_change-choiceoutcome_sigSession', 'fig');
    % saveas(gcf, 'MLR-change_choiceoutcome_sigSession','svg');
    %
    
    
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
    reg_cr_all.coeff= all_coeff_future;
    
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
    
    %
    % figure;
    % xAxis = 1:8;
    % bar(xAxis(2:8),reg_cr_future.coeff_bootave(2:8),'FaceColor',[0.7,0.7,0.7])
    %
    % hold on
    % erneg = reg_cr_future.coeff_bootave(2:8)-reg_cr_future.bootlow(2:8);
    % erpos =reg_cr_future.boothigh(2:8)-reg_cr_future.coeff_bootave(2:8);
    % er = errorbar(xAxis(2:8),reg_cr_future.coeff_bootave(2:8),erneg,erpos);
    %
    % er.Color = [0 0 0];
    % er.LineStyle = 'none';
    % set(gca,'xticklabel',{'C(n)','C(n+1)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)'})
    % hold off
    %
    % xtickangle(45)
    % ylabel('Coefficients (a.u.)');
    % title('Coefficient for pupil change - choice and reward');
    xtitle='Time from cue (s)';
    tlabel={'c(n+1)','c(n)','c(n-1)','c(n-2)','r(n+1)','r(n)', 'r(n-1)','r(n-2)',...
        'c(n+1)*r(n+1)','c(n)*r(n)','c(n-1)*r(n-1)','c(n-2)*r(n-2)','Reward Rate','Cumulative Reward'};
    pvalThresh=0.01;
    MP_plot_regrcoef_fluo(reg_cr_all,pvalThresh,tlabel,xtitle);
    print(gcf,'-dpng',['MLR-choiceoutcome_future_averageSession_subject_',animalList{aa}]);    %png format
    saveas(gcf, ['MLR-choiceoutcome_future_averageSession_subject_',animalList{aa}], 'fig');
    saveas(gcf, ['MLR-choiceoutcome_future_averageSession_subject_',animalList{aa}],'svg');
    
    % plot the figure as number of session that is significant
    reg_sig.coeff = all_coeff_future;
    reg_sig.pval = all_pval_future;
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
    print(gcf,'-dpng',['MLR-choiceoutcome_future_sigCell_subject_',animalList{aa}]);    %png format
    saveas(gcf, ['MLR-choiceoutcome_future_sigCell_subject_',animalList{aa}], 'fig');
    saveas(gcf, ['MLR-choiceoutcome_future_sigCell_subject_',animalList{aa}],'svg');
    
end

%%
close all

end