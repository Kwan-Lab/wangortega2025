function MP_GRAB_MLR_comp(dataIndex_ACh,dataIndex_NE)

%% ACH
nFiles = size(dataIndex_ACh,1);

subject_mask = [];
animalList = unique(dataIndex_ACh.Animal);

all_coeff_ACh = [];
all_pval_ACh = [];

riseTOutcome_ACh = [];
tau_ACh = [];
for ii = 1:nFiles
    savematpath = fullfile(dataIndex_ACh.BehPath{ii},'analysis-fluo');
    % load behavior files
     fn_beh = dir(fullfile(dataIndex_ACh.BehPath{ii},'beh_cut.mat'));
    
     saveRegName = fullfile(savematpath,'regCR_norm.mat');  % regression for fluo change
   
    if exist(saveRegName)
        load(saveRegName)
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex_ACh.Animal{ii}));
        
        % load choice and reward regression
%          reg_all.regr_time = reg_cr_change.regr_time;
%         reg_all.numPredictor = reg_cr_change.numPredictor;
%         reg_all.nback = reg_cr_change.nback;
%         reg_all.interaction = reg_cr_change.interaction;
%           all_coeff = cat(3,all_coeff, reg_cr_change.coeff);
%         all_pval = cat(3,all_pval, reg_cr_change.pval);
        
    
    % load the MLR with C(n+1)
        for rr = 1:length(reg_cr) 
            all_coeff_ACh = cat(3,all_coeff_ACh, reg_cr{rr}.coeff);
            all_pval_ACh = cat(3, all_pval_ACh, reg_cr{rr}.pval);
        end
        reg_all.regr_time = reg_cr{1}.regr_time;
        reg_all.numPredictor = reg_cr{1}.numPredictor;
        reg_all.nback = reg_cr{1}.nback;
        reg_all.interaction = reg_cr{1}.interaction;
        
        % get the curve of outcome
        sigOutcome = sum(all_pval_ACh(:,7,:)<0.01,3)/size(all_pval_ACh,1);
        %figure;plot(sigOutcome);
        % find the maximam value in 0<3
        [~,tInd] = max(sigOutcome(reg_all.regr_time<=3 & reg_all.regr_time>0));
        tRange = reg_all.regr_time(reg_all.regr_time<=3 & reg_all.regr_time>0);
        realT = tRange(tInd);
        
         % get rise time
        riseTOutcome_ACh = [riseTOutcome_ACh,realT];
        
       
        % get decay time
        f = fit(reg_all.regr_time(reg_all.regr_time>realT),sigOutcome(reg_all.regr_time>realT),'exp1');
        tau_ACh = [tau_ACh,f.b];
    end
end
     

%% NE
nFiles = size(dataIndex_NE,1);

subject_mask = [];
animalList = unique(dataIndex_NE.Animal);

all_coeff_NE = [];
all_pval_NE = [];

riseTOutcome_NE = [];
tau_NE = [];
for ii = 1:nFiles
    savematpath = fullfile(dataIndex_NE.BehPath{ii},'analysis-fluo');
    % load behavior files
     fn_beh = dir(fullfile(dataIndex_NE.BehPath{ii},'beh_cut.mat'));
    
     saveRegName = fullfile(savematpath,'regCR_norm.mat');  % regression for fluo change
   
    if exist(saveRegName)
        load(saveRegName)
        % get subject mask
        subject_mask(end+1) = find(strcmp(animalList, dataIndex_NE.Animal{ii}));
        
        % load choice and reward regression
%          reg_all.regr_time = reg_cr_change.regr_time;
%         reg_all.numPredictor = reg_cr_change.numPredictor;
%         reg_all.nback = reg_cr_change.nback;
%         reg_all.interaction = reg_cr_change.interaction;
%           all_coeff = cat(3,all_coeff, reg_cr_change.coeff);
%         all_pval = cat(3,all_pval, reg_cr_change.pval);
        
    
    % load the MLR with C(n+1)
        for rr = 1:length(reg_cr) 
            all_coeff_NE = cat(3,all_coeff_NE, reg_cr{rr}.coeff);
            all_pval_NE = cat(3, all_pval_NE, reg_cr{rr}.pval);
        end
        reg_all.regr_time = reg_cr{1}.regr_time;
        reg_all.numPredictor = reg_cr{1}.numPredictor;
        reg_all.nback = reg_cr{1}.nback;
        reg_all.interaction = reg_cr{1}.interaction;
        
        % get the curve of outcome
        sigOutcome = sum(all_pval_NE(:,7,:)<0.01,3)/size(all_pval_NE,1);
        %figure;plot(sigOutcome);
        % find the maximam value in 0<3
        [~,tInd] = max(sigOutcome(reg_all.regr_time<=3 & reg_all.regr_time>0));
        tRange = reg_all.regr_time(reg_all.regr_time<=3 & reg_all.regr_time>0);
        realT = tRange(tInd);
        
         % get rise time
        riseTOutcome_NE = [riseTOutcome_NE,realT];
        
       
        % get decay time
        f = fit(reg_all.regr_time(reg_all.regr_time>realT),sigOutcome(reg_all.regr_time>realT),'exp1');
        tau_NE = [tau_NE,f.b];
    end
end

[h,p] = ttest2(riseTOutcome_ACh,riseTOutcome_NE)
[h,p] = ttest2(tau_ACh,tau_NE)