% analysis the GRAB_NE and GRAB_ACh signal

clearvars;
close all;
setup_figprop;


root_path_NE = 'K:\GRAB2024\2PData\NE';
root_path_ACh = 'K:\GRAB2024\2PData\ACh';


%% matching pennies behavior

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('--------------------------------------------------open---------');

% Look for data files and create a database index

%% NE
logfilepath_NE = fullfile(root_path_NE,'data');
analysispath_NE = fullfile(root_path_NE,'analysis');
dataIndex_NE = makeDataIndex(logfilepath_NE, analysispath_NE);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex_NE = MP_GRAB_createBehMatFiles(dataIndex_NE,'NE');

% sort Index according to experiment data

dataIndex_NE = sortdataIndex(dataIndex_NE);

% Determine if each session fulfill performance criteria (Matching Pennies)
MP_determineBehCriteria(dataIndex_NE);

% ACh
logfilepath_ACh = fullfile(root_path_ACh,'data');
analysispath_ACh = fullfile(root_path_ACh,'analysis');
dataIndex_ACh = makeDataIndex(logfilepath_ACh, analysispath_ACh);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex_ACh = MP_GRAB_createBehMatFiles(dataIndex_ACh,'ACh');

% sort Index according to experiment data

dataIndex_ACh = sortdataIndex(dataIndex_ACh);

% Determine if each session fulfill performance criteria (Matching Pennies)
MP_determineBehCriteria(dataIndex_ACh);

%% behavior

% NE
nFiles_NE = size(dataIndex_NE);

for ii = 1:nFiles_NE
    savematpath_NE = dataIndex_NE.BehPath{ii};
    MP_GRAB_session(dataIndex_NE.BehPath{ii},dataIndex_NE.LogFileName{ii},savematpath_NE);
end
save_path_NE = fullfile(root_path_NE,'summary','figs_summary');


MP_GRAB_behaviorPerAnimal(dataIndex_NE,save_path_NE);

MP_GRAB_behaviorAll(dataIndex_NE, save_path_NE);

% ACh
nFiles_ACh = size(dataIndex_ACh,1);

for ii = 1:nFiles_ACh
    savematpath_ACh = dataIndex_ACh.BehPath{ii};
    MP_GRAB_session(dataIndex_ACh.BehPath{ii},dataIndex_ACh.LogFileName{ii},savematpath_ACh);
end
save_path_ACh = fullfile(root_path_ACh,'summary','figs_summary');


MP_GRAB_behaviorPerAnimal(dataIndex_ACh,save_path_ACh);

MP_GRAB_behaviorAll(dataIndex_ACh, save_path_ACh);

% summary
MP_GRAB_behComp(dataIndex_ACh, dataIndex_NE, save_path_ACh)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Behavior - Model fitting

% NE
model_path_NE = fullfile(root_path_NE,'mat_models');


MP_GRAB_fittingPerAnimal(dataIndex_NE,model_path_NE);

% save the predicted latent variable into each behavior mat files
MP_GRAB_saveLatent(dataIndex_NE, model_path_NE);

%ACh

model_path_ACh = fullfile(root_path_ACh,'mat_models');


MP_GRAB_fittingPerAnimal(dataIndex_ACh,model_path_ACh);

% save the predicted latent variable into each behavior mat files
MP_GRAB_saveLatent(dataIndex_ACh, model_path_ACh);
%% fluorescent signal preprocessing

% NE
MP_GRAB_preprocess(dataIndex_NE);
%ACh
MP_GRAB_preprocess(dataIndex_ACh);


MP_GRAB_checkbaseline(dataIndex_NE);
MP_GRAB_checkbaseline(dataIndex_ACh);

% in 910-MP-1213, one img trigger was missing at trial 450 (two trials 450
% + 451) in one imaging file
% solve later

%% simple plots
MP_GRAB_simpleplots(dataIndex_NE);
MP_GRAB_simpleplots_average(dataIndex_NE);
MP_GRAB_simpleplots_summary(dataIndex_NE,save_path_NE);
MP_GRAB_simpleplots(dataIndex_ACh);
MP_GRAB_simpleplots_average(dataIndex_ACh);
MP_GRAB_simpleplots_summary(dataIndex_ACh,save_path_ACh );



%% regression: observable variable, choice, outcome, average reward, cumulative reward

% NE
model_path_NE = fullfile(root_path_NE,'mat_models');

save_path_fluo_NE = fullfile(root_path_NE,'summary','figs_summary_fluo');

% ACh
model_path_ACh = fullfile(root_path_ACh,'mat_models');

save_path_fluo_ACh = fullfile(root_path_ACh,'summary','figs_summary_fluo');
%

%% Linear regression with choice and reward
%% for linear regression - load averagve grid intensity, and its corresponding linear regression R2


MP_GRAB_MLR(dataIndex_NE);


MP_GRAB_MLR_acrossAnimals(dataIndex_NE,save_path_fluo_NE);
MP_GRAB_MLR_acrossSessions(dataIndex_NE, save_path_fluo_NE)



% ACh
MP_GRAB_MLR(dataIndex_ACh);


MP_GRAB_MLR_acrossAnimals(dataIndex_ACh,save_path_fluo_ACh);
MP_GRAB_MLR_acrossSessions(dataIndex_ACh, save_path_fluo_ACh);


%% latent variable, ROIs sparsely modulated by latent variables
%% no clear conclusions drawn

%NE
MP_GRABRL_MLR(dataIndex_NE)
MP_GRABRL_MLR_acrossSessions(dataIndex_NE, save_path_fluo_NE)

MP_GRABRL_RPE_MLR(dataIndex_NE);
MP_GRAB_RPE_MLR_acrossAnimals(dataIndex_NE, save_path_fluo_NE)
MP_GRABRL_RPE_MLR_acrossSessions(dataIndex_NE, save_path_fluo_NE)

% ACh
MP_GRABRL_MLR(dataIndex_ACh)
MP_GRABRL_MLR_acrossSessions(dataIndex_ACh, save_path_fluo_ACh)

MP_GRABRL_RPE_MLR(dataIndex_ACh);
MP_GRAB_RPE_MLR_acrossAnimals(dataIndex_ACh, save_path_fluo_ACh)
MP_GRABRL_RPE_MLR_acrossSessions(dataIndex_ACh, save_path_fluo_ACh)

%% spatial & temporal patterns

% NE
% MP_GRAB_selectivitySpatial(dataIndex_NE);
% MP_GRAB_selectivityCorr(dataIndex_NE);
% MP_GRAB_selectivitySummary(dataIndex_NE, save_path_fluo_NE);

MP_GRAB_temporalCorr(dataIndex_NE);
save_path_mat_NE = fullfile(root_path_NE,'summary','data_summary');
MP_GRAB_temporalCorrSummary(dataIndex_NE,save_path_fluo_NE,save_path_mat_NE)
MP_GRAB_clusterEval(dataIndex_NE, save_path_fluo_NE, save_path_mat_NE)

% temporal correlation (same grid, correlation of different variables)

% ACh
% MP_GRAB_selectivitySpatial(dataIndex_ACh);
% MP_GRAB_selectivityCorr(dataIndex_ACh);
% MP_GRAB_selectivitySummary(dataIndex_ACh, save_path_fluo_ACh);


MP_GRAB_temporalCorr(dataIndex_ACh);
save_path_mat_ACh = fullfile(root_path_ACh,'summary','data_summary');
MP_GRAB_temporalCorrSummary(dataIndex_ACh,save_path_fluo_ACh, save_path_mat_ACh)
MP_GRAB_clusterEval(dataIndex_ACh, save_path_fluo_ACh, save_path_mat_ACh)

% compare NE and ACh result
MP_GRAB_tempComp(save_path_mat_NE, save_path_mat_ACh,save_path_fluo_ACh)

