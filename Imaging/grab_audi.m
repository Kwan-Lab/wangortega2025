% analyze grab_ne data with auditory cue
clearvars;
close all;
setup_figprop;

root_path = 'V:\HongliWang\GRAB\GRAB_ctrl';

audi_logfilepath = fullfile(root_path,'data','audi');
spon_logfilepath = fullfile(root_path,'data','spon');

audi_analysispath = fullfile(root_path,'analysis\audi');
spon_analysispath = fullfile(root_path,'analysis\spon');
audi_dataIndex = makeDataIndex_grabctrl(audi_logfilepath, audi_analysispath);
audi_dataIndex = audi_createBehMatFiles(audi_dataIndex);

spon_dataIndex = makeDataIndex_grabctrl(spon_logfilepath, spon_analysispath);
spon_dataIndex = audi_createBehMatFiles(spon_dataIndex);

%% pupil

createPupilFiles(spon_dataIndex);

% plot pupil-tone PSTH

%audi_pupil_plots(dataIndex([14],:));
%% fluo
MP_GRAB_preprocess(dataIndex);
audi_GRAB_plots(dataIndex);

%% check correlation
audi_GRAB_correlation(dataIndex([14],:));

% plot pupil PSTH in auditory stimulus trials
audi_pupil_plots(dataIndex([25],:));

%% how to find the patchness in the video?
% autocorrealtion?
GRAB_pattern(dataIndex(13,:));

%% auditory-evoked response

audi_GRAB_pupil_plots(audi_dataIndex);

savefigpath = fullfile(root_path,'summary','audi');
if ~exist(savefigpath)
    mkdir(savefigpath)
end
audi_GRAB_summary(audi_dataIndex,savefigpath);


%% spontaneous recordings

% spontaneous, GRAB_Ach, before MP. 5,8,10,12

% coherence and cross spectrum
fluo_pupil_plots(spon_dataIndex);

savefigpath = fullfile(root_path,'summary','spon');
if ~exist(savefigpath)
    mkdir(savefigpath)
end
fluo_pupil_GRAB_summary(spon_dataIndex, savefigpath);
% coherence
coherence_plots(dataIndex,savefigpath);

% cross-correlation of pupil and GRAB signals
fluo_pupil_GRAB_correlation(spon_dataIndex);
% summary

fluo_pupil_corrSummary(spon_dataIndex,savefigpath);
