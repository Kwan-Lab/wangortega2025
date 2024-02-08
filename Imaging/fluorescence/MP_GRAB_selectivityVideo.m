function MP_GRAB_selectivityVideo(dataIndex);

% calculate choice and reward selectivity
% choice_s = (df/f_L-df/f_R)/(df/f_L+df/f_R);
% reward_s = (df/f_R-df/f_UR)/df/f_R+df/f_UP);

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
        
        
        load(fullfile(fn_fluo.folder,fn_fluo.name));
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        saveSelName = fullfile(savematpath,'select_norm.mat');  % random forest
        savePrevSelName = fullfile(savematpath,'prev_select_norm.mat'); 
        if exist(saveSelName)
            
            
            
            load(saveSelName);
            colorRange=[-0.3 0.3];
            tlabel = 'Choice';
            savefigpath = savematpath;
            recordingsite = dataIndex.RecordingSite{ii};
             video_selectivity(choicesel,tlabel,colorRange,savefigpath,recordingsite)
             
            tlabel = 'Outcome';
             video_selectivity(outcomesel,tlabel,colorRange,savefigpath,recordingsite)
             
             load( savePrevSelName)
              tlabel = 'Prev choice';
             video_selectivity(choicesel,tlabel,colorRange,savefigpath,recordingsite)
              tlabel = 'Prev outcome';
             video_selectivity(outcomesel,tlabel,colorRange,savefigpath,recordingsite)
            
        end
    
    close all;
    end
end