function MP_GRAB_selectivity(dataIndex);

% calculate choice and reward selectivity
% choice_s = (df/f_L-df/f_R)/(df/f_L+df/f_R);
% reward_s = (df/f_R-df/f_UR)/df/f_R+df/f_UP);

nFiles = size(dataIndex,1);




%% go through every session, running RF separately
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
        saveSelName = fullfile(savematpath,'select_nonorm.mat');  % random forest
        
        if ~exist(saveSelName)
            params=[];
            params.trigTime = trialData.cueTimes;
            params.xtitle = 'Time from cue (s)';
            params.window = [-3:0.1:5];
            params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
            params.CI = 0.95;  %confidence interval
            params.minNumTrial = 50;
        
                       
%             if isfield(cells,'normdFF')
%                 sel_t= cells.t; sel_dFF = cells.normdFF;
%             else
                 sel_t= cells.t; sel_dFF = cells.dFF;
%             end
            
            fieldname1={'left'}; trialMask1 = getMask(trials,fieldname1);
            fieldname2={'right'}; trialMask2 = getMask(trials,fieldname2);
            fieldname3={'reward'}; trialMask3 = getMask(trials,fieldname3);
             fieldname4={'noreward'}; trialMask4 = getMask(trials,fieldname4);
             tic
        parfor tt = 1:length(sel_dFF)
            % get PSTH
            
            %choice
            
            psth_left = get_psth( sel_dFF{tt}, sel_t, params.trigTime(trialMask1), strjoin(fieldname1), params);
            
            psth_right = get_psth( sel_dFF{tt}, sel_t, params.trigTime(trialMask2), strjoin(fieldname2), params);
               
            % reward
            
            psth_reward= get_psth( sel_dFF{tt}, sel_t, params.trigTime(trialMask3), strjoin(fieldname3), params);
           
            psth_noreward = get_psth( sel_dFF{tt}, sel_t, params.trigTime(trialMask4), strjoin(fieldname4), params);
                
            choicesel{tt} = calc_selectivity(psth_left,psth_right);
            outcomesel{tt} = calc_selectivity(psth_reward,psth_noreward);
        end
        toc
        else
            load(saveSelName)
        end
        tlabel1 = 'Choice selectivity';
        xtitle = 'Time from cue{s}';
        colorRange=[-0.3 0.3];
        plot_selectivity(choicesel,tlabel1,xtitle,colorRange,savefluofigpath)
        
        tlabel1 = 'Outcome selectivity';
        xtitle = 'Time from cue{s}';
        colorRange=[-0.3 0.3];
        plot_selectivity(outcomesel,tlabel1,xtitle,colorRange,savefluofigpath)
        
        save(saveSelName,...
        'choicesel','outcomesel');
%     save(saveSelName,...
%         'choicesel');
    close all;
    end
end