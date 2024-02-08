function MP_GRAB_selectivity_px(dataIndex);

% calculate choice and reward selectivity
% choice_s = (df/f_L-df/f_R)/(df/f_L+df/f_R);
% reward_s = (df/f_R-df/f_UR)/df/f_R+df/f_UP);

nFiles = size(dataIndex,1);




%% go through every session, running RF separately
for ii = 1:1 % only calculate px choice selecvitity first
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    
    % load fluorescent files
    fn_fluo = "V:\HongliWang\GRAB_Ach\Ach_MP\910-MP-1204\mat";
    if isdir(fn_fluo) == 1
        
        savefluofigpath = fullfile(dataIndex.BehPath{ii},'figs-fluo');
        if ~exist(savefluofigpath,'dir')
            mkdir(savefluofigpath);
        end
        
        cd(savefluofigpath);
        
      
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        saveSelName = fullfile(savematpath,'select_px.mat');  % random forest
        
        %% load mat files, calculate df/f first
          stacks = dir(fullfile(fn_fluo,'*.mat'));
          % get number of pixels
         pic=load(fullfile(stacks(1).folder,stacks(1).name));
         
         % set ROIMask to every pixel
      
        f=cell(size(pic.stack,1)*size(pic.stack,2),1); n=cell(size(pic.stack,1)*size(pic.stack,1),1);
            for j=1:numel(stacks)
                disp(['Loading reg image file ' stacks(j).name]);
              
                pic=load(fullfile(stacks(j).folder,stacks(1).name));
                [nY nX nZ]=size(pic.stack);
                
                for k=1:size(pic,1)*size(pic,1)
                    tempf=[]; tempn=[];
                    for i=1:1:nZ
                        %get sum of pixel values within the ROI
                        tempf(i)=sum(sum(pic.stack(:,:,i).*uint16(shifted_cellmask(:,:,k))));
                        tempn(i)=sum(sum(pic.stack(:,:,i).*uint16(shifted_neuropilmask(:,:,k))));
                    end
                    if sum(sum(shifted_cellmask(:,:,k)))>0     %if there are pixels belonging the the ROI
                        if j==1 %if this is the first reg image, then start new variables
                            f{k}=tempf/sum(sum(shifted_cellmask(:,:,k)));   %per-pixel fluorescence
                            n{k}=tempn/sum(sum(shifted_neuropilmask(:,:,k)));   %per-pixel fluorescence
                        else
                            f{k}=[f{k} tempf/sum(sum(shifted_cellmask(:,:,k)))];   %per-pixel fluorescence
                            n{k}=[n{k} tempn/sum(sum(shifted_neuropilmask(:,:,k)))];   %per-pixel fluorescence
                        end
                    else %if the ROI is outside of the imaging field of view
                        f{k}=nan(size(tempf));
                        n{k}=nan(size(tempn));
                    end
                end
                clear pic;
            end
    end
        
        if ~exist(saveSelName)
            params=[];
            params.trigTime = trialData.cueTimes;
            params.xtitle = 'Time from cue (s)';
            params.window = [-3:0.1:5];
            params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
            params.CI = 0.95;  %confidence interval
            params.minNumTrial = 50;
        
                       
            if isfield(cells,'normdFF')
                sel_t= cells.t; sel_dFF = cells.normdFF;
            else
                 sel_t= cells.t; sel_dFF = cells.dFF;
            end
            
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
    
    close all;
    end
end