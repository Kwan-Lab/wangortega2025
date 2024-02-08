function MP_GRAB_simpleplots(dataIndex)

% load behavior and dff data
% make some simple plots

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
   fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));

    savefigpath = fullfile(fn_beh.folder,'figs-fluo','PSTH');
    if ~exist(savefigpath,'dir')
        mkdir(savefigpath);
    end
    
    if ~exist(fullfile(savefigpath,'cell196.png')) 
        
         
        load(fullfile(fn_beh.folder,fn_beh.name));
        % load dFF files
        load(fullfile(fn_beh.folder,'dff.mat'));
    
    
        %% Plot dF/F of all the cells
%         MP_plot_allCells( cells, trialData ,savefigpath);
        
        close all;
        %% check correaltion between different cells
%         autoCorr = zeros(length(cells.dFF));
%         for uu = 1:length(cells.dFF)
%             for vv = 1:length(cells.dFF)
%                 autoCorr(uu,vv) = corr(cells.dFF{uu}, cells.dFF{vv});
%             end
%         end
%         figure;
%         heatmap(autoCorr);
%         print(gcf,'-dpng',fullfile(savefigpath,'corr-cells'));    %png format
%         saveas(gcf, fullfile(savefigpath,'corr-cells'), 'fig');
        
        %% Plot cue-aligned dF/F for each cell
        % Fig. 3d in paper came from 140605 data set, cell 8 10 37 74
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-3:0.1:5];
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = 0.95;  %confidence interval
        params.minNumTrial = 50;
        %for r=1:numel(cells.dFF)
            % for NE: j = 710
            % for ACh: j = 33
            %j=710;  % NE
             j = 33;
            psth_output=[];
            for k=1:2
                fieldname=[];
                if k==1 %panel 1
                    fieldname{1}={'reward','contra'};
                    fieldname{2}={'noreward','contra'};
                else
                    fieldname{1}={'reward','ipsi'};
                    fieldname{2}={'noreward','ipsi'};
                   
                end
                for kk=1:numel(fieldname)
                    trialMask = getMask_GRAB(trials,fieldname{kk},dataIndex.RecordingSite{ii});
                 
                    psth_panel{k}.sig{kk} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
               
                end
            end
            
            if cells.isRedCell{j}
                tlabel = ['RED Cell ' int2str(j)];
            else
                tlabel = ['Cell ' int2str(j)];
            end
            plot_psth(psth_panel,tlabel,params.xtitle);
            print(gcf,'-dpng',fullfile(savefigpath,['cell' int2str(j)]));
            saveas(gcf, fullfile(savefigpath,['cell' int2str(j)]), 'fig');
            saveas(gcf, fullfile(savefigpath,['cell' int2str(j)]),'svg');
            close;
        end
        
        %% normalized df/f
%         if isfield(cells,'normdFF')
%             for j=1:numel(cells.dFF)
%             psth_output=[];
%             for k=1:2
%                 fieldname=[];
%                 if k==1 %panel 1
%                     fieldname{1}={'left','reward'};
%                     fieldname{2}={'left','noreward'};
%                 elseif k==2 %panel 2
%                     fieldname{1}={'right','reward'};
%                     fieldname{2}={'right','noreward'};
%                 end
%                 for kk=1:numel(fieldname)
%                     trialMask = getMask(trials,fieldname{kk});
%                     psth_panel(k).sig{kk} = get_psth( cells.normdFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
%                 end
%             end
%             
%             if cells.isRedCell{j}
%                 tlabel = ['RED Cell ' int2str(j)];
%             else
%                 tlabel = ['Cell ' int2str(j)];
%             end
%             plot_psth(psth_panel,tlabel,params.xtitle);
%             print(gcf,'-dpng',fullfile(savefigpath,['normcell' int2str(j)]));
%             close;
%             end
%         end
    end
    % plot fluorescent signal align to left/right lick
%     params_left = [];
%     params_right = [];
%     params_left.window = [-3:0.1:5];
%     params_right.window = [-3:0.1:5];
%     params_left.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
%     params_right.numBootstrapRepeat = 1000;
%     params_right.CI = 0.95;  %confidence interval
%     params_left.CI = 0.95;
%     params.xtitle = 'Time from lick (s)';
%     responseTimes = trialData.cueTimes + trialData.rt;
%     
%         % only include licks without outcome and reward (missed trials and
%         % in intertrial interval?
%     params_left.trigTime = setdiff(sessionData.lickTimes{1}, responseTimes);
%     params_right.trigTime = setdiff(sessionData.lickTimes{2}, responseTimes);
%     fieldname{1}={'left lick'};
%     fieldname{2}={'right lick'};
%      for j=1:numel(cells.dFF)
%         psth_output=[];
%        
%         psth_lick(1).sig{1} = get_psth( cells.dFF{j}, cells.t, params_left.trigTime, strjoin(fieldname{1}), params_left );
%         psth_lick(1).sig{2} = get_psth( cells.dFF{j}, cells.t, params_right.trigTime, strjoin(fieldname{2}), params_right );
%         
%         if cells.isRedCell{j}
%             tlabel = ['RED Cell' int2str(j)];
%         else
%             tlabel = ['Cell' int2str(j)];
%         end
%         plot_psth(psth_lick,tlabel,params.xtitle);
%         print(gcf,'-dpng',fullfile(savefigpath,['cell-lick' int2str(j)]));
%         close;
%     end
        
end