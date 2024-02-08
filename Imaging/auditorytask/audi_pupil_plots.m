function audi_pupil_plots(dataIndex)
% simple plots

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii},['*',date(1:6),'*_pup.mat']);
    fn_pup = dir(pup_name);
    if length(fn_pup) == 1
        load(fullfile(fn_pup.folder,fn_pup.name));
        
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
        %% Plot pupil for all trials
        
        % MP_plot_pupil( pupil, trialData );
        
        %% Plot cue-aligned pupil
        
        % plot a single trace, create a gif
        
       
        % tCue = trialData.cueTimes(1):trialData.cueTimes(2);
        
      
        % plot the whole session
        figure;

        plot(pupil.t,pupil.dia,'black','LineWidth',1);
        set(gca,'box','off');
        %set(gca,'XColor', 'none', 'YColor','none');
        hold on;
        xlim([-100 800]);
        %plot([5000,5000],[3,2],'k');
        xlabel('Time (s)');
        ylabel('Pupil diameter (z)');
        print(gcf,'-dpng',['whole-trace']);    %png format
        saveas(gcf,['whole-trace'], 'fig');
        saveas(gcf, 'whole-trace','svg');

        %% plot trial averaged pupil trace
        
      
        %% plot pupil change against cue (psth plot)
        savepsthname = fullfile(fn_pup.folder,'psthPup.mat');
        
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from cue (s)';
        params.window = [-1:0.1:3];
        params.minNumTrial = 5;
        params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
        params.CI = 0.95;  %confidence interval
        psth_output=[];
        
     
                % find the shortest timeline
                                
                psth_panel(1).sig{1} = get_psth(pupil.dia, pupil.t, params.trigTime(2:end-1), 'Pupil change', params );
           
        
        
        
        tlabel = 'pupil';
        
        MP_plot_psth(psth_panel,tlabel,params.xtitle);
        print(gcf,'-dpng',['cell_choice' int2str(j)]);
        
%         % plot reward aligned to rewardtime
%         params=[];
%         params.trigTime = trialData.outcomeTimes;
%         params.xtitle = 'Time from outcome (s)';
%         params.window = [-3:0.1:5];
%         params.minNumTrial = 5;
%         params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
%         params.CI = 0.95;  %confidence interval
%         psth_output=[];
%         for k=1:2
%         if k==1 % panel 3
%                 fieldname{1}={'left','reward'}; col{1} = 'r';
%                 fieldname{2}={'left','noreward'}; col{2} = 'b';
%             elseif k==2
%                 fieldname{1}={'right','reward'}; col{1} = 'r';
%                 fieldname{2}={'right','noreward'}; col{2} = 'b';
%             end
%             for kk=1:numel(fieldname)
%                 trialMask = getMask(trials,fieldname{kk});
%                 trialMask = trialMask(trialData.cueTimes<(pupil.t(end)-8));
%                 psth_panel(k).sig{kk} = MP_get_psth(pupil.dia, pupil.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
%                 psth_panel(k).col{kk} = col{kk};
%             end
%         end
%         
%         
%         
%         tlabel = 'pupil';
%         
%         MP_plot_psth(psth_panel,tlabel,params.xtitle);
%         print(gcf,'-dpng',['cell_reward' int2str(j)]);
        
        % plot a figure to show the extremties influence the error bar
        % (bootstrap)
        % MP_plot_extreme(psth_panel, tlabel, params.xtitle);
        save(
        close all;
    end

end

end