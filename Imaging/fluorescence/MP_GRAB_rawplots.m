function MP_GRAB_rawplots(dataIndex)

% load behavior and dff data
% make some simple plots

% plot trials: first 2-6 trials
% middle 5 trials (50 percentile)
% last 5 trials (not miss)

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    % load dFF files
    load(fullfile(fn_beh.folder,'dff.mat'));
    
    
    savefigpath = fullfile(fn_beh.folder,'figs-fluo','rawplot');
    if ~exist(savefigpath,'dir')
        mkdir(savefigpath);
    end
    plotROIs = [45, 54, 105, 157, 166];
    if ~exist(fullfile(savefigpath,'cell196.png'))
        %% Plot dF/F of all the cells
        numTrials = length(trialData.cue);
        plottrials = [2:6, ceil(numTrials/2):ceil(numTrials/2)+4, numTrials-4:numTrials];
        
        for ii = 1:length(plottrials)
            savefigname = ['rawdFF-trial',num2str(plottrials(ii))];
            t = plottrials(ii);
            startTime = trialData.cueTimes(t)-3; endTime = trialData.cueTimes(t)+5;
            plotIndex = (cells.t<endTime & cells.t>=startTime);
            figure;
            
            for jj = 1:length(plotROIs)
                subplot(5,1,jj)
                plot(cells.t(plotIndex)-trialData.cueTimes(t), cells.dFF{jj}(plotIndex),'black');
                hold on;
                % plot all trial event (lick, reward)
                leftLicks = trialData.leftlickTimes{t};
                rightLicks = trialData.rightlickTimes{t};
                outcome = trialData.outcomeTimes(t)-trialData.cueTimes(t);
                if ~isempty(leftLicks)
                    for uu = 1:length(leftLicks)
                        hold on; plot([leftLicks(uu),leftLicks(uu)],[0.1,0.2],'r');  %lick left
                    end
                end
                if ~isempty(rightLicks)
                    for uu = 1:length(rightLicks)
                        hold on; plot([rightLicks(uu),rightLicks(uu)],[0.1,0.2],'b');  %lick left
                    end
                end
                if trials.reward(t) == 1  % if reward
                    hold on; plot([outcome,outcome],[0.15,0.25],'g');
                end
                % response is not included in the licktimes
                if trialData.response(t) == 2
                    hold on; plot([trialData.rt(t),trialData.rt(t)],[0.1,0.2],'r');  %lick left
                elseif trialData.response(t) == 3
                    hold on; plot([trialData.rt(t),trialData.rt(t)],[0.1,0.2],'b');  %lick right
                end
                set(gca,'box','off');
                xlim([-3 5]);
                ylabel(['ROI #',num2str(plotROIs(jj))]);
            end
            sgtitle(['Trial #', num2str(t)]);
            print(gcf,'-dpng',fullfile(savefigpath,savefigname));
            saveas(gcf, fullfile(savefigpath,savefigname), 'fig');
        end
        close all
    end
end
end
