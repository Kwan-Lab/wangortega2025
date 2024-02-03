clearvars;
close all;
setup_figprop;

power = '3mW';
region = [''];
root_path = append('C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\preference\', power);

% perference  behavior
%    OUTCOME.REWARDLEFT = 5;
%     OUTCOME.REWARDRIGHT = 6;
%     OUTCOME.REWARDMANUAL = 7;
%     OUTCOME.NOREWARDLEFT = 75;
%     OUTCOME.NOREWARDRIGHT = 76;
%     OUTCOME.MISS = 8;      %miss
%     
%     RULE.LStimRew_RRew = 41; %Left stim + reward ; right reward only
%     RULE.LRew_RStimRew = 42;

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('--------------------------------------------------open---------');

% Look for data files and create a database index

% create data indexes 

    logfilepath = fullfile(root_path,region,'data');
    analysispath = fullfile(root_path, region, 'analysis');
    dataIndex = makeDataIndex(logfilepath, analysispath);
    
    % Parse and analyze each logfile, save as .mat files, as needed
    dataIndex = createBehMatFiles(dataIndex);
    
    % sort Index according to experiment data
    
    dataIndex = sortdataIndex(dataIndex);

for r = 1:size(dataIndex, 1)
    LogFileName = dataIndex.LogFileName{r}; 
    load(fullfile(dataIndex.BehPath{r},[LogFileName(1:end-4),'_beh.mat']));
    logfile = dataIndex.LogFileName{r}(1:end-4);

    trialData.response = trialData.response(~(trialData.response == 0)); %remove misses
    side.subj = sessionData.subject;
    side.strain = dataIndex.Strain{r};
    side.stimRew = 0;
    side.rew = 0; 
    side.stay = 0;
    side.consec = 0;
    side.switch = 0;
    side.totalConsec = 0;
    side.rule = 0;
    side.nTrial = size(trialData.response,1);
    side.choiceWaterStim = 0;
    side.left = sum(trialData.response == 2);
    side.right = sum(trialData.response == 3);
    side.waterSwitch = 0;
    side.waterStimSwitch = 0; 
    side.nrSwitch = 0;
    side.rewSwitch = 0;

    for k = 1:size(trialData.response,1)-1

        if trialData.rule(k) == 41 %Left = stim + reward, Right = reward only 
            if trialData.response(k) == 2 %left
                side.stimRew = side.stimRew + 1;
            elseif trialData.response(k) == 3 %right
                side.rew = side.rew + 1;
            end
            if trialData.outcome(k) == 5 %reward
                if trialData.response(k+1) == 3
                    side.waterStimSwitch = side.waterStimSwitch + 1;
                    side.rewSwitch = side.rewSwitch + 1;
                end 
            elseif trialData.outcome(k) == 75 %no reward left
                 if trialData.response(k+1) == 3 %right
                    side.nrSwitch = side.nrSwitch + 1;
                end
            end
        elseif trialData.rule(k) == 42
            if trialData.response(k) == 3 %right
                side.stimRew = side.stimRew + 1;
            elseif trialData.response(k) == 2 %left
                side.rew = side.rew + 1;
            end
            if trialData.outcome(k) == 6
                if trialData.response(k+1) == 2 %left
                    side.waterSwitch = side.waterSwitch + 1;
                    side.rewSwitch = side.rewSwitch + 1;
                end
            elseif trialData.outcome(k) == 76 %no reward right
                if trialData.response(k+1) == 2 %left
                    side.nrSwitch = side.nrSwitch + 1;
                end
            end

        end

    end

        nTrials = 200;
    start = size(trialData.stimulationRegion,1) - nTrials; %analyze last nTrials  
    if start <= 0
        start = 1;
        nTrials = size(trialData.stimulationRegion,1);
    end 
    for m = start:(size(trialData.stimulationRegion,1)-1)

        if trialData.response(m) == trialData.response(m+1)
            side.consec = side.consec + 1; 
        elseif trialData.response(m) ~= trialData.response(m+1)
            side.switch = side.switch + 1;
            idx = cell2mat({side.switch});
            side.totalConsec(idx) = side.consec; %save n of consecutive choices
            side.rule(idx) = trialData.rule(m); %what is the rule for these consec choices
            side.consec = 0; %reset consecutive choice counter

        end

        if trialData.response(m) == 2 && trialData.rule(m) == 41 %left choice and left is water + stim
            side.choiceWaterStim = side.choiceWaterStim + 1;
        elseif trialData.response(m) == 3 && trialData.rule(m) == 42 %right choice and right is water + stim
            side.choiceWaterStim = side.choiceWaterStim + 1;
        end 
    end
    
    side.preferStim = side.choiceWaterStim / nTrials ; 
    preference{r} = side;
end 

    cd 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\preference\'
    filename = append('preference_', power, '.mat');
    save(filename,'preference')
%%
cd 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\preference\'

for tt = 1:2
    if tt == 1
        load preference_0mW.mat
    
    elseif tt == 2
            load preference_3mW.mat
    end 
    
    d = 1;
    c = 1; 
    temp = [];
    for j = 1:size(preference,2)
        if strcmp(preference{j}.strain, 'Dbh') == 1
            temp.nTrials.dbh (d) = preference{j}.nTrial;
            temp.plot.dbh (d) = preference{j}.preferStim;
            temp.pSwitches.dbh (d) = preference{j}.switch / preference{j}.nTrial;
            temp.nSwitches.dbh (d) = preference{j}.switch;
            temp.pRight.dbh (d) = preference{j}.right / (preference{j}.left + preference{j}.right);
            temp.pWaterSwitch.dbh (d) = preference{j}.waterSwitch / preference{j}.nTrial;
            temp.pWaterStimSwitch.dbh (d) = preference{j}.waterStimSwitch / preference{j}.nTrial;
            temp.pnrSwitch.dbh (d) = preference{j}.nrSwitch / preference{j}.nTrial;
            temp.pRewSwitch.dbh (d) = preference{j}.rewSwitch / preference{j}.nTrial;
                d = d+1;
        elseif strcmp(preference{j}.strain, 'ChAt') == 1
            temp.nTrials.chat (c) = preference{j}.nTrial;
            temp.plot.chat (c) = preference{j}.preferStim;
            temp.pSwitches.chat (c) = preference{j}.switch / preference{j}.nTrial;
            temp.nSwitches.chat (c) = preference{j}.switch;
            temp.pRight.chat (c) = preference{j}.right / (preference{j}.left + preference{j}.right);
            temp.pWaterSwitch.chat (c) = preference{j}.waterSwitch / preference{j}.nTrial;
            temp.pWaterStimSwitch.chat (c) = preference{j}.waterStimSwitch / preference{j}.nTrial;
            temp.pnrSwitch.chat (c) = preference{j}.nrSwitch / preference{j}.nTrial;
            temp.pRewSwitch.chat (c) = preference{j}.rewSwitch / preference{j}.nTrial;
                c = c+1;
        end 
    end
    dat{tt} = temp; 
end 
   
        zero = dat{1};
        three = dat{2};

        cc = [0.9290 0.6940 0.1250];
        dc = [.1 0.6 0.6];
        g = [.7 .7 .7]; 
        sz = 500;
        xoffset_d0 = 0.3*rand(1, size(zero.plot.dbh,2));
        xoffset_c0 = 0.3*rand(1, size(zero.plot.chat,2));
        d0 = size(zero.plot.dbh,2);
        c0 = size(zero.plot.chat,2);

        xoffset_d3 = 0.3*rand(1, size(three.plot.dbh,2));
        xoffset_c3 = 0.3*rand(1, size(three.plot.chat,2));
        d3 = size(three.plot.dbh,2);
        c3 = size(three.plot.chat,2);


setup_figprop
fig = figure ('position', [0 0 1500 600])

subplot(1,4,1)
hold on
scatter(ones(1,d3)+ xoffset_d3, three.plot.dbh, sz, dc)
scatter(1.1, mean(three.plot.dbh), sz, dc, 'filled')
errorbar(1.1, mean(three.plot.dbh, 'omitnan'), std(three.plot.dbh) / sqrt(size(three.plot.dbh,2))*2,...
                         'color', dc);
scatter(2.5*ones(1,c3)+ xoffset_c3, three.plot.chat, sz, cc)
scatter(2.6, mean(three.plot.chat), sz, cc, "filled")
errorbar(2.6, mean(three.plot.chat, 'omitnan'), std(three.plot.chat) / sqrt(size(three.plot.chat,2))*2,...
                         'color', cc);

scatter(ones(1,d0)*1.5 + xoffset_d0, zero.plot.dbh, sz, g)
scatter(1.6, mean(zero.plot.dbh), sz, g, 'filled')
errorbar(1.6, mean(zero.plot.dbh, 'omitnan'), std(zero.plot.dbh) / sqrt(size(zero.plot.dbh,2))*2,...
                         'color', g);
scatter(3*ones(1,c0)+ xoffset_c0, zero.plot.chat, sz, g)
scatter(3.1, mean(zero.plot.chat), sz, g, "filled")
errorbar(3.1, mean(zero.plot.chat, 'omitnan'), std(zero.plot.chat) / sqrt(size(zero.plot.chat,2))*2,...
                         'color', g);

xlim([.5 3.5])
ylim([0 1])
%title('p(water + stimulation)')
xticks([])
%xticklabels(['Dbh '; 'Chat'])
yticks([0 .5 1])
H = gca;
H.FontSize = 20;
H.LineWidth = 3;
yline(.5,"LineStyle","--")
title("Preference for side with Stimulation")

h = ttest2(zero.plot.chat, three.plot.chat);
if h == 1
       plot([2.5, 3], [.8, .8], '-k')
end 

h = ttest2(zero.plot.dbh, three.plot.dbh);
if h == 1
       plot([1, 1.5], [.8, .8], '-k')
end 

subplot(1,4,2)
hold on
scatter(ones(1,d3) + xoffset_d3, three.pSwitches.dbh, sz, dc)
scatter(1.1, mean(three.pSwitches.dbh), sz, dc, 'filled')
scatter(2.5*ones(2,c3)+ xoffset_c3, three.pSwitches.chat, sz, cc)
scatter(2.6, mean(three.pSwitches.chat), sz, cc, 'filled')

scatter(ones(1,d0)*1.5 + xoffset_d0, zero.pSwitches.dbh, sz, g)
scatter(1.6, mean(zero.pSwitches.dbh), sz, g, 'filled')
scatter(3*ones(2,c0)+ xoffset_c0, zero.pSwitches.chat, sz, g)
scatter(3.1, mean(zero.pSwitches.chat), sz, g, 'filled')

errorbar(1.1, mean(three.pSwitches.dbh, 'omitnan'), std(three.pSwitches.dbh) / sqrt(size(three.pSwitches.dbh,2))*2,...
                         'color', dc);
errorbar(2.6, mean(three.pSwitches.chat, 'omitnan'), std(three.pSwitches.chat) / sqrt(size(three.pSwitches.chat,2))*2,...
                         'color', cc);
errorbar(1.6, mean(zero.pSwitches.dbh, 'omitnan'), std(zero.pSwitches.dbh) / sqrt(size(zero.pSwitches.dbh,2))*2,...
                         'color', g);
errorbar(3.1, mean(zero.pSwitches.chat, 'omitnan'), std(zero.pSwitches.chat) / sqrt(size(zero.pSwitches.chat,2))*2,...
                         'color', g);

h = ttest2(zero.pSwitches.chat, three.pSwitches.chat);
if h == 1
       plot([2.5, 3], [.8, .8], '-k')
end 

h = ttest2(zero.pSwitches.dbh, three.pSwitches.dbh);
if h == 1
       plot([1, 1.5], [.8, .8], '-k')
end 


xlim([.5 3.5])
ylim([0 1])
xticks([])
%xticklabels(['Dbh '; 'Chat'])
yticks([0 .5 1])
H = gca;
H.FontSize = 20;
H.LineWidth = 3;
title('Probability of Switch')

subplot(1,4,3)
hold on
scatter(ones(1,d3)+ xoffset_d3, three.nSwitches.dbh, sz, dc)
scatter(1.1, mean(three.nSwitches.dbh), sz, dc, 'filled')
scatter(2.5*ones(2,c3)+ xoffset_c3, three.nSwitches.chat, sz, cc)
scatter(2.6, mean(three.nSwitches.chat), sz, cc, 'filled')

scatter(ones(1,d0)*1.5+ xoffset_d0, zero.nSwitches.dbh, sz, g)
scatter(1.6, mean(zero.nSwitches.dbh), sz, g, 'filled')
scatter(3*ones(2,c0)+ xoffset_c0, zero.nSwitches.chat, sz, g)
scatter(3.1, mean(zero.nSwitches.chat), sz, g, 'filled')

errorbar(1.1, mean(three.nSwitches.dbh, 'omitnan'), std(three.nSwitches.dbh) / sqrt(size(three.nSwitches.dbh,2))*2,...
                         'color', dc);
errorbar(2.6, mean(three.nSwitches.chat, 'omitnan'), std(three.nSwitches.chat) / sqrt(size(three.nSwitches.chat,2))*2,...
                         'color', cc);
errorbar(1.6, mean(zero.nSwitches.dbh, 'omitnan'), std(zero.nSwitches.dbh) / sqrt(size(zero.nSwitches.dbh,2))*2,...
                         'color', g);
errorbar(3.1, mean(zero.nSwitches.chat, 'omitnan'), std(zero.nSwitches.chat) / sqrt(size(zero.nSwitches.chat,2))*2,...
                         'color', g);

h = ttest2(zero.nSwitches.chat, three.nSwitches.chat);
if h == 1
       plot([2.5, 3], [85, 85], '-k')
end 

h = ttest2(zero.nSwitches.dbh, three.nSwitches.dbh);
if h == 1
       plot([1, 1.5], [85, 85], '-k')
end 

xlim([.5 3.5])
ylim([0 100])
xticks([])
%xticklabels(['Dbh '; 'Chat'])
yticks([0 50 100])
H = gca;
H.FontSize = 20;
H.LineWidth = 3;
title('Number of Switches')

subplot(1,4,4)
hold on
scatter(ones(1,d3)+ xoffset_d3, three.nTrials.dbh, sz, dc)
scatter(1.1, mean(three.nTrials.dbh), sz, dc, 'filled')
errorbar(1.1, mean(three.nTrials.dbh, 'omitnan'), std(three.nTrials.dbh) / sqrt(size(three.nTrials.dbh,2)),...
                         'color', dc);
scatter(2.5*ones(1,c3)+ xoffset_c3, three.nTrials.chat, sz, cc)
scatter(2.6, mean(three.nTrials.chat), sz, cc, "filled")
errorbar(2.6, mean(three.nTrials.chat, 'omitnan'), std(three.nTrials.chat) / sqrt(size(three.nTrials.chat,2)),...
                         'color', cc);

scatter(ones(1,d0)*1.5 + xoffset_d0, zero.nTrials.dbh, sz, g)
scatter(1.6, mean(zero.nTrials.dbh), sz, g, 'filled')
errorbar(1.6, mean(zero.nTrials.dbh, 'omitnan'), std(zero.nTrials.dbh) / sqrt(size(zero.nTrials.dbh,2)),...
                         'color', g);
scatter(3*ones(1,c0)+ xoffset_c0, zero.nTrials.chat, sz, g)
scatter(3.1, mean(zero.nTrials.chat), sz, g, "filled")
errorbar(3.1, mean(zero.nTrials.chat, 'omitnan'), std(zero.nTrials.chat) / sqrt(size(zero.nTrials.chat,2)),...
                         'color', g);

xlim([.5 3.5])
ylim([0 400])
%title('p(water + stimulation)')
xticks([])
%xticklabels(['Dbh '; 'Chat'])
yticks([0 100 200 300 400])
H = gca;
H.FontSize = 20;
H.LineWidth = 3;
yline(.5,"LineStyle","--")
title('Number of Responses')

h = ttest2(zero.nTrials.chat, three.nTrials.chat);
if h == 1
       plot([2.5, 3], [.8, .8], '-k')
end 

h = ttest2(zero.nTrials.dbh, three.nTrials.dbh);
if h == 1
       plot([1, 1.5], [.8, .8], '-k')
end 
cd C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\preference\3mW

        saveas(fig, 'preference', 'tif')
        saveas(fig, 'preference', 'fig')


        figure 
        subplot(1,4,1)
        hold on
        scatter(ones(1,d3)+ xoffset_d3, three.pWaterSwitch.dbh, sz, dc)
        scatter(1.1, mean(three.pWaterSwitch.dbh), sz, dc, "filled")
        errorbar(1.1, mean(three.pWaterSwitch.dbh, 'omitnan'), std(three.pWaterSwitch.dbh) / sqrt(size(three.pWaterSwitch.dbh,2)),...
                         'color', dc);
        scatter(2*ones(1,d0)+ xoffset_d0, zero.pWaterSwitch.dbh, sz, g)
        scatter(2.1, mean(zero.pWaterSwitch.dbh), sz, g, "filled")
        errorbar(2.1, mean(zero.pWaterSwitch.dbh, 'omitnan'), std(zero.pWaterSwitch.dbh) / sqrt(size(zero.pWaterSwitch.dbh,2)),...
                         'color', g);

        scatter(4*ones(1,c3)+ xoffset_c3, three.pWaterSwitch.chat, sz, cc)
        scatter(4.1, mean(three.pWaterSwitch.chat), sz, cc, "filled")
        errorbar(4.1, mean(three.pWaterSwitch.chat, 'omitnan'), std(three.pWaterSwitch.chat) / sqrt(size(three.pWaterSwitch.chat,2)),...
                         'color', cc);
        scatter(5*ones(1,c0)+ xoffset_c0, zero.pWaterSwitch.chat, sz, g)
        scatter(5.1, mean(zero.pWaterSwitch.chat), sz, g, "filled")
        errorbar(5.1, mean(zero.pWaterSwitch.chat, 'omitnan'), std(zero.pWaterSwitch.chat) / sqrt(size(zero.pWaterSwitch.chat,2)),...
                         'color', g);
        ylim([0 .2])
        xlim([0 6])
        yticks([ 0 .1 .2])
        xticks([])
        title('probability to switch after water reward')


        subplot(1,4,2)
        hold on
        scatter(ones(1,d3)+ xoffset_d3, three.pWaterStimSwitch.dbh, sz, dc)
        scatter(1.1, mean(three.pWaterStimSwitch.dbh), sz, dc, "filled")
        errorbar(1.1, mean(three.pWaterStimSwitch.dbh, 'omitnan'), std(three.pWaterStimSwitch.dbh) / sqrt(size(three.pWaterStimSwitch.dbh,2)),...
                         'color', dc);
        scatter(2*ones(1,d0)+ xoffset_d0, zero.pWaterStimSwitch.dbh, sz, g)
        scatter(2.1, mean(zero.pWaterStimSwitch.dbh), sz, g, "filled")
        errorbar(2.1, mean(zero.pWaterStimSwitch.dbh, 'omitnan'), std(zero.pWaterStimSwitch.dbh) / sqrt(size(zero.pWaterStimSwitch.dbh,2)),...
                         'color', g);
        scatter(4*ones(1,c3)+ xoffset_c3, three.pWaterStimSwitch.chat, sz, cc)
        scatter(4.1, mean(three.pWaterStimSwitch.chat), sz, cc, "filled")
        errorbar(4.1, mean(three.pWaterStimSwitch.chat, 'omitnan'), std(three.pWaterStimSwitch.chat) / sqrt(size(three.pWaterStimSwitch.chat,2)),...
                         'color', cc);
        scatter(5*ones(1,c0)+ xoffset_c0, zero.pWaterStimSwitch.chat, sz, g)
        scatter(5.1, mean(zero.pWaterStimSwitch.chat), sz, g, "filled")
        errorbar(5.1, mean(zero.pWaterStimSwitch.chat, 'omitnan'), std(zero.pWaterStimSwitch.chat) / sqrt(size(zero.pWaterStimSwitch.chat,2)),...
                         'color', g);

        ylim([0 .2])
        xlim([0 6])
        yticks([ 0 .1 .2])
        xticks([])
        title('probability to switch after water+stim reward')

        subplot(1,4,3)
        hold on
        scatter(ones(1,d3)+ xoffset_d3, three.pnrSwitch.dbh, sz, dc)
        scatter(1.1, mean(three.pnrSwitch.dbh), sz, dc, "filled")
        errorbar(1.1, mean(three.pnrSwitch.dbh, 'omitnan'), std(three.pnrSwitch.dbh) / sqrt(size(three.pnrSwitch.dbh,2)),...
                         'color', dc);
        
        scatter(2*ones(1,d0)+ xoffset_d0, zero.pnrSwitch.dbh, sz, g)
        scatter(2.1, mean(zero.pnrSwitch.dbh), sz, g, "filled")
        errorbar(2.1, mean(zero.pnrSwitch.dbh, 'omitnan'), std(zero.pnrSwitch.dbh) / sqrt(size(zero.pnrSwitch.dbh,2)),...
                         'color', g);


        scatter(4*ones(1,c3)+ xoffset_c3, three.pnrSwitch.chat, sz, cc)
        scatter(4.1, mean(three.pnrSwitch.chat), sz, cc, "filled")
        errorbar(4.1, mean(three.pnrSwitch.chat, 'omitnan'), std(three.pnrSwitch.chat) / sqrt(size(three.pnrSwitch.chat,2)),...
                         'color', cc);

        scatter(5*ones(1,c0)+ xoffset_c0, zero.pnrSwitch.chat, sz, g)
        scatter(5.1, mean(zero.pnrSwitch.chat), sz, g, "filled")
        errorbar(5.1, mean(zero.pnrSwitch.chat, 'omitnan'), std(zero.pnrSwitch.chat) / sqrt(size(zero.pnrSwitch.chat,2)),...
                         'color', g);

        ylim([0 .2])
        xlim([0 6])
        yticks([ 0 .1 .2])
        xticks([])
        title('probability to switch after no reward')

        subplot(1,4,4)
        hold on
        scatter(ones(1,d3)+ xoffset_d3, three.pRewSwitch.dbh, sz, dc)
        scatter(1.1, mean(three.pRewSwitch.dbh), sz, dc, "filled")
        errorbar(1.1, mean(three.pRewSwitch.dbh, 'omitnan'), std(three.pRewSwitch.dbh) / sqrt(size(three.pRewSwitch.dbh,2)),...
                         'color', dc);
        scatter(2*ones(1,d0)+ xoffset_d0, zero.pRewSwitch.dbh, sz, g)
        scatter(2.1, mean(zero.pRewSwitch.dbh), sz, g, "filled")
        errorbar(2.1, mean(zero.pRewSwitch.dbh, 'omitnan'), std(zero.pRewSwitch.dbh) / sqrt(size(zero.pRewSwitch.dbh,2)),...
                         'color', g);

        scatter(4*ones(1,c3)+ xoffset_c3, three.pRewSwitch.chat, sz, cc)
        scatter(4.1, mean(three.pRewSwitch.chat), sz, cc, "filled")
        errorbar(4.1, mean(three.pRewSwitch.chat, 'omitnan'), std(three.pRewSwitch.chat) / sqrt(size(three.pRewSwitch.chat,2)),...
                         'color', cc);
        scatter(5*ones(1,c0)+ xoffset_c0, zero.pRewSwitch.chat, sz, g)
        scatter(5.1, mean(zero.pRewSwitch.chat), sz, g, "filled")
        errorbar(5.1, mean(zero.pRewSwitch.chat, 'omitnan'), std(zero.pRewSwitch.chat) / sqrt(size(zero.pRewSwitch.chat,2)),...
                         'color', g);

        ylim([0 .2])
        xlim([0 6])
        yticks([ 0 .1 .2])
        xticks([])
        title('probability to switch after any reward')

% figure
% scatter(ones(1,size(three.pRight.dbh,2))+ xoffset_d3, three.pRight.dbh, sz, dc)
% hold on
% scatter(4*ones(1,size(three.pRight.chat,2))+ xoffset_c3, three.pRight.chat, sz, cc)
% 
% scatter(ones(1,size(zero.pRight.dbh,2))*1.5+ xoffset_d0, zero.pRight.dbh, sz, g)
% 
% scatter(4.5*ones(1,size(zero.pRight.chat,2))+ xoffset_c0, zero.pRight.chat, sz, g)
% 
% ylim([0 1])
% xlim([0 5])
% xticks([1 4])
% xticklabels(["Dbh"; "ChAt"])
% ylabel("Right Preference")
% 
% 
%         saveas(fig, 'sidePreference', 'tif')
%         saveas(fig, 'sidePreference', 'fig')
%  
