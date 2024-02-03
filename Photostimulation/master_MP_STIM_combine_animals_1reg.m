clearvars;
close all;
setup_figprop;
regionList = ['LM2'; 'RM2'; 'LV1'; 'RV1'];


sidePreference = 0; %1 to analyze side preference by session

for k = 1:size(regionList,1)
        currReg = regionList(k,:);
    root_path = 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\1Region';

    sz = [0 4];
        varTypes = ["string","double","double","double"];
        varNames = ["logFileName","Entropy","nTrials","right/total"];
        lowEntropyTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    
    % matching pennies behavior
    
    disp('-----------------------------------------------------------');
    disp('--- HERE WE GO ');
    disp('--------------------------------------------------open---------');
    
    % Look for data files and create a database index
    
    % create data indexes 
    
        logfilepath = fullfile(root_path,currReg,'data');
        analysispath = fullfile(root_path, currReg, 'analysis');
        dataIndex = makeDataIndex(logfilepath, analysispath);
        
        % Parse and analyze each logfile, save as .mat files, as needed
        dataIndex = MP_STIM_createBehMatFiles(dataIndex);
        
        % sort Index according to experiment data
        
        dataIndex = sortdataIndex(dataIndex);
        
        % Determine if each session fulfill performance criteria (Matching Pennies)
        MP_determineBehCriteria(dataIndex);
    %
        % start separate analysis for each strain
        clear sessionData stats trialData
    
        nFiles = size(dataIndex(:,:),1);
        for ii = 1:nFiles
            savematpath = dataIndex.BehPath{ii};
             [lowEntropyTable] = MP_GRAB_session(dataIndex.BehPath{ii},dataIndex.LogFileName{ii},savematpath, lowEntropyTable);
        end
        save_path = fullfile(root_path,currReg,'summary','figs_summary');
        
    
        MP_GRAB_behaviorPerAnimal(dataIndex,save_path);
        
        MP_GRAB_behaviorAll(dataIndex, save_path);
    %
        % Calculate stats by stimulation region - by session
        % Animal
        regions = {'STIM', 'CTRL'};
        plotChoice = 0; %make 1 to plot session data
        for i = 1:2 
            if i ==1
                dataIndex_curr = dataIndex(strcmp(dataIndex.Strain, 'Dbh'),:);
                animalList = unique(dataIndex_curr.Animal); 
                for j = 1:size(animalList,1) %change to loop through animal then sessions, make animals separate rows
                    animal = animalList(j); 
                    tempIndex = dataIndex_curr(strcmp(dataIndex_curr.Animal, animal),:);
                    stats_dbh{1,j} = MP_STIM_stats(analysispath,tempIndex,plotChoice, regions, currReg, sidePreference);    
                end
                clear sessionData stats trialData
            else
                dataIndex_curr = dataIndex(strcmp(dataIndex.Strain, 'ChAt'),:);
                animalList = unique(dataIndex_curr.Animal); 
                for j = 1:size(animalList,1) %change to loop through animal then sessions, make animals separate rows
                    animal = animalList(j); 
                    tempIndex = dataIndex_curr(strcmp(dataIndex_curr.Animal, animal),:);
                        stats_chat{1,j} = MP_STIM_stats(analysispath,tempIndex,plotChoice, regions, currReg, sidePreference);
                end
                clear sessionData stats trialData
            end
        end
        % start plotting
    
        variables = fieldnames(stats_chat{1, 1}(1).beh_output_WSLS{1, 1}); 
        %regions = {'LM2', 'RM2','LV1','RV1','CTRL'};
    for j = 1:2 %iterate through genotypes
        if j == 1
            dat = stats_chat;
        else 
            dat = stats_dbh;
        end
    
        mat = {}; 
            for i = 1:size(dat,2) %iterate through animals
                animalStruct = dat{1,i};
                for m = 1:size(regions,2) %iterate through regions
                    regionStruct = animalStruct(1,m);
                    cRegion = regions{m};
                    sessionsStruct = regionStruct.beh_output_WSLS;
                    for k = 1:size(variables,1)
                        cVar = variables{k};
                        for yy = 1:size(dat{1, i}(1).beh_output_WSLS,2) %iterate through sessions
                            temp.(cVar)(1,yy) = sessionsStruct{1,yy}.(cVar);
    %                         ntrial(1,yy) = sessionsStruct{1,yy}.nTrial;
    %                         nLwin(1,yy) = sessionsStruct{1,yy}.nLWin;
    %                         nRwin(1,yy) = sessionsStruct{1,yy}.nRwin;
    %                         nLLose(1,yy) = sessionsStruct{1,yy}.nLLose;
    %                         nRLose (1,yy) = sessionsStruct{1,yy}.nRLose;
                        end
    %                     sumSessions = sum(temp.(cVar));
    %                     mat = setfield(mat, cVar, cRegion, {1,i}, sumSessions);
    %                     clear sumSessions 
                    end
                        mat.stats.(cRegion)(1,i) = temp;
                        mat.nTrial.(cRegion){1,i} = temp.left + temp.right;
                        mat.left.(cRegion){1,i} = temp.left;
                        mat.right.(cRegion){1,i} = temp.right;
                        mat.miss.(cRegion){1,i} = temp.miss;
                        mat.nwin.(cRegion){1,i} = temp.nWin;
                        mat.pStay.(cRegion)(1,i) = sum(temp.pStay) / sum(temp.nTrial); %(sum(temp.LWinStay) + sum(temp.LLoseStay) + sum(temp.RWinStay) + sum(temp.RLoseStay)) / (sum(temp.nTrial));
                        mat.pSwitch.(cRegion)(1,i) = sum(temp.pSwitch) / sum(temp.nTrial); %(sum(temp.LWinSwitch) + sum(temp.LLoseSwitch) + sum(temp.RWinSwitch) + sum(temp.RLoseSwitch)) / (sum(temp.nTrial));
                        
                        mat.prewNMinusOneStay.(cRegion)(1,i) = sum(temp.rewNMinusOneStay) / sum(temp.nTrial);
                        mat.prewNMinusOneSwitch.(cRegion)(1,i) = sum(temp.rewNMinusOneSwitch) / sum(temp.nTrial);
                        mat.pnoRewNMinusOneStay.(cRegion)(1,i) = sum(temp.noRewNMinusOneStay) / sum(temp.nTrial);
                        mat.pnoRewNMinusOneSwitch.(cRegion)(1,i) = sum(temp.noRewNMinusOneSwitch) / sum(temp.nTrial);

                        mat.pLWinSwitch.(cRegion)(1,i) = sum(temp.LWinSwitch) / (sum(temp.nLWin));
                        mat.pLWinStay.(cRegion)(1,i) = sum(temp.LWinStay) / (sum(temp.nLWin));
                        mat.pLWinMiss.(cRegion)(1,i) = sum(temp.LWinMiss) / (sum(temp.nLWin));
                        mat.pRWinSwitch.(cRegion)(1,i) = sum(temp.RWinSwitch) / (sum(temp.nRWin));
                        mat.pRWinStay.(cRegion)(1,i) = sum(temp.RWinStay) / (sum(temp.nRWin));
                        mat.pRWinMiss.(cRegion)(1,i) = sum(temp.RWinMiss) / (sum(temp.nRWin));
                        
                        
                        mat.pLLoseSwitch.(cRegion)(1,i) = sum(temp.LLoseSwitch) / sum(temp.nLLose);
                        mat.pLLoseStay.(cRegion)(1,i) = sum(temp.LLoseStay) / sum(temp.nLLose);
                        mat.pLLoseMiss.(cRegion)(1,i) = sum(temp.LLoseMiss) / sum(temp.nLLose);
                        mat.pRLoseSwitch.(cRegion)(1,i) = sum(temp.RLoseSwitch) / sum(temp.nRLose);
                        mat.pRLoseStay.(cRegion)(1,i) = sum(temp.RLoseStay) / sum(temp.nRLose);
                        mat.pRLoseMiss.(cRegion)(1,i) = sum(temp.RLoseMiss) / sum(temp.nRLose);
                        
                        
                        mat.pWinStay.(cRegion)(1,i) = (sum(temp.LWinStay) + sum(temp.RWinStay)) / sum(temp.nWin);
                        mat.pLoseStay.(cRegion)(1,i) = (sum(temp.LLoseStay) + sum(temp.RLoseStay)) / sum(temp.nLose);
                        mat.pWinSwitch.(cRegion)(1,i) = (sum(temp.LWinSwitch) + sum(temp.RWinSwitch)) / sum(temp.nWin);
                        mat.pLoseSwitch.(cRegion)(1,i) = (sum(temp.LLoseSwitch) + sum(temp.RLoseSwitch)) / sum(temp.nLose);
                        mat.pWinMiss.(cRegion)(1,i) = (sum(temp.LWinMiss) + sum(temp.RWinMiss)) / sum(temp.nWin);
                        mat.pLoseMiss.(cRegion)(1,i) = (sum(temp.LLoseMiss) + sum(temp.RLoseMiss)) / sum(temp.nLose);
                        
                        mat.pMissNextTrial.(cRegion)(1,i) = sum(temp.missNextTrial) / sum(temp.nTrial); 
                        mat.pRewardNextTrial.(cRegion)(1,i) = sum(temp.rewardNextTrial) / sum(temp.nTrial); 
                        
                        mat.pStayCurrentTrial.(cRegion)(1,i) = sum(temp.stayCurrentTrial) / sum(temp.nTrial); 
                        mat.pSwitchCurrentTrial.(cRegion)(1,i) = sum(temp.switchCurrentTrial) / sum(temp.nTrial); 
                        
                        mat.pWSLS.(cRegion)(1,i) = sum(temp.WSLS) / sum(temp.nTrial); 
                end 
            end 
        if j == 1
            WSLS_chat = mat;
        else
            WSLS_dbh = mat;
        end
    end 
    
    cd 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\1Region\'
    filename = append('1Region_', currReg, '.mat');
    save(filename,'WSLS_chat','WSLS_dbh',"lowEntropyTable","dataIndex")
end

%%

%basic trial stats
% MP_STIM_plot_all(WSLS_chat.nTrial, WSLS_dbh.nTrial, 'Number of Trials', save_path)
% MP_STIM_plot_all(WSLS_chat.miss, WSLS_dbh.miss, 'Number of Misses', save_path)
% MP_STIM_plot_all(WSLS_chat.left, WSLS_dbh.left, 'Number of Left Choices', save_path)
% MP_STIM_plot_all(WSLS_chat.right, WSLS_dbh.right, 'Number of Right Choices', save_path)

%Win v Lose 
MP_STIM_plot_all(WSLS_chat.pStay, WSLS_dbh.pStay, 'Probability of Stay', save_path)
MP_STIM_plot_all(WSLS_chat.pSwitch, WSLS_dbh.pSwitch, 'Probability of Switch', save_path)
MP_STIM_plot_all(WSLS_chat.pWinStay, WSLS_dbh.pWinStay, 'Probability of Win-Stay', save_path)
MP_STIM_plot_all(WSLS_chat.pLoseStay, WSLS_dbh.pLoseStay, 'Probability of Lose-Stay', save_path)

MP_STIM_plot_all(WSLS_chat.pLWinStay, WSLS_dbh.pLWinStay, 'Probability of Left then Win-Stay', save_path)
MP_STIM_plot_all(WSLS_chat.pRWinStay, WSLS_dbh.pRWinStay, 'Probability of Right then Win-Stay', save_path)
MP_STIM_plot_all(WSLS_chat.pLLoseStay, WSLS_dbh.pLLoseStay, 'Probability of Left then Lose-Stay', save_path)
MP_STIM_plot_all(WSLS_chat.pRLoseStay, WSLS_dbh.pRLoseStay, 'Probaility of Right then Lose-Stay', save_path)


MP_STIM_plot_all(WSLS_chat.pStayCurrentTrial, WSLS_dbh.pStayCurrentTrial, 'Probability to Stay on Current Trial', save_path)
%MP_STIM_plot_all(WSLS_chat.pSwitchCurentTrial, WSLS_dbh.pSwitchCurentTrial, 'Probability to Switch on Current Trial, save_path')


MP_STIM_plot_all(WSLS_chat.pRewardNextTrial, WSLS_dbh.pRewardNextTrial, 'Probability of Reward Next Trial', save_path)
MP_STIM_plot_all(WSLS_chat.pMissNextTrial, WSLS_dbh.pMissNextTrial, 'Probability to Miss Next Trial', save_path)
MP_STIM_plot_all(WSLS_chat.pWSLS, WSLS_dbh.pWSLS, 'Probability to WSLS', save_path)

%%
clear all
cd 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\1Region\'

load('1Region_LM2.mat')
dat.chat.LM2 = WSLS_chat;
dat.dbh.LM2 = WSLS_dbh;
lowEntropyTable.Strain = dataIndex.Strain;
stats.LM2 = lowEntropyTable;

load('1Region_RM2.mat')
dat.chat.RM2 = WSLS_chat;
dat.dbh.RM2 = WSLS_dbh;
lowEntropyTable.Strain = dataIndex.Strain;
stats.RM2 = lowEntropyTable;

load('1Region_LV1.mat')
dat.chat.LV1 = WSLS_chat;
dat.dbh.LV1 = WSLS_dbh;
lowEntropyTable.Strain = dataIndex.Strain;
stats.LV1 = lowEntropyTable;

load('1Region_RV1.mat')
dat.chat.RV1 = WSLS_chat;
dat.dbh.RV1 = WSLS_dbh;
lowEntropyTable.Strain = dataIndex.Strain;
stats.RV1 = lowEntropyTable;

%% figure out each mouse's preferred side

strains = fields(dat);
field = fields(dat.chat.LM2);
regions = fields(dat.chat);
stims = fields(dat.chat.LM2.left);
sides = [{'left'};{'right'}];

%loop produces k x p (regions x # of mouse) structure for each combination
%of genotype and side
for j = 1:size(strains,1)
    strain = strains{j};

    for k = 1:size(regions,1)
        region = regions{k};

        for m = 1:size(sides,1)
            side = sides{m};
                
            stim = dat.(strain).(region).(side).STIM;
            ctrl = dat.(strain).(region).(side).CTRL;

            for p = 1:size(stim,2)
               currStim = stim{p};
               currCtrl = ctrl{p};

               pref.(strain).(side){k,p} = currStim + currCtrl; 
            end 
        end
    end 
end 

%% plot with distinct color for each region and shape for each animal
sidePreference = [];
for h = 1:size(strains,1) 
    figure %new figure for each strain 
    sz = 20;
    strain = strains{h};
    for j = 1:size(pref.(strain).left,2) %iterate through animals
                %colors = {'b';'g';'r';'m'};
                temp = [];
        for i = 1:size(pref.(strain).left,1) %iterate through regions
                %symbols = {'o'; 'x'; 'd'; 'v'; '*';'h';'s'};
                rand_offset = 0.3*rand(1, size(pref.(strain).left{i,j},2));
                colors = {'b';'g';'r';'m';'k';'c';'y'};
                hold on 
                subplot(3,3,j)
                %plot right/total choices by region and animal
                scatter((pref.(strain).right{i,j}) ./ (pref.(strain).right{i,j} + pref.(strain).left{i,j}),...
                        i*ones(1, size(pref.(strain).right{i,j},2))+rand_offset,...
                        sz,colors{i},'o') %symbols{j}

                hold on %plot mean
                scatter(mean(pref.(strain).right{i,j} ./ (pref.(strain).right{i,j} + pref.(strain).left{i,j})),...
                         i,...
                         sz*2, colors{i},'+')

                hold on %plot error bars
                errorbar(mean(pref.(strain).right{i,j} ./ (pref.(strain).right{i,j} + pref.(strain).left{i,j})),...
                            i,...
                         std((pref.(strain).right{i,j}) ./ (pref.(strain).right{i,j} + pref.(strain).left{i,j})) / sqrt(size(pref.(strain).right{i,j},2))*2,...
                         colors{i}, 'horizontal');

                hold on %plot right/total by animal (all regions combined)
                
                scatter((pref.(strain).right{i,j}) ./ (pref.(strain).right{i,j} + pref.(strain).left{i,j}),...
                        5*ones(1, size(pref.(strain).right{i,j},2))+rand_offset,...
                        sz,'k','o') %symbols{j}
                temp = [temp,(pref.(strain).right{i,j}) ./ (pref.(strain).right{i,j} + pref.(strain).left{i,j})];

                xticks([0 .5 1])
                yticks([1 2 3 4 5])
                yticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1' ;'ALL'])
                xlim([0 1])
                ylim([0 6])
                hold on
                xline(0.5, '-k')

        end 
        hold on %plot mean of all regions
        scatter(mean(temp), 5, sz*2, 'k', '+')

        hold on %plot SE of all regions
        errorbar(mean(temp), 5,...
                         std(temp) / sqrt(size(temp,2)),...
                         'k', 'horizontal');
        %legend

        if mean(temp) > .5
            sidePreference.(strain){j} = 'right';
        elseif mean(temp) < .5 
            sidePreference.(strain){j} = 'left';
        else 
            sidePreference.(strain){j} = 'none';
        end
    end 
    sgtitle(strain)

end 

%% if mouse has left preference, switch left and right sided variables so all mice have same prefered side (right)

strains = fields(dat);
field = fields(dat.chat.LM2);
regions = fields(dat.chat);

for t = 1:size(strains,1) %for each strain
    strain = strains{t};
    for r = 1:size(sidePreference.(strain),2) %for each animal
        if strcmp(sidePreference.(strain){r},'left') == 1 %if animal has left preference
            for s = 1:size(regions,1) %for each region
                region = regions{s};
                for u = 1:size(field,1) %for each variable
                    var = field{u};
                    if strcmp(var(2), 'R') == 1
                        if strcmp(var, 'pRewardNextTrial')
                        else
                            tempR = var; %save R sided variable name
                            var(2) = 'L'; %create L sided variable name
                            tempVarR = dat.(strain).(region).(tempR); %save R sided data (non-prefered side)
                            dat.(strain).(region).(tempR) = dat.(strain).(region).(var); %overwrite prefered side onto Right
                            dat.(strain).(region).(var) = tempVarR; %overwrite non-prefered side date onto left
                        end
                    end
                end
            end 
        end 
    end 
end 


%% plot WSLS with stim and no stim on the same plot
% including stats for p<.01 and p<.001 in two way anova (region and
    %stim/no stim as factors)

%choose plotting format 
    % 1 = raw data
    % 2 = percent change from no stim trials
    % 3 = percent change from average of control regions (LV1 and RV1)
   setup_figprop;     
        d = 1; 

        %set p value threshold
        threshold = 0.05;

        %test for significance or not?
        testForSig = 0;



strains = fields(dat);
field = fields(dat.chat.LM2);
regions = fields(dat.chat);

for i = 1:size(field,1)
    name = field{i};
    if strcmp(name, 'stats') == 1 || strcmp(name, 'nTrial') == 1 || strcmp(name, 'left') == 1 ...
            || strcmp(name, 'right') == 1 || strcmp(name, 'miss') == 1 || strcmp(name, 'nwin') == 1
    
    else
        fig = figure('position',[0 0 1500 800])
        stats_dbh = []; %empty array for data used in anova
        stats_chat = [];
        stim_dbh = []; %empty array for STIM or CTRL category
        stim_chat = [];
        region_dbh = []; %empty array for region categories
        region_chat = [];

        for r = 1:size(regions,1)
            region = regions{r};
            chaty = ones(1, size(dat.chat.(region).(name).CTRL,2));
            rand_chat = 0.3*rand(1, size(dat.chat.(region).(name).CTRL,2));
            dbhy = ones(1, size(dat.dbh.(region).(name).CTRL,2));
            rand_dbh = 0.3*rand(1, size(dat.dbh.(region).(name).CTRL,2));
            sz = 100;

            c = [0.9290 0.6940 0.1250];
            if d == 1
                fileAppend = 'pt2_pt6Yaxis_stats_';
                offset = .5;
                grey = [.7 .7 .7]; 
                subplot(1,2,2)
                hold on
                c = [0.9290 0.6940 0.1250];
                scatter((2*r*chaty)+rand_chat, dat.chat.(region).(name).STIM, sz, c)

                hold on
                scatter(2*r, mean(dat.chat.(region).(name).STIM, 'omitnan'), (sz*2),c, 'filled')
                errorbar(2*r, mean(dat.chat.(region).(name).STIM, 'omitnan'), std(dat.chat.(region).(name).STIM, 'omitnan') / sqrt(size(dat.chat.(region).(name).STIM,2)),...
                         'Color',c)

                hold on
                scatter((2*r*chaty)+rand_chat+offset, dat.chat.(region).(name).CTRL, sz, grey)

                hold on
                scatter((2*r)+offset, mean(dat.chat.(region).(name).CTRL, 'omitnan'), (sz*2),grey, 'filled')
                errorbar((2*r)+offset, mean(dat.chat.(region).(name).CTRL, 'omitnan'), std(dat.chat.(region).(name).CTRL, 'omitnan') / sqrt(size(dat.chat.(region).(name).CTRL,2)),...
                         'color', grey)

                ax=gca;               
                ax.TickLength = [0 0];
                ax.XLim = [0 10];
                ax.YLim = [.17 .6];
                    ax.YTick = [.2 .4 .6];
                    ax.YTickLabel = ["20", "40", "60"];
                ax.XTick = [2 4 6 8];
                %ylabel('pTrials')
                ax.XTickLabel = ['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'];
                %title('ChAt')
                ax.XTickLabelRotation = 90;
                ax.FontSize  = 25;
                ax.TitleFontSizeMultiplier = 2;

                subplot(1,2,1)
                c = [.1 0.6 0.6];
                hold on
                scatter((2*r*dbhy)+rand_dbh, dat.dbh.(region).(name).STIM, sz,c)

                hold on
                scatter((2*r), mean(dat.dbh.(region).(name).STIM, 'omitnan'), (sz*2),c, 'filled')
                errorbar(2*r, mean(dat.dbh.(region).(name).STIM, 'omitnan'), std(dat.dbh.(region).(name).STIM, 'omitnan') / sqrt(size(dat.dbh.(region).(name).STIM,2)),...
                         'color', c)
            
                hold on
                scatter((2*r*dbhy)+rand_dbh+offset, dat.dbh.(region).(name).CTRL, sz,grey)

                hold on
                scatter((2*r)+offset, mean(dat.dbh.(region).(name).CTRL, 'omitnan'), (sz*2),grey, 'filled')
                errorbar((2*r)+offset, mean(dat.dbh.(region).(name).CTRL, 'omitnan'), std(dat.dbh.(region).(name).CTRL, 'omitnan') / sqrt(size(dat.dbh.(region).(name).CTRL,2)),...
                         'color', grey)

                ax =gca;
                ax.TickLength = [0 0];
                ax.XLim = [0 10];
                ax.YLim = [.17 .6];
                    ax.YTick = [.2 .4 .6];
                    ax.YTickLabel = ["20", "40", "60"];
                ax.XTick = [2 4 6 8];
                %ylabel('pTrials')
                ax.XTickLabel = ['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'];
                %title('DBH')
                ax.XTickLabelRotation = 90;
                ax.FontSize  = 25;
                ax.TitleFontSizeMultiplier = 2;

                %add data to arrays for anova

                stats_dbh = [stats_dbh, dat.dbh.(region).(name).STIM, dat.dbh.(region).(name).CTRL];
                stim_dbh = [stim_dbh; repmat("STIM", size(dat.dbh.(region).(name).STIM,2),1); repmat("CTRL", size(dat.dbh.(region).(name).CTRL,2) ,1)];
                region_dbh = [region_dbh; repmat(string(region), size(dat.dbh.(region).(name).STIM,2),1); repmat(string(region), size(dat.dbh.(region).(name).CTRL,2),1)];

                stats_chat = [stats_chat, dat.chat.(region).(name).STIM, dat.chat.(region).(name).CTRL];
                stim_chat = [stim_chat; repmat("STIM", size(dat.chat.(region).(name).STIM,2),1); repmat("CTRL", size(dat.chat.(region).(name).CTRL,2) ,1)];
                region_chat = [region_chat; repmat(string(region), size(dat.chat.(region).(name).STIM,2),1); repmat(string(region), size(dat.chat.(region).(name).CTRL,2),1)];


    
            elseif d == 2
                low = -100;
                high = 100;
                fileAppend = 'pChangeNoStim_';
                subplot(2,2,1)
                hold on
                c = [0.9290 0.6940 0.1250];
                scatter(r*chaty+rand_chat, (dat.chat.(region).(name).STIM - dat.chat.(region).(name).CTRL) ./ dat.chat.(region).(name).CTRL * 100, sz, c)
                scatter(r, mean((dat.chat.(region).(name).STIM - dat.chat.(region).(name).CTRL) ./ dat.chat.(region).(name).CTRL) * 100, (sz*2),c, 'filled')

                ylim([low high])
                ylabel('%ChangeNoStim')
                xlim ([0 5])
                xticks([1 2 3 4])
                xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
                title('ChAt w/ stim')

                subplot(2,2,2)
                hold on
                scatter(r*chaty+rand_chat, (dat.chat.(region).(name).CTRL - dat.chat.(region).(name).CTRL) ./ dat.chat.(region).(name).CTRL * 100, sz, c)
                scatter(r, mean((dat.chat.(region).(name).CTRL - dat.chat.(region).(name).CTRL) ./ dat.chat.(region).(name).CTRL) * 100, (sz*2),c, 'filled')

                xlim ([0 5])
                ylim([low high])
                xticks([1 2 3 4])
                ylabel('%ChangeNoStim')
                xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
                title('ChAt no stim')

                subplot(2,2,3)
                c = [.1 0.6 0.6];
                hold on
                scatter(r*dbhy+rand_dbh, (dat.dbh.(region).(name).STIM - dat.dbh.(region).(name).CTRL) ./ dat.dbh.(region).(name).CTRL * 100, sz,c)
                scatter(r, mean((dat.dbh.(region).(name).STIM - dat.dbh.(region).(name).CTRL) ./ dat.dbh.(region).(name).CTRL) * 100, (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([low high])
                xticks([1 2 3 4])
                ylabel('%ChangeNoStim')
                xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
                title('DBH w/ stim')
            
                subplot(2,2,4)
                hold on
                scatter(r*dbhy+rand_dbh, (dat.dbh.(region).(name).CTRL - dat.dbh.(region).(name).CTRL) ./ dat.dbh.(region).(name).CTRL * 100, sz,c)
                scatter(r, mean((dat.dbh.(region).(name).CTRL - dat.dbh.(region).(name).CTRL) ./ dat.dbh.(region).(name).CTRL) * 100, (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([low high])
                xticks([1 2 3 4])
                ylabel('%ChangeNoStim')
                xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
                title('DBH no stim')
            elseif d == 3
                fileAppend = 'pChangeV1_';
                bottom = -100;
                top = 100;

                chat_v1 = dat.chat.LV1.(name).STIM + dat.chat.RV1.(name).STIM ./ 2;
                dbh_v1 = dat.dbh.LV1.(name).STIM + dat.dbh.RV1.(name).STIM ./ 2;

                chat_v1_ctrl = dat.chat.LV1.(name).CTRL + dat.chat.RV1.(name).CTRL ./ 2;
                dbh_v1_ctrl = dat.dbh.LV1.(name).CTRL + dat.dbh.RV1.(name).CTRL ./ 2;


                subplot(2,2,1)
                hold on
                c = [0.9290 0.6940 0.1250];
                scatter(r*chaty+rand_chat, (dat.chat.(region).(name).STIM - chat_v1) ./ chat_v1 * 100, sz, c)
                scatter(r, mean((dat.chat.(region).(name).STIM - chat_v1) ./ chat_v1) * 100, (sz*2),c, 'filled')

                ylim([bottom top])
                ylabel('%ChangeV1')
                xlim ([0 5])
                xticks([1 2 3 4])
                xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
                title('ChAt w/ stim')

                subplot(2,2,2)
                hold on
                scatter(r*chaty+rand_chat, (dat.chat.(region).(name).CTRL - chat_v1_ctrl) ./ chat_v1_ctrl * 100, sz, c)
                scatter(r, mean((dat.chat.(region).(name).CTRL - chat_v1_ctrl) ./ chat_v1_ctrl) * 100, (sz*2),c, 'filled')

                xlim ([0 5])
                ylim([bottom top])
                xticks([1 2 3 4])
                ylabel('%ChangeV1')
                xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
                title('ChAt no stim')

                subplot(2,2,3)
                c = [.1 0.6 0.6];
                hold on
                scatter(r*dbhy+rand_dbh, (dat.dbh.(region).(name).STIM - dbh_v1) ./ dbh_v1 * 100, sz,c)
                scatter(r, mean((dat.dbh.(region).(name).STIM - dbh_v1) ./ dbh_v1) * 100, (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([bottom top])
                xticks([1 2 3 4])
                ylabel('%ChangeV1')
                xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
                title('DBH w/ stim')
            
                subplot(2,2,4)
                hold on
                scatter(r*dbhy+rand_dbh, (dat.dbh.(region).(name).CTRL - dbh_v1_ctrl) ./ dbh_v1_ctrl * 100, sz,c)
                scatter(r, mean((dat.dbh.(region).(name).CTRL - dbh_v1_ctrl) ./ dbh_v1_ctrl) * 100, (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([bottom top])
                xticks([1 2 3 4])
                ylabel('%ChangeV1')
                xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
                title('DBH no stim')
            end 
            
        end
    if testForSig == 1
               %test chat stats with a 2-way anova
            areg=[];%empty arrays to keep track of significance
            breg=[];
            aint=[];
            bint=[];
            aa=[];
            factors_chat = {stim_chat' region_chat'};
            [P1,T,STATS,TERMS] = anovan(stats_chat, factors_chat,'varnames',{'stim','region'}, 'display', 'off', 'model','interaction');

                if P1(3) < threshold %if interaction is significant 
                    [results,~,~,gnames] = multcompare(STATS, "display", "off", "dimension", [1 2]);
                    chatPostHocInteraction = array2table(results,"VariableNames", ...
                        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                    chatPostHocInteraction.("Group A")=gnames(chatPostHocInteraction.("Group A"));
                    chatPostHocInteraction.("Group B")=gnames(chatPostHocInteraction.("Group B"))
                    name

                     n = 0; % off set for p<.05 y value 
            
                       %add stats to graph for interaction
                        for m = 1:size(chatPostHocInteraction,1)
                            p = chatPostHocInteraction.("P-value")(m);
                             if strcmp(chatPostHocInteraction.("Group A"){m}(6:9), 'STIM') == 1
                                offsetA = 0; %no offset if stim category
                            else
                                offsetA = .5; %offset if CTRL
                            end 
                            
                            if strcmp(chatPostHocInteraction.("Group B"){m}(6:9), 'STIM') == 1
                                offsetB = 0; %no offset if stim category
                            else
                                offsetB = .5; %offset if CTRL
                            end
        
                            if strcmp(chatPostHocInteraction.("Group A"){m}(18:20), 'LM2') == 1
                                xcoordA = 2; %Y coordinate for region
                            elseif strcmp(chatPostHocInteraction.("Group A"){m}(18:20), 'RM2') == 1
                                xcoordA = 4; %Y coordinate for region
                            elseif strcmp(chatPostHocInteraction.("Group A"){m}(18:20), 'LV1') == 1
                                xcoordA = 6; %Y coordinate for region
                            elseif strcmp(chatPostHocInteraction.("Group A"){m}(18:20), 'RV1') == 1
                                xcoordA = 8; %Y coordinate for region
                            end 
                            
                            if strcmp(chatPostHocInteraction.("Group B"){m}(18:20), 'LM2') == 1
                                xcoordB = 2; %Y coordinate for region
                            elseif strcmp(chatPostHocInteraction.("Group B"){m}(18:20), 'RM2') == 1
                                xcoordB = 4; %Y coordinate for region
                            elseif strcmp(chatPostHocInteraction.("Group B"){m}(18:20), 'LV1') == 1
                                xcoordB = 6; %Y coordinate for region
                            elseif strcmp(chatPostHocInteraction.("Group B"){m}(18:20), 'RV1') == 1
                                xcoordB = 8; %Y coordinate for region
                            end 
        
        
%                             if p <= .001 
        
%                                     subplot(1,2,2)
%                                     hold on
%                                     y1 = .75 + n; 
%                                     aint = plot([xcoordA+offsetA, xcoordB+offsetB], [y1, y1], '-k')
%                                     plot(xcoordA+offsetA, y1, '+k')
%                                     plot(xcoordB+offsetB, y1, '+k')
%                                     
%                                     n = n + .01; %y offset for multiple p value plots
        
                            if  p <= threshold
                                    subplot(1,2,2)
                                    hold on
                                    y2 = .75 + n; 
                                   bint = plot([xcoordA+offsetA, xcoordB+offsetB],[y2, y2], '--k')
                                   hold on
                                    plot(xcoordA+offsetA, y2, '+k')
                                    hold on
                                    plot(xcoordB+offsetB, y2, '+k')
                                    n = n + .01; %y offset for multiple p value plots;
                            end
                            clear xcoordA xcoordB offsetA offsetB
                        end 
                end

                chatPostHocRegion = table();
                if P1(2) < threshold %if main effect of region is significant 
                    [results,~,~,gnames] = multcompare(STATS, "display", "off", "dimension", 2);
                    chatPostHocRegion = array2table(results,"VariableNames", ...
                        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                    chatPostHocRegion.("Group A")=gnames(chatPostHocRegion.("Group A"));
                    chatPostHocRegion.("Group B")=gnames(chatPostHocRegion.("Group B"))
                    name

                     n = 0; % off set for p<.05 y value 
            
                       %add stats to graph for interaction
                        for m = 1:size(chatPostHocRegion,1)
                            p = chatPostHocRegion.("P-value")(m);
        
                            if strcmp(chatPostHocRegion.("Group A"){m}(end-2:end), 'LM2') == 1
                                xcoordA = 2.25; %Y coordinate for region
                            elseif strcmp(chatPostHocRegion.("Group A"){m}(end-2:end), 'RM2') == 1
                                xcoordA = 4.25; %Y coordinate for region
                            elseif strcmp(chatPostHocRegion.("Group A"){m}(end-2:end), 'LV1') == 1
                                xcoordA = 6.25; %Y coordinate for region
                            elseif strcmp(chatPostHocRegion.("Group A"){m}(end-2:end), 'RV1') == 1
                                xcoordA = 8.25; %Y coordinate for region
                            end 
                            
                            if strcmp(chatPostHocRegion.("Group B"){m}(end-2:end), 'LM2') == 1
                                xcoordB = 2.25; %Y coordinate for region
                            elseif strcmp(chatPostHocRegion.("Group B"){m}(end-2:end), 'RM2') == 1
                                xcoordB = 4.25; %Y coordinate for region
                            elseif strcmp(chatPostHocRegion.("Group B"){m}(end-2:end), 'LV1') == 1
                                xcoordB = 6.25; %Y coordinate for region
                            elseif strcmp(chatPostHocRegion.("Group B"){m}(end-2:end), 'RV1') == 1
                                xcoordB = 8.25; %Y coordinate for region
                            end 

%                              if p <= .001 
        
%                                     subplot(1,2,1)
%                                     hold on
%                                     y1 = .8 + n; 
%                                     areg = plot([xcoordA, xcoordB], [y1, y1], '-b')
%                                     plot(xcoordA, y1, '+b')
%                                     plot(xcoordB, y1, '+b')
%                                     
%                                     n = n + .01; %y offset for multiple p value plots
        
                            if  p <= threshold
                                    subplot(1,2,1)
                                    hold on
                                    y2 = .8 + n; 
                                    breg = plot([xcoordA, xcoordB],[y2, y2], '--b')
                                    hold on
                                    plot(xcoordA, y2, '+b')
                                    hold on
                                    plot(xcoordB, y2, '+b')
                                    n = n + .01; %y offset for multiple p value plots;
                            end
                            clear xcoordA xcoordB
                         end
                end

                if P1(1) < threshold %if main effect of stimulation is significant
                        aa = plot(10,0, '--r'); 
                end 
  

        
           
                subplot(1,2,1)
                thresholdName = strcat("p<",num2str(threshold));
                legendMarkers = [];
                legendLabels = [];
                if ~isempty(aa)
                    legendMarkers = [legendMarkers,aa];
                    legendLabels = [legendLabels, "main effect of stimulation"];
                end
                if ~isempty(breg)
                    legendMarkers = [legendMarkers; breg];
                    tempMark = strcat("main effect region ", thresholdName);
                    legendLabels = [legendLabels; tempMark];
                end
%                 if ~isempty(areg)
%                     legendMarkers = [legendMarkers;areg];
%                     tempMark = strcat("main effect region ", 'p<.001');
%                     legendLabels = [legendLabels; tempMark];
%                 end
%                 if ~isempty(aint)
%                     legendMarkers = [legendMarkers; aint];
%                     tempMark = strcat("interaction ", "p<.001");
%                     legendLabels = [legendLabels; tempMark];
%                 end
                if ~isempty(bint)
                    legendMarkers = [legendMarkers; bint];
                    tempMark = strcat("interaction ", thresholdName);
                    legendLabels = [legendLabels; tempMark];
                end
         
                if ~isempty(aa) || ~isempty(areg) || ~isempty(breg) || ~isempty(aint) || ~isempty(bint)
                        legend(legendMarkers, legendLabels, 'Location','best')
                end
% 
%                 if ~isempty(areg) && ~isempty(breg) && ~isempty(aa) && ~isempty
%                     legend([areg,breg,aa],'p<.001', thresholdName,'main effect STIM', 'Location','best')
%                     ~isempty(areg) && ~isempty(breg) && ~isempty(aa)
%                     legend([areg,breg,aa],'p<.001', thresholdName,'main effect STIM', 'Location','best')
%                 elseif ~isempty(a) && ~isempty(b) 
%                     legend([a,b],'p<.001', thresholdName, 'Location','best')
%                 elseif ~isempty(a)
%                     legend(a,'p<.001','Location','best')
%                 elseif ~isempty(b)
%                     legend(b, thresholdName,'Location','best')
%                 end


        %test dbh stats with a 2-way anova

            areg=[];%empty arrays to keep track of significance
            breg=[];
            aint=[];
            bint=[];
            aa=[];
            factors_dbh = {stim_dbh' region_dbh'};
            [P2,T2,STATS2,TERMS2] = anovan(stats_dbh, factors_dbh,'varnames',{'stim','region'}, 'display', 'off', 'model','interaction');
                if P2(3) < threshold %if interaction is significant 
                    [results2,~,~,gnames2] = multcompare(STATS2, "display", "off", "dimension", [1 2]);
                    dbhPostHocInteraction = array2table(results2,"VariableNames", ...
                        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                    dbhPostHocInteraction.("Group A")=gnames2(dbhPostHocInteraction.("Group A"));
                    dbhPostHocInteraction.("Group B")=gnames2(dbhPostHocInteraction.("Group B"))
                    name

                     n = 0; % off set for p<.05 y value 
            
                       %add stats to graph for interaction
                        for m = 1:size(dbhPostHocInteraction,1)
                            p = dbhPostHocInteraction.("P-value")(m);
                             if strcmp(dbhPostHocInteraction.("Group A"){m}(6:9), 'STIM') == 1
                                offsetA = 0; %no offset if stim category
                            else
                                offsetA = .5; %offset if CTRL
                            end 
                            
                            if strcmp(dbhPostHocInteraction.("Group B"){m}(6:9), 'STIM') == 1
                                offsetB = 0; %no offset if stim category
                            else
                                offsetB = .5; %offset if CTRL
                            end
        
                            if strcmp(dbhPostHocInteraction.("Group A"){m}(18:20), 'LM2') == 1
                                xcoordA = 2; %Y coordinate for region
                            elseif strcmp(dbhPostHocInteraction.("Group A"){m}(18:20), 'RM2') == 1
                                xcoordA = 4; %Y coordinate for region
                            elseif strcmp(dbhPostHocInteraction.("Group A"){m}(18:20), 'LV1') == 1
                                xcoordA = 6; %Y coordinate for region
                            elseif strcmp(dbhPostHocInteraction.("Group A"){m}(18:20), 'RV1') == 1
                                xcoordA = 8; %Y coordinate for region
                            end 
                            
                            if strcmp(dbhPostHocInteraction.("Group B"){m}(18:20), 'LM2') == 1
                                xcoordB = 2; %Y coordinate for region
                            elseif strcmp(dbhPostHocInteraction.("Group B"){m}(18:20), 'RM2') == 1
                                xcoordB = 4; %Y coordinate for region
                            elseif strcmp(dbhPostHocInteraction.("Group B"){m}(18:20), 'LV1') == 1
                                xcoordB = 6; %Y coordinate for region
                            elseif strcmp(dbhPostHocInteraction.("Group B"){m}(18:20), 'RV1') == 1
                                xcoordB = 8; %Y coordinate for region
                            end 
        
        
%                             if p <= .001 
%         
%                                     subplot(1,2,2)
%                                     hold on
%                                     y1 = .75 + n; 
%                                     areg = plot([xcoordA+offsetA, xcoordB+offsetB], [y1, y1], '-k')
%                                     plot(xcoordA+offsetA, y1, '+k')
%                                     plot(xcoordB+offsetB, y1, '+k')
                                    
%                                     n = n + .01; %y offset for multiple p value plots
        
                            if  p <= threshold
                                    subplot(1,2,2)
                                    hold on
                                    y2 = .75 + n; 
                                    breg = plot([xcoordA+offsetA, xcoordB+offsetB],[y2, y2], '--b')
                                    hold on
                                    plot(xcoordA+offsetA, y2, '+b')
                                    hold on
                                    plot(xcoordB+offsetB, y2, '+b')
                   
                                    n = n + .01; %y offset for multiple p value plots;
                            end
                            clear xcoordA xcoordB offsetA offsetB
                        end 
                end
     
                if P2(2) < threshold %if main effect of region is significant 
                    [results2,~,~,gnames2] = multcompare(STATS2, "display", "off", "dimension", 2);
                    dbhPostHocRegion = array2table(results2,"VariableNames", ...
                        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                    dbhPostHocRegion.("Group A")=gnames2(dbhPostHocRegion.("Group A"));
                    dbhPostHocRegion.("Group B")=gnames2(dbhPostHocRegion.("Group B"))
                    name

                     n = 0; % off set for p<.05 y value 
            
                       %add stats to graph for interaction
                        for m = 1:size(dbhPostHocRegion,1)
                            p = dbhPostHocRegion.("P-value")(m);
        
                            if strcmp(dbhPostHocRegion.("Group A"){m}(end-2:end), 'LM2') == 1
                                xcoordA = 2.25; %Y coordinate for region
                            elseif strcmp(dbhPostHocRegion.("Group A"){m}(end-2:end), 'RM2') == 1
                                xcoordA = 4.25; %Y coordinate for region
                            elseif strcmp(dbhPostHocRegion.("Group A"){m}(end-2:end), 'LV1') == 1
                                xcoordA = 6.25; %Y coordinate for region
                            elseif strcmp(dbhPostHocRegion.("Group A"){m}(end-2:end), 'RV1') == 1
                                xcoordA = 8.25; %Y coordinate for region
                            end 
                            
                            if strcmp(dbhPostHocRegion.("Group B"){m}(end-2:end), 'LM2') == 1
                                xcoordB = 2.25; %Y coordinate for region
                            elseif strcmp(dbhPostHocRegion.("Group B"){m}(end-2:end), 'RM2') == 1
                                xcoordB = 4.25; %Y coordinate for region
                            elseif strcmp(dbhPostHocRegion.("Group B"){m}(end-2:end), 'LV1') == 1
                                xcoordB = 6.25; %Y coordinate for region
                            elseif strcmp(dbhPostHocRegion.("Group B"){m}(end-2:end), 'RV1') == 1
                                xcoordB = 8.25; %Y coordinate for region
                            end 
% 
%                              if p <= .001 
        
%                                     subplot(1,2,2)
%                                     hold on
%                                     y1 = .8 + n; 
%                                     aint = plot([xcoordA, xcoordB], [y1, y1], '-b')
%                                     plot(xcoordA, y1, '+b')
%                                     plot(xcoordB, y1, '+b')
%                                     
%                                     n = n + .01; %y offset for multiple p value plots
        
                           if  p <= threshold
                                    subplot(1,2,2)
                                    hold on
                                    y2 = .8 + n; 
                                   bint = plot([xcoordA, xcoordB],[y2, y2], '--k')
                                   hold on
                                    plot(xcoordA, y2, '+k')
                                    hold on
                                    plot(xcoordB, y2, '+k')
                                    n = n + .01; %y offset for multiple p value plots;
                            end
                            clear xcoordA xcoordB
                         end
                end

                if P2(1) < threshold %if main effect of stimulation is significant
                        aa = plot(10,0, '--r'); 
                end 
  

                subplot(1,2,2)
                thresholdName = strcat("p<",num2str(threshold));
                legendMarkers = [];
                legendLabels = [];
                if ~isempty(aa)
                    legendMarkers = [legendMarkers,aa];
                    legendLabels = [legendLabels, "main effect of stimulation"];
                end
                if ~isempty(breg)
                    legendMarkers = [legendMarkers; breg];
                    tempMark = strcat("main effect region ", thresholdName);
                    legendLabels = [legendLabels; tempMark];
                end
%                 if ~isempty(areg)
%                     legendMarkers = [legendMarkers;areg];
%                     tempMark = strcat("main effect region ", 'p<.001');
%                     legendLabels = [legendLabels; tempMark];
%                 end
%                 if ~isempty(aint)
%                     legendMarkers = [legendMarkers; aint];
%                     tempMark = strcat("interaction ", "p<.001");
%                     legendLabels = [legendLabels; tempMark];
%                 end
                if ~isempty(bint)
                    legendMarkers = [legendMarkers; bint];
                    tempMark = strcat("interaction ", thresholdName);
                    legendLabels = [legendLabels; tempMark];
                end
         
                if ~isempty(aa) || ~isempty(areg) || ~isempty(breg) || ~isempty(aint) || ~isempty(bint)
                        legend(legendMarkers, legendLabels, 'Location','best')
                end
        
           
%                 thresholdName = strcat('p<',num2str(threshold));
%                 if ~isempty(a) && ~isempty(b) && ~isempty(c)
%                     legend([a,b,aa],'p<.001', thresholdName,'main effect STIM', 'Location','best')
%                 elseif ~isempty(a) && ~isempty(b) 
%                     legend([a,b],'p<.001', thresholdName, 'Location','best')
%                 elseif ~isempty(a)
%                     legend(a,'p<.001','Location','best')
%                 elseif ~isempty(b)
%                     legend(b, thresholdName,'Location','best')
%                 elseif ~isempty(aa)
%                     legend(aa,'main effect stim','Location','best')
%                 end
            tablename = append ('2anovaTableRegion', name, '_dbh.csv');
            writecell(T2, tablename)

             tablename = append ('2anovaTableRegion', name, '_chat.csv');
            writecell(T, tablename)

            tablename = append ('2anovaTableRegion_postHoc', name, '_dbh.csv');
            writetable(dbhPostHocRegion, tablename)

                 tablename = append ('2anovaTableRegion_postHoc', name, '_chat.csv');
                writetable(chatPostHocRegion, tablename)
    end 
            sg = sgtitle(name);
            sg.FontSize = 36;
            
%             ylim([0 1])
    
            filename = strcat('2regions_', fileAppend, name);
            saveas(fig, filename, 'tif')
            saveas(fig, filename, 'fig')
    end 

end 

%%
strains = fields(dat);
field = fields(dat.chat.LM2);
regions = fields(dat.chat);

for r = 1:size(strains, 1)
        strain = strains{r};
    for i = 1:size(regions,1)
        region = regions{i};
        nTrial.(region){r} = cell2mat(dat.(strain).(region).nTrial.CTRL) + cell2mat(dat.(strain).(region).nTrial.STIM);
        percentStim.(region){r} = cell2mat(dat.(strain).(region).nTrial.STIM) ./ cell2mat(dat.(strain).(region).nTrial.CTRL);
        nleft.(region){r} = cell2mat(dat.(strain).(region).left.STIM) + cell2mat(dat.(strain).(region).left.CTRL);
        nright.(region){r} = cell2mat(dat.(strain).(region).right.STIM) + cell2mat(dat.(strain).(region).right.CTRL);
        nmiss.(region){r} = cell2mat(dat.(strain).(region).miss.STIM) + cell2mat(dat.(strain).(region).miss.CTRL);
        nReward.(region){r} = cell2mat(dat.(strain).(region).nwin.STIM) + cell2mat(dat.(strain).(region).nwin.CTRL);
        rewardRate.(region){r} = nReward.(region){r} ./ nTrial.(region){r};
        pLeft.(region){r} = nleft.(region){r} ./ nTrial.(region){r};
        pRight.(region){r} = nright.(region){r} ./ nTrial.(region){r};
    end 
end


%% plot WSLS data as a function of region and genotype with stim and no stim as separate plots
%choose plotting format 
    % 1 = raw data
    % 2 = percent change from no stim trials
    % 3 = percent change from average of control regions (LV1 and RV1)
        
        d = 1; 



strains = fields(dat);
field = fields(dat.chat.LM2);
regions = fields(dat.chat);

for i = 1:size(field,1)
    name = field{i};
    if strcmp(name, 'stats') == 1 || strcmp(name, 'nTrial') == 1 || strcmp(name, 'left') == 1 ...
            || strcmp(name, 'right') == 1 || strcmp(name, 'miss') == 1 || strcmp(name, 'nwin') == 1
    
    else
        fig = figure;
        for r = 1:size(regions,1)
            region = regions{r};
            chaty = ones(1, size(dat.chat.(region).(name).CTRL,2));
            rand_chat = 0.3*rand(1, size(dat.chat.(region).(name).CTRL,2));
            dbhy = ones(1, size(dat.dbh.(region).(name).CTRL,2));
            rand_dbh = 0.3*rand(1, size(dat.dbh.(region).(name).CTRL,2));
            sz = 100;
        
            subplot(2,2,1)
            hold on
            c = [0.9290 0.6940 0.1250];
            if d == 1
                fileAppend = 'rawData_';
                subplot(2,2,1)
                hold on
                c = [0.9290 0.6940 0.1250];
                scatter(r*chaty+rand_chat, dat.chat.(region).(name).STIM, sz, c)
                scatter(r, mean(dat.chat.(region).(name).STIM), (sz*2),c, 'filled')

                ylim([0 1])
                ylabel('pTrials')
                xlim ([0 5])
                xticks([1 2 3 4])
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('ChAt w/ stim')

                subplot(2,2,2)
                hold on
                scatter(r*chaty+rand_chat, dat.chat.(region).(name).CTRL, sz, c)
                scatter(r, mean(dat.chat.(region).(name).CTRL), (sz*2),c, 'filled')

                xlim ([0 5])
                ylim([0 1])
                xticks([1 2 3 4])
                ylabel('pTrials')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('ChAt no stim')

                subplot(2,2,3)
                c = [.1 0.6 0.6];
                hold on
                scatter(r*dbhy+rand_dbh, dat.dbh.(region).(name).STIM, sz,c)
               
                scatter(r, mean(dat.dbh.(region).(name).STIM), (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([0 1])
                xticks([1 2 3 4])
                ylabel('pTrials')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('DBH w/ stim')
            
                subplot(2,2,4)
                hold on
                scatter(r*dbhy+rand_dbh, dat.dbh.(region).(name).CTRL, sz,c)
                scatter(r, mean(dat.dbh.(region).(name).CTRL), (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([0 1])
                xticks([1 2 3 4])
                ylabel('pTrials')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('DBH no stim')



    
            elseif d == 2
                low = -100;
                high = 100;
                fileAppend = 'pChangeNoStim_';
                subplot(2,2,1)
                hold on
                c = [0.9290 0.6940 0.1250];
                scatter(r*chaty+rand_chat, (dat.chat.(region).(name).STIM - dat.chat.(region).(name).CTRL) ./ dat.chat.(region).(name).CTRL * 100, sz, c)
                scatter(r, mean((dat.chat.(region).(name).STIM - dat.chat.(region).(name).CTRL) ./ dat.chat.(region).(name).CTRL) * 100, (sz*2),c, 'filled')

                ylim([low high])
                ylabel('%ChangeNoStim')
                xlim ([0 5])
                xticks([1 2 3 4])
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('ChAt w/ stim')

                subplot(2,2,2)
                hold on
                scatter(r*chaty+rand_chat, (dat.chat.(region).(name).CTRL - dat.chat.(region).(name).CTRL) ./ dat.chat.(region).(name).CTRL * 100, sz, c)
                scatter(r, mean((dat.chat.(region).(name).CTRL - dat.chat.(region).(name).CTRL) ./ dat.chat.(region).(name).CTRL) * 100, (sz*2),c, 'filled')

                xlim ([0 5])
                ylim([low high])
                xticks([1 2 3 4])
                ylabel('%ChangeNoStim')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('ChAt no stim')

                subplot(2,2,3)
                c = [.1 0.6 0.6];
                hold on
                scatter(r*dbhy+rand_dbh, (dat.dbh.(region).(name).STIM - dat.dbh.(region).(name).CTRL) ./ dat.dbh.(region).(name).CTRL * 100, sz,c)
                scatter(r, mean((dat.dbh.(region).(name).STIM - dat.dbh.(region).(name).CTRL) ./ dat.dbh.(region).(name).CTRL) * 100, (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([low high])
                xticks([1 2 3 4])
                ylabel('%ChangeNoStim')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('DBH w/ stim')
            
                subplot(2,2,4)
                hold on
                scatter(r*dbhy+rand_dbh, (dat.dbh.(region).(name).CTRL - dat.dbh.(region).(name).CTRL) ./ dat.dbh.(region).(name).CTRL * 100, sz,c)
                scatter(r, mean((dat.dbh.(region).(name).CTRL - dat.dbh.(region).(name).CTRL) ./ dat.dbh.(region).(name).CTRL) * 100, (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([low high])
                xticks([1 2 3 4])
                ylabel('%ChangeNoStim')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('DBH no stim')
            elseif d == 3
                fileAppend = 'pChangeV1_';
                bottom = -100;
                top = 100;

                chat_v1 = dat.chat.LV1.(name).STIM + dat.chat.RV1.(name).STIM ./ 2;
                dbh_v1 = dat.dbh.LV1.(name).STIM + dat.dbh.RV1.(name).STIM ./ 2;

                chat_v1_ctrl = dat.chat.LV1.(name).CTRL + dat.chat.RV1.(name).CTRL ./ 2;
                dbh_v1_ctrl = dat.dbh.LV1.(name).CTRL + dat.dbh.RV1.(name).CTRL ./ 2;


                subplot(2,2,1)
                hold on
                c = [0.9290 0.6940 0.1250];
                scatter(r*chaty+rand_chat, (dat.chat.(region).(name).STIM - chat_v1) ./ chat_v1 * 100, sz, c)
                scatter(r, mean((dat.chat.(region).(name).STIM - chat_v1) ./ chat_v1) * 100, (sz*2),c, 'filled')

                ylim([bottom top])
                ylabel('%ChangeV1')
                xlim ([0 5])
                xticks([1 2 3 4])
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('ChAt w/ stim')

                subplot(2,2,2)
                hold on
                scatter(r*chaty+rand_chat, (dat.chat.(region).(name).CTRL - chat_v1_ctrl) ./ chat_v1_ctrl * 100, sz, c)
                scatter(r, mean((dat.chat.(region).(name).CTRL - chat_v1_ctrl) ./ chat_v1_ctrl) * 100, (sz*2),c, 'filled')

                xlim ([0 5])
                ylim([bottom top])
                xticks([1 2 3 4])
                ylabel('%ChangeV1')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('ChAt no stim')

                subplot(2,2,3)
                c = [.1 0.6 0.6];
                hold on
                scatter(r*dbhy+rand_dbh, (dat.dbh.(region).(name).STIM - dbh_v1) ./ dbh_v1 * 100, sz,c)
                scatter(r, mean((dat.dbh.(region).(name).STIM - dbh_v1) ./ dbh_v1) * 100, (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([bottom top])
                xticks([1 2 3 4])
                ylabel('%ChangeV1')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('DBH w/ stim')
            
                subplot(2,2,4)
                hold on
                scatter(r*dbhy+rand_dbh, (dat.dbh.(region).(name).CTRL - dbh_v1_ctrl) ./ dbh_v1_ctrl * 100, sz,c)
                scatter(r, mean((dat.dbh.(region).(name).CTRL - dbh_v1_ctrl) ./ dbh_v1_ctrl) * 100, (sz*2),c, 'filled')
                xlim ([0 5])
                ylim([bottom top])
                xticks([1 2 3 4])
                ylabel('%ChangeV1')
                xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1'])
                title('DBH no stim')
            end 
        end
    
        
        
            
            sg = sgtitle(name);
            sg.FontSize = 36;
%             ylim([0 1])
    
            filename = strcat('regions_', fileAppend, name);
            saveas(fig, filename, 'tif')
            saveas(fig, filename, 'fig')
    end 

end 

 %%
strains = fields(dat);
field = fields(dat.chat.LM2);
regions = fields(dat.chat);

for r = 1:size(strains, 1)
        strain = strains{r};
    for i = 1:size(regions,1)
        region = regions{i};
        nTrial.(region){r} = cell2mat(dat.(strain).(region).nTrial.CTRL) + cell2mat(dat.(strain).(region).nTrial.STIM);
        percentStim.(region){r} = cell2mat(dat.(strain).(region).nTrial.STIM) ./ cell2mat(dat.(strain).(region).nTrial.CTRL);
        nleft.(region){r} = cell2mat(dat.(strain).(region).left.STIM) + cell2mat(dat.(strain).(region).left.CTRL);
        nright.(region){r} = cell2mat(dat.(strain).(region).right.STIM) + cell2mat(dat.(strain).(region).right.CTRL);
        nmiss.(region){r} = cell2mat(dat.(strain).(region).miss.STIM) + cell2mat(dat.(strain).(region).miss.CTRL);
        nReward.(region){r} = cell2mat(dat.(strain).(region).nwin.STIM) + cell2mat(dat.(strain).(region).nwin.CTRL);
        rewardRate.(region){r} = nReward.(region){r} ./ nTrial.(region){r};
        pLeft.(region){r} = nleft.(region){r} ./ nTrial.(region){r};
        pRight.(region){r} = nright.(region){r} ./ nTrial.(region){r};
    end 
end


%%

%are 10% of trials stimulated
    %make plot using precentStim to show 10% of trials are stimulated
    regions = fields(percentStim);
    figure('position',[0 0 1500 800])
        for i = 1:size(regions,1)
            region = regions{i};
            dx = i*ones(1,size(percentStim.(region){2},2)) + 0.3*rand(1, size(percentStim.(region){2},2));
            cx = (i+4)*ones(1,size(percentStim.(region){1},2)) + 0.3*rand(1, size(percentStim.(region){1},2));

                dc = [.1 0.6 0.6];
                cc = [0.9290 0.6940 0.1250];


%             subplot(2,3,1)
%             hold on
%             scatter(cx, percentStim.(region){1},200,cc)
%             scatter(dx, percentStim.(region){2},200,dc)
% 
%             xticks([1 2 3 4 5 6 7 8])
%             xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1';'LM2'; 'RM2'; 'LV1'; 'RV1'])
%             xlim([0 10])
%             ylim([0 .5])
%             title('pTrials stimulated')

             subplot(1,3,2)
             hold on
            scatter(cx+1, nTrial.(region){1},100,cc)
            scatter(dx, nTrial.(region){2},100,dc)
            scatter(i+5, mean(nTrial.(region){1}),200, cc, 'filled')
            scatter(i, mean(nTrial.(region){2}),200,dc,'filled')

            errorbar(i+5, mean( nTrial.(region){1}, 'omitnan'), std(nTrial.(region){1},'omitnan') / sqrt(size( nTrial.(region){1},2))*2,...
                             'color', cc);
            errorbar(i, mean(nTrial.(region){2}, 'omitnan'), std(nTrial.(region){2},'omitnan') / sqrt(size(nTrial.(region){2},2))*2,...
                             'color', dc);

       
                                    ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTickLabelRotation = 90;
            xticks([1 2 3 4  6 7 8 9])
            xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R';'M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
            xlim([0 10])
            ylim([0 1000])
            yticks([0 200 400 600 800 1000])
            title('Number of Trials')

            subplot(1,3,3)
            hold on
            scatter(cx+1, rewardRate.(region){1},100,cc)
            scatter(i+5, mean(rewardRate.(region){1}),200, cc, 'filled')
            scatter(dx, rewardRate.(region){2},100,dc)
            scatter(i, mean(rewardRate.(region){2}),200,dc,'filled')


            errorbar(i+5, mean(rewardRate.(region){1}, 'omitnan'), std(rewardRate.(region){1},'omitnan') / sqrt(size(rewardRate.(region){1},2))*2,...
                             'color', cc);
            errorbar(i, mean(rewardRate.(region){2}, 'omitnan'), std(rewardRate.(region){2},'omitnan') / sqrt(size(rewardRate.(region){2},2))*2,...
                             'color', dc);

                                    ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTickLabelRotation = 90;
            xticks([1 2 3 4  6 7 8 9])
            xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R';'M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
            xlim([0 10])
            ylim([0 .6])
            yticks([0 .5])
            yticklabels(["0%", "50%" "100%"])
            title('Reward Rate')

% 
%             subplot(1,3,3)
%             hold on
%         
%             scatter(cx, nmiss.(region){1},100,cc)
%             scatter(dx, nmiss.(region){2},100,dc)
%             scatter(i, mean(nmiss.(region){1}),200,cc,'filled')
%             scatter(i+4, mean(nmiss.(region){2}),200,dc,'filled')
%             xticks([1 2 3 4 5 6 7 8])
%             xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1';'LM2'; 'RM2'; 'LV1'; 'RV1'])
%             xlim([0 10])
%             ylim([0 100])
%             yticks([0 50 100])
%             title('Number of Misses')
    
% 
%             subplot(2,3,5)
%             hold on
%             scatter(cx, pLeft.(region){1},200,cc)
%             scatter(dx, pLeft.(region){2},200,dc)
%             xticks([1 2 3 4 5 6 7 8])
%             xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1';'LM2'; 'RM2'; 'LV1'; 'RV1'])
%             xlim([0 10])
%             ylim([0 1])
%             title('pLeft Choices')
%     
%             subplot(2,3,6)
%             hold on
%             scatter(cx, pRight.(region){1},200,cc)
%             scatter(dx, pRight.(region){2},200,dc)
%             xticks([1 2 3 4 5 6 7 8])
%             xticklabels(['LM2'; 'RM2'; 'LV1'; 'RV1';'LM2'; 'RM2'; 'LV1'; 'RV1'])
%             xlim([0 10])
%             ylim([0 1])
%             title('pRight choices')
        end 

%         figure
        for i = 1:size(regions,1)
            subplot(1,3,1)
            region = regions{i}; 
            hold on 
            dbhIdx = ismember(string(stats.(region).Strain), "Dbh");
            chatIdx = ismember(string(stats.(region).Strain), "ChAt");

            cx = (i+4)*ones(1,size(stats.(region).Entropy(chatIdx),1)) + 0.3*rand(1, size(stats.(region).Entropy(chatIdx),1));
            dx = (i)*ones(1,size(stats.(region).Entropy(dbhIdx),1)) + 0.3*rand(1, size(stats.(region).Entropy(dbhIdx),1));

            scatter(cx+1, stats.(region).Entropy(chatIdx),100,cc)
            scatter(dx, stats.(region).Entropy(dbhIdx),100,dc)
            scatter(i+5, mean(stats.(region).Entropy(chatIdx), 'omitnan'),200, cc, 'filled')
            scatter(i, mean(stats.(region).Entropy(dbhIdx), 'omitnan'),200,dc,'filled')

            errorbar(i+5, mean(stats.(region).Entropy(chatIdx), 'omitnan'), std(stats.(region).Entropy(chatIdx),'omitnan') / sqrt(size(stats.(region).Entropy(chatIdx),2))*2,...
                             'color', cc);
            errorbar(i, mean(stats.(region).Entropy(dbhIdx), 'omitnan'), std(stats.(region).Entropy(dbhIdx),'omitnan') / sqrt(size(stats.(region).Entropy(dbhIdx),2))*2,...
                             'color', dc);

            
                        ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTickLabelRotation = 90;
            xticks([1 2 3 4  6 7 8 9])
            xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R';'M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
            xlim([0 10])
            ylim([0 3.5])
            yticks([0 1 2 3])
            title('Entropy')


        end 

        %%
clearvars;
close all;
setup_figprop;
regionList = ['LM2'; 'RM2'; 'LV1'; 'RV1'];

fullIndex = [];

for k = 1:size(regionList,1)
        currReg = regionList(k,:);
        root_path = 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\1Region';

        logfilepath = fullfile(root_path,currReg,'data');
        analysispath = fullfile(root_path, currReg, 'analysis');
        dataIndex = makeDataIndex(logfilepath, analysispath);
        dataIndex = MP_STIM_createBehMatFiles(dataIndex);
        animalList = unique(dataIndex.Animal);

       for jj=1:size(dataIndex,1)
            name = dataIndex.LogFileName{jj};
            name = name(1:end-4);
            %load(fullfile(animalData.BehPath{jj}, 'beh_cut.mat'));
            name = strcat(name, '_beh_cut.mat');
            load(fullfile(dataIndex.BehPath{jj}, name));
            dataIndex.nTrials(jj) = size(trialData.outcome,1);
            dataIndex.entropy(jj) = entro;
            dataIndex.region{jj} = currReg;
            dataIndex.rewardRate(jj) = (sum(trialData.outcome == 100) + sum(trialData.outcome == 111)) / (dataIndex.nTrials(jj) - sum(trialData.outcome == 77));
       end
       fullIndex.(currReg) = dataIndex;
        clear dataIndex
end 
       
%%
regionList = ["LM2"; "RM2"; "LV1"; "RV1"];
    figure('position',[0 0 1500 800])

                trials.dbh = [];
            trials.chat = [];
            rewardRate.dbh = [];
            rewardRate.chat = [];
            entropy.dbh = [];
            entropy.chat = [];
            grouping.dbh = [];
            grouping.chat = [];
        for i = 1:size(regionList,1)
            region = regionList(i,:);
            dbhTrials = fullIndex.(region).nTrials(fullIndex.(region).Strain == "Dbh");
            chatTrials = fullIndex.(region).nTrials(fullIndex.(region).Strain == "ChAt");
            dbhRR = fullIndex.(region).rewardRate(fullIndex.(region).Strain == "Dbh");
            chatRR = fullIndex.(region).rewardRate(fullIndex.(region).Strain == "ChAt");
            dbhE = fullIndex.(region).entropy(fullIndex.(region).Strain == "Dbh");
            chatE = fullIndex.(region).entropy(fullIndex.(region).Strain == "ChAt");
            
            trials.dbh = [trials.dbh; fullIndex.(region).nTrials(fullIndex.(region).Strain == "Dbh")];
            trials.chat = [trials.chat; fullIndex.(region).nTrials(fullIndex.(region).Strain == "ChAt")];
            rewardRate.dbh = [rewardRate.dbh; fullIndex.(region).rewardRate(fullIndex.(region).Strain == "Dbh")];
            rewardRate.chat = [rewardRate.chat; fullIndex.(region).rewardRate(fullIndex.(region).Strain == "ChAt")];
            entropy.dbh = [entropy.dbh; fullIndex.(region).entropy(fullIndex.(region).Strain == "Dbh")];
            entropy.chat = [entropy.chat; fullIndex.(region).entropy(fullIndex.(region).Strain == "ChAt")];
            grouping.dbh = [grouping.dbh; repmat(region, [size(fullIndex.(region).nTrials(fullIndex.(region).Strain == "Dbh"),1),1])];
            grouping.chat = [grouping.chat; repmat(region, [size(fullIndex.(region).nTrials(fullIndex.(region).Strain == "ChAt"),1),1])];
            dx = i*ones(1,size(dbhTrials,2)) + 0.3*rand(1, size(dbhTrials,2));
            cx = (i+4)*ones(1,size(chatTrials,2)) + 0.3*rand(1, size(chatTrials,2));

                dc = [.1 0.6 0.6];
                cc = [0.9290 0.6940 0.1250];

             subplot(1,3,2)
             hold on
            scatter(cx+1, chatTrials,100,cc)
            scatter(dx, dbhTrials,100,dc)
            scatter(i+5, mean(chatTrials, 'omitnan'),200, cc, 'filled')
            scatter(i, mean(dbhTrials, 'omitnan'),200,dc,'filled')

            errorbar(i+5, mean(chatTrials, 'omitnan'), std(chatTrials,'omitnan') / sqrt(size(chatTrials,2)),...
                             'color', cc, 'LineWidth', 1.5);
            errorbar(i, mean(dbhTrials, 'omitnan'), std(dbhTrials,'omitnan') / sqrt(size(dbhTrials,2)),...
                             'color', dc, 'LineWidth', 1.5);

       
                                    ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTickLabelRotation = 90;
                        ax.LineWidth = 5;
            xticks([1 2 3 4  6 7 8 9])
            xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R';'M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
            xlim([0 10])
            ylim([0 1000])
            yticks([0 200 400 600 800 1000])
            title('Number of Trials')

            subplot(1,3,3)
            hold on
            scatter(cx+1, chatRR,100,cc)
            scatter(i+5, mean(chatRR, 'omitnan'),200, cc, 'filled')
            scatter(dx, dbhRR,100,dc)
            scatter(i, mean(dbhRR, 'omitnan'),200,dc,'filled')


            errorbar(i+5, mean(chatRR, 'omitnan'), std(chatRR,'omitnan') / sqrt(size(chatRR,2)),...
                             'color', cc, 'LineWidth', 1.5);
            errorbar(i, mean(dbhRR, 'omitnan'), std(dbhRR,'omitnan') / sqrt(size(dbhRR,2)),...
                             'color', dc, 'LineWidth', 1.5);

                        ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTickLabelRotation = 90;
                        ax.LineWidth = 5;
            xticks([1 2 3 4  6 7 8 9])
            xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R';'M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
            xlim([0 10])
            ylim([0 .6])
            yticks([0 .5])
            yticklabels(["0%", "50%" "100%"])
            title('Reward Rate')

            subplot(1,3,1)
            hold on 

            scatter(cx+1, chatE,100,cc)
            scatter(dx, dbhE,100,dc)
            scatter(i+5, mean(chatE, 'omitnan'),200, cc, 'filled')
            scatter(i, mean(dbhE, 'omitnan'),200,dc,'filled')

            errorbar(i+5, mean(chatE, 'omitnan'), std(chatE,'omitnan') / sqrt(size(chatE,2)),...
                             'color', cc, 'LineWidth', 1.5);
            errorbar(i, mean(dbhE, 'omitnan'), std(dbhE,'omitnan') / sqrt(size(dbhE,2)),...
                             'color', dc, 'LineWidth', 1.5);

            
                        ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTickLabelRotation = 90;
                        ax.LineWidth = 5;
            xticks([1 2 3 4  6 7 8 9])
            xticklabels(['M2-L'; 'M2-R'; 'V1-L'; 'V1-R';'M2-L'; 'M2-R'; 'V1-L'; 'V1-R'])
            xlim([0 10])
            ylim([0 3.5])
            yticks([0 1 2 3])
            title('Entropy')


        end 

[p,tbl,stats] = anova1(trials.dbh, grouping.dbh);
            writecell(tbl, 'nTrials_dbh_regions.csv')
[p,tbl,stats] = anova1(trials.chat, grouping.chat);
    writecell(tbl, 'nTrials_chat_regions.csv')

[p,tbl,stats] = anova1(rewardRate.dbh, grouping.dbh);
    writecell(tbl, 'rewardRate_dbh_regions.csv')
[p,tbl,stats] = anova1(rewardRate.chat, grouping.chat);
    writecell(tbl, 'rewardRate_chat_regions.csv')
[p,tbl,stats] = anova1(entropy.dbh, grouping.dbh);
    writecell(tbl, 'entropy_dbh_regions.csv')
[p,tbl,stats] = anova1(entropy.chat, grouping.chat);
    writecell(tbl, 'entropy_chat_regions.csv')