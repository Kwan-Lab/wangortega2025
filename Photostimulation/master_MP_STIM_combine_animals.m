clearvars;
close all;
setup_figprop;
powerList = ['0.0mW'; '1.5mW'; '3.0mW'];

sidePreference = 0; %1 to analyze based on side preference per session
for k = 1:size(powerList,1)
            currPow = powerList(k,:);
    
    
    root_path = 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim';
    
    % matching pennies behavior
    
    disp('-----------------------------------------------------------');
    disp('--- HERE WE GO ');
    disp('--------------------------------------------------open---------');
    
    % Look for data files and create a database index
    
    % create data indexes 
    
        logfilepath = fullfile(root_path,currPow,'data');
        analysispath = fullfile(root_path,currPow,'analysis');
        dataIndex = makeDataIndex(logfilepath, analysispath);
        
        % Parse and analyze each logfile, save as .mat files, as needed
        dataIndex = MP_STIM_createBehMatFiles(dataIndex);
        
        % sort Index according to experiment data
        
        dataIndex = sortdataIndex(dataIndex);
        
        % Determine if each session fulfill performance criteria (Matching Pennies)
        MP_determineBehCriteria(dataIndex);
    
        % start separate analysis for each strain
        clear sessionData stats trialData
        nFiles = size(dataIndex(:,:),1);
        sz = [0 4];
        varTypes = ["string","double","double","double"];
        varNames = ["logFileName","Entropy","nTrials","right/total"];
        lowEntropyTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    
        for ii = 1:nFiles
            savematpath = dataIndex.BehPath{ii};
             [lowEntropyTable] = MP_STIM_session(dataIndex.BehPath{ii},dataIndex.LogFileName{ii},savematpath, lowEntropyTable);
        end
        save_path = fullfile(root_path,'summary','figs_summary');
        
        
        MP_GRAB_behaviorPerAnimal(dataIndex,save_path);
        
        MP_GRAB_behaviorAll(dataIndex, save_path);
   %
        % Calculate stats by stimulation region - by session
        % Animal
        regions = {'LM2','RM2','LV1','RV1','CTRL'};
        plotChoice = 0; %make 1 to plot session data
        for i = 1:2 
            if i ==1
                dataIndex_curr = dataIndex(strcmp(dataIndex.Strain, 'Dbh'),:);
                animalList = unique(dataIndex_curr.Animal); 
                for j = 1:size(animalList,1) %change to loop through animal then sessions, make animals separate rows
                    animal = animalList(j); 
                    tempIndex = dataIndex_curr(strcmp(dataIndex_curr.Animal, animal),:);
                    stats_dbh{1,j} = MP_STIM_stats(analysispath,tempIndex,plotChoice,regions,'',sidePreference);    
                end
                clear sessionData stats trialData
                    save_path = fullfile(root_path,'summary','figs_summary');
            else
                dataIndex_curr = dataIndex(strcmp(dataIndex.Strain, 'ChAt'),:);
                animalList = unique(dataIndex_curr.Animal); 
                for j = 1:size(animalList,1) %change to loop through animal then sessions, make animals separate rows
                    animal = animalList(j); 
                    tempIndex = dataIndex_curr(strcmp(dataIndex_curr.Animal, animal),:);
                        stats_chat{1,j} = MP_STIM_stats(analysispath,tempIndex,plotChoice,regions,'', sidePreference);
                end
                clear sessionData stats trialData
                    save_path = fullfile(root_path,'summary','figs_summary');
            end
        end
        % start plotting
    
        variables = fieldnames(stats_chat{1, 1}(1).beh_output_WSLS{1, 1}); 
        %regions = {'LM2', 'CTRL'};
    
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
    cd 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim';
    filename = append('powers_', currPow, '.mat');
    save(filename,'WSLS_chat','WSLS_dbh','lowEntropyTable','dataIndex')
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
cd 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim';
load("powers_0.0mW.mat")
dat.chat.zero = WSLS_chat;
dat.dbh.zero = WSLS_dbh;
lowEntropyTable.Strain = dataIndex.Strain;
stats.zero = lowEntropyTable;

load("powers_1.5mW.mat")
dat.chat.onefive = WSLS_chat;
dat.dbh.onefive = WSLS_dbh;
lowEntropyTable.Strain = dataIndex.Strain;
stats.onefive = lowEntropyTable;

load("powers_3.0mW.mat")
dat.chat.three = WSLS_chat;
dat.dbh.three = WSLS_dbh;
lowEntropyTable.Strain = dataIndex.Strain;
stats.three = lowEntropyTable;

%% test mixed effect linear regression 
field = fields(dat.chat.zero);
powers = fields(dat.chat);
strains = fields(dat);
regions = fields(dat.chat.zero.pWSLS);
cc = [0.9290 0.6940 0.1250]; %chat RBG color
dc = [.1 0.6 0.6]; %dbh RBG color
sz = 10;

for i = 1:size(field,1)
     name = field{i};
        %fig = figure;
        %f = sgtitle(name);
        %f.FontSize = 25;

                table_dbh = table('size', [0 5], 'VariableTypes', {'double','double','string','string','string'}, 'VariableNames',{'ID',name, 'genotype','power','region'});
                table_chat = table('size', [0 5], 'VariableTypes', {'double','double','string','string','string'}, 'VariableNames',{'ID',name, 'genotype','power','region'});
               
                chatStr = '';
                dbhStr = '';
                
        for k = 1:size(powers,1)
                power = powers{k};

                for j = 1:size(regions,1)
                    region = regions{j};
                    ndbh = size(dat.dbh.(power).(name).(region),2);
                    nchat = size(dat.chat.(power).(name).(region),2);
                    temp_dbh = table('size', [ndbh 2], 'VariableTypes', {'double','double'}, 'VariableNames',{'ID',name});
                    temp_chat = table('size', [nchat 2], 'VariableTypes', {'double','double'}, 'VariableNames',{'ID',name});

    
                    temp_chat.ID(:,1) = 1:nchat;
                    temp_dbh.ID(:,1) = 10*(1:ndbh);
    
                    temp_chat.genotype(1:nchat,1) = {'chat'};
                    temp_dbh.genotype(1:ndbh,1) = {'dbh'};

                    temp_chat.(name)(1:nchat,1) = dat.chat.(power).(name).(region);
                    temp_dbh.(name)(1:ndbh,1) = dat.dbh.(power).(name).(region);
                    
                    temp_chat.power(:,1) = {power};
                    temp_chat.region(:,1) = {region};
                    temp_dbh.power(:,1) = {power};
                    temp_dbh.region(:,1) = {region};

                    table_dbh = [temp_dbh; table_dbh];
                    table_chat = [temp_chat; table_chat];
                end
                 

        end
        modelEQ = append(name, ' ~ 1 + power + region + (1|ID)');
       
        fitlme(table_chat, modelEQ)
        fitlme(table_dbh, modelEQ)
        name



                
%                 if ranovatbl_chat.pValue(1) <= 0.05 %if there is a sig difference in time
%                    % do a post hoc test
%                         [H1,P1] = ttest(table_chat.zero, table_chat.onefive, 'Alpha', 0.05);
%                         [H2,P2] = ttest(table_chat.zero, table_chat.three, 'Alpha', 0.05);
%                         [H3,P3] = ttest(table_chat.three, table_chat.onefive, 'Alpha', 0.05);
%                         P1,P2,P3
% 
% %                         if P1 <= 0.05
% %                             subplot(3,2,k)
% %                             plot([1,2],[1,1], 'k')
% %                             %plot([1.5],[.95],'k*')
% %                         end
% %                         
% %                         if P2 <= 0.05
% %                             subplot(3,2,k)
% %                             plot([1,3],[.9,.9], 'k')
% %                             %plot([1.5],[.95],'k*')
% %                         end
% % 
% %                         if P3 <= 0.05
% %                             subplot(3,2,k)
% %                             plot([2,3],[.95,.95], 'k')
% %                         end
% 
%                 end
% 
%                 if ranovatbl_dbh.pValue(1) <= 0.05 %if there is a sig difference in this region for this variable
%                    % do a post hoc test
%                         [H1,P1] = ttest(table_dbh.zero, table_dbh.onefive, 'Alpha', 0.05);
%                         [H2,P2] = ttest(table_dbh.zero, table_dbh.three, 'Alpha', 0.05);
%                         [H3,P3] = ttest(table_dbh.three, table_dbh.onefive, 'Alpha', 0.05);
%                         P1,P2,P3
% %                         if P1 <= 0.05
% %                             subplot(3,2,k)
% %                             plot([6,7],[1,1], 'k')
% %                             %plot([1.5],[.95],'k*')
% %                         end
% %                         
% %                         if P2 <= 0.05
% %                             subplot(3,2,k)
% %                             plot([6,8],[.9,.9], 'k')
% %                             %plot([1.5],[.95],'k*')
% %                         end
% % 
% %                         if P3 <= 0.05
% %                             subplot(3,2,k)
% %                             plot([7,8],[.95,.95], 'k')
% %                         end





    fig = figure
        for j = 1:size(powers,1)
            power = powers{j};

            chatx = [1*ones(1,size(dat.chat.(power).pWSLS.LM2,2));2*ones(1,size(dat.chat.(power).pWSLS.LM2,2));...
            3*ones(1,size(dat.chat.(power).pWSLS.LM2,2))];
            dbhx = [6*ones(1,size(dat.dbh.(power).pWSLS.LM2,2));7*ones(1,size(dat.dbh.(power).pWSLS.LM2,2));...
            8*ones(1,size(dat.dbh.(power).pWSLS.LM2,2))];

             rand_chat = 0.3*rand(1, size(dat.chat.(power).pWSLS.LM2,2));
             rand_dbh = 0.3*rand(1, size(dat.dbh.(power).pWSLS.LM2,2));
            for k = 1:size(regions,1)
                region = regions{k};
                subplot (3,2,k)
                hold on
                chat_ind = plot(chatx(j,:)+rand_chat, dat.chat.(power).(name).(region),'o',...
                                    'Color',cc, ...
                                    'MarkerSize',sz,...
                                    'LineWidth', .2);
                dbh_ind = plot(dbhx(j,:)+rand_dbh, dat.dbh.(power).(name).(region),'o',...
                                    'Color', dc, ...
                                    'MarkerSize',sz,...
                                    'LineWidth', .2);
                dbh_mean = plot(dbhx(j,1), mean(dat.dbh.(power).(name).(region)),'.',...
                                    'Color', dc, ...
                                    'MarkerSize',40,...
                                    'LineWidth', .2);
                chat_mean = plot(chatx(j,1), mean(dat.chat.(power).(name).(region)),'.',...
                                    'Color', cc, ...
                                    'MarkerSize',40,...
                                    'LineWidth', .2);
                ax = gca;
                ax.TickLength = [0 0];
                ax.YLim = [0 1];
                ax.XLim = [0 9];
                ax.XTick = [1 2 3 6 7 8];
                ax.YTick = [0 1];
                ax.XTickLabel = ['0.0'; '1.5'; '3.0'; '0.0'; '1.5';'3.0'];
                ax.XTickLabelRotation = 45;
                ax.TitleFontSizeMultiplier = 0.8;
                title(region) 
            end
           
        end
        lg = legend([chat_mean chat_ind dbh_mean dbh_ind], 'ChAt Mean', 'ChAt Individual', 'Dbh Mean', 'Dbh Individual');
        lg.Position (1:2) = [.7 .1];
        sgtitle(name)

                

        filename = strcat('powers_', name);
        saveas(fig, filename, 'tif')
        saveas(fig, filename, 'fig')
end 

%% test anova
field = fields(dat.chat.zero);
powers = fields(dat.chat);
strains = fields(dat);
regions = fields(dat.chat.zero.pWSLS);
cc = [0.9290 0.6940 0.1250]; %chat RBG color
dc = [.1 0.6 0.6]; %dbh RBG color
sz = 100;

pval = .05; %set p value threshold
testForSig = 0; %test for significance? 

for i = 1:size(field,1)
     name = field{i};
     if strcmp(name, "stats") == 1 || strcmp(name, "nTrial") == 1 || strcmp(name, "left") == 1 || strcmp(name, "right") == 1 || strcmp(name, "miss") == 1 || strcmp(name, "nwin") == 1 
     else
            %fig = figure;
            %f = sgtitle(name);
            %f.FontSize = 25;
    
                    table_dbh = table('size', [0 5], 'VariableTypes', {'double','double','string','string','string'}, 'VariableNames',{'ID',name, 'genotype','power','region'});
                    table_chat = table('size', [0 5], 'VariableTypes', {'double','double','string','string','string'}, 'VariableNames',{'ID',name, 'genotype','power','region'});
                   
                    chatStr = '';
                    dbhStr = '';
                    
            for k = 1:size(powers,1)
                    power = powers{k};
    
                    for j = 1:size(regions,1)
                        region = regions{j};
                        ndbh = size(dat.dbh.(power).(name).(region),2);
                        nchat = size(dat.chat.(power).(name).(region),2);
                        temp_dbh = table('size', [ndbh 2], 'VariableTypes', {'double','double'}, 'VariableNames',{'ID',name});
                        temp_chat = table('size', [nchat 2], 'VariableTypes', {'double','double'}, 'VariableNames',{'ID',name});
    
        
                        temp_chat.ID(:,1) = 1:nchat;
                        temp_dbh.ID(:,1) = 10*(1:ndbh);
        
                        temp_chat.genotype(1:nchat,1) = {'chat'};
                        temp_dbh.genotype(1:ndbh,1) = {'dbh'};
    
                        temp_chat.(name)(1:nchat,1) = dat.chat.(power).(name).(region);
                        temp_dbh.(name)(1:ndbh,1) = dat.dbh.(power).(name).(region);
                        
                        temp_chat.power(:,1) = {power};
                        temp_chat.region(:,1) = {region};
                        temp_dbh.power(:,1) = {power};
                        temp_dbh.region(:,1) = {region};
    
                        table_dbh = [temp_dbh; table_dbh];
                        table_chat = [temp_chat; table_chat];
                    end
                     
    
            end
            
            
    %                 if ranovatbl_chat.pValue(1) <= 0.05 %if there is a sig difference in time
    %                    % do a post hoc test
    %                         [H1,P1] = ttest(table_chat.zero, table_chat.onefive, 'Alpha', 0.05);
    %                         [H2,P2] = ttest(table_chat.zero, table_chat.three, 'Alpha', 0.05);
    %                         [H3,P3] = ttest(table_chat.three, table_chat.onefive, 'Alpha', 0.05);
    %                         P1,P2,P3
    % 
    % %                         if P1 <= 0.05
    % %                             subplot(3,2,k)
    % %                             plot([1,2],[1,1], 'k')
    % %                             %plot([1.5],[.95],'k*')
    % %                         end
    % %                         
    % %                         if P2 <= 0.05
    % %                             subplot(3,2,k)
    % %                             plot([1,3],[.9,.9], 'k')
    % %                             %plot([1.5],[.95],'k*')
    % %                         end
    % % 
    % %                         if P3 <= 0.05
    % %                             subplot(3,2,k)
    % %                             plot([2,3],[.95,.95], 'k')
    % %                         end
    % 
    %                 end
    % 
    %                 if ranovatbl_dbh.pValue(1) <= 0.05 %if there is a sig difference in this region for this variable
    %                    % do a post hoc test
    %                         [H1,P1] = ttest(table_dbh.zero, table_dbh.onefive, 'Alpha', 0.05);
    %                         [H2,P2] = ttest(table_dbh.zero, table_dbh.three, 'Alpha', 0.05);
    %                         [H3,P3] = ttest(table_dbh.three, table_dbh.onefive, 'Alpha', 0.05);
    %                         P1,P2,P3
    % %                         if P1 <= 0.05
    % %                             subplot(3,2,k)
    % %                             plot([6,7],[1,1], 'k')
    % %                             %plot([1.5],[.95],'k*')
    % %                         end
    % %                         
    % %                         if P2 <= 0.05
    % %                             subplot(3,2,k)
    % %                             plot([6,8],[.9,.9], 'k')
    % %                             %plot([1.5],[.95],'k*')
    % %                         end
    % % 
    % %                         if P3 <= 0.05
    % %                             subplot(3,2,k)
    % %                             plot([7,8],[.95,.95], 'k')
    % %                         end
    
    
    
    
    setup_figprop
        fig = figure('Position',[0 0 2500 800])
        for rr = 1:size(strains,1)
            strain = strains{rr};
            if strcmp(strain, 'chat') == 1
                sp = 2;
            elseif strcmp(strain, 'dbh') == 1
                sp = 1;
            end
    
            for jj = 1:size(powers,1)
                power = powers{jj};
                
                if jj == 1
                        poweroffset = 0;
                elseif jj == 2
                    poweroffset = 8;
                else
                    poweroffset = 16;
                end 
                xcoord = [1*ones(1,size(dat.(strain).(power).pWSLS.LM2,2));2*ones(1,size(dat.(strain).(power).pWSLS.LM2,2));...
                3*ones(1,size(dat.(strain).(power).pWSLS.LM2,2)); 4*ones(1,size(dat.(strain).(power).pWSLS.LM2,2));...
                5*ones(1,size(dat.(strain).(power).pWSLS.LM2,2));];
              
                 xoffset = 0.3*rand(1, size(dat.(strain).(power).pWSLS.LM2,2));
    
                for k = 1:size(regions,1)
                    if k == 5
                        c = [.7 .7 .7];
                    elseif strcmp(strain, 'chat') == 1
                        c = cc;
                    else 
                        c = dc;
                    end 
                    region = regions{k};
                    subplot (1,2,sp)
                    hold on 
                    individual = scatter(2*(xcoord(k,:))+xoffset+2*poweroffset, dat.(strain).(power).(name).(region),sz,c);
                                     ;
                    hold on
                    group_mean = scatter(2*xcoord(k,1)+2*poweroffset, mean(dat.(strain).(power).(name).(region), 'omitnan'),2*sz,c,'filled');
    
                    hold on
                    sem = errorbar(2*xcoord(k,1)+2*poweroffset, mean(dat.(strain).(power).(name).(region), 'omitnan'), std(dat.(strain).(power).(name).(region), 'omitnan') / sqrt(size(dat.(strain).(power).(name).(region),2)),...
                             'color', c);
                end
               
            end
                    ax = gca;
                    ax.TickLength = [0 0];
                    ax.XLim = [0 50];
                    ax.XTick = [2:2:10, 18:2:26, 34:2:42];
                    ax.YLim = [.17 .6];
                    ax.YTick = [.2 .4 .6];
                    ax.YTickLabel = ["20", "40", "60"];
                    ax.XTickLabel = ["M2-L"; "M2-R"; "V1-L"; "V1-R"; "CTRL"; "M2-L"; "M2-R"; "V1-L"; "V1-R"; "CTRL";...
                        "M2-L"; "M2-R"; "V1-L"; "V1-R"; "CTRL"];
                    ax.XTickLabelRotation = 90;
                    ax.FontSize  = 25;
                    ax.TitleFontSizeMultiplier = 2;
                    %title(strain) 
        end 
    %         lg = legend([chat_mean chat_ind dbh_mean dbh_ind], 'ChAt Mean', 'ChAt Individual', 'Dbh Mean', 'Dbh Individual');
    %         lg.Position (1:2) = [.7 .1];
    
        if testForSig == 1
            factors_chat = {table_chat.power table_chat.region};
            [P1,T1,STATS1,TERMS1] = anovan(table_chat.(name), factors_chat,'varnames',{'power','region'}, 'display', 'off', 'model','interaction');
            name
            T1
    
            factors_dbh = {table_dbh.power table_dbh.region};
            [P2,T2,STATS2,TERMS2] = anovan(table_dbh.(name), factors_dbh,'varnames',{'power','region'}, 'display', 'off', 'model','interaction');
            T2
    
            if P1(1) < .05 %main effect of power
                n = 0;
                 [results,~,~,gnames] = multcompare(STATS1, "display", "off", "dimension", 1);
                        chatPostHocPower = array2table(results,"VariableNames", ...
                            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                        chatPostHocPower.("Group A")=gnames(chatPostHocPower.("Group A"));
                        chatPostHocPower.("Group B")=gnames(chatPostHocPower.("Group B"));
                        name
                        xcoordA = [18.5 18.5 10.5];
                        xcoordB = [10.5 2.5 2.5];
                    for m = 1:size(chatPostHocPower,1)
                        p = chatPostHocPower.("P-value")(m);
                        if p < .05 
                                        subplot(1,2,1)
                                        hold on
                                        areg = plot([xcoordA(m), xcoordB(m)], [.8+n, .8+n], '-b')
                                        plot(xcoordA(m), .8+n, '+b')
                                        plot(xcoordB(m), .8+n, '+b')
    
                                        n = n+ .01;
                        end 
                    end 
              
            elseif P1(2) < .05 %main effect region
                ugh = 'main effect region for'
                name
                 [results,~,~,gnames] = multcompare(STATS1, "display", "off", "dimension", 2);
                        chatPostHocRegion = array2table(results,"VariableNames", ...
                            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                        chatPostHocRegion.("Group A")=gnames(chatPostHocRegion.("Group A"));
                        chatPostHocRegion.("Group B")=gnames(chatPostHocRegion.("Group B"));
            elseif P1(3) < .05 %interaction
                ugh = 'interaction for'
                name
                 [results,~,~,gnames] = multcompare(STATS1, "display", "off", "dimension", [1 2]);
                        chatPostHocInteraction = array2table(results,"VariableNames", ...
                            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                        chatPostHocInteraction.("Group A")=gnames(chatPostHocInteraction.("Group A"));
                        chatPostHocInteraction.("Group B")=gnames(chatPostHocInteraction.("Group B"));
            end 
    
            if P2(1) < .05 %main effect of power
                 n = 0; % off set for p<.05 y value 
                 [results,~,~,gnames] = multcompare(STATS2, "display", "off", "dimension", 1);
                        dbhPostHocPower = array2table(results,"VariableNames", ...
                            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                        dbhPostHocPower.("Group A")=gnames(dbhPostHocPower.("Group A"));
                        dbhPostHocPower.("Group B")=gnames(dbhPostHocPower.("Group B"));
                        name
                        xcoordA = [18.5 18.5 10.5];
                        xcoordB = [10.5 2.5 2.5];
                    for m = 1:size(dbhPostHocPower,1)
                        p = dbhPostHocPower.("P-value")(m);
                        if p < .05 
                                        subplot(1,2,2)
                                        hold on
                                        areg = plot([xcoordA(m), xcoordB(m)], [.8+n , .8+n], '-b')
                                        plot(xcoordA(m), .8+n, '+b')
                                        plot(xcoordB(m), .8+n, '+b')
                                        n = n +.01; 
                        end 
                    end 
              
            elseif P2(2) < .05 %main effect region
                ugh = 'main effect region for'
                name
                 [results,~,~,gnames] = multcompare(STATS2, "display", "off", "dimension", 2);
                        dbhPostHocRegion = array2table(results,"VariableNames", ...
                            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                        dbhPostHocRegion.("Group A")=gnames(dbhPostHocRegion.("Group A"));
                        dbhPostHocRegion.("Group B")=gnames(dbhPostHocRegion.("Group B"));
                        
              
            elseif P2(3) < .05 %interaction
                ugh = 'interaction for'
                name
                 [results,~,~,gnames] = multcompare(STATS2, "display", "off", "dimension", [1 2]);
                        dbhPostHocInteraction = array2table(results,"VariableNames", ...
                            ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
                        dbhPostHocInteraction.("Group A")=gnames(dbhPostHocInteraction.("Group A"));
                        dbhPostHocInteraction.("Group B")=gnames(dbhPostHocInteraction.("Group B"));
            end 
        end 
    
            sgtitle(name)
    
                    
    
            filename = strcat('2powers_figure_version_smallAxis_', name);
            saveas(fig, filename, 'tif')
            saveas(fig, filename, 'fig')
            
%            tablename = append ('2anovaTablePowers', name, '_dbh.csv');
%             writecell(T2, tablename)
% 
%              tablename = append ('2anovaTablePowers', name, '_chat.csv');
%             writecell(T1, tablename)
% 
%             tablename = append ('2anovaTablePowers_postHoc', name, '_dbh.csv');
%             writetable(dbhPostHocPower, tablename)
% 
%              tablename = append ('2anovaTablePowers_postHoc', name, '_chat.csv');
%             writetable(chatPostHocPower, tablename)
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
        fig = figure
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
powers = fields(dat.chat);
name = fields(dat.chat.zero);
regions = fields(dat.chat.zero.nTrial);
nTrial = [];
nleft = [];
nright = [];


for r = 1:size(strains, 1)
        strain = strains{r};
        for k = 1:size(powers,1)
            power = powers{k};
                for i = 1:size(regions,1)
                    for p = 1:size(dat.(strain).(power).nTrial.CTRL,2)
                        region = regions{i};
                        nTrial.(power).(strain)(i,p) = sum(dat.(strain).(power).nTrial.(region){p});
                        nleft.(power).(strain)(i,p) = sum(dat.(strain).(power).left.(region){p});
                        nright.(power).(strain)(i,p) = sum(dat.(strain).(power).right.(region){p});
                        nmiss.(power).(strain)(i,p) = sum(dat.(strain).(power).miss.(region){p});
                        nReward.(power).(strain)(i,p) = sum(dat.(strain).(power).nwin.(region){p});
                        %pLeft.(power).(strain)(i,p) = nleft.(power).(region)(r,p) / nTrial.(power).(region)(r,p);
                        %pRight.(power).(strain)(i,p) = nright.(power).(region)(r,p) / nTrial.(power).(region)(r,p);
                        %rewardRate.(power).(strain)(i,p)) = nReward.(power).(region)(r,p) / nTrial.(power).(region)(r,p);
                    end 
                end 

                rewardRate.(power).(strain) = sum(nReward.(power).(strain),1) ./ sum(nTrial.(power).(strain),1);
        end 
end

%%
%are 10% of trials stimulated
    %make plot using precentStim to show 10% of trials are stimulated
    powers = fields(nTrial);
    regions = fields(nTrial.zero);
    idx = 0;
    anova_entropy.chat = [];
    anova_entropy.dbh = [];
    figure
        for i = 1:size(powers,1)
            power = powers{i};
             for k = 1:size(regions,1)
                 region = regions{k};
                 idx = idx + 1;
                        dx = idx*ones(1,size(nTrial.(power).(region)(2,:),2)-1) + 0.3*rand(1, size(nTrial.(power).(region)(2,:),2)-1);
                        cx = (idx+15)*ones(1,size(nTrial.(power).(region)(1,:),2)) + 0.3*rand(1, size(nTrial.(power).(region)(1,:),2));
            
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
            
%                          subplot(1,3,1)
%                          hold on
%                          chatx = nTrial.(power).LM2(1,:) + nTrial.(power).RM2(1,:) + nTrial.(power).LV1(1,:) + nTrial.(power).RV1(1,:) + nTrial.(power).CTRL(1,:);
%                          dbhx = nTrial.(power).LM2(2,1:6) + nTrial.(power).RM2(2,1:6) + nTrial.(power).LV1(2,1:6) + nTrial.(power).RV1(2,1:6) + nTrial.(power).CTRL(2,1:6);
%                         scatter(cx, chatx,100,cc)
%                         scatter(dx, dbhx,100,dc)
%                         scatter(idx, mean(chatx),200, cc, 'filled')
%                         scatter(idx+15, mean(dbhx),200,dc,'filled')
%             
%                    
%                         xticks([1:30])
%                         xticklabels(["0mW"; "1.5mW"; "3mW"; "0mW"; "1.5mW"; "3mW";"0mW"; "1.5mW"; "3mW"])
%                         xlim([0 8])
% %                         ylim([0 1400])
% %                         yticks([0 400 800 1200])
%                         title('Number of Trials')
%             
%                         subplot(1,3,2)
%                         hold on
%                         scatter(cx, rewardRate.(power).(region)(1,:),100,cc)
%                         scatter(idx, mean(rewardRate.(power).(region)(1,:)),200, cc, 'filled')
%                         scatter(dx, rewardRate.(power).(region)(2,1:6),100,dc)
%                         scatter(idx+15, mean(rewardRate.(power).(region)(2,1:6)),200,dc,'filled')
%             
%                         xticks([1:30])
%                         xticklabels(["M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";"M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";"M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";...
%                             "M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";"M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";"M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL"])
%                         xlim([0 31])
%                         ylim([0 .6])
%                         yticks([0 .5 ])
%                         yticklabels(["0%", "50%"])
%                         title('Reward Rate')
% %             
% %             
%                         subplot(1,3,3)
%                         hold on
%                     
%                         scatter(cx, nmiss.(power).(region)(1,:),100,cc)
%                         scatter(dx, nmiss.(power).(region)(2,1:6),100,dc)
%                         scatter(idx, mean(nmiss.(power).(region)(1,:)),200,cc,'filled')
%                         scatter(idx+15, mean(nmiss.(power).(region)(2,1:6)),200,dc,'filled')
%                        xticks([1:30])
%                         xticklabels(["M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";"M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";"M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";...
%                             "M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";"M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL";"M2-L"; "M2-R"; "V1-L"; "V1-R";"CTRL"])
%                         xlim([0 31])
%                         ylim([0 350])
%                         yticks([0 100 200 300])
%                         title('Number of Misses')
                
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
        end 

                figure
                setup_figprop;
        for i = 1:size(powers,1)
            subplot(1,3,1)
            power = powers{i}; 
            hold on 
            dbhIdx = ismember(string(stats.(power).Strain), "Dbh");
            chatIdx = ismember(string(stats.(power).Strain), "ChAt");

            dx = i*ones(1,size(stats.(power).Entropy(dbhIdx),1)) + 0.3*rand(1, size(stats.(power).Entropy(dbhIdx),1));
            cx = (i+4)*ones(1,size(stats.(power).Entropy(chatIdx),1)) + 0.3*rand(1, size(stats.(power).Entropy(chatIdx),1));
            
            anova_entropy.chat = [anova_entropy.chat; stats.(power).Entropy(chatIdx)];
            anova_entropy.dbh = [anova_entropy.dbh; stats.(power).Entropy(dbhIdx)];
            %powerList = 
            scatter(cx, stats.(power).Entropy(chatIdx),100,cc)
            scatter(dx, stats.(power).Entropy(dbhIdx),100,dc)
            scatter(i+4, mean(stats.(power).Entropy(chatIdx), 'omitnan'),200, cc, 'filled')
            scatter(i, mean(stats.(power).Entropy(dbhIdx), 'omitnan'),200,dc,'filled')

            errorbar(i+4, mean(stats.(power).Entropy(chatIdx), 'omitnan'), std(stats.(power).Entropy(chatIdx),'omitnan') / sqrt(size(stats.(power).Entropy(chatIdx),2)),...
                             'color', cc);
            errorbar(i, mean(stats.(power).Entropy(dbhIdx), 'omitnan'), std(stats.(power).Entropy(dbhIdx),'omitnan') / sqrt(size(stats.(power).Entropy(dbhIdx),2)),...
                             'color', dc);


             ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTick = [1 2 3 5 6 7];
                       ax.XTickLabel = (["0mW"; "1.5mW"; "3mW"; "0mW"; "1.5mW"; "3mW"]);
                       ax.XLim = [0 8];
                        ax.YLim = [0 3.5];
                        ax.YTick = [0:3];
                        ax.XTickLabelRotation = 90;
            title('Entropy')

            subplot(1,3,2)
            hold on
            scatter(cx, stats.(power).nTrials(chatIdx),100,cc)
            scatter(dx, stats.(power).nTrials(dbhIdx),100,dc)
            scatter(i+4, mean(stats.(power).nTrials(chatIdx), 'omitnan'),200, cc, 'filled')
            scatter(i, mean(stats.(power).nTrials(dbhIdx), 'omitnan'),200,dc,'filled')

           errorbar(i+4, mean(stats.(power).nTrials(chatIdx), 'omitnan'), std(stats.(power).nTrials(chatIdx),'omitnan') / sqrt(size(stats.(power).nTrials(chatIdx),2)),...
                             'color', cc);
            errorbar(i, mean(stats.(power).nTrials(dbhIdx), 'omitnan'), std(stats.(power).nTrials(dbhIdx),'omitnan') / sqrt(size(stats.(power).nTrials(dbhIdx),2)),...
                             'color', dc);
            
            ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTick = [1 2 3 5 6 7];
                       ax.XTickLabel = (["0mW"; "1.5mW"; "3mW"; "0mW"; "1.5mW"; "3mW"]);
                       ax.XLim = [0 8];
                        ax.YLim = [0 750];
                        ax.YTick = [100 200 300 400 500 600 700];
                        ax.XTickLabelRotation = 90;
            title('Number of Trials')

            subplot(1,3,3)
            hold on
            scatter(cx(1:7), rewardRate.(power).chat,100,cc)
            scatter(i+4, mean(rewardRate.(power).chat),200, cc, 'filled')
            scatter(dx(1:6),rewardRate.(power).dbh,100,dc)
            scatter(i, mean(rewardRate.(power).dbh),200,dc,'filled')

            errorbar(i+4, mean(rewardRate.(power).chat, 'omitnan'), std(rewardRate.(power).chat,'omitnan') / sqrt(size(rewardRate.(power).chat,2)),...
                             'color', cc);
            errorbar(i, mean(rewardRate.(power).dbh, 'omitnan'), std(rewardRate.(power).dbh,'omitnan') / sqrt(size(rewardRate.(power).dbh,2)),...
                             'color', dc);
            
                        ax = gca; 
                        ax.TickLength = [0 0];
                        ax.XTick = [1 2 3 5 6 7];
                       ax.XTickLabel = (["0mW"; "1.5mW"; "3mW"; "0mW"; "1.5mW"; "3mW"]);
                       ax.XLim = [0 8];
                        ax.YLim = [0 .6];
                        ax.YTick = [0 .25 .5];
                        ax.YTickLabel = (["0%"; "25%"; "50%"]);
                        ax.XTickLabelRotation = 90;
                        title('Reward Rate')



        end 


                %%
clearvars;
close all;
setup_figprop;
powerList = ['0.0mW'; '1.5mW'; '3.0mW'];
fullIndex = [];

for k = 1:size(powerList,1)
        currPow = powerList(k,:);
        if k == 1
            power = 'zero';
        elseif k == 2
            power = 'onefive';
        elseif k == 3
            power = 'three';
        end
        root_path = 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim';


        logfilepath = fullfile(root_path,currPow,'data');
        analysispath = fullfile(root_path, currPow, 'analysis');
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
            dataIndex.power{jj} = power;
            dataIndex.rewardRate(jj) = (sum(trialData.outcome == 100) + sum(trialData.outcome == 111)) / (dataIndex.nTrials(jj) - sum(trialData.outcome == 77));
       end
       fullIndex.(power) = dataIndex;
        clear dataIndex
end 
       
%%
powerList = ["zero"; "onefive"; "three"];
            trials.dbh = [];
            trials.chat = [];
            rewardRate.dbh = [];
            rewardRate.chat = [];
            entropy.dbh = [];
            entropy.chat = [];
            grouping.dbh = [];
            grouping.chat = [];

            setup_figprop;


    figure('position',[0 0 1500 800])
        for i = 1:size(powerList,1)
            region = powerList(i,:);
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
                        ax.LineWidth=5;
            xticks([1 2 3 6 7 8])
            xticklabels(["0mW"; "1.5mW"; "3mW";"0mW"; "1.5mW"; "3mW"])
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
                        ax.LineWidth=5;
            xticks([1 2 3 6 7 8])
            xticklabels(["0mW"; "1.5mW"; "3mW";"0mW"; "1.5mW"; "3mW"])
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
                        ax.LineWidth=5;
            xticks([1 2 3 6 7 8])
            xticklabels(["0mW"; "1.5mW"; "3mW";"0mW"; "1.5mW"; "3mW"])
            xlim([0 10])
            ylim([0 3.5])
            yticks([0 1 2 3])
            title('Entropy')


        end 


[p,tbl,stats] = anova1(trials.dbh, grouping.dbh);
            writecell(tbl, 'nTrials_dbh_powers.csv')
[p,tbl,stats] = anova1(trials.chat, grouping.chat);
    writecell(tbl, 'nTrials_chat_powers.csv')

[p,tbl,stats] = anova1(rewardRate.dbh, grouping.dbh);
    writecell(tbl, 'rewardRate_dbh_powers.csv')
[p,tbl,stats] = anova1(rewardRate.chat, grouping.chat);
    writecell(tbl, 'rewardRate_chat_powers.csv')
[p,tbl,stats] = anova1(entropy.dbh, grouping.dbh);
    writecell(tbl, 'entropy_dbh_powers.csv')
[p,tbl,stats] = anova1(entropy.chat, grouping.chat);
    writecell(tbl, 'entropy_chat_powers.csv')
