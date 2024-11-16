% analysis the GRAB_NE signal

clearvars;
close all;
setup_figprop;


root_path_NE = 'C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\22522\1.5mW';
%root_path_ACh = 'V:\HongliWang\GRAB\ACh_analysis_784';
%root_path = 'K:\Ach_784';


%% matching pennies behavior

disp('-----------------------------------------------------------');
disp('--- HERE WE GO ');
disp('--------------------------------------------------open---------');

% Look for data files and create a database index

%% NE
logfilepath_NE = fullfile(root_path_NE,'data');
analysispath_NE = fullfile(root_path_NE,'analysis');
dataIndex = makeDataIndex(logfilepath_NE, analysispath_NE);

% Parse and analyze each logfile, save as .mat files, as needed
dataIndex = MP_GRAB_createBehMatFiles(dataIndex,'NE');

% sort Index according to experiment data

dataIndex = sortdataIndex(dataIndex);

% Determine if each session fulfill performance criteria (Matching Pennies)
MP_determineBehCriteria(dataIndex);

% % ACh
% logfilepath_ACh = fullfile(root_path_ACh,'data');
% analysispath_ACh = fullfile(root_path_ACh,'analysis');
% dataIndex_ACh = makeDataIndex(logfilepath_ACh, analysispath_ACh);
% 
% % Parse and analyze each logfile, save as .mat files, as needed
% dataIndex_ACh = MP_GRAB_createBehMatFiles(dataIndex_ACh,'ACh');
% 
% % sort Index according to experiment data
% 
% dataIndex_ACh = sortdataIndex(dataIndex_ACh);

% % Determine if each session fulfill performance criteria (Matching Pennies)
% MP_determineBehCriteria(dataIndex_ACh);
%% WSLS
save_path = analysispath_NE; 
for j = 1: size(dataIndex,1)
    load(fullfile(dataIndex.BehPath{j},[dataIndex.LogFileName{j}(1:end-4),'_beh.mat']));
%    %% analysis of behavioral performance
%         % plot choice behavior - whole sessions
%         n_plot = 100*ceil(numel(stats.c)/100); %plot up to the nearest 100 trials
%         % if n_plot > 1000   %for computer simulations, #trials can be too high to plot effectively
%         %     n_plot = 1000;
%         % end
%         name = cellstr(dataIndex.Animal{j}); 
%         plot_session_task(stats,n_plot,name)
%         jj = num2str(j); 
%         filename = append('session_', jj); 
%         print(gcf,'-dpng',fullfile(save_path,filename));
%         saveas(gcf, fullfile(save_path,filename), 'fig');
    % Add stimulation info
    % stimulated side
    stats.stRegion = trialData.stimulationRegion-1000;
    
     %stimulation: yes=1; no=0;
    stats.st=nan(numel(trialData.stimulationRegion),1);
    stats.st = stats.stRegion;
    stats.st (stats.stRegion==1) =1;
    stats.st (stats.stRegion==2) =1;
    stats.st (stats.stRegion==3) =1;
    stats.st (stats.stRegion==4) =1;
   
  
%           %% plot lick timing around cue by stimulation type
%         plot_session_beh_vert(trialData, trials, stats, 'ALM stimulation', [-2 6]);
%         date = num2str(dataIndex.DateNumber(j));
%         fileName = append('sessionVert-', dataIndex.Animal(j),'-',date); 
%         fileName = char(fileName); 
%         %saveas(gcf, fullfile(save_path, fileName), 'fig');
%         print(gcf,'-dpng',fullfile(save_path, fileName)); 
    % now divide the blocks for stimulated vs non-stimulated
    %nTrials = numel(stats.c);
    %statsOriginal = stats; % stats will be used below function
    %trialsOriginal = trials;
    %trialDataOriginal = trialData;
    
    
    for kk = 1:5 % region 1 (ALM L), 2 (ALM R), 3 (V1 L), 4 (V1 R), 5 (NONE) 
       
            if kk == 5 %no stim
                regionInd = 0;
            else 
                regionInd = kk;
            end
        
         dat(kk).beh_output_WSLS{j}=WSLS_byRegion(trialData, stats, regionInd);
         dat(kk).beh_output_stim{j}=stats;
    end
end

%% calc WSLS w/o region
 save_path = analysispath_NE; 
for j = 1: size(dataIndex,1)
    load(fullfile(dataIndex.BehPath{j},[dataIndex.LogFileName{j}(1:end-4),'_beh.mat']));
    beh_output_WSLS{j}=WSLS(trialData);
    %cat.beh_output_stim{j}=stats;
end 

for tt= 1:size(dataIndex,1)
    wStay = beh_output_WSLS{1, tt}.RWinStay  + beh_output_WSLS{1, tt}.LWinStay;  
    wSwitch = beh_output_WSLS{1, tt}.RWinSwitch  + beh_output_WSLS{1, tt}.LWinSwitch;  
    lStay = beh_output_WSLS{1, tt}.RLoseStay  + beh_output_WSLS{1, tt}.LLoseStay;
    lSwitch = beh_output_WSLS{1, 1}.RLoseSwitch  + beh_output_WSLS{1, tt}.LLoseSwitch;
    
    
     
    y = [wStay;wSwitch; lStay; lSwitch];
    x = categorical({'wStay', 'wSwitch', 'lStay', 'lSwitch'});
    
    figure
    h = bar(x,y);
    title('Pre stim')
end

% plot WSLS for all sessions separated by region and side (left vs right)
for yy = 1:size(dataIndex,1)
    wStayR(yy) = beh_output_WSLS{1, yy}.RWinStay;
    wStayL(yy)=  beh_output_WSLS{1, yy}.LWinStay;  
    
    wSwitchR(yy) = beh_output_WSLS{1, yy}.RWinSwitch;
    wSwitchL(yy) = beh_output_WSLS{1, yy}.LWinSwitch;  
   
    lStayR(yy) = beh_output_WSLS{1, yy}.RLoseStay;
    lStayL(yy) = beh_output_WSLS{1, yy}.LLoseStay;
   
    lSwitchR(yy) = beh_output_WSLS{1, yy}.RLoseSwitch;
    lSwitchL(yy) = beh_output_WSLS{1, yy}.LLoseSwitch;
end

    wStayR = sum(wStayR);
    wStayL = sum (wStayL);


   wSwitchR = sum(wSwitchR); 
    wSwitchL = sum(wSwitchL); 
    
   lStayR = sum(lStayR);
    lStayL = sum(lStayL);
  
   lSwitchR = sum(lSwitchR);
    lSwitchL = sum(lSwitchL);
 

y = [wStayL; wStayR;
    wSwitchL; wSwitchR;
    lStayL; lStayR;
    lSwitchL; lSwitchR];

x = categorical({'wStayL', 'wStayR', 'wSwitchL', 'wSwitchR', 'lStayL', 'lStayR', 'lSwitchL', 'lSwitchR'});

figure
h = bar(x,y);
title('PreStim')
ylabel('Raw count of occurances')

% plot WSLS data converted to proportion (# of occurances / number of trials with that region stimulated)



for yy = 1:size(dataIndex,1)  
    trials(yy) = beh_output_WSLS{1,yy}.nTrial;
    wStayR(yy) = beh_output_WSLS{1, yy}.RWinStay;
    wStayL(yy)=  beh_output_WSLS{1, yy}.LWinStay;  
 
    
    wSwitchR(yy) = beh_output_WSLS{1, yy}.RWinSwitch;
    wSwitchL(yy) = beh_output_WSLS{1, yy}.LWinSwitch;  
 
    lStayR(yy) = beh_output_WSLS{1, yy}.RLoseStay;
    lStayL(yy) = beh_output_WSLS{1, yy}.LLoseStay;
   
    
    lSwitchR(yy) = beh_output_WSLS{1, yy}.RLoseSwitch;
    lSwitchL(yy) = beh_output_WSLS{1, yy}.LLoseSwitch;
end

    prop_wStayR = sum(wStayR)/sum(trials);
    prop_wStayL = sum (wStayL)/sum(trials);
 

    prop_wSwitchR = sum(wSwitchR/sum(trials)); 
    prop_wSwitchL = sum(wSwitchL)/sum(trials); 
 

    prop_lStayR = sum(lStayR)/sum(trials);
    prop_lStayL = sum(lStayL)/sum(trials);
 
    
    prop_lSwitchR = sum(lSwitchR)/sum(trials);
    prop_lSwitchL = sum(lSwitchL)/sum(trials);
  

y = [prop_wStayL; prop_wStayR;
    prop_wSwitchL; prop_wSwitchR;
    prop_lStayL; prop_lStayR;
    prop_lSwitchL; prop_lSwitchR];

x = categorical({'wStayL', 'wStayR', 'wSwitchL', 'wSwitchR', 'lStayL', 'lStayR', 'lSwitchL', 'lSwitchR'});

figure
h = bar(x,y);
title('Pre stim')
ylabel('# occurances / # trials')

% plot WSLS data converted to proportion (# of occurances / number of WSLS trial type (ex leftWin or rightLose)
numSessions = size(dataIndex,1);
leftLose = zeros(numSessions, 5);
leftWin = zeros(numSessions, 5);
rightWin = zeros(numSessions, 5);
rightLose = zeros(numSessions, 5);
for dd = 1:numSessions
    leftLose(dd) = beh_output_WSLS{1, dd}.leftLose;
    leftWin(dd) =  beh_output_WSLS{1, dd}.leftWin;
    rightWin(dd) = beh_output_WSLS{1, dd}.rightWin;
    rightLose(dd) = beh_output_WSLS{1, dd}.rightLose;
end

leftLose = sum(leftLose);
leftWin = sum(leftWin);
rightWin = sum(rightWin);
rightLose = sum(rightLose);
% LM2 = 1;
% RM2 = 2;
% LV1 = 3;
% RV1 = 4;
% NS = 5;

    WSLS_wStayR = sum(wStayR)/rightWin(1);
    WSLS_wStayL = sum (wStayL)/leftWin(1);

    WSLS_wSwitchR = sum(wSwitchR)/rightWin(1);
    WSLS_wSwitchL = sum(wSwitchL)/leftWin(1);
   
    WSLS_lStayR = sum(lStayR)/rightLose(1);
    WSLS_lStayL = sum(lStayL)/leftLose(1);
   
    
    WSLS_lSwitchR = sum(lSwitchR)/rightLose(1);
    WSLS_lSwitchL = sum(lSwitchL)/leftLose(1);

y = [ WSLS_wStayL;  WSLS_wStayR;
     WSLS_wSwitchL;  WSLS_wSwitchR;
     WSLS_lStayL;  WSLS_lStayR;
    WSLS_lSwitchL;  WSLS_lSwitchR];

x = categorical({'wStayL', 'wStayR', 'wSwitchL', 'wSwitchR', 'lStayL', 'lStayR', 'lSwitchL', 'lSwitchR'});

figure
h = bar(x,y);
title('Pre stim')
ylabel('# occurances / # similar trials (leftLose)')


% plot combined L&R effects
    ipsi_wStayL = (sum(wStayL)+sum(wStayR))/(leftWin(1)+rightWin(1));

    ipsi_wSwitchL = (sum(wSwitchL)+sum(wSwitchR))/(leftWin(1)+rightWin(1));

    ipsi_lStayL = (sum(lStayL)+sum(lStayR))/(leftLose(1)+rightLose(1));
    
    ipsi_lSwitchL = (sum(lSwitchL)+sum(lSwitchR))/(leftLose(1)+rightLose(1));

y = [ ipsi_wStayL;
     ipsi_wSwitchL;
     ipsi_lStayL;
    ipsi_lSwitchL];

x = categorical({'wStay','wSwitch', 'lStay', 'lSwitch'});

figure
h = bar(x,y);
title('Collapsed L/R changes')
ylabel('# occurances / # similar trials (leftLose)')
%%
% NE
clear sessionData stats trialData
nFiles_NE = size(dataIndex(:,:),1);

for ii = 1:nFiles_NE
    savematpath_NE = dataIndex.BehPath{ii};
    MP_GRAB_session(dataIndex.BehPath{ii},dataIndex.LogFileName{ii},savematpath_NE);
end
save_path_NE = fullfile(root_path_NE,'summary','figs_summary');


MP_STIM_behaviorPerAnimal(dataIndex,save_path_NE);

MP_STIM_behaviorAll(dataIndex, save_path_NE);
%%
% NE
model_path = fullfile(root_path_NE,'mat_models\NS');


MP_GRAB_fittingPerAnimal(dataIndex,model_path);
%MP_STIM_fittingPerSession(dataIndex,model_path);

% save the predicted latent variable into each behavior mat files
MP_GRAB_saveLatent(dataIndex, model_path);




%% plot WSLS data
%plot each day WSLS by region bar graphs

for tt= 1:size(dataIndex,1)
LM2_wStay = dat(1).beh_output_WSLS{1, tt}.RWinStay  + dat(1).beh_output_WSLS{1, tt}.LWinStay;  
RM2_wStay = dat(2).beh_output_WSLS{1, tt}.RWinStay  + dat(2).beh_output_WSLS{1, tt}.LWinStay;
LV1_wStay = dat(3).beh_output_WSLS{1, tt}.RWinStay  + dat(3).beh_output_WSLS{1, tt}.LWinStay;
RV1_wStay = dat(4).beh_output_WSLS{1, tt}.RWinStay  + dat(4).beh_output_WSLS{1, tt}.LWinStay;

LM2_wSwitch = dat(1).beh_output_WSLS{1, tt}.RWinSwitch  + dat(1).beh_output_WSLS{1, tt}.LWinSwitch;  
RM2_wSwitch = dat(2).beh_output_WSLS{1, tt}.RWinSwitch  + dat(2).beh_output_WSLS{1, tt}.LWinSwitch;
LV1_wSwitch = dat(3).beh_output_WSLS{1, tt}.RWinSwitch  + dat(3).beh_output_WSLS{1, tt}.LWinSwitch;
RV1_wSwitch = dat(4).beh_output_WSLS{1, tt}.RWinSwitch  + dat(4).beh_output_WSLS{1, tt}.LWinSwitch;

LM2_lStay = dat(1).beh_output_WSLS{1, tt}.RLoseStay  + dat(1).beh_output_WSLS{1, tt}.LLoseStay;
RM2_lStay = dat(2).beh_output_WSLS{1, tt}.RLoseStay  + dat(2).beh_output_WSLS{1, tt}.LLoseStay;
LV1_lStay = dat(3).beh_output_WSLS{1, tt}.RLoseStay  + dat(3).beh_output_WSLS{1, tt}.LLoseStay;
RV1_lStay = dat(4).beh_output_WSLS{1, tt}.RLoseStay  + dat(4).beh_output_WSLS{1, tt}.LLoseStay;

LM2_lSwitch = dat(1).beh_output_WSLS{1, 1}.RLoseSwitch  + dat(1).beh_output_WSLS{1, tt}.LLoseSwitch;
RM2_lSwitch = dat(2).beh_output_WSLS{1, 1}.RLoseSwitch  + dat(2).beh_output_WSLS{1, tt}.LLoseSwitch;
LV1_lSwitch = dat(3).beh_output_WSLS{1, 1}.RLoseSwitch  + dat(3).beh_output_WSLS{1, tt}.LLoseSwitch;
RV1_lSwitch= dat(4).beh_output_WSLS{1, 1}.RLoseSwitch  + dat(4).beh_output_WSLS{1, tt}.LLoseSwitch;


 
y = [LM2_wStay RM2_wStay LV1_wStay RV1_wStay; LM2_wSwitch RM2_wSwitch LV1_wSwitch RV1_wSwitch; LM2_lStay RM2_lStay LV1_lStay RV1_lStay; LM2_lSwitch RM2_lSwitch LV1_lSwitch RV1_lSwitch];

x = categorical({'wStay', 'wSwitch', 'lStay', 'lSwitch'});

figure
h = bar(x,y);
legend('LM2', 'RM2', 'LV1', 'RV1') 
title(dataIndex.LogFileName(tt))
end


% plot WSLS for all sessions separated by region and side (left vs right)
for yy = 1:size(dataIndex,1)
    LM2_wStayR(yy) = dat(1).beh_output_WSLS{1, yy}.RWinStay;
    LM2_wStayL(yy)=  dat(1).beh_output_WSLS{1, yy}.LWinStay;  
    RM2_wStayR(yy) = dat(2).beh_output_WSLS{1, yy}.RWinStay; 
    RM2_wStayL(yy) = dat(2).beh_output_WSLS{1, yy}.LWinStay;
    LV1_wStayR(yy) = dat(3).beh_output_WSLS{1, yy}.RWinStay; 
    LV1_wStayL(yy) = dat(3).beh_output_WSLS{1, yy}.LWinStay;
    RV1_wStayR(yy) = dat(4).beh_output_WSLS{1, yy}.RWinStay;
    RV1_wStayL(yy) = dat(4).beh_output_WSLS{1, yy}.LWinStay;
    
    LM2_wSwitchR(yy) = dat(1).beh_output_WSLS{1, yy}.RWinSwitch;
    LM2_wSwitchL(yy) = dat(1).beh_output_WSLS{1, yy}.LWinSwitch;  
    RM2_wSwitchR(yy) = dat(2).beh_output_WSLS{1, yy}.RWinSwitch;
    RM2_wSwitchL(yy) = dat(2).beh_output_WSLS{1, yy}.LWinSwitch;
    LV1_wSwitchR(yy) = dat(3).beh_output_WSLS{1, yy}.RWinSwitch;
    LV1_wSwitchL(yy) = dat(3).beh_output_WSLS{1, yy}.LWinSwitch;
    RV1_wSwitchR(yy) = dat(4).beh_output_WSLS{1, yy}.RWinSwitch;
    RV1_wSwitchL(yy) = dat(4).beh_output_WSLS{1, yy}.LWinSwitch;
    
    LM2_lStayR(yy) = dat(1).beh_output_WSLS{1, yy}.RLoseStay;
    LM2_lStayL(yy) = dat(1).beh_output_WSLS{1, yy}.LLoseStay;
    RM2_lStayR(yy) = dat(2).beh_output_WSLS{1, yy}.RLoseStay;
    RM2_lStayL(yy) = dat(2).beh_output_WSLS{1, yy}.LLoseStay;
    LV1_lStayR(yy) = dat(3).beh_output_WSLS{1, yy}.RLoseStay;
    LV1_lStayL(yy) = dat(3).beh_output_WSLS{1, yy}.LLoseStay;
    RV1_lStayR(yy) = dat(4).beh_output_WSLS{1, yy}.RLoseStay;
    RV1_lStayL(yy) = dat(4).beh_output_WSLS{1, yy}.LLoseStay;
    
    LM2_lSwitchR(yy) = dat(1).beh_output_WSLS{1, yy}.RLoseSwitch;
    LM2_lSwitchL(yy) = dat(1).beh_output_WSLS{1, yy}.LLoseSwitch;
    RM2_lSwitchR(yy) = dat(2).beh_output_WSLS{1, yy}.RLoseSwitch;
    RM2_lSwitchL(yy) = dat(2).beh_output_WSLS{1, yy}.LLoseSwitch;
    LV1_lSwitchR(yy) = dat(3).beh_output_WSLS{1, yy}.RLoseSwitch;
    LV1_lSwitchL(yy) = dat(3).beh_output_WSLS{1, yy}.LLoseSwitch;
    RV1_lSwitchR(yy)= dat(4).beh_output_WSLS{1, yy}.RLoseSwitch;
    RV1_lSwitchL(yy) = dat(4).beh_output_WSLS{1, yy}.LLoseSwitch;
end

    LM2_wStayR = sum(LM2_wStayR);
    LM2_wStayL = sum (LM2_wStayL);
    RM2_wStayR = sum(RM2_wStayR);
    RM2_wStayL = sum(RM2_wStayL);
    LV1_wStayR = sum(LV1_wStayR); 
    LV1_wStayL = sum(LV1_wStayL);
    RV1_wStayR = sum(RV1_wStayR);
    RV1_wStayL = sum(RV1_wStayL); 


    LM2_wSwitchR = sum(LM2_wSwitchR); 
    LM2_wSwitchL = sum(LM2_wSwitchL); 
    RM2_wSwitchR = sum(RM2_wSwitchR); 
    RM2_wSwitchL = sum(RM2_wSwitchL); 
    LV1_wSwitchR = sum(LV1_wSwitchR);
    LV1_wSwitchL = sum(LV1_wSwitchL); 
    RV1_wSwitchR = sum(RV1_wSwitchR);
    RV1_wSwitchL = sum(RV1_wSwitchL);
    
    LM2_lStayR = sum(LM2_lStayR);
    LM2_lStayL = sum(LM2_lStayL);
    RM2_lStayR = sum(RM2_lStayR);
    RM2_lStayL = sum(RM2_lStayL);
    LV1_lStayR = sum(LV1_lStayR);
    LV1_lStayL = sum(LV1_lStayL);
    RV1_lStayR = sum(RV1_lStayR);
    RV1_lStayL = sum(RV1_lStayL);
    
    LM2_lSwitchR = sum(LM2_lSwitchR);
    LM2_lSwitchL = sum(LM2_lSwitchL);
    RM2_lSwitchR = sum(RM2_lSwitchR);
    RM2_lSwitchL = sum(RM2_lSwitchL);
    LV1_lSwitchR = sum(LV1_lSwitchR);
    LV1_lSwitchL = sum(LV1_lSwitchL);
    RV1_lSwitchR = sum(RV1_lSwitchR);
    RV1_lSwitchL = sum(RV1_lSwitchL);

y = [LM2_wStayL RM2_wStayL LV1_wStayL RV1_wStayL; LM2_wStayR RM2_wStayR LV1_wStayR RV1_wStayR;
    LM2_wSwitchL RM2_wSwitchL LV1_wSwitchL RV1_wSwitchL; LM2_wSwitchR RM2_wSwitchR LV1_wSwitchR RV1_wSwitchR;
    LM2_lStayL RM2_lStayL LV1_lStayL RV1_lStayL; LM2_lStayR RM2_lStayR LV1_lStayR RV1_lStayR;
    LM2_lSwitchL RM2_lSwitchL LV1_lSwitchL RV1_lSwitchL; LM2_lSwitchR RM2_lSwitchR LV1_lSwitchR RV1_lSwitchR];

x = categorical({'wStayL', 'wStayR', 'wSwitchL', 'wSwitchR', 'lStayL', 'lStayR', 'lSwitchL', 'lSwitchR'});

figure
h = bar(x,y);
legend('LM2', 'RM2', 'LV1', 'RV1') 
title(dataIndex.LogFileName(tt))
ylabel('Raw count of occurances')

% plot WSLS data converted to proportion (# of occurances / number of trials with that region stimulated)



for yy = 1:size(dataIndex,1)
    LM2_trials(yy) = sum(dat(1).beh_output_stim{1,yy}.stRegion==1);
    RM2_trials(yy) = sum(dat(2).beh_output_stim{1,yy}.stRegion==2);
    LV1_trials(yy) = sum(dat(3).beh_output_stim{1,yy}.stRegion==3);
    RV1_trials(yy) = sum(dat(4).beh_output_stim{1,yy}.stRegion==4);
    Ctrl_trials(yy) = sum(dat(5).beh_output_stim{1,yy}.stRegion==0);

    LM2_wStayR(yy) = dat(1).beh_output_WSLS{1, yy}.RWinStay;
    LM2_wStayL(yy)=  dat(1).beh_output_WSLS{1, yy}.LWinStay;  
    RM2_wStayR(yy) = dat(2).beh_output_WSLS{1, yy}.RWinStay; 
    RM2_wStayL(yy) = dat(2).beh_output_WSLS{1, yy}.LWinStay;
    LV1_wStayR(yy) = dat(3).beh_output_WSLS{1, yy}.RWinStay; 
    LV1_wStayL(yy) = dat(3).beh_output_WSLS{1, yy}.LWinStay;
    RV1_wStayR(yy) = dat(4).beh_output_WSLS{1, yy}.RWinStay;
    RV1_wStayL(yy) = dat(4).beh_output_WSLS{1, yy}.LWinStay;
    Ctrl_wStayR(yy) = dat(5).beh_output_WSLS{1, yy}.RWinStay;
    Ctrl_wStayL(yy) = dat(5).beh_output_WSLS{1, yy}.LWinStay;
    
    LM2_wSwitchR(yy) = dat(1).beh_output_WSLS{1, yy}.RWinSwitch;
    LM2_wSwitchL(yy) = dat(1).beh_output_WSLS{1, yy}.LWinSwitch;  
    RM2_wSwitchR(yy) = dat(2).beh_output_WSLS{1, yy}.RWinSwitch;
    RM2_wSwitchL(yy) = dat(2).beh_output_WSLS{1, yy}.LWinSwitch;
    LV1_wSwitchR(yy) = dat(3).beh_output_WSLS{1, yy}.RWinSwitch;
    LV1_wSwitchL(yy) = dat(3).beh_output_WSLS{1, yy}.LWinSwitch;
    RV1_wSwitchR(yy) = dat(4).beh_output_WSLS{1, yy}.RWinSwitch;
    RV1_wSwitchL(yy) = dat(4).beh_output_WSLS{1, yy}.LWinSwitch;
    Ctrl_wSwitchR(yy) = dat(5).beh_output_WSLS{1, yy}.RWinSwitch;
    Ctrl_wSwitchL(yy) = dat(5).beh_output_WSLS{1, yy}.LWinSwitch;
    
    LM2_lStayR(yy) = dat(1).beh_output_WSLS{1, yy}.RLoseStay;
    LM2_lStayL(yy) = dat(1).beh_output_WSLS{1, yy}.LLoseStay;
    RM2_lStayR(yy) = dat(2).beh_output_WSLS{1, yy}.RLoseStay;
    RM2_lStayL(yy) = dat(2).beh_output_WSLS{1, yy}.LLoseStay;
    LV1_lStayR(yy) = dat(3).beh_output_WSLS{1, yy}.RLoseStay;
    LV1_lStayL(yy) = dat(3).beh_output_WSLS{1, yy}.LLoseStay;
    RV1_lStayR(yy) = dat(4).beh_output_WSLS{1, yy}.RLoseStay;
    RV1_lStayL(yy) = dat(4).beh_output_WSLS{1, yy}.LLoseStay;
    Ctrl_lStayR(yy) = dat(5).beh_output_WSLS{1, yy}.RLoseStay;
    Ctrl_lStayL(yy) = dat(5).beh_output_WSLS{1, yy}.LLoseStay;
    
    LM2_lSwitchR(yy) = dat(1).beh_output_WSLS{1, yy}.RLoseSwitch;
    LM2_lSwitchL(yy) = dat(1).beh_output_WSLS{1, yy}.LLoseSwitch;
    RM2_lSwitchR(yy) = dat(2).beh_output_WSLS{1, yy}.RLoseSwitch;
    RM2_lSwitchL(yy) = dat(2).beh_output_WSLS{1, yy}.LLoseSwitch;
    LV1_lSwitchR(yy) = dat(3).beh_output_WSLS{1, yy}.RLoseSwitch;
    LV1_lSwitchL(yy) = dat(3).beh_output_WSLS{1, yy}.LLoseSwitch;
    RV1_lSwitchR(yy)= dat(4).beh_output_WSLS{1, yy}.RLoseSwitch;
    RV1_lSwitchL(yy) = dat(4).beh_output_WSLS{1, yy}.LLoseSwitch;
    Ctrl_lSwitchR(yy)= dat(5).beh_output_WSLS{1, yy}.RLoseSwitch;
    Ctrl_lSwitchL(yy) = dat(5).beh_output_WSLS{1, yy}.LLoseSwitch;
end

    prop_LM2_wStayR = sum(LM2_wStayR)/sum(LM2_trials);
    prop_LM2_wStayL = sum (LM2_wStayL)/sum(LM2_trials);
    prop_RM2_wStayR = sum(RM2_wStayR)/sum(RM2_trials);
    prop_RM2_wStayL = sum(RM2_wStayL)/sum(RM2_trials);
    prop_LV1_wStayR = sum(LV1_wStayR)/sum(LV1_trials); 
    prop_LV1_wStayL = sum(LV1_wStayL)/sum(LV1_trials);
    prop_RV1_wStayR = sum(RV1_wStayR)/sum(RV1_trials);
    prop_RV1_wStayL = sum(RV1_wStayL)/sum(RV1_trials); 
    prop_Ctrl_wStayR = sum(Ctrl_wStayR)/sum(Ctrl_trials);
    prop_Ctrl_wStayL = sum(Ctrl_wStayL)/sum(Ctrl_trials);

    prop_LM2_wSwitchR = sum(LM2_wSwitchR/sum(LM2_trials)); 
    prop_LM2_wSwitchL = sum(LM2_wSwitchL)/sum(LM2_trials); 
    prop_RM2_wSwitchR = sum(RM2_wSwitchR)/sum(RM2_trials); 
    prop_RM2_wSwitchL = sum(RM2_wSwitchL)/sum(RM2_trials); 
    prop_LV1_wSwitchR = sum(LV1_wSwitchR)/sum(LV1_trials);
    prop_LV1_wSwitchL = sum(LV1_wSwitchL)/sum(LV1_trials); 
    prop_RV1_wSwitchR = sum(RV1_wSwitchR)/sum(RV1_trials);
    prop_RV1_wSwitchL = sum(RV1_wSwitchL)/sum(RV1_trials);
    prop_Ctrl_wSwitchR = sum(Ctrl_wSwitchR)/sum(Ctrl_trials);
    prop_Ctrl_wSwitchL = sum(Ctrl_wSwitchL)/sum(Ctrl_trials);

    prop_LM2_lStayR = sum(LM2_lStayR)/sum(LM2_trials);
    prop_LM2_lStayL = sum(LM2_lStayL)/sum(LM2_trials);
    prop_RM2_lStayR = sum(RM2_lStayR)/sum(RM2_trials);
    prop_RM2_lStayL = sum(RM2_lStayL)/sum(RM2_trials);
    prop_LV1_lStayR = sum(LV1_lStayR)/sum(LV1_trials);
    prop_LV1_lStayL = sum(LV1_lStayL)/sum(LV1_trials);
    prop_RV1_lStayR = sum(RV1_lStayR)/sum(RV1_trials);
    prop_RV1_lStayL = sum(RV1_lStayL)/sum(RV1_trials);
    prop_Ctrl_lStayR = sum(Ctrl_lStayR)/sum(Ctrl_trials);
    prop_Ctrl_lStayL = sum(Ctrl_lStayL)/sum(Ctrl_trials);
    
    prop_LM2_lSwitchR = sum(LM2_lSwitchR)/sum(LM2_trials);
    prop_LM2_lSwitchL = sum(LM2_lSwitchL)/sum(LM2_trials);
    prop_RM2_lSwitchR = sum(RM2_lSwitchR)/sum(RM2_trials);
    prop_RM2_lSwitchL = sum(RM2_lSwitchL)/sum(RM2_trials);
    prop_LV1_lSwitchR = sum(LV1_lSwitchR)/sum(LV1_trials);
    prop_LV1_lSwitchL = sum(LV1_lSwitchL)/sum(LV1_trials);
    prop_RV1_lSwitchR = sum(RV1_lSwitchR)/sum(RV1_trials);
    prop_RV1_lSwitchL = sum(RV1_lSwitchL)/sum(RV1_trials);
    prop_Ctrl_lSwitchR = sum(Ctrl_lSwitchR)/sum(Ctrl_trials);
    prop_Ctrl_lSwitchL = sum(Ctrl_lSwitchL)/sum(Ctrl_trials);

y = [prop_LM2_wStayL prop_RM2_wStayL prop_LV1_wStayL prop_RV1_wStayL prop_Ctrl_wStayL; prop_LM2_wStayR prop_RM2_wStayR prop_LV1_wStayR prop_RV1_wStayR prop_Ctrl_wStayR;
    prop_LM2_wSwitchL prop_RM2_wSwitchL prop_LV1_wSwitchL prop_RV1_wSwitchL prop_Ctrl_wSwitchL; prop_LM2_wSwitchR prop_RM2_wSwitchR prop_LV1_wSwitchR prop_RV1_wSwitchR prop_Ctrl_wSwitchR;
    prop_LM2_lStayL prop_RM2_lStayL prop_LV1_lStayL prop_RV1_lStayL prop_Ctrl_lStayL; prop_LM2_lStayR prop_RM2_lStayR prop_LV1_lStayR prop_RV1_lStayR prop_Ctrl_lStayR;
    prop_LM2_lSwitchL prop_RM2_lSwitchL prop_LV1_lSwitchL prop_RV1_lSwitchL prop_Ctrl_lSwitchL; prop_LM2_lSwitchR prop_RM2_lSwitchR prop_LV1_lSwitchR prop_RV1_lSwitchR prop_Ctrl_lSwitchR];

x = categorical({'wStayL', 'wStayR', 'wSwitchL', 'wSwitchR', 'lStayL', 'lStayR', 'lSwitchL', 'lSwitchR'});

figure
h = bar(x,y);
legend('LM2', 'RM2', 'LV1', 'RV1', 'Ctrl') 
title('22522 DbhxAi32 5 sessions combined')
ylabel('# occurances / # trials stimulated per region')

% plot WSLS data converted to proportion (# of occurances / number of WSLS trial type (ex leftWin or rightLose)
numSessions = size(dataIndex,1);
leftLose = zeros(numSessions, 5);
leftWin = zeros(numSessions, 5);
rightWin = zeros(numSessions, 5);
rightLose = zeros(numSessions, 5);
for dd = 1:numSessions
    for ss = 1:5
    leftLose(dd,ss) = dat(ss).beh_output_WSLS{1, dd}.leftLose;
    leftWin(dd,ss) =  dat(ss).beh_output_WSLS{1, dd}.leftWin;
    rightWin(dd,ss) = dat(ss).beh_output_WSLS{1, dd}.rightWin;
    rightLose(dd,ss) =dat(ss).beh_output_WSLS{1, dd}.rightLose;
    end
end

leftLose = sum(leftLose,1);
leftWin = sum(leftWin,1);
rightWin = sum(rightWin,1);
rightLose = sum(rightLose,1);
% LM2 = 1;
% RM2 = 2;
% LV1 = 3;
% RV1 = 4;
% NS = 5;

    WSLS_LM2_wStayR = sum(LM2_wStayR)/rightWin(1);
    WSLS_LM2_wStayL = sum (LM2_wStayL)/leftWin(1);
    WSLS_RM2_wStayR = sum(RM2_wStayR)/rightWin(2);
    WSLS_RM2_wStayL = sum(RM2_wStayL)/leftWin(2);
    WSLS_LV1_wStayR = sum(LV1_wStayR)/rightWin(3); 
    WSLS_LV1_wStayL = sum(LV1_wStayL)/leftWin(3);
    WSLS_RV1_wStayR = sum(RV1_wStayR)/rightWin(4);
    WSLS_RV1_wStayL = sum(RV1_wStayL)/leftWin(4); 
    WSLS_Ctrl_wStayR = sum(Ctrl_wStayR)/rightWin(5);
    WSLS_Ctrl_wStayL = sum(Ctrl_wStayL)/leftWin(5);

    WSLS_LM2_wSwitchR = sum(LM2_wSwitchR)/rightWin(1);
    WSLS_LM2_wSwitchL = sum(LM2_wSwitchL)/leftWin(1);
    WSLS_RM2_wSwitchR = sum(RM2_wSwitchR)/rightWin(2); 
    WSLS_RM2_wSwitchL = sum(RM2_wSwitchL)/leftWin(2); 
    WSLS_LV1_wSwitchR = sum(LV1_wSwitchR)/rightWin(3);
    WSLS_LV1_wSwitchL = sum(LV1_wSwitchL)/leftWin(3);
    WSLS_RV1_wSwitchR = sum(RV1_wSwitchR)/rightWin(4);
    WSLS_RV1_wSwitchL = sum(RV1_wSwitchL)/leftWin(4);
    WSLS_Ctrl_wSwitchR = sum(Ctrl_wSwitchR)/rightWin(5);
    WSLS_Ctrl_wSwitchL = sum(Ctrl_wSwitchL)/leftWin(5);

    WSLS_LM2_lStayR = sum(LM2_lStayR)/rightLose(1);
    WSLS_LM2_lStayL = sum(LM2_lStayL)/leftLose(1);
    WSLS_RM2_lStayR = sum(RM2_lStayR)/rightLose(2);
    WSLS_RM2_lStayL = sum(RM2_lStayL)/leftLose(2);
    WSLS_LV1_lStayR = sum(LV1_lStayR)/rightLose(3);
    WSLS_LV1_lStayL = sum(LV1_lStayL)/leftLose(3);
    WSLS_RV1_lStayR = sum(RV1_lStayR)/rightLose(4);
    WSLS_RV1_lStayL = sum(RV1_lStayL)/leftLose(4);
    WSLS_Ctrl_lStayR = sum(Ctrl_lStayR)/rightLose(5);
    WSLS_Ctrl_lStayL = sum(Ctrl_lStayL)/leftLose(5);
    
    WSLS_LM2_lSwitchR = sum(LM2_lSwitchR)/rightLose(1);
    WSLS_LM2_lSwitchL = sum(LM2_lSwitchL)/leftLose(1);
    WSLS_RM2_lSwitchR = sum(RM2_lSwitchR)/rightLose(2);
    WSLS_RM2_lSwitchL = sum(RM2_lSwitchL)/leftLose(2);
    WSLS_LV1_lSwitchR = sum(LV1_lSwitchR)/rightLose(3);
    WSLS_LV1_lSwitchL = sum(LV1_lSwitchL)/leftLose(3);
    WSLS_RV1_lSwitchR = sum(RV1_lSwitchR)/rightLose(4);
    WSLS_RV1_lSwitchL = sum(RV1_lSwitchL)/leftLose(4);
    WSLS_Ctrl_lSwitchR = sum(Ctrl_lSwitchR)/rightLose(5);
    WSLS_Ctrl_lSwitchL = sum(Ctrl_lSwitchL)/leftLose(5);

y = [ WSLS_LM2_wStayL  WSLS_RM2_wStayL  WSLS_LV1_wStayL  WSLS_RV1_wStayL  WSLS_Ctrl_wStayL;  WSLS_LM2_wStayR  WSLS_RM2_wStayR  WSLS_LV1_wStayR  WSLS_RV1_wStayR  WSLS_Ctrl_wStayR;
     WSLS_LM2_wSwitchL  WSLS_RM2_wSwitchL  WSLS_LV1_wSwitchL  WSLS_RV1_wSwitchL  WSLS_Ctrl_wSwitchL;  WSLS_LM2_wSwitchR  WSLS_RM2_wSwitchR  WSLS_LV1_wSwitchR  WSLS_RV1_wSwitchR  WSLS_Ctrl_wSwitchR;
     WSLS_LM2_lStayL WSLS_RM2_lStayL  WSLS_LV1_lStayL  WSLS_RV1_lStayL  WSLS_Ctrl_lStayL;  WSLS_LM2_lStayR  WSLS_RM2_lStayR  WSLS_LV1_lStayR  WSLS_RV1_lStayR  WSLS_Ctrl_lStayR;
    WSLS_LM2_lSwitchL  WSLS_RM2_lSwitchL  WSLS_LV1_lSwitchL  WSLS_RV1_lSwitchL  WSLS_Ctrl_lSwitchL;  WSLS_LM2_lSwitchR  WSLS_RM2_lSwitchR WSLS_LV1_lSwitchR  WSLS_RV1_lSwitchR  WSLS_Ctrl_lSwitchR];

x = categorical({'wStayL', 'wStayR', 'wSwitchL', 'wSwitchR', 'lStayL', 'lStayR', 'lSwitchL', 'lSwitchR'});

figure
h = bar(x,y);
legend('LM2', 'RM2', 'LV1', 'RV1', 'Ctrl') 
title('22522 DbhxAi32 1 session 0.0mW')
ylabel('# occurances / # similar trials (leftLose) per region')
% plot ipsilateral effects
    ipsi_LM2_wStayL = (sum(LM2_wStayL)+sum(RM2_wStayR))/(leftWin(1)+rightWin(2));
    ipsi_LV1_wStayL = (sum(LV1_wStayL)+sum(RV1_wStayR))/(leftWin(3)+rightWin(4));

    ipsi_LM2_wSwitchL = (sum(LM2_wSwitchL)+sum(RM2_wSwitchR))/(leftWin(1)+rightWin(2));
    ipsi_LV1_wSwitchL = (sum(LV1_wSwitchL)+sum(RV1_wSwitchR))/(leftWin(3)+rightWin(4));

    ipsi_LM2_lStayL = (sum(LM2_lStayL)+sum(RM2_lStayR))/(leftLose(1)+rightLose(2));
    ipsi_LV1_lStayL = (sum(LV1_lStayL)+sum(RV1_lStayR))/(leftLose(3)+rightLose(4));
    
    ipsi_LM2_lSwitchL = (sum(LM2_lSwitchL)+sum(RM2_lSwitchR))/(leftLose(1)+rightLose(2));
    ipsi_LV1_lSwitchL = (sum(LV1_lSwitchL)+sum(RV1_lSwitchR))/(leftLose(3)+rightLose(4));

y = [ ipsi_LM2_wSwitchL ipsi_LV1_wSwitchL WSLS_Ctrl_wSwitchL;
    ipsi_LM2_lSwitchL  ipsi_LV1_lSwitchL  WSLS_Ctrl_lSwitchL ];

x = categorical({'wSwitch','lSwitch'});

figure
h = bar(x,y);
legend('M2', 'V1', 'Ctrl') 
title('ipsilateral changes')
ylabel('# occurances / # similar trials (leftLose) per region')

% plot contralateral effects
    contra_M2_wStay = (sum(LM2_wStayR)+sum(RM2_wStayL))/(rightWin(1)+leftWin(2));
    contra_V1_wStay = (sum(LV1_wStayR)+sum(RV1_wStayL))/(rightWin(3)+leftWin(4));

    contra_M2_wSwitch = (sum(LM2_wSwitchR)+sum(RM2_wSwitchL))/(rightWin(1)+leftWin(2));
    contra_V1_wSwitch = (sum(LV1_wSwitchR)+sum(RV1_wSwitchL))/(rightWin(3)+leftWin(4));

    contra_M2_lStay = (sum(LM2_lStayR)+sum(RM2_lStayL))/(rightLose(1)+leftLose(2));
    contra_V1_lStay = (sum(LV1_lStayR)+sum(RV1_lStayL))/(rightLose(3)+leftLose(4));
    
    contra_M2_lSwitch = (sum(LM2_lSwitchR)+sum(RM2_lSwitchL))/(rightLose(1)+leftLose(2));
    contra_V1_lSwitch = (sum(LV1_lSwitchR)+sum(RV1_lSwitchL))/(rightLose(3)+leftLose(4));

y = [contra_M2_wSwitch contra_V1_wSwitch WSLS_Ctrl_wSwitchL;
    contra_M2_lSwitch  contra_V1_lSwitch  WSLS_Ctrl_lSwitchL ];

%x = categorical({'wStay','wSwitch', 'lStay', 'lSwitch'});
x = categorical({'wSwitch','lSwitch'});

figure
h = bar(x,y);
legend('M2', 'V1', 'Ctrl') 
title('contralateral changes')
ylabel('# occurances / # similar trials (leftLose) per region')
% plot % change from no stim

    LM2_wStayR_NS = (prop_LM2_wStayR-prop_Ctrl_wStayR)/prop_Ctrl_wStayR*100;
    LM2_wStayL_NS = (prop_LM2_wStayL-prop_Ctrl_wStayL)/prop_Ctrl_wStayL*100;
    RM2_wStayR_NS = (prop_RM2_wStayR-prop_Ctrl_wStayR)/prop_Ctrl_wStayR*100;
    RM2_wStayL_NS = (prop_RM2_wStayL-prop_Ctrl_wStayL)/prop_Ctrl_wStayL*100;
    LV1_wStayR_NS = (prop_LV1_wStayR-prop_Ctrl_wStayR)/prop_Ctrl_wStayR*100;
    LV1_wStayL_NS = (prop_LV1_wStayL-prop_Ctrl_wStayL)/prop_Ctrl_wStayL*100;
    RV1_wStayR_NS = (prop_RV1_wStayR-prop_Ctrl_wStayR)/prop_Ctrl_wStayR*100;
    RV1_wStayL_NS = (prop_RV1_wStayL-prop_Ctrl_wStayL)/prop_Ctrl_wStayL*100;


    LM2_wSwitchR_NS = (prop_LM2_wSwitchR-prop_Ctrl_wSwitchR)/prop_Ctrl_wSwitchR*100; 
    LM2_wSwitchL_NS = (prop_LM2_wSwitchL-prop_Ctrl_wSwitchL)/prop_Ctrl_wSwitchL*100;
    RM2_wSwitchR_NS = (prop_RM2_wSwitchR-prop_Ctrl_wSwitchR)/prop_Ctrl_wSwitchR*100;
    RM2_wSwitchL_NS = (prop_RM2_wSwitchL-prop_Ctrl_wSwitchL)/prop_Ctrl_wSwitchL*100; 
    LV1_wSwitchR_NS = (prop_LV1_wSwitchR-prop_Ctrl_wSwitchR)/prop_Ctrl_wSwitchR*100;
    LV1_wSwitchL_NS = (prop_LV1_wSwitchL-prop_Ctrl_wSwitchL)/prop_Ctrl_wSwitchL*100;
    RV1_wSwitchR_NS = (prop_RV1_wSwitchR-prop_Ctrl_wSwitchR)/prop_Ctrl_wSwitchR*100;
    RV1_wSwitchL_NS = (prop_RV1_wSwitchL-prop_Ctrl_wSwitchL)/prop_Ctrl_wSwitchL*100;


    LM2_lStayR_NS = (prop_LM2_lStayR-prop_Ctrl_lStayR)/prop_Ctrl_lStayR*100;
    LM2_lStayL_NS = (prop_LM2_lStayL-prop_Ctrl_lStayL)/prop_Ctrl_lStayL*100;
    RM2_lStayR_NS = (prop_RM2_lStayR-prop_Ctrl_lStayR)/prop_Ctrl_lStayR*100;
    RM2_lStayL_NS = (prop_RM2_lStayL-prop_Ctrl_lStayL)/prop_Ctrl_lStayL*100;
    LV1_lStayR_NS = (prop_LV1_lStayR-prop_Ctrl_lStayR)/prop_Ctrl_lStayR*100;
    LV1_lStayL_NS = (prop_LV1_lStayL-prop_Ctrl_lStayL)/prop_Ctrl_lStayL*100;
    RV1_lStayR_NS = (prop_RV1_lStayR-prop_Ctrl_lStayR)/prop_Ctrl_lStayR*100;
    RV1_lStayL_NS = (prop_RV1_lStayL-prop_Ctrl_lStayL)/prop_Ctrl_lStayL*100;
    
    LM2_lSwitchR_NS = (prop_LM2_lSwitchR-prop_Ctrl_lSwitchR)/prop_Ctrl_lSwitchR*100;
    LM2_lSwitchL_NS = (prop_LM2_lSwitchL-prop_Ctrl_lSwitchL)/prop_Ctrl_lSwitchL*100;
    RM2_lSwitchR_NS = (prop_RM2_lSwitchR-prop_Ctrl_lSwitchR)/prop_Ctrl_lSwitchR*100;
    RM2_lSwitchL_NS = (prop_RM2_lSwitchL-prop_Ctrl_lSwitchL)/prop_Ctrl_lSwitchL*100;
    LV1_lSwitchR_NS = (prop_LV1_lSwitchR-prop_Ctrl_lSwitchR)/prop_Ctrl_lSwitchR*100;
    LV1_lSwitchL_NS = (prop_LV1_lSwitchL-prop_Ctrl_lSwitchL)/prop_Ctrl_lSwitchL*100;
    RV1_lSwitchR_NS = (prop_RV1_lSwitchR-prop_Ctrl_lSwitchR)/prop_Ctrl_lSwitchR*100;
    RV1_lSwitchL_NS = (prop_RV1_lSwitchL-prop_Ctrl_lSwitchL)/prop_Ctrl_lSwitchL*100;


y = [LM2_wStayL_NS RM2_wStayL_NS LV1_wStayL_NS RV1_wStayL_NS; LM2_wStayR_NS RM2_wStayR_NS LV1_wStayR_NS RV1_wStayR_NS;
    LM2_wSwitchL_NS RM2_wSwitchL_NS LV1_wSwitchL_NS RV1_wSwitchL_NS; LM2_wSwitchR_NS RM2_wSwitchR_NS LV1_wSwitchR_NS RV1_wSwitchR_NS;
    LM2_lStayL_NS RM2_lStayL_NS LV1_lStayL_NS RV1_lStayL_NS; LM2_lStayR_NS RM2_lStayR_NS LV1_lStayR_NS RV1_lStayR_NS;
    LM2_lSwitchL_NS RM2_lSwitchL_NS LV1_lSwitchL_NS RV1_lSwitchL_NS; LM2_lSwitchR_NS RM2_lSwitchR_NS LV1_lSwitchR_NS RV1_lSwitchR_NS];

x = categorical({'wStayL', 'wStayR', 'wSwitchL', 'wSwitchR', 'lStayL', 'lStayR', 'lSwitchL', 'lSwitchR'});

figure
h = bar(x,y);
legend('LM2', 'RM2', 'LV1', 'RV1') 
title('22522 DbhxAi32 5 sessions combined')
ylabel('% change from no stim')



% plot % change from average of V1 control regions
    V1_trials = sum(LV1_trials)+sum(RV1_trials);
    prop_V1_lStayR = (sum(LV1_lStayR)+sum(RV1_lStayR))/V1_trials;
    prop_V1_lStayL = (sum(LV1_lStayL)+sum(RV1_lStayL))/V1_trials;
    prop_V1_wStayR = (sum(LV1_wStayR)+sum(RV1_wStayR))/V1_trials;
    prop_V1_wStayL = (sum(LV1_wStayL)+sum(RV1_wStayL))/V1_trials;
    prop_V1_lSwitchR = (sum(LV1_lSwitchR)+sum(RV1_lSwitchR))/V1_trials;
    prop_V1_lSwitchL = (sum(LV1_lSwitchL)+sum(RV1_lSwitchL))/V1_trials;
    prop_V1_wSwitchR = (sum(LV1_wSwitchR)+sum(RV1_wSwitchR))/V1_trials;
    prop_V1_wSwitchL = (sum(LV1_wSwitchL)+sum(RV1_wSwitchL))/V1_trials;



    LM2_wStayR_V1 = (prop_LM2_wStayR-prop_V1_wStayR)/prop_V1_wStayR*100;
    LM2_wStayL_V1 = (prop_LM2_wStayL-prop_V1_wStayL)/prop_V1_wStayL*100;
    RM2_wStayR_V1 = (prop_RM2_wStayR-prop_V1_wStayR)/prop_V1_wStayR*100;
    RM2_wStayL_V1 = (prop_RM2_wStayL-prop_V1_wStayL)/prop_V1_wStayL*100;
    LV1_wStayR_V1 = (prop_LV1_wStayR-prop_V1_wStayR)/prop_V1_wStayR*100;
    LV1_wStayL_V1 = (prop_LV1_wStayL-prop_V1_wStayL)/prop_V1_wStayL*100;
    RV1_wStayR_V1 = (prop_RV1_wStayR-prop_V1_wStayR)/prop_V1_wStayR*100;
    RV1_wStayL_V1 = (prop_RV1_wStayL-prop_V1_wStayL)/prop_V1_wStayL*100;


    LM2_wSwitchR_V1 = (prop_LM2_wSwitchR-prop_V1_wSwitchR)/prop_V1_wSwitchR*100; 
    LM2_wSwitchL_V1 = (prop_LM2_wSwitchL-prop_V1_wSwitchL)/prop_V1_wSwitchL*100; 
    RM2_wSwitchR_V1 = (prop_RM2_wSwitchR-prop_V1_wSwitchR)/prop_V1_wSwitchR*100; 
    RM2_wSwitchL_V1 = (prop_RM2_wSwitchL-prop_V1_wSwitchL)/prop_V1_wSwitchL*100;  
    LV1_wSwitchR_V1 = (prop_LV1_wSwitchR-prop_V1_wSwitchR)/prop_V1_wSwitchR*100; 
    LV1_wSwitchL_V1 = (prop_LV1_wSwitchL-prop_V1_wSwitchL)/prop_V1_wSwitchL*100; 
    RV1_wSwitchR_V1 = (prop_RV1_wSwitchR-prop_V1_wSwitchR)/prop_V1_wSwitchR*100; 
    RV1_wSwitchL_V1 = (prop_RV1_wSwitchL-prop_V1_wSwitchL)/prop_V1_wSwitchL*100; 


    LM2_lStayR_V1 = (prop_LM2_lStayR-prop_V1_lStayR)/prop_V1_lStayR*100;
    LM2_lStayL_V1 = (prop_LM2_lStayL-prop_V1_lStayL)/prop_V1_lStayL*100;
    RM2_lStayR_V1 = (prop_RM2_lStayR-prop_V1_lStayR)/prop_V1_lStayR*100;
    RM2_lStayL_V1 = (prop_RM2_lStayL-prop_V1_lStayL)/prop_V1_lStayL*100;
    LV1_lStayR_V1 = (prop_LV1_lStayR-prop_V1_lStayR)/prop_V1_lStayR*100;
    LV1_lStayL_V1 = (prop_LV1_lStayL-prop_V1_lStayL)/prop_V1_lStayL*100;
    RV1_lStayR_V1 = (prop_RV1_lStayR-prop_V1_lStayR)/prop_V1_lStayR*100;
    RV1_lStayL_V1 = (prop_RV1_lStayL-prop_V1_lStayL)/prop_V1_lStayL*100;
    
    LM2_lSwitchR_V1 = (prop_LM2_lSwitchR-prop_V1_lSwitchR)/prop_V1_lSwitchR*100; 
    LM2_lSwitchL_V1 = (prop_LM2_lSwitchL-prop_V1_lSwitchL)/prop_V1_lSwitchL*100;
    RM2_lSwitchR_V1 = (prop_RM2_lSwitchR-prop_V1_lSwitchR)/prop_V1_lSwitchR*100;
    RM2_lSwitchL_V1 = (prop_RM2_lSwitchL-prop_V1_lSwitchL)/prop_V1_lSwitchL*100;
    LV1_lSwitchR_V1 = (prop_LV1_lSwitchR-prop_V1_lSwitchR)/prop_V1_lSwitchR*100;
    LV1_lSwitchL_V1 = (prop_LV1_lSwitchL-prop_V1_lSwitchL)/prop_V1_lSwitchL*100;
    RV1_lSwitchR_V1 = (prop_RV1_lSwitchR-prop_V1_lSwitchR)/prop_V1_lSwitchR*100;
    RV1_lSwitchL_V1 = (prop_RV1_lSwitchL-prop_V1_lSwitchL)/prop_V1_lSwitchL*100;


y = [LM2_wStayL_V1 RM2_wStayL_V1 LV1_wStayR_V1 RV1_wStayL_V1; LM2_wStayR_V1 RM2_wStayR_V1 LV1_wStayR_V1 RV1_wStayR_V1;
    LM2_wSwitchL_V1 RM2_wSwitchL_V1 LV1_wSwitchL_V1 RV1_wSwitchL_V1; LM2_wSwitchR_V1 RM2_wSwitchR_V1 LV1_wSwitchR_V1 RV1_wSwitchR_V1;
    LM2_lStayL_V1 RM2_lStayL_V1 LV1_lStayL_V1 RV1_lStayL_V1; LM2_lStayR_V1 RM2_lStayR_V1 LV1_lStayR_V1 RV1_lStayR_V1;
    LM2_lSwitchL_V1 RM2_lSwitchL_V1 LV1_lSwitchL_V1 RV1_lSwitchL_V1; LM2_lSwitchR_V1 RM2_lSwitchR_V1 LV1_lSwitchR_V1 RV1_lSwitchR_V1];

x = categorical({'wStayL', 'wStayR', 'wSwitchL', 'wSwitchR', 'lStayL', 'lStayR', 'lSwitchL', 'lSwitchR'});

figure
h = bar(x,y);
legend('LM2', 'RM2', 'LV1', 'RV1') 
title('22522 DbhxAi32 5 sessions combined')
ylabel('% change from V1')

%% plot model fitting data
clearvars;
close all;
setup_figprop;

powerList = ["pre", "1.0mW","1.5mW","2.5mW"];
sz = 25; %size of scatter dots
colors = ["b", "r", "g", "y"];

% powerList = ["pre","1.5mW"];
% sz = 25; %size of scatter dots
% colors = ["b", "r"];

figure
for o = 1:size(powerList,2)
    c = colors(o);
    currPower = powerList(o);
    root_path = fullfile('C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\22522\', currPower);
    logfilepath = fullfile(root_path,'data');
    analysispath = fullfile(root_path,'analysis');
    dataIndex = makeDataIndex(logfilepath, analysispath);
    nSessions = size(dataIndex,1);
    model = table(...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1)...
        );

    model.Properties.VariableNames = {...
        'prob_WSLS',...
        'Q_RPE_a',...
        'Q_RPE_b',...
        'DQ_RPE_ar',...
        'DQ_RPE_b',...
        'DQ_RPE_an',...
        'FQ_RPE_ar',...
        'FQ_RPE_b',...
        'FQ_RPE_CK_a',...
        'FQ_RPE_CK_b',...
        'FQ_RPE_CK_tau',...
        'FQ_RPE_CK_phi'...
        };

    for l = 1:size(dataIndex,1)
        load(fullfile(dataIndex.BehPath{l},[dataIndex.LogFileName{l}(1:end-4),'_beh_cut.mat']));
        model.Session {l} = (dataIndex.LogFileName{l}(1:end-4));
        model.prob_WSLS {l} = fitpar.WSLS;
        model.Q_RPE_a {l} = fitpar.Q_RPE(1);
        model.Q_RPE_b {l} = fitpar.Q_RPE(2);
        model.DQ_RPE_ar {l} = fitpar.DQ_RPE(1);
        model.DQ_RPE_b {l} = fitpar.DQ_RPE(2);
        model.DQ_RPE_an {l} = fitpar.DQ_RPE(3);
        model.FQ_RPE_ar {l} = fitpar.FQ_RPE(1);
        model.FQ_RPE_b {l} = fitpar.FQ_RPE(2);
        model.FQ_RPE_CK_a {l} = fitpar.FQ_RPE_CK(1);
        model.FQ_RPE_CK_b {l} = fitpar.FQ_RPE_CK(2);
        model.FQ_RPE_CK_tau {l} = fitpar.FQ_RPE_CK(3);
        model.FQ_RPE_CK_phi {l} = fitpar.FQ_RPE_CK(4);
    end

    for f = 1:12
        name = string(model.Properties.VariableNames (f));
        hold on
        subplot(3,4,f)
        x = zeros(1,nSessions);
        for w = 1:nSessions
            temp = getfield(model, name);
            x (w) = temp {w,1};
        end
        y = ones(1,nSessions);
        scatter(y,x,sz,c,'filled');
        title(name)
        top = max(x);
        ax = [0 top];
        ylim(ax)
    end
end
%% plot basic stats by session
clearvars;
close all;
setup_figprop;
animal = '22522\'; 
powerList = ["pre", "1.0mW","1.5mW","2.5mW"];
sz = 25; %size of scatter dots
colors = ["bl", "r", "g", "b"];
top = 0;
bottom = 1000;
figure
for n = 1:size(powerList,2)
    c = colors(n);
    currPower = powerList(n);
    root_path = fullfile('C:\Users\Cornell\Documents\MATLAB\MP_GRABNE-main\MP_Stim\', animal, currPower);
    logfilepath = fullfile(root_path,'data');
    analysispath = fullfile(root_path,'analysis');
    dataIndex = makeDataIndex(logfilepath, analysispath);

    nSessions = size(dataIndex,1);
    session = table(...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1),...
        cell(nSessions,1)...
        );

    session.Properties.VariableNames = {...
        'ntrial',...
        'nmiss',...
        'nreward',...
        'rewardRate'...
        };

    for j = 1:nSessions
        load(fullfile(dataIndex.BehPath{j},[dataIndex.LogFileName{j}(1:end-4),'_beh_cut.mat']));
        session.ntrial {j} = sessionData.nTrials;
        session.nmiss {j} = sum((trialData.outcome == 77)); 
        session.nreward {j} = sum((trialData.outcome == 100)) + sum((trialData.outcome == 111));
        session.rewardRate {j} = (session.nreward {j} / (session.ntrial {j} - session.nmiss {j}))*100;
    end
    for k = 1:4
        hold on
        name = string(session.Properties.VariableNames (k));
        subplot(2,2,k)
        temp = cell2mat(getfield(session, name));
        scatter(n(ones(1,nSessions)), temp, sz, c,'filled')
        if k == 1
            ay = [0 700];
            yticks([200 400 600])
        end
        if k == 2
            ay = [0 20];
        end 
        if k == 3
            ay = [0 400];
            yticks([100 200 300])
        end
        if k == 4
            ay = [0 60];
            yticks([10 20 30 40 50])
        end
        ylim(ay)
        xlim([0 5])
        xticks([1 2 3 4])
        xticklabels(powerList)
        title(name)
    end 
end 
