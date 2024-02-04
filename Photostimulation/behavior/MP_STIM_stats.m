function dat = MP_STIM_stats (analysispath, dataIndex, plotChoice, regions, region, sidePreference)
save_path = analysispath; 
for j = 1: size(dataIndex,1)
    load(fullfile(dataIndex.BehPath{j},[dataIndex.LogFileName{j}(1:end-4),'_beh.mat']));

    if plotChoice == 1
   % analysis of behavioral performance
        % plot choice behavior - whole sessions
        n_plot = 100*ceil(numel(stats.c)/100); %plot up to the nearest 100 trials
        % if n_plot > 1000   %for computer simulations, #trials can be too high to plot effectively
        %     n_plot = 1000;
        % end
        name = cellstr(dataIndex.Animal{j}); 
        plot_session_task(stats,n_plot,name)
        jj = num2str(j); 
        filename = append('session_', jj); 
        print(gcf,'-dpng',fullfile(save_path,filename));
        saveas(gcf, fullfile(save_path,filename), 'fig');
    end 

    if sidePreference == 1 %if you want side preference analyzed per session

        if sum(trialData.response ==2) == 0 || sum(trialData.response == 3) == 0
            %if no choices on one side
        
        
        elseif sum(trialData.response == 2) / (sum(trialData.response == 3) + sum(trialData.response == 2)) > .5

                %     OUTCOME.REWARDLEFT = 100;
                %     OUTCOME.REWARDRIGHT = 111;
                %     OUTCOME.NOREWARDLEFT = 101;
                %     OUTCOME.NOREWARDRIGHT = 110;
                %     OUTCOME.REWARDMANUAL = 10;
                %     OUTCOME.MISS = 77;      %miss
                %     RESP.LEFT=2;
                %     RESP.RIGHT=3;

                %switch left and right response values
                maskLeft = [];
                maskLeft = (trialData.response == 2)*3; %change positions for left choices to '3' for right choice
                
                maskRight = [];
                maskRight = (trialData.response == 3)*2; 

                maskLeftWin = [];
                maskLeftWin = (trialData.outcome == 100)*111; %change left reward to right reward

                maskLeftLose = [];
                maskLeftLose = (trialData.outcome == 101)*110; %change left noreward to right noreward

                maskRightWin =[]; 
                maskRightWin = (trialData.outcome == 111)*100;

                maskRightLose =[]; 
                maskRightLose = (trialData.outcome == 110)*101;

                maskMiss = [];
                maskMiss = (trialData.outcome == 77)*77; %keep miss 

                maskManual = [];
                maskManual = (trialData.outcome == 10)*10; %keep manual rewards

                trialData.response = [];
                trialData.response = maskLeft + maskRight;
                
                trialData.outcome = [];
                trialData.outcome = maskManual + maskMiss + maskRightLose ...
                    + maskRightWin + maskLeftLose + maskLeftWin;

                
                
        end
    end 
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

    stats.Animal = dataIndex.Animal; 
    stats.Strain = dataIndex.Strain; 
       
    
    for kk = 1:size(regions,2) % region 1 (ALM L), 2 (ALM R), 3 (V1 L), 4 (V1 R), 5 (NONE) 
       if size(regions,2) == 2 %catch for one region stimulated
           if kk == 1
               if strcmp(region, 'LM2') == 1
                   regionInd = 1;
               elseif strcmp(region, 'RM2') == 1
                   regionInd = 2;
               elseif strcmp(region, 'LV1') == 1
                   regionInd = 3;
               elseif strcmp(region, 'RV1') == 1
                   regionInd = 4;
               end
           elseif kk == 2
               regionInd = 0;
           end
       else
            if kk == size(regions,2) %no stim
                regionInd = 0;
            else 
                regionInd = kk;
            end
       end 
        
         dat(kk).beh_output_WSLS{j}=WSLS_byRegion(trialData, stats, regionInd, stats.st);
         dat(kk).beh_output_stim{j}=stats;
    end
end
