
function MP_GRAB_behComp(dataIndex_ACh, dataIndex_NE, save_path_ACh)

logInd = 1;
nFiles_NE = size(dataIndex_NE,1);
nFiles_ACh = size(dataIndex_ACh,1);

for ii = 1:nFiles_NE
     
    
        load(fullfile(dataIndex_NE.BehPath{ii},'beh_cut.mat'));


        % calculate ITI time for every session
        iti_time = zeros(1, length(trialData.cueTimes)-1);
        
        for tt=1:length(trialData.cueTimes)-1
            iti_time(tt) = trialData.cueTimes(tt+1) - trialData.outcomeTimes(tt);
        end
        iti_trueTime_NE{ii} = iti_time;
        %iti_num(logInd) = length(trialData.itiTimes);
        trial_num_NE(ii) = length(trialData.cueTimes);
        
        lick_trType_array_NE{ii}=lick_trType;
    
        lregRCUC_array_NE{ii}=lregRCUC_output;
        lregCRInt_array_NE{ii}=lregCRInt_output;
    
        iti_array_NE{ii}=iti_trType;
        respTime_array_NE{ii}=respTime_trType;
        trueRespTime_NE{ii} = trialData.rt;
        choiceBySession_NE{ii} = stats;
        
        nTrial_array_NE(ii)=sum(stats.c(:,1)==-1)+sum(stats.c(:,1)==1);
        entro_array_NE(ii)=entro;
        rrate_array_NE(ii)=sum(stats.r==1)/(sum(stats.r==1)+sum(stats.r==0));
        
      
    end

    close all


for ii = 1:nFiles_ACh
     
    
        load(fullfile(dataIndex_ACh.BehPath{ii},'beh_cut.mat'));


        % calculate ITI time for every session
        iti_time = zeros(1, length(trialData.cueTimes)-1);
        
        for tt=1:length(trialData.cueTimes)-1
            iti_time(tt) = trialData.cueTimes(tt+1) - trialData.outcomeTimes(tt);
        end
        iti_trueTime_ACh{ii} = iti_time;
        %iti_num(logInd) = length(trialData.itiTimes);
        trial_num_ACh(ii) = length(trialData.cueTimes);
        
        lick_trType_array_ACh{ii}=lick_trType;
    
        lregRCUC_array_ACh{ii}=lregRCUC_output;
        lregCRInt_array_ACh{ii}=lregCRInt_output;
    
        iti_array_ACh{ii}=iti_trType;
        respTime_array_ACh{ii}=respTime_trType;
        trueRespTime_ACh{ii} = trialData.rt;
        choiceBySession_ACh{ii} = stats;
       
        nTrial_array_ACh(ii)=sum(stats.c(:,1)==-1)+sum(stats.c(:,1)==1);
        entro_array_ACh(ii)=entro;
        rrate_array_ACh(ii)=sum(stats.r==1)/(sum(stats.r==1)+sum(stats.r==0));
        
      
    end


%% plot the behavior summary of ACh and NE group
Group = [ones(1,nFiles_NE),2*ones(1,nFiles_ACh)];
figure;

% average trials 
subplot(2,3,1); hold on;
boxplot([nTrial_array_NE,nTrial_array_ACh]',Group,'PlotStyle','compact');
%boxplot(nTrial_array,'Colors','k','Notch','off');
ylim([0 1000]); 
ylabel(['Trials performed']);
set(gca,'box','off') 
[h,p] = ttest2(nTrial_array_NE,nTrial_array_ACh)

% entropy
subplot(2,3,2); hold on;
boxplot([entro_array_NE,entro_array_ACh]',Group','PlotStyle','compact');

plot([-1 3],[3 3],'k--','LineWidth',2);
ylim([2 3.1]); 
ylabel('Entropy (bits)');
set(gca,'box','off') 
[h,p] = ttest2(entro_array_NE,entro_array_ACh)

% reward rate
subplot(2,3,3); hold on;
%plot(rand(1,numel(rrate_array)),100*rrate_array,'k^','MarkerSize',15);
%boxplot(rrate_array*100,'Colors','k','Symbol','k+','Notch','off','Labels',{[char(956),'=',num2str((mean(rrate_array))*100),'%']});
boxplot([rrate_array_NE*100,rrate_array_ACh*100]',Group','PlotStyle','compact');
plot([-1 3],[50 50],'k--','LineWidth',2);
ylim([20 65]);
ylabel('Reward rate (%)');
set(gca,'box','off') 
[h,p]=ttest2(rrate_array_NE*100,rrate_array_ACh*100)

[h,p]=ttest(rrate_array_NE*100,50)
[h,p]=ttest(rrate_array_ACh*100,50)

print(gcf,'-dpng',fullfile(save_path_ACh,'summary-All'));    %png format
saveas(gcf,fullfile(save_path_ACh,'summary-All'), 'fig');
saveas(gcf, fullfile(save_path_ACh,'summary-All'),'svg');


