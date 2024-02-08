function audi_GRAB_summary(dataIndex,savefigpath)

% summarize all subjects
% Ach right before MP recording: 14, 19,22,25
% NE tone response
% load PSTH results, calculate the mean (and STE)
Ind = 1:size(dataIndex,1);
ACh = Ind(strcmp(dataIndex.GRAB,'ACh'));
NE = Ind(strcmp(dataIndex.GRAB,'NE'));

nFiles = size(dataIndex,1);

% rise time: from cue to 90% amplitude
% decay time, from peak to exponential decay constant
PSTH_ave_ACh = [];
peakTime_ACh = [];
decayTau_ACh = [];
PSTH_ave_AChPup = [];
PSTH_ave_NE = [];
PSTH_ave_NEPup= [];
peakTime_NE= [];
decayTau_NE = [];

for ii = 1:length(ACh)
     fn_beh = dir(fullfile(dataIndex.BehPath{ACh(ii)},'*beh.mat'));
      date = num2str(dataIndex.DateNumber(ACh(ii)));
    fn_pup = dir(fullfile(dataIndex.BehPath{ACh(ii)},['*',date(1:6),'*_pup.mat']));
    savepsthname = fullfile(fn_beh.folder,'psthMat.mat');
      savepupname = fullfile(fn_pup.folder,'psthPup.mat');
    load(savepsthname);
    load(savepupname);
    psth_ACh.t = psth_panel(1).sig{1}.t;
    psth_AChPup.t = psthpupil_panel(1).sig{1}.t;
    sesACh = []; tPeak = []; tDecay = [];
    for jj = 1:length(psth_panel)
        PSTH_ave_ACh = [PSTH_ave_ACh;psth_panel(jj).sig{1}.signal'];
        % fit for rise time and decay time
        sesACh = [sesACh;psth_panel(jj).sig{1}.signal'];
        [~,ind] = max(psth_panel(jj).sig{1}.signal);
        tPeak = [tPeak,psth_ACh.t(ind)];
        f = fit(psth_ACh.t(ind:end),psth_panel(jj).sig{1}.signal(ind:end),'exp1');
        tDecay = [tDecay, f.b];
    end
    %avecurve = nanmean(sesACh,1);
    
    % find peak time and decau tau
    %[~,ind] = max(avecurve);
    peakTime_ACh = [peakTime_ACh,nanmedian(tPeak)];

    % decay tau
    decayTau_ACh = [decayTau_ACh, nanmedian(tDecay)];

    %PSTH_ave_AChPup = [PSTH_ave_AChPup;psthpupil_panel(1).sig{1}.signal'];
end

for ii = 1:length(NE)
     fn_beh = dir(fullfile(dataIndex.BehPath{NE(ii)},'*beh.mat'));
      date = num2str(dataIndex.DateNumber(NE(ii)));
    fn_pup = dir(fullfile(dataIndex.BehPath{NE(ii)},['*',date(1:6),'*_pup.mat']));
    savepsthname = fullfile(fn_beh.folder,'psthMat.mat');
      %savepupname = fullfile(fn_pup.folder,'psthPup.mat');
    load(savepsthname)
      %load(savepupname);
    psth_NE.t = psth_panel(1).sig{1}.t;
     %psth_NEPup.t = psthpupil_panel(1).sig{1}.t;
     sesNE = [];tPeak=[];tDecay = [];
    for jj = 1:length(psth_panel)
        PSTH_ave_NE = [PSTH_ave_NE;psth_panel(jj).sig{1}.signal'];
        sesNE = [sesNE;psth_panel(jj).sig{1}.signal'];
  
        [~,ind] = max(psth_panel(jj).sig{1}.signal);
        tPeak = [tPeak,psth_NE.t(ind)];
        if ind<80
            f = fit(psth_NE.t(ind:end),psth_panel(jj).sig{1}.signal(ind:end),'exp1');
            tDecay = [tDecay, f.b];
        else
            tDecay = [tDecay,NaN];
        end
    end
    %avecurve = nanmean(sesACh,1);
    
    % find peak time and decau tau
    %[~,ind] = max(avecurve);
    peakTime_NE = [peakTime_NE, nanmedian(tPeak)];

    % decay tau
    decayTau_NE = [decayTau_NE, nanmedian(tDecay)];
end

% bootstrap to get the average
psth_ACh.coeff= PSTH_ave_ACh;
psth_NE.coeff= PSTH_ave_NE;
%psth_AChPup.coeff = PSTH_ave_AChPup;
%psth_NEPup.coeff = PSTH_ave_NEPup;
% use bootstrp to get coefficient
psth_ACh = getBootstrp(psth_ACh, 0, 0.05);
psth_NE = getBootstrp(psth_NE, 0, 0.05);
%psth_AChPup = getBootstrp(psth_AChPup, 0, 0.05);
%psth_NEPup = getBootstrp(psth_NEPup, 0, 0.05);
% normalize baseline
baseACh = nanmean(psth_ACh.coeff_bootave(psth_ACh.t<0));
baseNE = nanmean(psth_NE.coeff_bootave(psth_NE.t<0));
%baseAChPup = nanmean(psth_AChPup.coeff_bootave(psth_AChPup.t<0));
%baseNEPup = nanmean(psth_NEPup.coeff_bootave(psth_NEPup.t<0));

figure;
plot(psth_ACh.t,psth_ACh.coeff_bootave-baseACh,'Color',[255 189 53]/255)
hold on;
plot(psth_NE.t,psth_NE.coeff_bootave-baseNE,'Color',[63,167,150]/255)
hold on;
errorshade(psth_ACh.t,psth_ACh.bootlow-baseACh,psth_ACh.boothigh-baseACh,[255 189 53]/255,0.5);
hold on;
errorshade(psth_NE.t,psth_NE.bootlow-baseNE,psth_NE.boothigh-baseNE,[63,167,150]/255,0.5);
hold on;
plot(psth_NE.t,psth_NE.coeff_bootave-baseNE,'Color',[63,167,150]/255)
hold on;
plot(psth_ACh.t,psth_ACh.coeff_bootave-baseACh,'Color',[255,189,53]/255)

% plot pupil
% hold on;plot(psth_AchPup.t, psth_AchPup.coeff_bootave-baseAchPup);
% hold on;plot(psth_NEPup.t, psth_NEPup.coeff_bootave-baseNEPup);
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('df/f');
legend('Ach','NE');
ax1 = gca; 


title('GRAB tone response');
print(gcf,'-dpng',fullfile(savefigpath,'GRAB-tone'));    %png format
saveas(gcf, fullfile(savefigpath,'GRAB-tone'), 'fig');
saveas(gcf, fullfile(savefigpath,'GRAB-tone'),'svg');

%% box plot for peak time and decay time

groupInd = [1*ones(1,length(NE)),2*ones(1,length(ACh))];
figure;
subplot(1,2,1)
boxplot([peakTime_NE,peakTime_ACh],groupInd,'PlotStyle','compact');
set(gca,'box','off')
ylabel('PeakTime');
subplot(1,2,2)
boxplot([decayTau_NE,decayTau_ACh],groupInd,'PlotStyle','compact')
set(gca,'box','off')
ylabel('DecayTime');
print(gcf,'-dpng',fullfile(savefigpath,'GRAB-tone-stat'));    %png format
saveas(gcf, fullfile(savefigpath,'GRAB-tone-stat'), 'fig');
saveas(gcf, fullfile(savefigpath,'GRAB-tone-stat'),'svg');
p = ranksum(peakTime_ACh, peakTime_NE)
p = ranksum(decayTau_ACh, decayTau_NE)


