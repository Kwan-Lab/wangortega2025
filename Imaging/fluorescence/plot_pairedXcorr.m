function plot_pairedXcorr(xCorrBoot,regCoeff,sigVar,t,ind1,ind2,varName,save_path_fluo)

% plot cross-correlation results with regression coefficient of each
% varibles
figure;
subplot(3,1,1)
plot(t,xCorrBoot{ind1,ind2}.coeff_bootave);
gray = [0.7 0.7 0.7];
hold on;
errorshade(t, xCorrBoot{ind1,ind2}.boothigh, xCorrBoot{ind1,ind2}.bootlow,gray);
hold on;
plot(t,xCorrBoot{ind1,ind2}.coeff_bootave,'k.-');
set(gca,'box','off');
title(['Cross-correlation of ',varName{ind1},' and ',varName{ind2}]);

subplot(3,1,2)
plot(t,regCoeff{ind1}.coeff_bootave,'Color',[0 0.4470 0.7410]);
hold on;plot(t,regCoeff{ind2}.coeff_bootave,'Color',[0.8500 0.3250 0.0980]);
legend(varName{ind1},varName{ind2});
errorshade(t, regCoeff{ind1}.boothigh, regCoeff{ind1}.bootlow,[0 0.4470 0.7410],0.5);
errorshade(t, regCoeff{ind2}.boothigh, regCoeff{ind2}.bootlow,[0.8500 0.3250 0.0980],0.5);
plot(t,regCoeff{ind1}.coeff_bootave,'Color',[0 0.4470 0.7410]);
hold on;plot(t,regCoeff{ind2}.coeff_bootave,'Color',[0.8500 0.3250 0.0980]);
set(gca,'box','off');
legend('box','off');
title('Regression coefficient');



% get A:sig/notsig; B:notsig/sig; C:sig/sig
SNS = squeeze(sum(sum(sigVar(:,:,:,:,ind1),1),2)/(size(sigVar,1)*size(sigVar,2)));
NSS = squeeze(sum(sum(sigVar(:,:,:,:,ind2),1),2)/(size(sigVar,1)*size(sigVar,2)));
SS = squeeze(sum(sum((sigVar(:,:,:,:,ind1)&sigVar(:,:,:,:,ind2)),1),2)/(size(sigVar,1)*size(sigVar,2)));

%bootstrap
sns.coeff = SNS'; sns = getBootstrp(sns,0,0.05);
   nss.coeff = NSS'; nss = getBootstrp(nss,0,0.05);
ss.coeff = SS'; ss = getBootstrp(ss,0,0.05);

subplot(3,1,3)

plot(t,sns.coeff_bootave,'Color',[0 0.4470 0.7410]);
hold on;plot(t,nss.coeff_bootave,'Color',[0.8500 0.3250 0.0980]);
hold on;plot(t,ss.coeff_bootave,'k-');

legend(varName{ind1},varName{ind2},'Overlap');
legend('box','off');
errorshade(t, ss.boothigh, ss.bootlow,[0.7 0.7 0.7]);
errorshade(t, sns.boothigh, sns.bootlow,[0 0.4470 0.7410],0.5);

errorshade(t, nss.boothigh, nss.bootlow,[0.8500 0.3250 0.0980],0.5);

% errorshade(t, regCoeff{ind2}.boothigh, regCoeff{ind2}.bootlow,[0.8500 0.3250 0.0980],0.5);
plot(t,ss.coeff_bootave,'k-');
plot(t,sns.coeff_bootave,'Color',[0 0.4470 0.7410]);
hold on;plot(t,nss.coeff_bootave,'Color',[0.8500 0.3250 0.0980]);
set(gca,'box','off');
title('Fraction of significant grids');
%hold on;plot(t,regCoeff{ind2}.coeff_bootave,'Color',[0.8500 0.3250 0.0980]);
       

xlabel('Time from cue');

   print(gcf,'-dpng',fullfile(save_path_fluo,['Cross-correlation of ',varName{ind1},' and ',varName{ind2}]));
    saveas(gcf, fullfile(save_path_fluo,['Cross-correlation of ',varName{ind1},' and ',varName{ind2}]), 'fig');
    saveas(gcf, fullfile(save_path_fluo,['Cross-correlation of ',varName{ind1},' and ',varName{ind2}]), 'svg');