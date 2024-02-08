function plot_xcorrCoeff(input1,input2,input3,thresh,save_path_fluo,label)

input2 = getBootstrp(input2, 0, 0.05);
input3 = getBootstrp(input3, 0, 0.05);

figure;
subplot(3,1,1)
plot(input1.t,sum(sum(input1.coeff<thresh,1),3)/(size(input1.coeff,1)*size(input1.coeff,3)),'-k');
ylabel('Significant Grids');
title(['Fraction of significant grids (',label,')']);
set(gca,'box','off')

subplot(3,1,2)
plot(input2.t,input2.coeff_bootave);
hold on;
gray=[0.7 0.7 0.7];
errorshade(input2.t,input2.bootlow,input2.boothigh,gray);
hold on;
plot(input2.t,input2.coeff_bootave,'k-');
set(gca,'box','off');
ylabel('Coefficient');
title(['Cross-correlation coefficient (',label,')']);
%ylim([-0.2,1]);

subplot(3,1,3)
plot(input3.t,input3.coeff_bootave);
hold on;
gray=[0.7 0.7 0.7];
errorshade(input3.t,input3.bootlow,input3.boothigh,gray);
hold on;
plot(input3.t,input3.coeff_bootave,'k-');
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('Offset');
title(['Offset of maximum correlation coefficient (',label,')']);


print(gcf,'-dpng',fullfile(save_path_fluo,[label,' xcorrCoeff']));
saveas(gcf, fullfile(save_path_fluo,[label,' xcorrCoeff']), 'fig');
%bootstrap the results and plot it

% plot histogram to show the distribution of the offsets
% input 3 is offset
% making a histogram
maxValue = max(input3.coeff(:));
minValue = 0;
binStep = 1;
numBin = ceil(maxValue)/binStep;
offsetHist = zeros(numBin,size(input3.coeff,2));
matHist3 = zeros(2,size(input3.coeff,1)*size(input3. coeff,2));
tHist3 = [];
vHist3 = [];
for tt = 1:size(input3.coeff,2)
    tHist3 = [tHist3;ones(size(input3.coeff,1),1)*input3.t(tt)];
    vHist3 = [vHist3;input3.coeff(:,tt)];
end
matHist = [tHist3';vHist3']';
edges = {0:1:maxValue input3.t};
figure;hist3([vHist3';tHist3']','Edges',edges,'CDataMode','auto','FaceColor','interp')
ylabel('Time from cue(s)');
xlabel('Offsets')
title(['Distribions of offsets - ',label]);
colorbar;


print(gcf,'-dpng',fullfile(save_path_fluo,[label,' xcorrCoeff-distribution']));
saveas(gcf, fullfile(save_path_fluo,[label,' xcorrCoeff-distribution']), 'fig');

% individual animals
% tHist3 = [];
% vHist3 = [];
% for tt = 1:size(input3.coeff,2)
%     tHist3 = [tHist3;ones(10,1)*input3.t(tt)];
%     vHist3 = [vHist3;input3.coeff(21:30,tt)];
% end
% matHist = [tHist3';vHist3']';
% edges = {0:1:maxValue input3.t};
% figure;hist3([vHist3';tHist3']','Edges',edges,'CDataMode','auto','FaceColor','interp')
% % figure;plot(mean(offsetHist,2))
% % for hh = 1:numBin
%     startV = 0+binStep*(hh-1);
%     endV = binStep*hh;
%     offsetHist(hh,tt) = sum(input3.coeff(:,tt)>=startV & input3.coeff(:,tt)<endV);
% end
% end
% figure;surf(offsetHist);

% figure;
% for ii=1:size(input.coeff,1)
%     hold on;plot(input.t,input.coeff(ii,:));
% end