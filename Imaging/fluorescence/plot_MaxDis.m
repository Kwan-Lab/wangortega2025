function plot_MaxDis(input,save_path_fluo,label)

%bootstrap the results and plot it
input = getBootstrp(input, 0, 0.05);

figure;plot(input.t,input.coeff_bootave);
hold on;
gray=[0.7 0.7 0.7];
errorshade(input.t,input.bootlow,input.boothigh,gray);
hold on;
plot(input.t,input.coeff_bootave,'k-');
set(gca,'box','off');


% figure;
% for ii=1:size(input.coeff,1)
%     hold on;plot(input.t,input.coeff(ii,:));
% end