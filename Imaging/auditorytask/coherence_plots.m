function coherence_plots(dataIndex,savefigpath)
% plot spontaneous pupil & fluorescent dynamics

nFiles = size(dataIndex,1);
cxyAll_NE.coeff = [];
cxyAll_Ach.coeff = [];

pxyAll_Ach.coeff = [];
pxyAll_NE.coeff = [];

Ach = [4:14];
NE = [1:3,15:18];

%% load Ach
for ii = 1:length(Ach)
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{Ach(ii)},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(Ach(ii)));
    pup_name = fullfile(dataIndex.BehPath{Ach(ii)},['*',date(1:6),'*_pup.mat']);
    dff_name = fullfile(fn_beh.folder,'dff.mat');
    
    fn_pup = dir(pup_name);
    fn_dff = dir(dff_name);
   
    
    if length(fn_pup) == 1 & length(fn_dff) == 1
         savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        % coherence analysis
       
       
        
        savematname='coherence.mat';
        cohfilename = fullfile(savematpath,savematname);
        load(cohfilename);

            cxyAll_Ach.coeff = [cxyAll_Ach.coeff; cxy.coeff];
            pxyAll_Ach.coeff = [pxyAll_Ach.coeff; pxy.coeff];
    end
        
end

%% load NE

for ii = 1:length(NE)
    
    % load behavior files
     fn_beh = dir(fullfile(dataIndex.BehPath{NE(ii)},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(NE(ii)));
    pup_name = fullfile(dataIndex.BehPath{NE(ii)},['*',date(1:6),'*_pup.mat']);
    dff_name = fullfile(fn_beh.folder,'dff.mat');
    
    fn_pup = dir(pup_name);
    fn_dff = dir(dff_name);
   
    
    if length(fn_pup) == 1 & length(fn_dff) == 1
         savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        % coherence analysis
       
       
        
        savematname='coherence.mat';
        cohfilename = fullfile(savematpath,savematname);
        load(cohfilename);

            cxyAll_NE.coeff = [cxyAll_NE.coeff; cxy.coeff];
            pxyAll_NE.coeff = [pxyAll_NE.coeff; pxy.coeff];

    end
        
end

%% plot
% plot total of 588 grids
cxyAll_NE = getBootstrp(cxyAll_NE, 0, 0.05);
cxyAll_Ach = getBootstrp(cxyAll_Ach, 0, 0.05);        

pxyAll_NE = getBootstrp(pxyAll_NE, 0, 0.05);
pxyAll_Ach = getBootstrp(pxyAll_Ach,0,0.05);

c_NE = [255 189 53]/255;
c_Ach = [63,167,150]/255;
figure;



%patch([fc fliplr(fc)], [cxy.boothigh  fliplr(cxy.bootlow)], [0.7 0.7 0.7])
plot(fc(2:end),cxyAll_NE.coeff_bootave(2:end),'-','Color',c_NE);
hold on;
plot(fc(2:end),cxyAll_Ach.coeff_bootave(2:end),'-','Color',c_Ach);
legend({'NE','Ach'});
legend boxoff
hold on;
errorshade(fc(2:end),cxyAll_NE.bootlow(2:end),cxyAll_NE.boothigh(2:end),c_NE,0.2);
hold on;
plot(fc(2:end),cxyAll_NE.coeff_bootave(2:end),'-','Color',c_NE);

hold on;
errorshade(fc(2:end),cxyAll_Ach.bootlow(2:end),cxyAll_Ach.boothigh(2:end),c_Ach,0.2);
hold on;
plot(fc(2:end),cxyAll_Ach.coeff_bootave(2:end),'-','Color',c_Ach);

xlim([0.001, 1]);
set(gca,'box','off')
set(gca, 'XScale', 'log')

xlabel('Frequency (Hz)');
ylabel('Coherence with pupil');

print(gcf,'-dpng',fullfile(savefigpath,['coh_spon_pupil_NE_Ach']));
saveas(gcf, fullfile(savefigpath,['coh_spon_pupil_NE_Ach']), 'fig');
saveas(gcf, fullfile(savefigpath,['coh_spon_pupil_NE_Ach']), 'svg');

figure;
plot(F(2:end),pxyAll_NE.coeff_bootave(2:end),'-','Color',c_NE);
hold on;
plot(F(2:end),pxyAll_Ach.coeff_bootave(2:end),'-','Color',c_Ach);
legend({'NE','Ach'});
legend boxoff
hold on;
errorshade(F(2:end),pxyAll_Ach.bootlow(2:end),pxyAll_Ach.boothigh(2:end),c_Ach,0.2);
hold on;
plot(F(2:end),pxyAll_Ach.coeff_bootave(2:end),'-','Color',c_Ach);

hold on;
errorshade(F(2:end),pxyAll_NE.bootlow(2:end),pxyAll_NE.boothigh(2:end),c_NE,0.2);
hold on;
plot(F(2:end),pxyAll_NE.coeff_bootave(2:end),'-','Color',c_NE);

xlim([0.001, 1]);
set(gca,'box','off')
set(gca, 'XScale', 'log')
xlabel('Frequency');
ylabel('Lag (\times\pi rad)')
print(gcf,'-dpng',fullfile(savefigpath,['xspec_spon_pupil_NE_Ach']));
saveas(gcf, fullfile(savefigpath,['xspec_spon_pupil_NE_Ach']), 'fig');
saveas(gcf, fullfile(savefigpath,['xspec_spon_pupil_NE_Ach']), 'svg');

end