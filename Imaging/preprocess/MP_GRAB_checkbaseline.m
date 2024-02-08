function MP_GRAB_checkbaseline(dataIndex)

% simple plots, check alignment

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % check if dff is already computed
    fn_dff = dir(fullfile(fn_beh.folder, 'dff.mat'));
    
    % check to see that imaging and behavior timing is sync'ed properly
    savefluofigpath = fullfile(fn_beh.folder,'figs_analysis');
    if ~exist(savefluofigpath)
        mkdir(savefluofigpath)
    end
    
    if ~isempty(fn_dff)
     load(fullfile(fn_dff.folder,fn_dff.name));
    
    tic
    if isfield(cells,'normdFF')
        display('Movement corrected')
    else
        cells = MP_GRAB_normalizedFF(cells, trialData);
    end
    toc
    %figure;histogram(cells.dFF{1});
    
    % need find the change point, and substract the average value difference of
    % the two groups from the signals with higher baseline
    % there are roughly 200 change points in one session, that's start with
    % this
    
    % use the grand mean? 
    %edges = 0:0.01:0.25;
    %figure;histogram(aveBaseBC,edges)
%     
  
    
    % two peaks
      savefluofigpath = fullfile(fn_beh.folder,'figs_analysis');
    aveBaseline = zeros(1, length(trialData.cue));
    for tt = 1:length(trialData.cue)
        aveBaseline(tt) = nanmean(cells.dFF{1}(cells.t>trialData.cueTimes(tt)-1 & cells.t<trialData.cueTimes(tt)));
    end
    figure; histogram(aveBaseline);
    savefilename = 'baselinecheck';
    print(gcf,'-dpng',fullfile(savefluofigpath,savefilename));    %png format
    saveas(gcf,fullfile(savefluofigpath,savefilename), 'fig');
    close;
    % detect change point? too sensitive sometimes, try to align byhel cue
    % onset
    save(fullfile(fn_dff.folder,fn_dff.name),'cells','-v7.3');
    else
        display('Fluorescent file not found');
    end
end
