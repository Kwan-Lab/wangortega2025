function fluo_pupil_GRAB_correlation(dataIndex)

%% calculate correlation between fluorescent and pupillometry with time lag
nFiles = size(dataIndex,1);


 colors=cbrewer('div','RdBu',256);
colors=flipud(colors);
for ii = 1:nFiles
    
    % load behavior files
     fn_beh = dir(fullfile(dataIndex.BehPath{ii},'*beh.mat'));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % load pupil files
    date = num2str(dataIndex.DateNumber(ii));
    pup_name = fullfile(dataIndex.BehPath{ii},['*',date(1:6),'*_pup.mat']);
    dff_name = fullfile(fn_beh.folder,'dff.mat');
    
    fn_pup = dir(pup_name);
    fn_dff = dir(dff_name);
   
    
    if length(fn_pup) == 1 & length(fn_dff) == 1
        load(fullfile(fn_pup.folder,fn_pup.name));
        load(fullfile(fn_dff.folder,fn_dff.name));
        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-pupil');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        savepupilfigpath = fullfile(dataIndex.BehPath{ii},[date(1:6),'_figs-pupil']);
        if ~exist(savepupilfigpath,'dir')
            mkdir(savepupilfigpath);
        end
        
        cd(savepupilfigpath);
        
        
           % interpolate pupil signal,delete NaNs first
           pupilZ = pupil.dia(~isnan(pupil.dia));
           pupilT = pupil.t(~isnan(pupil.dia));
           
            pupilInt = interp1(pupilT,pupilZ,cells.t);
            % calculate cross-correlation
            dt = mean(diff(cells.t));
            maxTime = 20; % only consider 5s lags
        for tt=1:length(cells.dFF) 
            [c{tt},lags] = xcorr(pupilInt(~isnan(pupilInt)),cells.dFF{tt}(~isnan(pupilInt)),'normalized');
            newlag{tt} = lags*dt;
            stem(newlag,c);
        end
        savefilename = fullfile(savematpath,'xcorr-pupil-fluo.mat');
        save(savefilename,'c','newlag');
    end
end