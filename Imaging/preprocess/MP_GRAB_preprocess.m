function MP_GRAB_preprocess(dataIndex)

% simple plots, check alignment 

nFiles = size(dataIndex,1);

for ii = 1:nFiles
    
    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh.mat']));
    load(fullfile(fn_beh.folder,fn_beh.name));
    
    % check if dff is already computed
    fn_dff = dir(fullfile(fn_beh.folder, 'dff.mat'));
    if isempty(fn_dff) 
        % load the fluorescence data extracted using cellROI, calculate dF/F
        stackPath = fullfile(dataIndex.LogFilePath{ii},'stack_info.mat');
        if exist(stackPath)
        stackInfo = load(stackPath);
        
        % get trigger delay times
        stackInfo.trigDelay = zeros(1, length(stackInfo.nFrames));
        for uu = 1:length(stackInfo.nFrames)
            stackInfo.trigDelay(uu) = stackInfo.scim_header(uu).acqTriggerTimestamps_sec(1);
        end
        frameTimeStamp = stackInfo.scim_header(1).frameTimestamps_sec;
        stackInfo.frameRate = 1/mean(diff(frameTimeStamp));
        stackInfo.roiDir = dataIndex.LogFilePath{ii};
        [ cells ] = calc_dFF( stackInfo, trialData.cueTimes );
        
        % check to see that imaging and behavior timing is sync'ed properly
        savefluofigpath = fullfile(fn_beh.folder,'figs_analysis');
        if ~exist(savefluofigpath)
            mkdir(savefluofigpath)
        end
        stackInfo.savFile_name = 'time align';
        check_imageFrameTimes ( cells, stackInfo, savefluofigpath );
        check_trialTimes(trialData, stackInfo, savefluofigpath);
        
        savematpath = fn_beh.folder;
        save(fullfile(savematpath,'dff.mat'),...
            'cells', '-v7.3');
        
        close all;
    else
        display('dFF already computed');
        end
    end
    end

end
