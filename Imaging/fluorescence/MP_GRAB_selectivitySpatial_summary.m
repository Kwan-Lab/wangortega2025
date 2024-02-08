function MP_GRAB_selectivitySpatial_summary(dataIndex, save_path_fluo);

% spatial properties of choice/outcome selectivity
%% probabily need more grids - pixel selectivity
colors=cbrewer('div','RdBu',256);
colors=flipud(colors);
colorRange = [-1 1];

nFiles = size(dataIndex,1);

subject_mask = [];
animalList = unique(dataIndex.Animal);

for aa = 1:numel(animalList)
    sessionInclude = [];
    for tt = 1:nFiles
        if dataIndex.Animal{tt}==animalList{aa}
            sessionInclude = [sessionInclude, tt];
        end
    end
    % initialize matrix
    posNum = zeros(length(sessionInclude),1);
    negNum = zeros(length(sessionInclude),1);
    posArea = [];
    negArea = [];

    for ii = 1:numel(sessionInclude)
        % load behavior files
        fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));
        load(fullfile(fn_beh.folder,fn_beh.name));
    
        % load fluorescent files
        fn_fluo = dir(fullfile(dataIndex.BehPath{ii},'dff.mat'));

        % make folders to save analysis and plots
        savematpath = fullfile(dataIndex.BehPath{ii},'analysis-fluo');
        if ~exist(savematpath,'dir')
            mkdir(savematpath);
        end
        saveSelName = fullfile(savematpath,'select_norm.mat');  
         savemaskpath = fullfile(savematpath,'choiceselMask.mat');
        %savePrevSelName = fullfile(savematpath,'prev_select_norm.mat'); 

        % load the areas according to animal and pos/neg selectivity
        if exist(savemaskpath) 
            load(savemaskpath);
            if isempty(Mask1)
                negNum(ii) =0;
            else
                negNum(ii) = size(Mask1,3);
                for nn = 1:size(Mask1,3)
                    negArea = [negArea, sum(sum(Mask1(:,:,nn)))];
                end
            end
            if isempty(Mask2)
                posNum(ii) =0;
            else
                posNum(ii) = size(Mask2,3);
                for nn = 1:size(Mask2,3)
                    posArea = [posArea, sum(sum(Mask2(:,:,nn)))];
                end
            end
        end
    end
    posNumTotal{aa} = posNum; negNumTotal{aa} = negNum;
    posAreaTotal{aa} = posArea; negAreaTotal{aa} = negArea;    
end
           
%% monte carlo simulation of
                % find the center of neg/pos areas
                % model the area by average density distribution, with
                % error
                % calculate the distribution of center distance
                % generate infinite plane by center distance and
                % distribution
               % poisson point process: number of points follows poisson
               % distribution


end