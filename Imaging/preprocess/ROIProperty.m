% get the mean and variance of the grid ROI intensity
%% to do:
%% if grid_ROI, check if mean intensity is above 55% (might be changed later)
%% batch process

clearvars;
tic;
 
%% create directory structure for loading and saving files
saveroot = "F:\GRAB_analysis\analysis";
%default values
default_scim_ver = '5';
default_root_dir = 'F:\GRABiCorre\893-02162021';
%default_scim_ver = '5';
%default_data_dir = '/Users/alexkwan/Desktop/ongoing data analysis/ROI extraction/testdata_resonant/';
default_reg_subdir = 'mat';
default_grid_roi = '0';  
default_batch_process = '0';

%ask user
prompt = {'ScanImage verison (3 or 5):','Root directory (for saving):','Subdirectory with registered .tiff images (to be analyzed):','Grid ROIs?','Batch Process?'};
dlg_title = 'Load ROIs and extract signals from registered images';
num_lines = 1;
defaultans = {default_scim_ver,default_root_dir,default_reg_subdir, default_grid_roi, default_batch_process};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
scim_ver=str2double(answer{1});
root_dir=answer{2};
reg_subdir=answer{3};
grid_roi = str2num(answer{4});
batch_process = str2num(answer{5});


if batch_process
    temp = dir(root_dir); %Edit to specify data directories
    data_dirs = {temp.name};
    temp = ~(strcmp(data_dirs,'.') | strcmp(data_dirs,'..'));
    data_dirs = data_dirs(temp); %remove '.' and '..' from directory list
    disp('Directories for movement correction:');
    disp(data_dirs');
else
    data_dirs{1} = root_dir;
    disp('Directories for movement correction:');
    disp(data_dirs);
end    
    % get all sub directories
    
% h=msgbox('Select the directory with reference ROI .mat files');
% uiwait(h);
% save_dir=uigetdir(data_dir,'Select the directory with reference ROI .mat files');

%% load the ROI masks
for jj = 1:length(data_dirs)
    
    subject = data_dirs{jj}(1:3);
    savedir = [data_dirs{jj}(1:3),'_',data_dirs{jj}(5:12)];
    savematpath = fullfile(saveroot,subject,savedir);
    if batch_process
        curr_dir = fullfile(root_dir,data_dirs{jj});
    else
        curr_dir = root_dir;
    end
    % get the directory with ROI masks
    temp = dir(fullfile(curr_dir,'*ROI*'));
    save_dir = temp([temp.isdir]==1).name;
    
    cd(fullfile(curr_dir,save_dir));
    roifiles=dir('cell*.mat');  % there is another roiData.mat file in cellROI2.0
    
    cellmask=[];
    for k=1:numel(roifiles)
        load(roifiles(k).name);
        
        cellmask(:,:,k)=bw;
        if grid_roi  % no need to calculate neuropil signal for grid rois
            neuropilmask(:,:,k) = 0;
        else
            neuropilmask(:,:,k)=subtractmask;
        end
        %isred(k)=isRedCell;   % no red channels
    end
    clear bw subtractmask;
    
    
    %% load the mean projection, calculate the ROI intensity distribution
    cd(curr_dir);
    meanProjPath = dir('*Mean*.tif');
    meanProj = loadtiffseq([],meanProjPath.name);
    
    % get the mean intensity in ROI
    meanIntensity = zeros(1, size(cellmask,3));
    for uu = 1:size(cellmask,3)
        meanIntensity(uu) = mean(meanProj(cellmask(:,:,uu)==1));
    end
    savemeanpath = fullfile(savematpath,'ROIMeanIntensity.mat');
    save(savemeanpath, 'meanIntensity');
    %% get the variance of th rois from mat files
    
    % remove roi with lower intensity later
    
    mat_subdir = 'mat';
    
    if isdir(mat_subdir)
        cd(mat_subdir);
        
            stacks = dir('*.mat');
            f=cell(numel(roifiles),1); n=cell(numel(roifiles),1);
            for j=1:numel(stacks)
                disp(['Loading reg image file ' stacks(j).name]);
                cd(fullfile(curr_dir,mat_subdir));
                pic=load(stacks(j).name);
                [nY nX nZ]=size(pic.stack);
                
                parfor k=1:numel(roifiles)
                    tempf=[]; tempn=[];
                    for i=1:1:nZ
                        %get sum of pixel values within the ROI
                        tempf(i)=sum(sum(pic.stack(:,:,i).*uint16(cellmask(:,:,k))));
                        tempn(i)=sum(sum(pic.stack(:,:,i).*uint16(neuropilmask(:,:,k))));
                    end
                    if sum(sum(shifted_cellmask(:,:,k)))>0     %if there are pixels belonging the the ROI
                        if j==1 %if this is the first reg image, then start new variables
                            f{k}=tempf/sum(sum(cellmask(:,:,k)));   %per-pixel fluorescence
                            n{k}=tempn/sum(sum(shifted_neuropilmask(:,:,k)));   %per-pixel fluorescence
                        else
                            f{k}=[f{k} tempf/sum(sum(shifted_cellmask(:,:,k)))];   %per-pixel fluorescence
                            n{k}=[n{k} tempn/sum(sum(shifted_neuropilmask(:,:,k)))];   %per-pixel fluorescence
                        end
                    else %if the ROI is outside of the imaging field of view
                        f{k}=nan(size(tempf));
                        n{k}=nan(size(tempn));
                    end
                end
                clear pic;
            end
        
    else
        cd(reg_subdir);
        stacks=dir('*.tif');
        if ~exist('ROI','dir')
            f=[]; n=[];
            for j=1:numel(stacks)
                disp(['Loading reg image file ' stacks(j).name]);
                cd(fullfile(curr_dir,reg_subdir));
                warning('off','all');   %scim_openTif generates warning
                pic=loadtiffseq([],stacks(j).name);
                warning('on','all');   %scim_openTif generates warning
                [nY nX nZ]=size(pic);
                
                parfor k=1:numel(roifiles)
                    tempf=[]; tempn=[];
                    for i=1:1:nZ
                        %get sum of pixel values within the ROI
                        tempf(i)=sum(sum(pic(:,:,i).*uint16(shifted_cellmask(:,:,k))));
                        tempn(i)=sum(sum(pic(:,:,i).*uint16(shifted_neuropilmask(:,:,k))));
                    end
                    if sum(sum(shifted_cellmask(:,:,k)))>0     %if there are pixels belonging the the ROI
                        if j==1 %if this is the first reg image, then start new variables
                            f{k}=tempf/sum(sum(shifted_cellmask(:,:,k)));   %per-pixel fluorescence
                            n{k}=tempn/sum(sum(shifted_neuropilmask(:,:,k)));   %per-pixel fluorescence
                        else
                            f{k}=[f{k} tempf/sum(sum(shifted_cellmask(:,:,k)))];   %per-pixel fluorescence
                            n{k}=[n{k} tempn/sum(sum(shifted_neuropilmask(:,:,k)))];   %per-pixel fluorescence
                        end
                    else %if the ROI is outside of the imaging field of view
                        f{k}=nan(size(tempf));
                        n{k}=nan(size(tempn));
                    end
                end
                clear pic;
            end
        end
    end
    
    %% save the extracted signals
    cd(curr_dir);
    
    
    if isdir(mat_subdir)
        cd(mat_subdir);
    else
        cd(reg_subdir);
    end
    
    if ~exist('ROI','dir')
        [~,~,~]=mkdir('ROI');
        cd('ROI');
        for k=1:numel(roifiles)
            cellf=f{k};
            neuropilf=n{k};
            bw=shifted_cellmask(:,:,k);
            subtractmask=shifted_neuropilmask(:,:,k);
            %isRedCell=isred(k);
            
            temp = sprintf('%03d',k);
            save(strcat('cell',temp,'.mat'),'cellf','neuropilf','bw','subtractmask');
        end
    end
    disp(['Processed ' int2str(numel(roifiles)) ' ROIs --- All done!']);
    close all;
   
end
%
figure;
subplot(2,1,1);
plot(cellf);
axis tight; title('downsampled');
subplot(2,1,2);
plot(f{numel(roifiles)});
axis tight; title('extracted');


