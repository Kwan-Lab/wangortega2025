function MP_GRAB_simpleplots_average(dataIndex)

% load behavior and dff data
% make some simple plots

nFiles = size(dataIndex,1);

for ii = 1:nFiles

    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));

    savefigpath = fullfile(fn_beh.folder,'figs-fluo','PSTH');
    if ~exist(savefigpath,'dir')
        mkdir(savefigpath);
    end



    load(fullfile(fn_beh.folder,fn_beh.name));
    % load dFF files
    load(fullfile(fn_beh.folder,'dff.mat'));


    %% Plot dF/F of all the cells
    %         MP_plot_allCells( cells, trialData ,savefigpath);

    close all;
    %% check correaltion between different cells
    %         autoCorr = zeros(length(cells.dFF));
    %         for uu = 1:length(cells.dFF)
    %             for vv = 1:length(cells.dFF)
    %                 autoCorr(uu,vv) = corr(cells.dFF{uu}, cells.dFF{vv});
    %             end
    %         end
    %         figure;
    %         heatmap(autoCorr);
    %         print(gcf,'-dpng',fullfile(savefigpath,'corr-cells'));    %png format
    %         saveas(gcf, fullfile(savefigpath,'corr-cells'), 'fig');

    %% Plot cue-aligned dF/F for each cell
    % Fig. 3d in paper came from 140605 data set, cell 8 10 37 74
    params=[];
    params.trigTime = trialData.cueTimes;
    params.xtitle = 'Time from cue (s)';
    params.window = [-3:0.1:5];
    params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
    params.CI = 0.95;  %confidence interval
    params.minNumTrial = 50;
    psth_output = cell(1,2);
    
    % average across ROIs
    for k=1:3
        fieldname=[];
        if k==1 %panel 1
            fieldname{1}={'reward','contra'};
            fieldname{2}={'noreward','contra'};
            fName{1} = fieldname;
        elseif k==2
            fieldname{1}={'reward','ipsi'};
            fieldname{2}={'noreward','ipsi'};
            fName{2} = fieldname;
        elseif k==3
            fieldname{1} = {'miss early'};
            fieldname{2} = {'miss late'};
            fName{3} = fieldname;
        end
        
            for ff = 1:length(fieldname)
                psth_output{k}{ff} = [];
            end
        
        tic

        % for NE: j = 710
        % for ACh: j = 33
        %j=710;  % N
         
        for kk=1:numel(fieldname)
            trialMask = getMask_GRAB(trials,fieldname{kk},dataIndex.RecordingSite{ii});
            results = zeros(numel(cells.dFF),length(params.window)-1);
            for j=1:numel(cells.dFF)
                psth_panel= get_singleTrace( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
                 results(j,:) = psth_panel.signal;
            end
            psth_output{k}{kk} = results;
        end
        toc
    end


save(fullfile(savefigpath,'PSTH.mat'),'psth_output');
plot_psth_average(psth_output,fName,params);
print(gcf,'-dpng',fullfile(savefigpath,'average PSTH'));
saveas(gcf, fullfile(savefigpath,'average PSTH'), 'fig');
saveas(gcf, fullfile(savefigpath,'average PSTH'),'svg');
close;

end
end


