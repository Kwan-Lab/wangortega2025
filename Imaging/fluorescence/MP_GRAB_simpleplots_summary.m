function MP_GRAB_simpleplots_summary(dataIndex, savesumpath)

% load behavior and dff data
% make some simple plots

nFiles = size(dataIndex,1);
psth_output_sum = cell(1,3);

params=[];
    
    params.xtitle = 'Time from cue (s)';
    params.window = [-3:0.1:5];
    params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
    params.CI = 0.95;  %confidence interval
    params.minNumTrial = 50;

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
            fieldname{1}={'miss early'};
            fieldname{2}={'miss late'};
            fName{3} = fieldname;
        end
  end
for ii = 1:nFiles

    % load behavior files
    fn_beh = dir(fullfile(dataIndex.BehPath{ii},'beh_cut.mat'));

    savefigpath = fullfile(fn_beh.folder,'figs-fluo','PSTH');
    savematpath = fullfile(savefigpath,'PSTH.mat');



    load(savematpath);
    % load dFF file
    for k = 1:length(psth_output_sum)
        for kk = 1:length(psth_output{k})
            if ii == 1
                psth_output_sum{k}{kk} = [];
            end
            psth_output_sum{k}{kk} = [psth_output_sum{k}{kk};nanmean(psth_output{k}{kk},1)];
        end
    end

end
plot_psth_average(psth_output_sum,fName,params);


print(gcf,'-dpng',fullfile(savesumpath,'summary PSTH'));
saveas(gcf, fullfile(savesumpath,'summary PSTH'), 'fig');
saveas(gcf, fullfile(savesumpath,'summary PSTH'),'svg');
close;

end



