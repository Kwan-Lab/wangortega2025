function MP_STIM_plot_all (input1, input2, tlabel, savepath)

%plot value of interest (yaxis) for each stimulation condition (x axis)

%ChAt x axis values (+1 to value pulled from logfiles)
%LM2 = 2, RM2 = 3, LV1 = 4, RV1 = 5, Ctrl = 1

%Dbh x axis values (+ 7  to value pulled from logfiles )
%LM2 = 9, RM2 = 10, LV1 = 11, RV1 = 12, Ctrl = 8

%get values for plotting and separate for color coding green (chat),
%magenta (dbh), black (ctrl)
y1 = cell2mat(struct2cell(input1));
%y1 = struct2cell(input1);
n1 = size(y1,2);
y1_test = y1(1:end-1,:);
y_ctrl = y1(end,:);
x1_test = [1*ones(1,n1);2*ones(1,n1);3*ones(1,n1);4*ones(1,n1)];
x_ctrl = 5*ones(1,n1);


y2 = cell2mat(struct2cell(input2));
%y2 = struct2cell(input2);
n2 = size(y2,2);
y2_test = y2(1:end-1,:);
y_ctrl = horzcat(y_ctrl, y2(end,:));
x2_test = [8*ones(1,n2);9*ones(1,n2);10*ones(1,n2);11*ones(1,n2)];
x_ctrl = horzcat(x_ctrl, 12*ones(1,n2));

figure

plot(x1_test,y1_test,'og','MarkerFaceColor','w', 'MarkerSize',3, 'LineStyle','none');
hold on;
plot(x2_test, y2_test,'om','MarkerFaceColor','w', 'MarkerSize',3,'LineStyle','none')
plot(x_ctrl, y_ctrl, 'ok', 'MarkerFaceColor','w', 'MarkerSize',3, 'LineStyle','none')
xlim([0 14])
xticks([0:14])
xlabel({'ChAT-cre x Ai32                  Dbh-cre x Ai32         '})
xticklabels({'', 'LM2', 'RM2', 'LV1', 'RV1', 'Ctrl','', '', 'LM2', 'RM2', 'LV1', 'RV1','Ctrl'})

h = gca;
h.XAxis.TickLength = [0,0];

if strcmp(tlabel, 'Number of Trials') == 1 || strcmp(tlabel, 'Number of Misses') == 1 || strcmp(tlabel, 'Number of Left Choices') == 1 || strcmp(tlabel, 'Number of Right Choices') == 1
    ylim auto
else 
    ylim([0 1])
end 

if n1 > 1

    y1_mean = mean(y1(1:4,:),2);
    y1_std = std(y1(1:4,:),'',2); 
    x_error = [1:4];
    errorbar(x_error,y1_mean,y1_std,'-go','MarkerSize',10,...
    'MarkerEdgeColor','g','MarkerFaceColor','g', 'LineStyle','none'); 
    
    ctrl_mean = mean(y1(end,:),2);
    ctrl_std = std(y1(end,:),'',2);
    x_ctrl = 5;
    errorbar(x_ctrl,ctrl_mean,ctrl_std,'-ko','MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','k', 'LineStyle','none'); 
    line([0, 5],[ctrl_mean,ctrl_mean], 'LineStyle', '--', 'Color', 'k')

end
if n2 >1 
    y2_mean = mean(y2(1:4,:),2);
    y2_std = std(y2(1:4,:),'',2); 
    x_error = [8:11];
    errorbar(x_error,y2_mean,y2_std,'-mo','MarkerSize',10,...
    'MarkerEdgeColor','m','MarkerFaceColor','m', 'LineStyle','none'); 
    
    ctrl_mean = mean(y2(end,:),2);
    ctrl_std = std(y2(end,:),'',2);
    x_ctrl = 12;
    errorbar(x_ctrl,ctrl_mean,ctrl_std,'-ko','MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','k', 'LineStyle','none'); 
    line([8, 12],[ctrl_mean,ctrl_mean], 'LineStyle', '--', 'Color', 'k')
end 

title(tlabel)

if ~exist(savepath)
    mkdir(savepath)
end
cd(savepath);
    print(gcf,'-dpng',tlabel);    %png format
    saveas(gcf, tlabel, 'fig');
end 
