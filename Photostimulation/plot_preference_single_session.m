function plot_preference_single_session(trialData)

% % plot_preference_single_session %
%PURPOSE:   Plot performance of a single session of preference paradigm
%AUTHORS:   HK Ortega 231003
%
%INPUT ARGUMENTS
%   trialData:  stats of the game including choice, outcome, and rule

%%
ugh = figure ('position', [0 0 2000 500]);
setup_figprop
nTrial = size(trialData.response,1);
hold on

    for r = 1:nTrial
        choiceColor = 'w';
        outcomeColor = 'w';
        outcomeLoc = 1.1;
        if trialData.response (r) == 2 %left choice
            choiceColor = 'r';
            choiceLoc = [.5,1];
            outcomeLoc = [.7/2, .3/2];
        elseif trialData.response(r) == 3 %right choice
            choiceColor = 'b';
            choiceLoc = [1,1.5];
            outcomeLoc = [3.3/2, 3.7/2];
        end 
        if trialData.rule(r) == 41 %left water + stim
            if trialData.outcome(r) == 5 %water + stim
                outcomeColor = [0.4940 0.1840 0.5560];
            elseif trialData.outcome(r) == 6 %water only
                outcomeColor = 'k';
            end 

        elseif trialData.rule(r)== 42 %right water + stim
            if trialData.outcome(r) == 6
                outcomeColor = [0.4940 0.1840 0.5560];
            elseif trialData.outcome(r) == 5
                outcomeColor = 'k';
            end
        end 
            plot([r,r],choiceLoc,choiceColor,'LineWidth',2)
            plot([r,r],outcomeLoc,'Color', outcomeColor,'LineWidth',2)
            xticks([0 100 200 300 400 500])
            xlim([0 nTrial])
            xlabel("Trials")
            ylim([0 2])
            yticks([.3/2 .5 1.5 3.7/2])
            set(gca, 'FontSize',20)
            H=gca;
            H.LineWidth=2;
            
            yticklabels(["Outcome", "Left", "Right", "Outcome"])
      
    end
end