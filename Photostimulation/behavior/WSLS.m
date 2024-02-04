function output=WSLS(trialData)
% % beh_performance %
%PURPOSE:   Calculate basic performance metrics for a session
%AUTHORS:   AC Kwan 191210
%
%INPUT ARGUMENTS
%   stats:  stats of the task
%
%OUTPUT ARGUMENTS
%   output:     numbers used to plot figure
%
% To plot the output, use plot_behperf().

%%
%remove misses

trialData.response = trialData.response(~(trialData.outcome == 77));
trialData.outcome  = trialData.outcome(~(trialData.outcome == 77));

ntrial=size(trialData.response,1);
output.left = 0;
output.right = 0;
output.miss = 0;
output.nTrial = 0;
output.leftWin = 0;
output.rightWin = 0;
output.leftLose = 0;
output.rightLose = 0;
output.LWinStay = 0;
output.LWinSwitch = 0;
output.LLoseStay = 0;
output.LLoseSwitch = 0;
output.RWinStay = 0;
output.RWinSwitch = 0;
output.RLoseStay = 0;
output.RLoseSwitch = 0;

%     OUTCOME.REWARDLEFT = 100;
%     OUTCOME.REWARDRIGHT = 111;
%     OUTCOME.NOREWARDLEFT = 101;
%     OUTCOME.NOREWARDRIGHT = 110;
%     OUTCOME.REWARDMANUAL = 10;
%     OUTCOME.MISS = 77;      %miss
%     RESP.LEFT=2;
%     RESP.RIGHT=3;

for k=1:ntrial
     output.nTrial = output.nTrial + 1; %add 1 to number of trials
     if trialData.response (k) == 2 
         output.left = output.left + 1; %add 1 to number of left choices
     end
     if trialData.response(k) == 3
         output.right = output.right +1;%add 1 to number of right choices
     end
     if trialData.response(k) == 0
         output.miss = output.miss + 1;
     end 
     if k == ntrial %catch for last trial where there is no next trial

         elseif trialData.response(k) == 2 && trialData.outcome(k) == 100 %left choice and rewarded
             output.leftWin = output.leftWin + 1;
             if trialData.response(k+1) == 2
                output.LWinStay = output.LWinStay + 1;
             elseif trialData.response(k+1) == 3
                output.LWinSwitch = output.LWinSwitch + 1;
             end
        elseif trialData.response(k) == 2 && trialData.outcome(k) == 110 %left choice and not rewarded
             output.leftLose = output.leftLose + 1;
             if trialData.response(k+1) == 2
                output.LLoseStay = output.LLoseStay + 1;
             elseif trialData.response(k+1) == 3
                output.LLoseSwitch = output.LLoseSwitch + 1;
             end
        elseif trialData.response(k) == 3 && trialData.outcome(k) == 111 %right choice and rewarded
             output.rightWin = output.rightWin + 1;
             if trialData.response(k+1) == 3
                output.RWinStay = output.RWinStay + 1;
             elseif trialData.response(k+1) == 2
                output.RWinSwitch = output.RWinSwitch + 1;
             end
        elseif trialData.response(k) == 3 && trialData.outcome(k) == 101 %right choice and not rewarded
             output.rightLose = output.rightLose + 1;
             if trialData.response(k+1) == 3
                output.RLoseStay = output.RLoseStay + 1;
             elseif trialData.response(k+1) == 2
                output.RLoseSwitch = output.RLoseSwitch + 1;
             end
      end
 end
end 


