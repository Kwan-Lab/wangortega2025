function output=WSLS_byRegion_includeMisses(trialData, stats, regionInd)
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

% trialData.response = trialData.response(~(trialData.outcome == 77));
% stats.stRegion = stats.stRegion(~(trialData.outcome == 77));
% trialData.outcome  = trialData.outcome(~(trialData.outcome == 77));

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
output.LWinMiss = 0;
output.LLoseStay = 0;
output.LLoseSwitch = 0;
output.LLoseMiss = 0;
output.RWinStay = 0;
output.RWinSwitch = 0;
output.RWinMiss = 0;
output.RLoseStay = 0;
output.RLoseSwitch = 0;
output.RLoseMiss = 0;
output.rewardNextTrial = 0; 
output.missNextTrial = 0; 
output.switchCurrentTrial = 0;
output.stayCurrentTrial = 0;
output.WSLS = 0; 

%     OUTCOME.REWARDLEFT = 100;
%     OUTCOME.REWARDRIGHT = 111;
%     OUTCOME.NOREWARDLEFT = 101;
%     OUTCOME.NOREWARDRIGHT = 110;
%     OUTCOME.REWARDMANUAL = 10;
%     OUTCOME.MISS = 77;      %miss
%     RESP.LEFT=2;
%     RESP.RIGHT=3;

for k=1:ntrial
  if regionInd ==  stats.stRegion(k) %if region of trial k is region of interest
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

     if k > 1 %catch for starting after trial 1 (no n-1) 
         if trialData.response(k) == 3 && trialData.response(k-1) == 3
             output.stayCurrentTrial = output.stayCurrentTrial + 1;
         end 
         if trialData.response(k) == 2 && trialData.response(k-1) == 2
             output.stayCurrentTrial = output.stayCurrentTrial + 1;
         end 
         if trialData.response(k) == 2 && trialData.response(k-1) == 3
             output.switchCurrentTrial = output.switchCurrentTrial + 1;
         end 
         if trialData.response(k) == 3 && trialData.response(k-1) == 2
             output.switchCurrentTrial = output.switchCurrentTrial + 1;
         end 
     end 
     if k < ntrial %catch for last trial where there is no next trial
        if trialData.outcome(k+1) == 100 || trialData.outcome(k+1) == 111
         output.rewardNextTrial = output.rewardNextTrial + 1;
        end

        if trialData.response(k+1) == 0 
         output.missNextTrial = output.missNextTrial + 1;
        end
        
        if trialData.response(k) == 2 && trialData.outcome(k) == 100 %left choice and rewarded
             output.leftWin = output.leftWin + 1;
             if trialData.response(k+1) == 2
                output.LWinStay = output.LWinStay + 1;
                output.WSLS = output.WSLS + 1;
             elseif trialData.response(k+1) == 3
                output.LWinSwitch = output.LWinSwitch + 1;
             elseif trialData.response(k+1) == 0
                 output.LWinMiss = output.LWinMiss + 1;
             end
        end

        if trialData.response(k) == 2 && trialData.outcome(k) == 110 %left choice and not rewarded
             output.leftLose = output.leftLose + 1;
             if trialData.response(k+1) == 2
                output.LLoseStay = output.LLoseStay + 1;
             elseif trialData.response(k+1) == 3
                output.LLoseSwitch = output.LLoseSwitch + 1;
                output.WSLS = output.WSLS + 1;
             elseif trialData.response(k+1) == 0
                 output.LLoseMiss = output.LLoseMiss + 1;
             end
        end

        if trialData.response(k) == 3 && trialData.outcome(k) == 111 %right choice and rewarded
             output.rightWin = output.rightWin + 1;
             if trialData.response(k+1) == 3
                output.RWinStay = output.RWinStay + 1;
                output.WSLS = output.WSLS + 1;
             elseif trialData.response(k+1) == 2
                output.RWinSwitch = output.RWinSwitch + 1;
             elseif trialData.response(k+1) == 0
                output.RWinMiss = output.RWinMiss +1; 
             end
        end

        if trialData.response(k) == 3 && trialData.outcome(k) == 101 %right choice and not rewarded
             output.rightLose = output.rightLose + 1;
             if trialData.response(k+1) == 3
                output.RLoseStay = output.RLoseStay + 1;
             elseif trialData.response(k+1) == 2
                output.RLoseSwitch = output.RLoseSwitch + 1;
                output.WSLS = output.WSLS + 1;
             elseif trialData.response(k+1) == 0
                 output.RLoseMiss = output.RLoseMiss + 1;
             end
        end
      end
   end
end 
% calculated stats
output.pStay = (output.LWinStay + output.LLoseStay + output.RWinStay + output.RLoseStay) / (output.LWinSwitch + output.LLoseSwitch + output.RWinSwitch + output.RLoseSwitch + output.LWinStay + output.LLoseStay + output.RWinStay + output.RLoseStay);
output.pSwitch = (output.LWinSwitch + output.LLoseSwitch + output.RWinSwitch + output.RLoseSwitch) / (output.LWinSwitch + output.LLoseSwitch + output.RWinSwitch + output.RLoseSwitch + output.LWinStay + output.LLoseStay + output.RWinStay + output.RLoseStay);

output.pLWinSwitch = output.LWinSwitch / (output.LWinSwitch + output.LWinStay + output.LWinMiss);
output.pLWinStay = output.LWinStay / (output.LWinSwitch + output.LWinStay + output.LWinMiss);
output.pLWinMiss = output.LWinMiss / (output.LWinSwitch + output.LWinStay + output.LWinMiss);
output.pRWinSwitch = output.RWinSwitch / (output.RWinSwitch + output.RWinStay + output.RWinMiss);
output.pRWinStay = output.RWinStay / (output.RWinSwitch + output.RWinStay + output.RWinMiss);
output.pRWinMiss = output.RWinMiss / (output.RWinSwitch + output.RWinStay + output.RWinMiss);


output.pLLoseSwitch = output.LLoseSwitch / (output.LLoseSwitch + output.LLoseStay + output.LLoseMiss);
output.pLLoseStay = output.LLoseStay / (output.LLoseSwitch + output.LLoseStay + output.LLoseMiss);
output.pLLoseMiss = output.LLoseMiss / (output.LLoseSwitch + output.LLoseStay + output.LLoseMiss);
output.pRLoseSwitch = output.RLoseSwitch / (output.RLoseSwitch + output.RLoseStay + output.RLoseMiss);
output.pRLoseStay = output.RLoseStay / (output.RLoseSwitch + output.RLoseStay + output.RLoseMiss);
output.pRLoseMiss = output.RLoseMiss / (output.RLoseSwitch + output.RLoseStay + output.RLoseMiss);


output.pWinStay = (output.LWinStay + output.RWinStay) / (output.LWinSwitch + output.LWinStay + output.RWinSwitch + output.RWinStay + output.LWinMiss + output.RWinMiss);
output.pLoseStay = (output.LLoseStay + output.RLoseStay) / (output.LLoseSwitch + output.LLoseStay + output.RLoseSwitch + output.RLoseStay + output.LLoseMiss + output.RLoseMiss);
output.pWinSwitch = (output.LWinSwitch + output.RWinSwitch) / (output.LWinSwitch + output.LWinStay + output.RWinSwitch + output.RWinStay + output.LWinMiss + output.RWinMiss);
output.pLoseSwitch = (output.LLoseSwitch + output.RLoseSwitch) / (output.LLoseSwitch + output.LLoseStay + output.RLoseSwitch + output.RLoseStay + output.LLoseMiss + output.RLoseMiss);
output.pWinMiss = (output.LWinMiss + output.RWinMiss) / (output.LWinSwitch + output.LWinStay + output.RWinSwitch + output.RWinStay + output.LWinMiss + output.RWinMiss);
output.pLoseMiss = (output.LLoseMiss + output.RLoseMiss) / (output.LLoseSwitch + output.LLoseStay + output.RLoseSwitch + output.RLoseStay + output.LLoseMiss + output.RLoseMiss);


output.pMissNextTrial = output.missNextTrial / output.nTrial; 
output.pRewardNextTrial = output.rewardNextTrial / output.nTrial; 

output.pStayCurrentTrial = output.stayCurrentTrial / output.nTrial; 
output.pSwitchCurrentTrial = output.switchCurrentTrial / output.nTrial; 

output.pWSLS = output.WSLS / output.nTrial; 
end


