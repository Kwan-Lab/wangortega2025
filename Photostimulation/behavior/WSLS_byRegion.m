function output=WSLS_byRegion(trialData, stats, regionInd, st)
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

%%itializie variables
output.pStay = 0;
output.pSwitch = 0;
output.pSwitch = 0;
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
output.rewNMinusOneStay = 0;
output.rewNMinusOneSwitch = 0;
output.noRewNMinusOneStay = 0;
output.noRewNMinusOneSwitch = 0;

%calculate miss-related variables first
ntrialMiss=size(trialData.response,1);
for k=1:ntrialMiss
  if regionInd ==  st(k) %if region of trial k is region of interest
      if trialData.response(k) == 0
         output.miss = output.miss + 1;
      end
      if k < ntrialMiss
           if trialData.response(k+1) == 0 
                output.missNextTrial = output.missNextTrial + 1;
           end

            if trialData.outcome(k) == 100 && trialData.response(k+1) == 0 %left choice and rewarded
                output.LWinMiss = output.LWinMiss +1;
            end 
            if trialData.outcome(k) == 101 && trialData.response(k+1) == 0  %left choice and not rewarded
                output.LLoseMiss = output.LLoseMiss + 1; 
            end
            if trialData.outcome(k) == 111 && trialData.response(k+1) == 0 %right choice and rewarded
                output.RWinMiss = output.RWinMiss + 1;
            end
            if trialData.outcome(k) == 110 && trialData.response(k+1) == 0 %right choice and not rewarded
                output.RLoseMiss = output.RLoseMiss + 1; 
            end
      end 
  end
end


%remove misses
trialData.response = trialData.response(~(trialData.outcome == 77));
stats.stRegion = stats.stRegion(~(trialData.outcome == 77));
trialData.outcome  = trialData.outcome(~(trialData.outcome == 77));

ntrial=size(trialData.response,1);

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

     if k > 1 %catch for starting after trial 1 (no n-1) 
         if trialData.outcome(k-1) == 100 %trial before stim left reward
             if trialData.outcome(k) == 100 || trialData.outcome(k) == 101 %trial of stim left
                 output.rewNMinusOneStay = output.rewNMinusOneStay + 1; 
             elseif trialData.outcome(k) == 110 || trialData.outcome(k) == 111 %trial of stim right
                 output.rewNMinusOneSwitch = output.rewNMinusOneSwitch + 1;
             end 
         elseif trialData.outcome(k-1) == 110 %trial before stim right no reward
             if trialData.outcome(k) == 100 || trialData.outcome(k) == 101 %trial of stim left
                 output.noRewNMinusOneSwitch = output.noRewNMinusOneSwitch + 1; 
             elseif trialData.outcome(k) == 110 || trialData.outcome(k) == 111 %trial of stim right
                 output.noRewNMinusOneStay = output.noRewNMinusOneStay + 1;
             end 
         elseif trialData.outcome(k-1) == 101 %trial before sitm left no reward
              if trialData.outcome(k) == 100 || trialData.outcome(k) == 101 %trial of stim left
                 output.noRewNMinusOneStay = output.noRewNMinusOneStay + 1; 
             elseif trialData.outcome(k) == 110 || trialData.outcome(k) == 111 %trial of stim right
                 output.noRewNMinusOneSwitch = output.noRewNMinusOneSwitch + 1;
              end 
         elseif trialData.outcome(k-1) == 111 %trial before sitm right reward
              if trialData.outcome(k) == 100 || trialData.outcome(k) == 101 %trial of stim left
                 output.rewNMinusOneSwitch = output.rewNMinusOneSwitch + 1; 
             elseif trialData.outcome(k) == 110 || trialData.outcome(k) == 111 %trial of stim right
                 output.rewNMinusOneStay = output.rewNMinusOneStay + 1;
             end 

         end 
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
        if trialData.response(k) == 3 && trialData.response(k+1) == 3
            output.pStay = output.pStay + 1;
        elseif trialData.response(k) == 2 && trialData.response(k+1) == 2
            output.pStay = output.pStay + 1;
        end

        if trialData.response(k) == 2 && trialData.response(k+1) == 3
            output.pSwitch = output.pSwitch + 1;
        elseif trialData.response(k) == 3 && trialData.response(k+1) == 2
            output.pSwitch = output.pSwitch + 1;
        end
        
        if trialData.outcome(k) == 100 %left choice and rewarded
             output.leftWin = output.leftWin + 1;
             if trialData.response(k+1) == 2
                output.LWinStay = output.LWinStay + 1;
                output.WSLS = output.WSLS + 1;
             elseif trialData.response(k+1) == 3
                output.LWinSwitch = output.LWinSwitch + 1;
             end
        end

        if trialData.outcome(k) == 101 %left choice and not rewarded
             output.leftLose = output.leftLose + 1;
             if trialData.response(k+1) == 2
                output.LLoseStay = output.LLoseStay + 1;
             elseif trialData.response(k+1) == 3
                output.LLoseSwitch = output.LLoseSwitch + 1;
                output.WSLS = output.WSLS + 1;
             end
        end

        if  trialData.outcome(k) == 111 %right choice and rewarded
             output.rightWin = output.rightWin + 1;
             if trialData.response(k+1) == 3
                output.RWinStay = output.RWinStay + 1;
                output.WSLS = output.WSLS + 1;
             elseif trialData.response(k+1) == 2
                output.RWinSwitch = output.RWinSwitch + 1;
             end
        end

        if trialData.outcome(k) == 110 %right choice and not rewarded
             output.rightLose = output.rightLose + 1;
             if trialData.response(k+1) == 3
                output.RLoseStay = output.RLoseStay + 1;
             elseif trialData.response(k+1) == 2
                output.RLoseSwitch = output.RLoseSwitch + 1;
                output.WSLS = output.WSLS + 1;
             end
        end
      end
   end
end 
% calculated stats w/ misses
output.nLWin = output.LWinSwitch + output.LWinStay;
output.nRWin = output.RWinSwitch + output.RWinStay;
output.nWin = output.nLWin + output.nRWin;
output.nLLose = output.LLoseSwitch + output.LLoseStay;
output.nRLose = output.RLoseSwitch + output.RLoseStay;
output.nLose = output.nLLose + output.nRLose;
end


