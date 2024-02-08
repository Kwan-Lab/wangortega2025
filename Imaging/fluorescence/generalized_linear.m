function [ output ] = generalized_linear ( signal, t, event, eventTime, trialSubset, params)
% % linear_regr %
%% change the regression equation to Sul et al.2011
%PURPOSE:  random forest regressor
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   signal:         time-series signal (e.g., dF/F)
%   t:              time points corresponding to the signal
%   event:          table of different predictors (categorical or continuous)
%   eventTime:      the event times
%   trialSubset:    the subset of trials to investigate, all else set to NaN
%   params.window:  the time window around which to align signal
%   params.nback:   consider events from up to this # trials back


%OUTPUT ARGUMENTS
%   output:         structure containing regression analysis results
%
% To plot the output, use MP_plot_regr().

%% determine levels in Predictors
% replace undefined with inf to work with function unique()

%% interpolate signal

window=params.window;

% interpolate the signal to a finer time scale
interdt=0.01;
intert=[t(1):interdt:t(end)]';

% t might contain similar time point
intersig=interp1(t,signal,intert);

% align signal to the event
% use window slightly wider than the regression, so regression analysis
% won't run into the boundaries of this variable
[sigbyTrial, tbyTrial]=align_signal(intert,intersig,eventTime,[window(1)-1 window(end)+1]);

%% the factors in the linear regression
output.numPredictor=size(event,2);  %number of predictor
output.nback = params.nback;
output.interaction = params.interaction;
%% parameters of the regression
step_dur = nanmean(diff(window));
step_size = step_dur;       %if step size is step duration, then doing this in non-overlapping windows
output.regr_time=[window(1)+step_dur/2:step_size:window(end)-step_dur/2]';
        
%% random forest with moving window 
% to see which behavioral events explain signal variability

warning('off','MATLAB:singularMatrix');
warning('off','stats:pvaluedw:ExactUnavailable');

% train and test trials
numTrial = size(eventTime,1);

%which trials are picked for test/train
numTest = round(0.2*numTrial);
drawNum = randsample(numTrial,numTrial,'false'); %each time draw another set without replacement
testTrials = drawNum(1:numTest);
trainTrials = drawNum(numTest+1:end);

%take out the first and last events, so don't have to deal with boundary
%issues when identifying imaging frames
testTrials = testTrials(testTrials~=1 & testTrials~=numTrial);
trainTrials = trainTrials(trainTrials~=1 & trainTrials~=numTrial);


%identify the fit that has the smallest mean least-square error
%[~,idx]=min(mean(sqerr,2));

%using that specific lambda parameter, fit again using training +
%validation set
for jj=1:numel(output.regr_time)
    
    idx1=sum(tbyTrial<=output.regr_time(jj)-step_dur/2);     %signal corresponding to current time step
    idx2=sum(tbyTrial<(output.regr_time(jj)+step_dur/2));
    tempsig=squeeze(nanmean(sigbyTrial(idx1:idx2,:),1));
    
    
    % run multiple times to get a list of coefficients,
    % then use bootstrap to determine the CI 
    
    numIter = 100;
    coeff = zeros(numIter,size(event,2));
    for nn = 1:numIter
        MDl=fitrlinear(event(trainTrials,:)',tempsig(trainTrials)','ObservationsIn','columns','Learner','leastsquares','Solver','sgd','Regularization','ridge','Lambda',params.lambda,'FitBias',true);
        coeff(nn,:) = MDl.Beta;
    end
    
    % bootstrap to get mean and confidence interval
   greg.coeff = coeff;
   greg = getBootstrp(greg, 0, 0.001);
    % save the fitted coefficient
%      figure;plot(greg.coeff_bootave)
%       hold on;plot(greg.boothigh);
%        hold on;plot(greg.bootlow);
%        
        for kk=1:size(event,2)
            output.coeff(jj,kk)=greg.coeff_bootave(kk);        %coefficient
            output.pval(jj,kk)=signtest(greg.coeff(:,kk),0);   %pvalue for coefficient
            %output.rsquare(jj,kk) = stats.rsquare;
            %output.adjrsquare(jj,kk) = stats.adjrsquare;
%             if kk == 1
%                 output.SumSq(jj,kk) = tbl.SumSq(end);  % last one is error term
%             else
%                 output.SumSq(jj,kk) = tbl.SumSq(kk-1);
%             end
        end
%evaluate fit using test set
%     ytest_fit=predict(MDl,event(testIdx,:));
% 
%     output_fit.ytest = nan(size(tempsig));
%     output_fit.ytest(testIdx) = signal(testIdx);  %the test signal, measured
%     output_fit.ytest_fit = nan(size(signal));
%     output_fit.ytest_fit(testIdx) = ytest_fit;    %the test signal, predicted

end
warning('on','MATLAB:singularMatrix');
warning('on','stats:pvaluedw:ExactUnavailable');

end
