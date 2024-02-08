function [ output ] = generalized_linear_lambda ( signal, t, event, eventTime, trialSubset, params,savefluopath)
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

testLambda = [0.05:0.05:0.6];
%testLambda=10.^(-(5:-0.1:0.1));
numRepeat = 5; %testCC(ll,tt,jj)
 testCC=zeros(numel(testLambda),numRepeat,numel(output.regr_time));
for jj=1:numel(output.regr_time)
    %display(jj);
    idx1=sum(tbyTrial<=output.regr_time(jj)-step_dur/2);     %signal corresponding to current time step
    idx2=sum(tbyTrial<(output.regr_time(jj)+step_dur/2));
    tempsig=squeeze(nanmean(sigbyTrial(idx1:idx2,:),1));
    
    for tt = 1:numRepeat  %do 5 repeats for each lambda value
        beta = zeros(14,numel(testLambda));
        %among the train trials, use 80% to train and 20% to test the lambda
        numTrainTrial = numel(trainTrials);
        numVal = round(0.2*numTrainTrial);
        
        drawNum = randsample(numTrainTrial,numTrainTrial,'false'); %each time draw another set without replacement
        
        trainvalTrials = trainTrials(drawNum(1:numVal));
        traintrainTrials = trainTrials(drawNum(numVal+1:end));
        
        %     traintrainIdx=false(1,size(signal,1));
        %     for j=1:numel(traintrainTrials)
        %         idx1=sum(eventTime(traintrainTrials(j),1)>=t);
        %         idx2=sum(eventTime(traintrainTrials(j),2)>=t);
        %         traintrainIdx(idx1:idx2)=true;
        %     end
        %     trainvalIdx=false(1,size(signal,1));
        %     for j=1:numel(trainvalTrials)
        %         idx1=sum(trigTime(trainvalTrials(j),1)>=t);
        %         idx2=sum(trigTime(trainvalTrials(j),2)>=t);
        %         trainvalIdx(idx1:idx2)=true;
        %     end
        %
        for ll=1:numel(testLambda) %test a range of regularization parameter lambda
            %L2-penalized linear least square regressoin (ridge regression)
            
            
            MDl=fitrlinear(event(traintrainTrials,:)',tempsig(traintrainTrials)','ObservationsIn','columns','Learner','leastsquares','Solver','sgd','Regularization','ridge','Lambda',testLambda(ll),'FitBias',true);
            yfit=predict(MDl,event(trainvalTrials,:));
            
            %computer corr. coefficient between predicted and measured activity
            tempCC=corrcoef(tempsig(trainvalTrials),yfit,'rows','pairwise');
            testCC(ll,tt,jj)=tempCC(1,2);
            %beta(:,ll) = MDl.Beta;
            %computer squared error between predicted and measured activity
            %note -- in my hands, this seems to identify the same index more consistently than CC
            %        sqerr(ll,jj) = sum((yfit - signal(trainvalIdx)).^2);
        end
        
    end
end


output=squeeze(mean(testCC,2));


warning('on','MATLAB:singularMatrix');
warning('on','stats:pvaluedw:ExactUnavailable');

end
