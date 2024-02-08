function [ output ] = bayesian_linear ( signal, t, event, eventTime, trialSubset, params)
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
for jj=1:numel(output.regr_time)
    idx1=sum(tbyTrial<=output.regr_time(jj)-step_dur/2);     %signal corresponding to current time step
    idx2=sum(tbyTrial<(output.regr_time(jj)+step_dur/2));
    tempsig=squeeze(nanmean(sigbyTrial(idx1:idx2,:),1));
    

        % bagged ensemble of regression tree
        p = size(event,2);
        prior = bayeslm(p);
        postMdl = estimate(prior,event,tempsig');
        plot(postMdl);
        
        for kk=1:(size(event,2)+1)
            output.coeff(jj,kk)=stats.beta(kk);        %coefficient
            output.pval(jj,kk)=stats.p(kk);   %pvalue for coefficient
            
        end
        
end
warning('on','MATLAB:singularMatrix');
warning('on','stats:pvaluedw:ExactUnavailable');

end
