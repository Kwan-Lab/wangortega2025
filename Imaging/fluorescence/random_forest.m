function [ output ] = random_forest ( signal, t, event, eventTime, trialSubset, params)
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
% temp_event = event;
% for ii = size(event,2)
%     temp =temp_event.(ii);
%     if iscategorical(temp)
%         temp(isundefined(temp)) =  categorical(inf);
%     else % if continuous
%         temp(isnan(temp)) = inf;
%     end
%    temp_event.(ii) = temp;
% end
% countLevels = @(x)numel(categories(categorical(x)));
% numLevels = varfun(countLevels,temp_event,'OutputFormat','uniform');

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
        t = templateTree('NumVariablesToSample','all',...
            'PredictorSelection','interaction-curvature','Surrogate','on');
        rng(1); % For reproducibility
        warning('off','all')
        Mdl = fitrensemble(event,tempsig,'Method','Bag','NumLearningCycles',100, ...
            'Learners',t,'NumBins',params.bin);

        % estimate R^2
        yHat = oobPredict(Mdl);
        rsquare = corr(Mdl.Y,yHat)^2;
        
         output.rsquare(jj)=rsquare;        % r-square
        % Estimate predictor importance values by permuting out-of-bag observations among the trees.
        try
            impOOB = oobPermutedPredictorImportance(Mdl);
            for kk=1:size(event,2)
               
                output.predImp(jj,kk) = impOOB(kk);
                %output.pval(jj,kk)=stats.tstat.pval(kk);   %pvalue for coefficient
                %output.rsquare(jj,kk) = stats.rsquare;
                %output.adjrsquare(jj,kk) = stats.adjrsquare;
%                 if kk == 1
%                     output.SumSq(jj,kk) = tbl.SumSq(end);  % last one is error term
%                 else
%                     output.SumSq(jj,kk) = tbl.SumSq(kk-1);
%                 end
            end
        catch
            for kk=1:size(event,2)
                output.predImp(jj,kk) = NaN;
            end
        end
        try
            [impGain,predAssociation] = predictorImportance(Mdl);
            for kk=1:size(event,2)
                output.impGain(jj,kk) = impGain(kk);
            end
            output.predAssociation(:,:,jj) = predAssociation;
        catch
            for kk=1:size(event,2)
                output.impGain(jj,kk) = NaN;
            end
            output.predAssociation(:,:,jj) = NaN;
        end
end
warning('on','MATLAB:singularMatrix');
warning('on','stats:pvaluedw:ExactUnavailable');

end
