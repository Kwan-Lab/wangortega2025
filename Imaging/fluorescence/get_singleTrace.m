function [ output ] = get_singleTrace ( signal, t, eventTime, psth_label, params)
% % get_psth %
%PURPOSE:   Given a time-series signal and event times, find average response signal
%           by binning
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   signal:         time-series signal (e.g., dF/F)
%   t:              time points corresponding to the signal
%   eventTime:      the event times
%   psth_label:     text saying what the event is about
%   params.window:   the time window around which to align signal
%   params.numBootstrapRepeat:  number of bootstrap repeats
%   params.CI:                  confidence interval
%
%OUTPUT ARGUMENTS
%   output:         structure containing peri-event time histogram
%
% To plot the output, use plot_psth().

%%

output.psth_label=psth_label;

window=params.window;
window_dt=nanmean(diff(window));
output.t=window(1:end-1)+window_dt/2; %use center of bin as the time associated with the average signal
output.t=output.t(:);

nCells = 1;
output.signal=nan(numel(output.t),nCells);

%% find average (by binning/histogram)
%find the segment of signal and its relative time, around each event
tempSig=[]; tempTime=[]; tempEvent=[];
nEvent = numel(eventTime);
valbyTrial = zeros(1,numel(eventTime));
for j=1:numel(eventTime)
    relTime=t-eventTime(j);     %time relative to the event
    
    tempSig=[tempSig; signal(relTime>=window(1) & relTime<=window(end))];
    tempTime=[tempTime; relTime(relTime>=window(1) & relTime<=window(end))];
    tempEvent=[tempEvent; j*ones(sum(relTime>=window(1) & relTime<=window(end)),1)];    %the event# associated with this signal segment

    %mean value for entire window for each trial
    output.valbyTrial(j)=nanmean(signal(relTime>=window(1) & relTime<=window(end)));
end

for j=1:numel(window)-1 %for each bin of time, what is the average signal?
    output.signal(j)=nanmean(tempSig(tempTime>=window(j) & tempTime<(window(j)+window_dt)));    
end


end