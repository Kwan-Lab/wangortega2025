function [ sessionData, trialData ] = audi_getSessionData( logData )
% % audi_getSessionData %
%PURPOSE: Retrieve session data for unexpected auditory response.
%AUTHORS: AC Kwan 170518
%
%INPUT ARGUMENTS
%   logdata:    Structure obtained with a call to MP_parseLogfile().
%   
%OUTPUT VARIABLES
%   sessionData:    Structure containing these fields:
%                   {subject, dateTime, nTrials, *lickTimes, *nSwitch}. 
%                   * lickTimes([1 2]):=[left right] lick times.
%   trialData:      Fields: 
%                   {startTimes, cueTimes, outcomeTimes, *cue, *response, *outcome}
%               *cue, response, and outcome for each trial contains
%               correspoinding eventCode from NBS Presentation

%COPY FROM LOGDATA
sessionData.subject = logData.subject;
sessionData.dateTime = logData.dateTime;

%SESSION DATA <<logData.header: 'Subject' 'Trial' 'Event Type' 'Code' 'Time'>>
TYPE = logData.values{3}; %Intersectional approach necessary, because values 2,3 were reused; 
CODE = logData.values{4}; %change in future Presentation scripts: with unique codes, only CODE would be needed to parse the logfile...

% 4:  start experiment (followed by 1 min spontaneous recording
% 5: start trial, 4kHz, 50 ms sound cue
% 100: 2 s intertrial interval 

% Get the event codes that are used in this log file, and then based on the set of codes used,
% make an educated guess on the correct event code set
tempidx=(strcmp(TYPE,'Nothing') | strcmp(TYPE,'Sound')); %do not consider RESPONSE or MANUAL
codeUsed = unique(CODE(tempidx));         %List the set of event codes found in this logfile


    SPON = 4; STIM = 5; OUTCOME = 100; END = 9;
    %interpulse is not used in later presentation files, 100 is actually an
    %outcome code


% Get rid of any CODE associated with the last, unfinished trial
lastCode=find(CODE==END,1,'last'); 
if ~isempty(lastCode)
    TYPE = TYPE(1:lastCode);
    CODE = CODE(1:lastCode);
end

% How many trials?

    sessionData.nTrials = sum(CODE==SPON | CODE == STIM);


%sessionData.nRules = numel(fieldnames(RULE));
%sessionData.rule_labels = fieldnames(RULE);

% Set up the time axis, and identify lick times
eventCodes = [SPON,STIM];  %looking for the first of any startexpt Events
time_0 = logData.values{5}(find(ismember(CODE,eventCodes),1,'first'));
time = logData.values{5}-time_0;   %time starts at first instance of startExpt
time = double(time)/10000;         %time as double in seconds

% Stimlus trials - which occurred and when?
cueCodes = [SPON,STIM]; %stimlus-associated codes as vector
trialData.cue =  CODE(ismember(CODE,cueCodes));
trialData.cueTimes = time(ismember(CODE,cueCodes));






end

