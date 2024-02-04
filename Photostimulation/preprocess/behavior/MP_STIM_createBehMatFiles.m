function newDataIndex = MP_STIM_createBehMatFiles(dataIndex)
% % MP_createBehMatFiles %
%PURPOSE:   Analyze each logfile specified in dataIndex and save the
%           results in a behavioral .mat file, for matching pennies task
%AUTHORS:   H Atilgan and AC Kwan 191127
%
%INPUT ARGUMENTS
%   dataIndex:  a table of the data files
%
%OUTPUT ARGUMENTS
%   newDataIndex:  a table of the data files, now with the BehFileCreated =
%                  true for those files that are processed
%

%% Create logfile-info related table

nFile = size(dataIndex,1);

behIndex = table(...
    cell(nFile,1),...
    cell(nFile,1),...
    NaN(nFile,1),...
    NaN(nFile,1),...
    cell(nFile,1)...
    );

behIndex.Properties.VariableNames = {...
    'Animal',...
    'Experiment',... Name of the scenario file ran by NBS Presentation
    'DateNumber',... Date/Time e.g, 1806291321 = 2018 June 29 13:21
    'BehCreated',...  Has the behavioral .mat file been created
    'Strain'...  add genotype
    };

%% parse and plot the logfiles specified in dataIndex

disp(['-----------------------------------------------------------']);
disp(['--- Detected: ' int2str(nFile) ' behavioral logfiles.']);
disp(['-----------------------------------------------------------']);

for ii = 1:nFile

    % Is the analysis file already created?
    fn = dir(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh.mat']));
    if size(fn,1)>0
        behIndex.BehCreated(ii) = 1;
    end

    if isnan(behIndex.BehCreated(ii))
        disp(['Parsing file #' int2str(ii) '.']);
        disp(['    ' dataIndex.LogFileName{ii}]);

        % extract data for each trial from raw logfile
        [ logData ] = MP_parseLogfile(dataIndex.LogFilePath{ii}, dataIndex.LogFileName{ii});

        logfileData.Animal = logData.subject;
        logfileData.Experiment = logData.scenario;
        yr=num2str(logData.dateTime{1}(9:10));
        mo=num2str(logData.dateTime{1}(1:2));
        day=num2str(logData.dateTime{1}(4:5));
        hr=num2str(logData.dateTime{2}(1:2));
        min=num2str(logData.dateTime{2}(4:5));
        logfileData.DateNumber=str2num([yr mo day hr min]);


        [ sessionData, trialData ] = MP_getSessionData(logData);
        trials = MP_getTrialMasks(trialData);
        %save behavioral .mat file
        save(fullfile(dataIndex.BehPath{ii}, [dataIndex.LogFileName{ii}(1:end-4),'_beh'])',...
            'sessionData','trialData','logfileData', 'trials');

        behIndex.BehCreated(ii) = 1;

    else  %else if the behavioral .mat file already created
        disp(['Loading file #' int2str(ii) '.']);
        disp(['    ' dataIndex.LogFileName{ii}]);

        load(fullfile(dataIndex.BehPath{ii},[dataIndex.LogFileName{ii}(1:end-4),'_beh.mat']));
    end

    Ind = strfind(dataIndex.LogFilePath{ii},filesep);


    behIndex.Animal(ii) = logfileData.Animal;
    behIndex.Experiment(ii) = logfileData.Experiment;
    behIndex.DateNumber(ii) = logfileData.DateNumber;

    %%
    clearvars -except nFile i dataIndex behIndex Mol

end

%add strain info
    subList = {'22522';'23514';'23513';'23517';'23518';'23509';'22526';'23526';'23524';'23519';'23508';'23538';'23527';'23536';'23541';'23529';'23534'};
    Strain = {'Dbh';'ChAt';'ChAt';'ChAt';'ChAt';'Dbh';'Dbh';'ChAt';'ChAt';'ChAt';'Dbh';'Dbh';'ChAt';'Dbh';'Dbh';'ChAt';'Dbh'};
for ii = 1:nFile
     Index = strcmp(subList, behIndex.Animal{ii});
     behIndex.Strain{ii} = Strain{Index};
end

%% Add the logfile-extracted information into the database index

newDataIndex = [dataIndex behIndex];

end
