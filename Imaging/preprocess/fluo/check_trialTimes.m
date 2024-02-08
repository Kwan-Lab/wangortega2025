function check_trialTimes(trialData, stackInfo, savefluofigpath)


% plot the trial length taken from behavior logfile and imaging data
trialTimes_Beh = diff(trialData.cueTimes);

trialTimes_Img = zeros(1, length(stackInfo.nFrames)-1);
for tt = 1:min(length(trialData.cueTimes),length(stackInfo.nFrames))-1
    trialTimes_Img(tt) = stackInfo.scim_header(tt+1).frameTimestamps_sec(1)-stackInfo.scim_header(tt).frameTimestamps_sec(1);
end

maxlen = min(length(trialTimes_Beh), length(trialTimes_Img)); 
figure;
scatter(trialTimes_Beh(1:maxlen), trialTimes_Img(1:maxlen), 'black', 'filled')
hold on; plot([0 30], [0 30], 'lineWidth',1);
axis square;
xlabel('Trial time according to behavior .log file');
ylabel('Trial time according to imaging file');
title('Value should not deviate from diagonal');

savefilename = fullfile(savefluofigpath, 'check-trialtiming');
print(gcf,'-dpng',savefilename);    %png format
saveas(gcf, savefilename, 'fig');

% plot the difference between behavior and imaging trial timing
figure;
histogram(trialTimes_Beh(1:maxlen)' - trialTimes_Img(1:maxlen));
ylabel('Time difference between behavior and imaging data (s)');
xlabel('Number of trials');
savefilename = fullfile(savefluofigpath, 'trialTimingDif');
print(gcf,'-dpng',savefilename);    %png format
saveas(gcf, savefilename, 'fig');
