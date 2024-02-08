% check results of random forest with different trees

 rf_cr_50 = load("F:\GRAB_analysis\analysis\891\891_01132021\analysis-fluo\RF_CR_50.mat");
 rf_cr_100 = load("F:\GRAB_analysis\analysis\891\891_01132021\analysis-fluo\RF_CR_100.mat");
rf_cr_200 = load("F:\GRAB_analysis\analysis\891\891_01132021\analysis-fluo\RF_CR_200.mat");

% check coefficient

predTime = rf_cr_50.rf_cr{1}.regr_time; numPred = size(rf_cr_50.rf_cr{1}.predImp,2);
numROIs = length(rf_cr_50.rf_cr);
% 50,100
corrCoeff1 = zeros(numROIs,numPred);
% 100,200
corrCoeff2 = zeros(numROIs,numPred);

for ii = 1:numROIs
    for jj = 1:numPred
        tempcorr = corrcoef(rf_cr_50.rf_cr{ii}.predImp(:,jj),rf_cr_100.rf_cr{ii}.predImp(:,jj));
        corrCoeff1(ii,jj) = tempcorr(1,2);
        tempcorr = corrcoef(rf_cr_100.rf_cr{ii}.predImp(:,jj),rf_cr_200.rf_cr{ii}.predImp(:,jj));
        corrCoeff2(ii,jj) = tempcorr(1,2);
    end
end

figure; 
for tt = 1:numPred
    subplot(2,5,tt)
    plot(corrCoeff1(:,tt));
    hold on; plot(corrCoeff2(:,tt));
end