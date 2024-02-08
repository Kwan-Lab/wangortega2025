function xcorrData = getRegSelxCorr(reg_cr1,reg_cr2,var1,var2,Ind1,Ind2, thresh, savefluofigpath)

%% get cross correlation between different variables
label=['Cross-correlation of ',var1,' and ',var2];

% get the corresponding regression results
regCoeff1 = [];regPval1 = [];
regCoeff2 = [];regPval2 = [];
for ii = 1:length(reg_cr1.reg_cr)
    regCoeff1 = [regCoeff1;reg_cr1.reg_cr{ii}.coeff(:,Ind1)'];
    regPval1 = [regPval1;reg_cr1.reg_cr{ii}.pval(:,Ind1)'];
end

for ii = 1:length(reg_cr2.reg_cr)
    regCoeff2 = [regCoeff2;reg_cr2.reg_cr{ii}.coeff(:,Ind2)'];
    regPval2 = [regPval2;reg_cr2.reg_cr{ii}.pval(:,Ind2)'];
end

Mat1 = selectRg2D(reg_cr1.reg_cr,Ind1,[]);
Mat2 = selectRg2D(reg_cr2.reg_cr,Ind2,[]);
% take in the data, identify the significant time period

% variable 1
Mat1.sig = Mat1.pval;
Mat1.sig(Mat1.pval<thresh) = 1; Mat1.sig(Mat1.pval>=thresh) = 0;

% set not significant coefficient to 0
Mat1.sigCoeff = Mat1.coeff;
Mat1.sigCoeff(Mat1.pval>=thresh) = 0;

% variable 2
Mat2.sig = Mat2.pval;
Mat2.sig(Mat2.pval<thresh) = 1; Mat2.sig(Mat2.pval>=thresh) = 0;

% set not significant coefficient to 0
Mat2.sigCoeff = Mat2.coeff;
Mat2.sigCoeff(Mat2.pval>=thresh) = 0;


% calculate frame-by-frame correlation, get the Index with maximum
% correlation achieved
sigCorr = zeros(size(Mat1.sig,1)*2-1,size(Mat1.sig,2)*2-1,size(Mat1.sig,3)-1);
sigMaxInd = zeros(2,size(Mat1.sig,3)-1);
sigMaxDis = zeros(1,size(Mat1.sig,3)-1);  % find the coordinate of the maximum crosscorrelation, calculate its distance from the center
sigMaxValue = zeros(1,size(Mat1.sig,3)-1);

centerCoor = [size(Mat1.sig,1), size(Mat1.sig,2)];

for ii = 1:size(Mat1.sig,3)
    mat1 = Mat1.sig(:,:,ii); mat2 = Mat2.sig(:,:,ii);
    
    % check if there is NaN in the data
    if sum(sum(isnan(mat1)))>0 | sum(sum(isnan(mat2)))>0
        display(['Warning! NaN in regression',savefluofigpath])
        mat1(isnan(mat1))=0;
        mat2(isnan(mat2))=0;
    end

    if all(mat1(:)==0) | all(mat2(:)==0)
        sigCorr(:,:,ii) = NaN;
        sigMaxInd(:,ii)=[NaN,NaN];
        sigMaxDis(ii) = NaN;
        sigMaxValue(ii) = NaN;
    else
        temp = normxcorr2(mat1,mat2);
        sigCorr(:,:,ii) = temp;
        [ssr,snd] = max(temp(:));
        sigMaxValue(ii) = ssr;
        [ij,ji] = ind2sub(size(temp),snd);
        sigMaxInd(:,ii)=[ij,ji];
        sigMaxDis(ii) = sqrt((ij-centerCoor(1))^2+(ji-centerCoor(2))^2);
    end
end

%% use significant coefficient only
coeffCorr = zeros(size(Mat1.coeff,1)*2-1,size(Mat1.coeff,2)*2-1,size(Mat1.coeff,3)-1);
coeffMaxInd = zeros(2,size(Mat1.coeff,3)-1);
coeffMaxDis = zeros(1,size(Mat1.coeff,3)-1);  % find the coordinate of the maximum crosscorrelation, calculate its distance from the center
coeffMaxValue = zeros(1, size(Mat1.coeff,3)-1);

for ii = 1:size(Mat1.coeff,3)
    mat1 = Mat1.sigCoeff(:,:,ii); mat2 = Mat2.sigCoeff(:,:,ii);
     if sum(sum(isnan(mat1)))>0 | sum(sum(isnan(mat2)))>0
        display(['Warning! NaN in regression',savefluofigpath])
        mat1(isnan(mat1))=0;
        mat2(isnan(mat2))=0;
     end

      if sum(sum(isinf(mat1)))>0 | sum(sum(isinf(mat2)))>0
        display(['Warning! inf in regression',savefluofigpath])
        mat1(isinf(mat1))=0;
        mat2(isinf(mat2))=0;
      end

    if all(mat1(:)==0) | all(mat2(:)==0)
        coeffCorr(:,:,ii) = NaN;
         coeffMaxInd(:,ii)=[NaN,NaN];
        coeffMaxDis(ii) = NaN;
        coeffMaxValue(ii) = NaN;
    else
        
        temp = normxcorr2(mat1,mat2);
        coeffCorr(:,:,ii) = temp;
        [ssr,snd] = max(temp(:));
        coeffMaxValue(ii) = ssr;
        [ij,ji] = ind2sub(size(temp),snd);
        coeffMaxInd(:,ii)=[ij,ji];
        coeffMaxDis(ii) = sqrt((ij-centerCoor(1))^2+(ji-centerCoor(2))^2);
    end
end

t= -2.95:0.1:4.95;

% plot the cross-correlation value of significance and coefficient
figure;
h1=subplot(4,1,1);
plot(reg_cr1.reg_cr{1}.regr_time,sum(regPval1<thresh,1)/size(regPval1,1));
hold on;
plot(reg_cr2.reg_cr{1}.regr_time,sum(regPval2<thresh,1)/size(regPval2,1));
legend(var1,var2);
legend('box','off');
ylabel('Significant Grids');
title(['Fraction of significant grids (',label,')']);
set(gca,'box','off')
h1.Position = [0.13,0.8094,0.7750,0.1];

h2=subplot(4,1,2);
plot(t,squeeze(sigCorr(28,28,:)),'-r');
hold on;plot(t,squeeze(coeffCorr(28,28,:)),'-k');
set(gca,'box','off');
ylabel('Coefficient');
title(['Cross-correlation coefficient (',label,')']);
legend({'Significance','Coefficient'});
legend('box','off');
ylim([-0.2,1]);
h2.Position = [0.13,0.5903,0.7750,0.1];

h3=subplot(4,1,3);
plot(t,sigMaxDis,'-r');
hold on;plot(t,coeffMaxDis,'-k');
set(gca,'box','off');
ylabel('Offset');
title(['Offset of maximum correlation coefficient (',label,')']);
legend({'Significance','Coefficient'});
legend('box','off');
h3.Position = [0.13,0.38,0.7750,0.1];

% plot the maximum cross correlation coefficient;
h4=subplot(4,1,4);
plot(t,sigMaxValue,'-r');
hold on;plot(t,coeffMaxValue,'-k');
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('Frame-by-frame cross-correlation');
title(['Max cross-correlation coefficient (',label,')']);
legend({'Significance','Coefficient'});
legend('box','off');
ylim([-0.2,1]);
h4.Position = [0.13,0.17,0.7750,0.1];


print(gcf,'-dpng',fullfile(savefluofigpath,[label,' xcorrCoeff']));
saveas(gcf, fullfile(savefluofigpath,[label,' xcorrCoeff']), 'fig');% plot the distance from center point of maximum crosscorrelation
saveas(gcf, fullfile(savefluofigpath,[label,' xcorrCoeff']), 'svg');

% save the data into a structure
xcorrData.sigCorr = sigCorr;
xcorrData.sigMaxDis = sigMaxDis;
xcorrData.sigMaxInd = sigMaxInd;
xcorrData.sigMaxValue = sigMaxValue;
xcorrData.regPval1 = regPval1;
xcorrData.regCoeff1 = regCoeff1;
xcorrData.regPval2 = regPval2;
xcorrData.regCoeff2 = regCoeff2;
xcorrData.coeffCorr = coeffCorr;
xcorrData.coeffMaxDis = coeffMaxDis;
xcorrData.coeffMaxInd = coeffMaxInd;
xcorrData.coeffMaxValue = coeffMaxValue;

close all;
