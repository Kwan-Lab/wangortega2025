function crossCorr = getRegTempCorr(reg_cr1,reg_cr2,var1,var2,Ind1,Ind2, thresh, savefluofigpath)

%% get cross correlation between different variables
timeSeriesLength = length(reg_cr1.reg_cr{1}.coeff(:,Ind1));
crossCorr = zeros(timeSeriesLength*2-1,2,length(reg_cr1.reg_cr));
label=['Cross-correlation of ',var1,' and ',var2];

% get the corresponding regression results
regCoeff1 = [];regPval1 = [];
regCoeff2 = [];regPval2 = [];
for ii = 1:length(reg_cr1.reg_cr)
    
    coeff1 = reg_cr1.reg_cr{ii}.coeff(:,Ind1);
    pval1 = reg_cr1.reg_cr{ii}.pval(:,Ind1);
    sigCoeff1 = coeff1; sigCoeff1(pval1>=thresh.alpha) = 0;
    
    coeff2 = reg_cr2.reg_cr{ii}.coeff(:,Ind2);
    pval2 = reg_cr2.reg_cr{ii}.pval(:,Ind2); 
    sigCoeff2 = coeff2; sigCoeff2(pval2>=thresh.alpha) = 0;

   
    if sum(isnan(sigCoeff1))>0

        sigCoeff1 = fillmissing(sigCoeff1,'linear');
    end
      if sum(isnan(sigCoeff2))>0
        sigCoeff2 = fillmissing(sigCoeff2,'linear');
    end
    [crossCorr(:,1,ii),crossCorr(:,2,ii)] = xcorr(sigCoeff1,sigCoeff2,'normalized');

end



