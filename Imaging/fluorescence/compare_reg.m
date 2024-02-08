function  compare_reg(reg_cr_ori,reg_cr_boot)

t=reg_cr_ori.regr_time;

nPre = size(reg_cr_ori.coeff,2)

% get bootstrap average
nBoot = length(reg_cr_boot);
coeffBoot = zeros(size(reg_cr_ori.coeff,1),size(reg_cr_ori.coeff,2),nBoot);
for bb = 1:nBoot
        coeffBoot(:,:,bb) = reg_cr_boot{bb}.coeff;
end



figure;
for pp = 2:nPre
subplot(4,4,pp-1)
plot(t, reg_cr_ori.coeff(:,pp));
aveCoeff = squeeze(nanmedian(coeffBoot(:,pp,:),3));
hold on; plot(t,aveCoeff);
end
