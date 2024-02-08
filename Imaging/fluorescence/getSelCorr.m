function xcorrData = getSelCorr(sel,recordingSite, label,savefluofigpath)

Mat = select2D(sel,recordingSite);

% take in the data, identify the significant time period

coeffCorr = zeros(size(Mat,1)*2-1,size(Mat,2)*2-1,size(Mat,3)-1);
coeffMaxInd = zeros(2,size(Mat,3)-1);
coeffMaxDis = zeros(1,size(Mat,3)-1);  % find the coordinate of the maximum crosscorrelation, calculate its distance from the center
centerCoor = [size(Mat,1), size(Mat,2)];
findMaxLim = 10;  % find maximum value within +- 10 grids
for ii = 1:size(Mat,3)-1
    mat1 = Mat(:,:,ii); mat2 = Mat(:,:,ii+1);
    
    if all(mat1(:)==0) & all(mat2(:)==0)
        coeffCorr(:,:,ii) = NaN;
         coeffMaxInd(:,ii)=[NaN,NaN];
        coeffMaxDis(ii) = NaN;
    else
        temp = normxcorr2(Mat(:,:,ii),Mat(:,:,ii+1));
        coeffCorr(:,:,ii) = temp;
        maxTemp = temp(centerCoor(1)-findMaxLim:centerCoor(1)+findMaxLim,centerCoor(2)-findMaxLim:centerCoor(2)+findMaxLim);
        [ssr,snd] = max(maxTemp(:));
        [ij,ji] = ind2sub(size(maxTemp),snd);
        coeffMaxInd(:,ii)=[ij+centerCoor(1)-findMaxLim-1,ji+centerCoor(2)-findMaxLim-1];
        coeffMaxDis(ii) = sqrt((ij-findMaxLim-1)^2+(ji-findMaxLim-1)^2);
    end
end

t= -2.95:0.1:4.85;

% plot the cross-correlation value of significance and coefficient
figure;
plot(t,squeeze(coeffCorr(28,28,:)),'-k');
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('Frame-by-frame cross-correlation');
title(['Cross-correlation coefficient (',label,')']);
print(gcf,'-dpng',fullfile(savefluofigpath,[label,' cross-correaltion']));
saveas(gcf, fullfile(savefluofigpath,[label,' cross-correaltion']), 'fig');


% plot the distance from center point of maximum crosscorrelation
figure;
plot(t,coeffMaxDis,'-k');
set(gca,'box','off');
xlabel('Time from cue(s)');
ylabel('Offset of maximum correlation coefficient');
title(['Offset of maximum correlation coefficient (',label,')']);
print(gcf,'-dpng',fullfile(savefluofigpath,[label,' Offset of maximum correlation coefficient']));
saveas(gcf, fullfile(savefluofigpath,[label,' Offset of maximum correlation coefficient']), 'fig');

% save the data into a structure
xcorrData.coeffCorr = coeffCorr;
xcorrData.coeffMaxDis = coeffMaxDis;
xcorrData.coeffMaxInd = coeffMaxInd;

close all;