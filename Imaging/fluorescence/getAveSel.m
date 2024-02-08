function aveSel = getAveSel(selMat, mask)

%% go through every masks, calculate corresponding average selecvitity as a function of time
aveSel = [];
if ~isempty(mask)
    for ii = 1:size(mask,3)
        tempMask = logical(mask(:,:,ii));
        tempSel = zeros(size(selMat,3),1);
        for tt = 1:size(selMat,3)
            tempselt = selMat(:,:,tt);
            ttempsel = tempselt(tempMask);
            tempSel(tt) = mean(ttempsel(:));
        end
        aveSel = [aveSel;tempSel'];
    end
end
