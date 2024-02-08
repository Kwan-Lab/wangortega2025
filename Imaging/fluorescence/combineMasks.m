function sumMask = combineMasks(maskList)

% masks are stored in a list, combine them in a 28*28 matrix
sumMask = zeros(28,28); 
    for ii = 1:size(maskList,3)
        sumMask(maskList(:,:,ii)==1) = 1;
    end

