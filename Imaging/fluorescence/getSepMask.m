function mask = getSepMask(L,minpx)

% input: L output of matlab function bwlabel
%        minpx: minimum area to keep
numObjs = max(L(:));

mask = [];
for ii = 1:numObjs
    tempMask = zeros(size(L,1),size(L,2));
    tempMask(L==ii) = 1;
    if sum(tempMask(:)) > minpx % only keep areas that larger than 10 px
        mask = cat(3, mask,tempMask);
    end
end
