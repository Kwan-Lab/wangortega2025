function [Mask] = getSigMask(sigMask, sigThresh)

%% use mergeMask to find areas with relatively constant pos/neg selectivity
%% need to get 1/0 selectivity later

area =  mean(sigMask,3);
BW1 = (area<sigThresh); %Get logical mask of pixels exceeding threshold
BW1 = ~bwareaopen(~BW1,10,4); %Remove small holes from pixel mask
AMask1 = bwareafilt(BW1,6,4);
L1 = bwlabel(AMask1);
%function to get separate masks of objects in L
Mask = getSepMask(L1,10);

