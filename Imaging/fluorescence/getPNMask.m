function [Mask1, Mask2] = getPNMask(mergeMask, Thresh)

%% use mergeMask to find areas with relatively constant pos/neg selectivity
%% need to get 1/0 selectivity later

area =  mean(mergeMask,3);
BW1 = (area<Thresh.neg); %Get logical mask of pixels exceeding threshold
BW1 = ~bwareaopen(~BW1,10,4); %Remove small holes from pixel mask
AMask1 = bwareafilt(BW1,6,4);
L1 = bwlabel(AMask1);
%function to get separate masks of objects in L
Mask1 = getSepMask(L1,10);

BW2 = (area>Thresh.pos); %Get logical mask of pixels exceeding threshold
BW2 = ~bwareaopen(~BW2,10,4); %Remove small holes from pixel mask
AMask2 = bwareafilt(BW2,6,4);
%CC = bwconncomp(AMask2);
L2 = bwlabel(AMask2);
Mask2 = getSepMask(L2,10);