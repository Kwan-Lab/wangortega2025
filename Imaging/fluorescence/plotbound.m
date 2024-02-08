function plotbound(mask)

%% plot the mask boundaries
if ~isempty(mask)
    for ii = 1:size(mask,3)
        [B,L]=bwboundaries(mask(:,:,ii),'noholes');
        hold on;
        plot(B{1}(:,2),B{1}(:,1),'r','LineWidth',1);
    end
end
