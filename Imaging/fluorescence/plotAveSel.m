function plotAveSel(t,aveSel,color)

if ~isempty(aveSel)
    %if size(aveSel,1) ~= 1
        for ii = 1:size(aveSel,1)
            hold on;
            plot(t,aveSel(ii,:),'Color',color);
        end
%     else
%         hold on;
%         plot(t,aveSel,color)
%     end
end