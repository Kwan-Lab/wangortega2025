function cells = MP_GRAB_normalizedFF(cells, trialData)

% input: cells structure after computing df/f
% ouput: same cells structure with normalized df/f and timemask
% normalize df/f based on average

% find consistant change point in the end

for ii = 1:numel(cells.dFF)
    %
    dffMask = zeros(1,length(cells.t));
    try
        [ipt,~] = findchangepts(cells.dFF{ii},'Statistic','mean','MaxNumChanges',200);
        % if the change point is close to cue onset, keep it . otherwise
        % discard.
        iptKept = [];
        trialInd = [];
        for cc = 1:numel(ipt)
            if min(abs(trialData.cueTimes-cells.t(ipt(cc))))< 0.3
                iptKept = [iptKept,ipt(cc)];
                [~,I] = min(abs(trialData.cueTimes-cells.t(ipt(cc))));
                trialInd = [trialInd,I];
            end
        end
%         figure;plot(cells.t,cells.dFF{ii});
%         yplot = ones(1,numel(iptKept))*0.3;

%         for tt = 1:length(iptKept)
%             if trialData.response(trialInd(tt)) == 2
%                 hold on;
%                 scatter(cells.t(iptKept(tt)),0.3,400,'r.');
%             elseif trialData.response(trialInd(tt)) == 3
%                 hold on;
%                 scatter(cells.t(iptKept(tt)),0.3,400,'b.');
%             end
%         end
        % check average baseline between changepoints
        aveBaseBC = zeros(1, length(iptKept)+1);
        for aa = 1:length(aveBaseBC)
            if aa == 1
                aveBaseBC(aa) = nanmean(cells.dFF{ii}(1:iptKept(aa)));
            elseif aa == length(aveBaseBC)
                aveBaseBC(aa) = nanmean(cells.dFF{ii}(iptKept(aa-1):end));
            else
                aveBaseBC(aa) = nanmean(cells.dFF{ii}(iptKept(aa-1):iptKept(aa)));
            end
        end
        sortAveBase = sort(aveBaseBC);
        [ipt,~] = findchangepts(sortAveBase,'Statistic','mean','MaxNumChanges',1);
        lowmean = nanmean(sortAveBase(1:ipt-1));
        highmean = nanmean(sortAveBase(ipt:end));

        % mark the timepoint, 0 for low baseline, 1 for high baseline

        for aa = 1:length(aveBaseBC)
            if abs(aveBaseBC(aa)-highmean) < abs(aveBaseBC(aa)-lowmean)
                % high baseline
                if aa == 1
                    dffMask(1:iptKept(aa)) = 1;
                elseif aa == length(aveBaseBC)
                    dffMask(iptKept(aa-1):end) = 1;
                else
                    dffMask(iptKept(aa-1):iptKept(aa)) = 1;
                end

            end
        end
        cells.dffMask{ii} = dffMask;
        cells.normdFF{ii} = cells.dFF{ii};
        cells.normdFF{ii}(cells.dffMask{ii}==1) = cells.normdFF{ii}(cells.dffMask{ii}==1)-(highmean-lowmean);

        % plot distributions of average baseline
%         figure;
%         edge = [0:0.01:0.3];
%         histogram(aveBaseBC,edge);
%         hold on; plot([lowmean lowmean],[0 100]);
%         hold on;plot([highmean highmean],[0 100]);
%         %
%         close all
    catch ME
        warning(ME.identifier, ME.message);
        cells.dffMask{ii} = dffMask;
        cells.normdFF{ii} = cells.dFF{ii};
        continue
    end

    %      figure;plot(cells.t,cells.dFF{ii});
    %     yplot = ones(1,numel(ipt))*0.3;
    %     hold on;scatter(cells.t(ipt),yplot,400,'k.');

end

