function output = select2D(input,recordingSite)

%% read in selectivity data in cells
% convert to 2D matrix

edgelength = sqrt(numel(input));
output = zeros(edgelength,edgelength,length(input{1}.signal));


for cc = 1:numel(input)
    if mod(cc,edgelength) == 0
        Ind2 = edgelength;
    else
        Ind2 = mod(cc,(edgelength));
    end
    if mod(cc,edgelength) == 0
        Ind1 = cc/edgelength;
    else
        Ind1 = floor(cc/edgelength)+1;
    end
    
    % change left-right to contra-ipsi
    
    if ~isempty(recordingSite) % if recordingSite is empty, the input is outcome selectivity
        % thus no need to convert to ipsi/contra
        if strcmp(recordingSite, 'left')
            output(Ind1, Ind2,:) = -input{cc}.signal;
        else
            output(Ind1, Ind2,:) = input{cc}.signal;
        end
    else
        output(Ind1, Ind2,:) = input{cc}.signal;
    end
end
