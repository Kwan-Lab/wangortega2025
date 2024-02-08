function output = selectRg2D(input, varInd, recordingSite)

%% read in selectivity data (from linear regression) in cells
% convert to 2D matrix
% varInd:  which predictor in the regression is needed

edgelength = sqrt(numel(input));
output.coeff = zeros(edgelength,edgelength,length(input{1}.regr_time));
output.pval = zeros(edgelength,edgelength,length(input{1}.regr_time));


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
    
    output.coeff(Ind1, Ind2,:) = input{cc}.coeff(:,varInd);
    output.pval(Ind1,Ind2,:) = input{cc}.pval(:,varInd);
end
