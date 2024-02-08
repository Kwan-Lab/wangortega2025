function [h, chi, p] = chi2ind(x, alpha)
    expected = (sum(x(1:end, :),2)*sum(x(:,1:end)))/sum(sum(x));
    OE = ((x - expected).^2)./expected;
    chi = sum(sum(OE,2));
    df = size(chi, 1)-1*size(chi,2);
    p = 1 - chi2cdf(chi,df);
    h = p>=alpha;
end