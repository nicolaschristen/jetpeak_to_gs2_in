%% Skewed Gaussian function
%
% Input : xpts -- points at which to evaluate the function
%         mpos -- position of the peak
%         width -- width of the peak
%         nrm -- overall scaling factor
%         skew -- skewness of the peak (< 0 => to the left,
%                 > 0 => to the right)
%
% Output: ypts -- value of the function at location xpts
%
function ypts = gauSkew(r, rmax, width, nrm, skew)

% Normal distribution
% Probability density function
nrmDist.pdf = @(x) 1/sqrt(2*pi)*exp(-1*x.^2/2);
% Cumulative distribution function
nrmDist.cdf = @(x) 0.5*(1+erf(x/sqrt(2)));

ypts =  nrm*2/width * nrmDist.pdf((xpts-mpos)/width) ...
    .* nrmDist.cdf(skew*(xpts-mpos)/width);

end
