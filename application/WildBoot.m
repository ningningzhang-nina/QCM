function WB = WildBoot(x, nboot)
% Generate Wild bootstrap samples 
% Input:          x     (nx1) Sample data. e.g. regression residuals
%                 nboot (1x1) Number of bootstrap replication
% Output:           WB  (n x nboot) wild bootstrap samples
% version 1:      YH, 23 June 2014, yujia.hu@gmail.com
% reference: Mammen, E. (1993), Bootstrap and Wild Bootstrap for High 
%               Dimensional Linear Models, Annals of Statistics 21
    if nargin == 1
        nboot = length(x);
    elseif nargin ~= 2
        error('Provide at least one argument')
    end
    if nboot<1
        error('Number of bootstrap replication has to be greater then zero')
    end
    n   = length(x);
    p   = binornd(1,(sqrt(5)+1)/(2*sqrt(5)),n,nboot);
    z   = (1 - sqrt(5))/2 * p + (1 + sqrt(5))/2 * (1 - p);
    WB  = z.*repmat(x,1,nboot);
end

