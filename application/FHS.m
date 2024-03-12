function VaR = FHS(returns,p)
T = size(returns,1);  
model = arima(AR=NaN,Distribution="t",Variance=egarch(1,1));
options = optimoptions(@fmincon,Display="off",Diagnostics="off", ...
    Algorithm="sqp",TolCon=1e-7);

fit = estimate(model,returns,Options=options);  % Fit the model
[residuals,variances] = infer(fit,returns);   % Infer residuals and variances
standardizedResiduals = residuals./sqrt(variances);
s = RandStream.getGlobalStream();
reset(s)

nTrials = 20000;                  % # of independent random trials
horizon = 1;                     % VaR forecast horizon

bootstrappedResiduals = standardizedResiduals(unidrnd(T,horizon,nTrials));
Y0 = returns(end);                            % Presample returns
Z0 = residuals(end)./sqrt(variances(end));    % Presample model standardized residuals
V0 = variances(end);                          % Presample variances

portfolioReturns = filter(fit,bootstrappedResiduals, ...
                         Y0=Y0,Z0=Z0,V0=V0);

%cumulativeReturns = sum(portfolioReturns);
VaR = quantile(portfolioReturns',p);