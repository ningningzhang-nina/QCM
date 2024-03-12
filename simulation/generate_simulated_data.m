function output=generate_simulated_data(alpha0,alpha1,beta,a0,a1,b1,lambda, lambda2, mu,n,ncut,plists,seeds)

%######################################################
%####### Simulate Diagonal MN-GARCH(1,1) model ########
%######################################################
sigma2 = zeros(2, n+ncut);
eps = zeros(1, n + ncut);
Y = zeros(1, n + ncut);

for t = 2:(n + ncut)
  sigma2(:,t) = alpha0 + alpha1 * eps(t-1)^2 + beta * sigma2(:,t-1);
  rng(seeds(t))
  components = binornd(1,lambda2) + 1;
  mus = mu;
  sds = sqrt(sigma2(:,t));
  eps(t) = normrnd(mus(components),sds(components));
  Y(t) = a0 + a1 * Y(t-1) + eps(t) + b1 * eps(t-1);
end

sigma2 = sigma2(:,(ncut):(n+ncut));
Y = Y((ncut):(n+ncut)); % save n+1 observations for later calculation
eps = eps((ncut):(n+ncut));

%######################################################
%############### Theoretical Moments ##################
%######################################################
mu_0 = zeros(1, n+1);
var_0 = zeros(1, n+1);
skew_0 = zeros(1, n+1);
kurt_0 = zeros(1, n+1);

for t = 2:(n+1)
  mu_0(t) = a0 + a1 * Y(t-1) + b1 * eps(t-1);
  var_0(t) = lambda(1) * (mu(1)^2 + sigma2(1,t)) + ...
      lambda(2) * (mu(2)^2 + sigma2(2,t)) - ...
      (lambda(1) * mu(1) + lambda(2) * mu(2))^2;
  
  skew_0(t) = (lambda(1) * (mu(1)^3 + 3 * sigma2(1,t) * mu(1)) + ...
                 lambda(2) * (mu(2)^3 + 3 * sigma2(2,t) * mu(2))) / ...
    (lambda(1) * (mu(1)^2 + sigma2(1,t)) + ...
       lambda(2) * (mu(2)^2 + sigma2(2,t)))^(3/2);
  
  kurt_0(t) = (lambda(1) * (mu(1)^4 + 6 * mu(1)^2 * sigma2(1,t) + 3 * sigma2(1,t)^2) + ...
                 lambda(2) * (mu(2)^4 + 6 * mu(2)^2 * sigma2(2,t) + 3 * sigma2(2,t)^2)) / ...
    (lambda(1) * (mu(1)^2 + sigma2(1,t)) + ...
       lambda(2) * (mu(2)^2 + sigma2(2,t)))^2;
end

mu_0 = mu_0(2:(n+1));
var_0 = var_0(2:(n+1));
skew_0 = skew_0(2:(n+1));
kurt_0 = kurt_0(2:(n+1));
Y = Y(2:(n+1));
eps = eps(2:(n+1));
sigma2 = sigma2(:,2:(n+1));

sigma_0 = sqrt(var_0);

%######################################################
%############### Theoretical Quantiles ################
%######################################################
% evaluate the function at the point x, where the components 
% of the mixture have weights w, means stored in u, and std deviations
% stored in s - all must have the same length.
Q = zeros(length(plists), n);
cnt = 0;
for p = plists
    cnt = cnt + 1;
    Q(cnt,:) = Q_theory(p, lambda, mu, sigma2, n, mu_0);
end

%######################################################
%#################### output ##########################
%######################################################
output.Y = Y;
output.mu_0 = mu_0;
output.sigma_0 = sigma_0;
output.skew_0 = skew_0;
output.kurt_0 = kurt_0;
output.Q = Q;
output.sigma2 = sigma2;
end
        
function y = cdf0(x,w,u,s) 
    y = w(1) * normcdf(x, u(1), s(1)) + w(2) * normcdf(x, u(2), s(2));
end
function y = cdf_inv(p,w,u,s)
    G = @(x) cdf0(x,w,u,s) - p; 
    y = fzero(G,0);
end
function Q = Q_theory(p, lambda, mu, sigma2, n, mu_0)
  Q = zeros(1, n);
  for t = 1:n
    Q(t) = cdf_inv(p, lambda, mu, sqrt(sigma2(:,t))) + mu_0(t);
  end
end