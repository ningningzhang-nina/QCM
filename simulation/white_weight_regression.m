function [w_theta_hat,estCov] = white_weight_regression(y,x_mat)
% y: n*1 array
% x_mat: n*p, p is the number of coefficients.

para0 = inv(x_mat' * x_mat) * x_mat' * y; %ols
res = y-x_mat*para0; %residuals
omega_hat = diag(res'.^2);   % white 
w = inv(omega_hat); %weighted matrix
w_theta_hat = inv(x_mat'*w*x_mat)*(x_mat'*w*y); % weighted ols
w_eta_hat = y-x_mat*w_theta_hat; %weighted residuals
estCov = inv(x_mat'*x_mat)*x_mat'*diag(w_eta_hat'.^2)*x_mat*inv(x_mat'*x_mat); % white Heteroscedasticity covariance estimators
%se = sqrt(diag(estCov)); % standard deviation
end

