function [moments] = estimate_moments(y,Q_hat,dqtest_in,pp)
% y: stock index returns
% Q_hat: estimated quantiles by 4 method of CAViaR
% dqtest_in: p values of in-sample for testing the validaity of estimates
% quantiles
%pp: quantile levels;
threshold = 0.1 ;
qq = [];
ppp=[];
for method=1:4
    index = find(dqtest_in(method,:)>threshold);
    qq = [qq;Q_hat{method}(index,:)];
    ppp = [ppp,pp(index)];
end
x_mat=zeros(length(ppp),4);
row = @(p) [1,norminv(p),(norminv(p)^2 - 1),(norminv(p)^3 - 3 * norminv(p))];
cnt=0;
for p = ppp
    cnt = cnt+1;
    x_mat(cnt,:)=row(p);
end
moments = zeros(3,length(y));
%new_moments = zeros(3,length(y));
%se_hat = zeros(3,length(y));
for i=1:length(y)
    y1 = qq(:,i);
    %[w_theta_hat,estCov,omega_hat] = white_weight_regression(y1,x_mat);
    w_theta_hat=inv(x_mat'*x_mat)*x_mat'*y1;
    moments(1,i)=w_theta_hat(1);
    moments(2,i)=w_theta_hat(2);
    moments(3,i)=w_theta_hat(3)/w_theta_hat(2)*6;
    moments(4,i)=w_theta_hat(4)/w_theta_hat(2)*24+3;
%     J = [2*moments(1,i) 0 0;...
%         -6*w_theta_hat(2)/(moments(1,i).^2) 6/moments(1,i) 0;...
%         -24*(w_theta_hat(3))/(moments(1,i).^2) 0 24/moments(1,i)];
%         cov_matrix = J*estCov*J';
%         se_hat(1,i) = sqrt(cov_matrix(1,1));
%         se_hat(2,i) = sqrt(cov_matrix(2,2));
%         se_hat(3,i) = sqrt(cov_matrix(3,3));
%    [parameters, LLF] = nonlinear_constraint_linear(y1,x_mat,omega_hat,w_theta_hat,[]) ;
%    new_moments(1,i) = parameters(1);
%    new_moments(2,i)=parameters(2)/parameters(1)*6;
%    new_moments(3,i)=parameters(3)/parameters(1)*24+3;
end

% test sigma
% sigma0 = moments(1,:);
% %size((sigma0').^2)
% %size((y).^2)
% [p2,~] = waldtest_sd((sigma0').^2,(y).^2);
% % test skewness and kurtosis
% skew0 = moments(2,:);
% kurt0 = moments(3,:);
% [p3,~,p4,~] = waldtest_skewness_kurtosis((skew0)',(y).^3./((sigma0').^3),...
%                  ((kurt0))',(y).^4./((sigma0').^4));
% 
% % p values of four wald test
% pvalues = [p2,p3,p4];
end