function [moments] = newest_estimate_moments(Q_hat,dqtest_in,pp)
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
moments = zeros(3,1);
y1 = qq(:,end);
w_theta_hat=inv(x_mat'*x_mat)*x_mat'*y1;
moments(1)=w_theta_hat(1);
moments(2)=w_theta_hat(2);
moments(3)=w_theta_hat(3)/w_theta_hat(2)*6;
moments(4)=w_theta_hat(4)/w_theta_hat(2)*24+3;

end