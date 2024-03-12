function [moments,num] = new_new_new_estimate_moments(y,Q_hat,dqtest_in,pp,threshold)
% y: stock index returns
% Q_hat: estimated quantiles by 4 method of CAViaR
% dqtest_in: p values of in-sample for testing the validaity of estimates
% quantiles
%pp: quantile levels;
%threshold = 0.1 ;
qq = [];
ppp=[];
for method=1:size(dqtest_in,1)
    if threshold==0
        index=1:99;
    else
        index = find(dqtest_in(method,:)>threshold);
    end
    qq = [qq;Q_hat{method}(index,:)];
    ppp = [ppp,pp(index)];
end
x_mat=zeros(length(ppp),4);
row=@(p)[1,norminv(p),norminv(p)^2-1,norminv(p)^3-3*norminv(p)];
% row = @(p) [norminv(p),1/4*(norminv(p)^3-3*norminv(p)),...
%     1/96*(5*norminv(p)^5-56*norminv(p)^3+75*norminv(p)),...
%    1/384*(3*norminv(p)^7-81*norminv(p)^5+417*norminv(p)^3-315*norminv(p))];
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
    x=inv(x_mat'*x_mat)*x_mat'*y1;
%     fun=@(x)y1-x_mat*[x(1);x(1)/(x(2)-2);x(1)/(x(2)-2)/(x(2)-2);x(1)/(x(2)-2)/(x(2)-2)/(x(2)-2)];
%     x0=[sqrt(1.2),6];
%     lb=[0,4.5];
%     ub=[100,100];
%     x = lsqnonlin(fun,x0,lb,ub)
%     [x(1)^2;0;6/(x(2)-4)+3]
%     moments(1:3,i)=[x(1)^2;0;6/(x(2)-4)+3];
moments(1:3,i)=[x(2)^2;x(3)/x(2)*6;x(4)/x(2)*24+3];
%     [mu,h,s,k]=approximate(y1,ppp);
%     moments_ap(1:4,i)=[mu;h;s;k];
%     % point_wise_checking
%     epsilon_hat = y1-x_mat*w_theta_hat;
%     [W,pvalue]=point_wise_checking(epsilon_hat,x_mat);
%     point_wise_pvalue(1,i)=pvalue;
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
num=length(y1);
% pvalues = test_moments_unbias(y,moments);
% % test sigma
% sigma0 = moments(2,:);
% %size((sigma0').^2)
% %size((y).^2)
% [p2,~] = waldtest_sd((sigma0').^2,(y).^2);
% % test skewness and kurtosis
% skew0 = moments(3,:);
% kurt0 = moments(4,:);
% [p3,~,p4,~] = waldtest_skewness_kurtosis((skew0)',(y).^3./((sigma0').^3),...
%                  ((kurt0))',(y).^4./((sigma0').^4));
% 
% % p values of four wald test
% pvalues = [p2,p3,p4];
end
% function [h,s,k]=newqcm(yit,n,ppp)
%     row = @(p) [(p-0.5),(p-0.5).^2,(p-0.5).^3,(p-0.5).^4,(p-0.5).^5,(p-0.5).^6,...
%     (p-0.5).^7,(p-0.5).^8,(p-0.5).^9,(p-0.5).^(10),(p-0.5).^(11),norminv(p),norminv(p)^2-1,norminv(p)^3-3*norminv(p)];
%     new_x_mat = zeros(length(ppp),14);
%     cnt = 0;
%     for p = ppp
%       cnt = cnt + 1;
%       new_x_mat(cnt,:) = row(p);
%     end
%     new_new_x_mat=new_x_mat(1:end,[1:3,4:(n+3)]);
%     mdl = fitlm(new_new_x_mat,yit)
%     betabar=inv(new_new_x_mat.' * new_new_x_mat) * new_new_x_mat.' * yit;
%     newbeta=betabar(end-2:end,1:end);
%     h=newbeta(1,1:end).^2;
%     s=newbeta(2,1:end)./newbeta(1,1:end)*6;
%     k=newbeta(3,1:end)./newbeta(1,1:end)*24+3;
% end