clc
clear all
replication_time=100;
for repeat=1:replication_time
    tic;
    try
mu=0;
%
omega=0.01;
alpha=0.08;
beta=0.85;
T=1500;

% 
sigma2t=zeros(1,T);
epsilont=zeros(1,T);
yt=zeros(1,T);
Q_0=zeros(99,T);
w_theta_hat=zeros(4,T);
row=@(p)[1,norminv(p),norminv(p)^2-1,norminv(p)^3-3*norminv(p)];
plist=0.01:0.01:0.99;
x_mat = zeros(length(plist),4);
cnt = 0;
for p = plist
  cnt = cnt + 1;
  x_mat(cnt,:) = row(p);
end

for t=2:T
    sigma2t(t)=omega+alpha*epsilont(t-1).^2+beta*sigma2t(t-1);
    rng((repeat-1)*T+t)
    kk=3+27*rand(1);
    rng((repeat-1)*T+t)
    epsilont(t)=sqrt(sigma2t(t)/(2*kk))*(chi2rnd(kk, 1, 1)-kk);
    yt(t)=epsilont(t);
    for k=1:99
        Q_0(k,t)=sqrt(sigma2t(t)/(2*kk))*(chi2inv(k*0.01,kk)-kk);
    end
    beta_hat=inv(x_mat'*x_mat)*x_mat'*Q_0(1:99,t);
    estimate_moments{repeat}(1:3,t)=[beta_hat(2)^2;beta_hat(3)/beta_hat(2)*6;beta_hat(4)/beta_hat(2)*24+3];
    true_moments{repeat}(1:3,t)=[sigma2t(t);sqrt(8/kk);12/kk+3];
    error0(repeat,t)=estimate_moments{repeat}(3,t)-true_moments{repeat}(3,t);
   
end
y(repeat,:)=yt(501:1500);
Q0{repeat}=Q_0(1:99,501:1500);
Y=yt(501:1500);

    pp = 0.01:0.01:0.99;
    cnt = 0;
    for ppp = pp
        cnt = cnt + 1;
        for method = 1:4
            tic;
            if ppp<0.5
                [~, outputs] = new_new_CAViaROptimisation(Y', method, ppp);
                Q_hat{method}(cnt,:) = -outputs.VaR(:,1);
                
                dqtest_in{repeat}(method,cnt) = outputs.DQinSample;
            else
                [~, outputs] = new_new_CAViaROptimisation(-Y', method, 1-ppp);
                Q_hat{method}(cnt,:) = outputs.VaR(:,1);
                dqtest_in{repeat}(method,cnt) = outputs.DQinSample;
            end
            toc;
            
        end
    end
    estimate_Q_hat1{repeat}=Q_hat{1};
    estimate_Q_hat2{repeat}=Q_hat{2};
    estimate_Q_hat3{repeat}=Q_hat{3};
    estimate_Q_hat4{repeat}=Q_hat{4};
    threshold=0.1;
    [moments,num] = new_new_new_estimate_moments(Y',Q_hat,dqtest_in{repeat},pp,threshold);
    estimate_moments_caviar{repeat}=moments;
%     estimate_moments_caviar_new{repeat,cou}=moments1;
    catch
        fprintf("no")
    end
    toc;
end 
save("new_DGP6_chisquare.mat")