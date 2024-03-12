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
    nu=5+15*rand(1);
    rng((repeat-1)*T+t)
    epsilont(t)=sqrt(sigma2t(t)*(nu-2)/nu)*trnd(nu);
    yt(t)=epsilont(t);
    for k=1:99
        Q_0(k,t)=sqrt(sigma2t(t)*(nu-2)/nu)*tinv(k*0.01,nu);
    end
    beta_hat=inv(x_mat'*x_mat)*x_mat'*Q_0(1:99,t);
    estimate_moments{repeat}(1:3,t)=[beta_hat(2)^2;beta_hat(3)/beta_hat(2)*6;beta_hat(4)/beta_hat(2)*24+3];
    true_moments{repeat}(1:3,t)=[sigma2t(t);0;6/(nu-4)+3];
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
Mdl=garch(1,1);
for rep=1:100
    EstMdl = estimate(Mdl,y(rep,:)');
    ht=infer(EstMdl,y(rep,:)');
    newQ=sqrt(ht).*quantile(y(rep,:)'./sqrt(ht),0.01:0.01:0.99);
    moments2=inv(x_mat'*x_mat)*x_mat'*newQ';
    sigma2hat2(rep,:)=moments2(2,:).^2;
    skewhat2(rep,:)=moments2(3,:)./moments2(2,:)*6;
    kurthat2(rep,:)=moments2(4,:)./moments2(2,:)*24+3;
    e21(rep,:)=sigma2hat2(rep,1:end)-true_moments{rep}(1,501:end);
    e22(rep,:)=skewhat2(rep,1:end)-true_moments{rep}(2,501:end);
    e23(rep,:)=kurthat2(rep,1:end)-true_moments{rep}(3,501:end);

end
muhat=readtable("new_dgp2_case5_muhat.csv");
sigma2hat=readtable("new_dgp2_case5_sigma2hat.csv");
skewhat=readtable("new_dgp2_case5_skewhat.csv");
kurthat=readtable("new_dgp2_case5_kurthat.csv");
Mdl = garch(1,1);
for rep=1:100%[1:85,87:100]
    EstMdl = estimate(Mdl,y(rep,:)');
    ht=infer(EstMdl,y(rep,:)');
    newskew(rep,:)=estimate_moments_caviar{rep}(2,:).*sqrt(estimate_moments_caviar{rep}(1,:))./sqrt(ht(:)');
    newkurt(rep,:)=(estimate_moments_caviar{rep}(3,:)-3).*sqrt(estimate_moments_caviar{rep}(1,:))./sqrt(ht')+3;
    e31(rep,:)=ht'-true_moments{rep}(1,501:end);
    e32(rep,:)=newskew(rep,:)-true_moments{rep}(2,501:end);
    e33(rep,:)=newkurt(rep,:)-true_moments{rep}(3,501:end);
    e11(rep,:)=estimate_moments{rep}(1,501:end)-true_moments{rep}(1,501:end);
    e12(rep,:)=estimate_moments{rep}(2,501:end)-true_moments{rep}(2,501:end);
    e13(rep,:)=estimate_moments{rep}(3,501:end)-true_moments{rep}(3,501:end);
    e41(rep,:)=estimate_moments_caviar{rep}(1,:)-true_moments{rep}(1,501:end);
    e42(rep,:)=estimate_moments_caviar{rep}(2,:)-true_moments{rep}(2,501:end);
    e43(rep,:)=estimate_moments_caviar{rep}(3,:)-true_moments{rep}(3,501:end);
    e51(rep,:)=sigma2hat{rep,1:end}-true_moments{rep}(1,501:end);
    e52(rep,:)=skewhat{rep,1:end}-true_moments{rep}(2,501:end);
    e53(rep,:)=kurthat{rep,1:end}-true_moments{rep}(3,501:end);
end
count=0;

for threshold=[0,0.1,0.3,0.5]
    count=count+1;
    for repeat=1:100%[1:85,87:100]
        Q_hat{1}=estimate_Q_hat{repeat,1};
        Q_hat{2}=estimate_Q_hat{repeat,2};
        Q_hat{3}=estimate_Q_hat{repeat,3};
        Q_hat{4}=estimate_Q_hat{repeat,4};
        Y=y(repeat,:);
        pp = 0.01:0.01:0.99;
        [moments,num] = new_new_new_estimate_moments(Y',Q_hat,dqtest_in{repeat},pp,threshold);
        estimate_moments_caviar_threshold{count,repeat}=moments;
        if threshold==0
            e91(repeat,:)=moments(1,:)-true_moments{repeat}(1,501:end);
            e92(repeat,:)=moments(2,:)-true_moments{repeat}(2,501:end);
            e93(repeat,:)=moments(3,:)-true_moments{repeat}(3,501:end);
        elseif threshold==0.1
            e101(repeat,:)=moments(1,:)-true_moments{repeat}(1,501:end);
            e102(repeat,:)=moments(2,:)-true_moments{repeat}(2,501:end);
            e103(repeat,:)=moments(3,:)-true_moments{repeat}(3,501:end);
        elseif threshold==0.3
            e111(repeat,:)=moments(1,:)-true_moments{repeat}(1,501:end);
            e112(repeat,:)=moments(2,:)-true_moments{repeat}(2,501:end);
            e113(repeat,:)=moments(3,:)-true_moments{repeat}(3,501:end);
        else
            e121(repeat,:)=moments(1,:)-true_moments{repeat}(1,501:end);
            e122(repeat,:)=moments(2,:)-true_moments{repeat}(2,501:end);
            e123(repeat,:)=moments(3,:)-true_moments{repeat}(3,501:end);
        end
    end
end
save("new_DGP2_studentt.mat")