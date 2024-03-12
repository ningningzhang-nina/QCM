clear all

cnts = 100; % set the times of realization
threshold=0.1;
counts=0;

for repeat =1:100
    counts=counts+1;
    tic;
    % generate simulated data by GARCH(1,1) model
    c=0.01;
    alpha=0.08;
    beta=0.85;
    n=1500;
    plist = 0.01:0.01:0.99;
    sigma2_t=zeros(1,n+1);
    epsilon_t=zeros(1,n+1);
    Q0=zeros(99,n+1);
    for t=2:(1500+1)
        counts=counts+1;
        sigma2_t(t)=c+alpha*epsilon_t(t-1).^2+beta*sigma2_t(t-1);
        rng(counts)
        epsilon_t(t)=normrnd(0,sqrt(sigma2_t(t)));
        for k=1:99
            Q0(k,t)=norminv(plist(k),0,sqrt(sigma2_t(t)));
        end
    end
    Y=epsilon_t(501:n); 
    Q0=Q0(1:99,501:n);
    yy(repeat,1:1000)=Y;
    Q_theory{repeat}=Q0;
    sigma2_theory(repeat,:)=sigma2_t(501:n);


    % ##########################################################################
    % ####################implied conditional moments###########################
    % ##########################################################################

    % use theoratical quantiles to estimate four conditional moments 
    row = @(p) [1,norminv(p),(norminv(p)^2 - 1),(norminv(p)^3 - 3 * norminv(p))];
    x_mat = zeros(length(plist),4);
    cnt = 0;
    for p = plist
      cnt = cnt + 1;
      x_mat(cnt,:) = row(p);
    end
    for t=1:length(Y)
        w_theta_hat = inv(x_mat.' * x_mat) * x_mat.' * Q0(:,t);

        mu_1(t) = w_theta_hat(1);
        sigma_1(t) = w_theta_hat(2);
        skew_1(t) = w_theta_hat(3)/w_theta_hat(2)*6;
        kurt_1(t) = w_theta_hat(4)/w_theta_hat(2)*24+3;
    end
    mu_hat{repeat}(1,:)=mu_1;
    sigma_hat{repeat}(1,:)=sigma_1;
    skew_hat{repeat}(1,:)=skew_1;
    kurt_hat{repeat}(1,:)=kurt_1;

    clear Q_hat
    pp = 0.01:0.01:0.99;
    cnt = 0;
    for ppp = pp
        cnt = cnt + 1;
        for method = 1:4

            if ppp<0.5
                [~, outputs] = new_new_CAViaROptimisation(Y', method, ppp);
                Q_hat{method}(cnt,:) = -outputs.VaR(:,1);

                dqtest_in{repeat}(method,cnt) = outputs.DQinSample;
            else
                [~, outputs] = new_new_CAViaROptimisation(-Y', method, 1-ppp);
                Q_hat{method}(cnt,:) = outputs.VaR(:,1);
                dqtest_in{repeat}(method,cnt) = outputs.DQinSample;
            end

        end
    end
    estimate_Q_hat{repeat,1}=Q_hat{1};
    estimate_Q_hat{repeat,2}=Q_hat{2};
    estimate_Q_hat{repeat,3}=Q_hat{3};
    estimate_Q_hat{repeat,4}=Q_hat{4};
    [moments,num] = new_new_new_estimate_moments(Y',Q_hat,dqtest_in{repeat},pp,threshold);
    estimate_moments{repeat}=moments;
    toc;
end
%%
Mdl=garch(1,1);
for rep=1:100
    EstMdl = estimate(Mdl,yy(rep,:)');
    ht=infer(EstMdl,yy(rep,:)');
    newQ=sqrt(ht).*quantile(yy(rep,:)'./sqrt(ht),0.01:0.01:0.99);
    moments2=inv(x_mat'*x_mat)*x_mat'*newQ';
    sigma2hat2(rep,:)=moments2(2,:).^2;
    skewhat2(rep,:)=moments2(3,:)./moments2(2,:)*6;
    kurthat2(rep,:)=moments2(4,:)./moments2(2,:)*24+3;
    e21(rep,:)=sigma2hat2(rep,:)-sigma2_theory(rep,:);
    e22(rep,:)=skewhat2(rep,:)-0;
    e23(rep,:)=kurthat2(rep,:)-3;

end
muhat=readtable("new_dgp1_case5_muhat.csv");
sigma2hat=readtable("new_dgp1_case5_sigma2hat.csv");
skewhat=readtable("new_dgp1_case5_skewhat.csv");
kurthat=readtable("new_dgp1_case5_kurthat.csv");
Mdl = garch(1,1);
for rep=1:100
    EstMdl = estimate(Mdl,yy(rep,:)');
    ht=infer(EstMdl,yy(rep,:)');
    newskew(rep,:)=estimate_moments{rep}(2,:).*sqrt(estimate_moments{rep}(1,:))./sqrt(ht');
    newkurt(rep,:)=(estimate_moments{rep}(3,:)-3).*sqrt(estimate_moments{rep}(1,:))./sqrt(ht')+3;
    e31(rep,:)=ht'-sigma2_theory(rep,:);
    e32(rep,:)=newskew(rep,:)-0;
    e33(rep,:)=newkurt(rep,:)-3;
    e11(rep,:)=sigma_hat{rep}(1,:).^2-sigma2_theory(rep,:);
    e12(rep,:)=skew_hat{rep}(1,:)-0;
    e13(rep,:)=kurt_hat{rep}(1,:)-3;
    e41(rep,:)=estimate_moments{rep}(1,:)-sigma2_theory(rep,:);
    e42(rep,:)=estimate_moments{rep}(2,:)-0;
    e43(rep,:)=estimate_moments{rep}(3,:)-3;
    e51(rep,:)=sigma2hat{rep,1:end}-sigma2_theory(rep,:);
    e52(rep,:)=skewhat{rep,1:end}-0;
    e53(rep,:)=kurthat{rep,1:end}-3;
end
%%
count=0;

for threshold=[0,0.1,0.3,0.5]
    count=count+1;
    for repeat=1:100
        Q_hat{1}=estimate_Q_hat{repeat,1};
        Q_hat{2}=estimate_Q_hat{repeat,2};
        Q_hat{3}=estimate_Q_hat{repeat,3};
        Q_hat{4}=estimate_Q_hat{repeat,4};
        Y=yy(repeat,:);
        pp = 0.01:0.01:0.99;
        [moments,num] = new_new_new_estimate_moments(Y',Q_hat,dqtest_in{repeat},pp,threshold);
        estimate_moments_caviar_threshold{count,repeat}=moments;
        if threshold==0
            e91(repeat,:)=moments(1,:)-sigma2_theory(repeat,:);
            e92(repeat,:)=moments(2,:)-0;
            e93(repeat,:)=moments(3,:)-3;
        elseif threshold==0.1
            e101(repeat,:)=moments(1,:)-sigma2_theory(repeat,:);
            e102(repeat,:)=moments(2,:)-0;
            e103(repeat,:)=moments(3,:)-3;
        elseif threshold==0.3
            e111(repeat,:)=moments(1,:)-sigma2_theory(repeat,:);
            e112(repeat,:)=moments(2,:)-0;
            e113(repeat,:)=moments(3,:)-3;
        else
            e121(repeat,:)=moments(1,:)-sigma2_theory(repeat,:);
            e122(repeat,:)=moments(2,:)-0;
            e123(repeat,:)=moments(3,:)-3;
        end
    end
end
save("new_DGP1_garch.mat")