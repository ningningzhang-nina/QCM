clear all
cnts = 100; % set the times of realization
threshold=0.1;
counts=0;
for repeat =1:cnts
    counts=counts+1;
    tic;
    % generate simulated data by GARCH(1,1)-snp model
    n=1500;
    c=0.01;alpha=0.08;beta=0.85;
    phi01=0.05;phi11=0.60;phi211=0.15;phi212=-0.02;phi31=0;
    phi02=0.04;phi12=0.50;phi221=0.2;phi222=-0.05;phi32=0;
    plist = 0.01:0.01:0.99;
    sigma2_t=zeros(1,n+1);epsilon_t=zeros(1,n+1);zt=zeros(1,n+1);xt=zeros(1,n+1);
    nu1=zeros(1,n+1);nu2=zeros(1,n+1);
    gamma1=zeros(1,n+1);gamma2=zeros(1,n+1);gamma3=zeros(1,n+1);gamma4=zeros(1,n+1);
    skew=zeros(1,n+1);kurt=zeros(1,n+1);
    at=zeros(1,n+1);bt=zeros(1,n+1);
    Q0=zeros(99,n+1);

    for t=2:1500
        counts=counts+1;
        sigma2_t(t)=c+alpha*epsilon_t(t-1).^2+beta*sigma2_t(t-1);
        nu1(t)=phi01+phi11*nu1(t-1)+phi211*(1+phi31*abs(zt(t-1)))*max(zt(t-1),0)+...
            phi212*(1+phi31*abs(zt(t-1)))*min(zt(t-1),0);
        nu2(t)=phi02+phi12*nu2(t-1)+phi221*(1+phi32*abs(zt(t-1)))*max(zt(t-1),0)+...
            phi222*(1+phi32*abs(zt(t-1)))*min(zt(t-1),0);
        
        %[nu1(t),nu2(t)]
        gamma1(t)=2*nu1(t)*(1+sqrt(2)*nu2(t))/(1+nu1(t)^2+nu2(t)^2);
        gamma2(t)=sqrt(2)*(nu1(t)^2+2*nu2(t)^2+sqrt(2)*nu2(t))/(1+nu1(t)^2+nu2(t)^2);
        gamma3(t)=2*sqrt(3)*nu1(t)*nu2(t)/(1+nu1(t)^2+nu2(t)^2);
        gamma4(t)=sqrt(6)*nu2(t)^2/(1+nu1(t)^2+nu2(t)^2);
        % dd=@(x)(normpdf(x).*(H(0,x)+gamma1(t)*H(1,x)+gamma2(t)*H(2,x)+gamma3(t)*H(3,x)+gamma4(t)*H(4,x)));
        % figure
        % plot(-4:0.1:4,dd(-4:0.1:4))
        mux1=gamma1(t);
        mux2=sqrt(2)*gamma2(t)+1;
        mux3=6*nu1(t)*(1+2*sqrt(2)*nu2(t))/(1+nu1(t)^2+nu2(t)^2);
        mux4=12*(nu1(t)^2+3*nu2(t)^2+sqrt(2)*nu2(t))/(1+nu1(t)^2+nu2(t)^2)+3;
        sigma2_x=mux2-mux1^2;
        bt(t)=1/sqrt(sigma2_x);
        at(t)=-bt(t)*mux1;
        skew(t)=at(t)^3+3*at(t)^2*bt(t)*mux1+3*at(t)*bt(t)^2*mux2+bt(t)^3*mux3;
        kurt(t)=at(t)^4+4*at(t)^3*bt(t)*mux1+6*at(t)^2*bt(t)^2*mux2+4*at(t)*bt(t)^3*mux3+bt(t)^4*mux4;
        rng(counts)
        randr=rand(1);
        f1=@(rt_star)normcdf(rt_star)-normpdf(rt_star).*(gamma1(t)/sqrt(1).*H(0,rt_star)+...
                gamma2(t)/sqrt(2).*H(1,rt_star)+gamma3(t)/sqrt(3).*H(2,rt_star)+...
                gamma4(t)/sqrt(4).*H(3,rt_star))-randr;
        % figure
        % plot(-4:0.1:4,f1(-4:0.1:4)+randr)
        s=fzero(f1,0);
        xt(t)=s;
        % syms rt_star
        % s=vpasolve(simplify(normcdf(rt_star)-normpdf(rt_star)*(gamma1(t)/sqrt(1)*H(0,rt_star)+...
        %         gamma2(t)/sqrt(2)*H(1,rt_star)+gamma3(t)/sqrt(3)*H(2,rt_star)+...
        %         gamma4(t)/sqrt(4)*H(3,rt_star))-randr==0));
        zt(t)=s*bt(t)+at(t);
        epsilon_t(t)=zt(t)*sqrt(sigma2_t(t));

        for k=1:99
            f2=@(rt_star)(normcdf(rt_star)-normpdf(rt_star)*(gamma1(t)/sqrt(1)*H(0,rt_star)+...
                gamma2(t)/sqrt(2)*H(1,rt_star)+gamma3(t)/sqrt(3)*H(2,rt_star)+...
                gamma4(t)/sqrt(4)*H(3,rt_star))-plist(k));
            s=fzero(f2,0);
            Q0(k,t)=s*sqrt(sigma2_t(t))*bt(t)+at(t)*sqrt(sigma2_t(t));
        end
    end
    Y=epsilon_t(501:n); 
    Q0=Q0(1:99,501:n);
    yy(repeat,1:1000)=Y;
    Q_theory{repeat}=Q0;
    sigma2_theory(repeat,:)=sigma2_t(501:n);
    skew_theory(repeat,:)=skew(501:n);
    kurt_theory(repeat,:)=kurt(501:n);
    

%    use theoratical quantiles to estimate four conditional moments 
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
    %clear Q_hat
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
    e22(rep,:)=skewhat2(rep,:)-skew_theory(rep,:);
    e23(rep,:)=kurthat2(rep,:)-kurt_theory(rep,:);

end
muhat=readtable("new_dgp3_case5_muhat.csv");
sigma2hat=readtable("new_dgp3_case5_sigma2hat.csv");
skewhat=readtable("new_dgp3_case5_skewhat.csv");
kurthat=readtable("new_dgp3_case5_kurthat.csv");
Mdl = garch(1,1);
for rep=1:100
    EstMdl = estimate(Mdl,yy(rep,:)');
    ht=infer(EstMdl,yy(rep,:)');
    newskew(rep,:)=estimate_moments{rep}(2,:).*sqrt(estimate_moments{rep}(1,:))./sqrt(ht');
    newkurt(rep,:)=(estimate_moments{rep}(3,:)-3).*sqrt(estimate_moments{rep}(1,:))./sqrt(ht')+3;
    e31(rep,:)=ht'-sigma2_theory(rep,:);
    e32(rep,:)=newskew(rep,:)-skew_theory(rep,:);
    e33(rep,:)=newkurt(rep,:)-kurt_theory(rep,:);
    e11(rep,:)=sigma_hat{rep}(1,:).^2-sigma2_theory(rep,:);
    e12(rep,:)=skew_hat{rep}(1,:)-skew_theory(rep,:);
    e13(rep,:)=kurt_hat{rep}(1,:)-kurt_theory(rep,:);
    e41(rep,:)=estimate_moments{rep}(1,:)-sigma2_theory(rep,:);
    e42(rep,:)=estimate_moments{rep}(2,:)-skew_theory(rep,:);
    e43(rep,:)=estimate_moments{rep}(3,:)-kurt_theory(rep,:);
    e51(rep,:)=sigma2hat{rep,1:end}-sigma2_theory(rep,:);
    e52(rep,:)=skewhat{rep,1:end}-skew_theory(rep,:);
    e53(rep,:)=kurthat{rep,1:end}-kurt_theory(rep,:);
end
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
            e92(repeat,:)=moments(2,:)-skew_theory(repeat,:);
            e93(repeat,:)=moments(3,:)-kurt_theory(repeat,:);
        elseif threshold==0.1
            e101(repeat,:)=moments(1,:)-sigma2_theory(repeat,:);
            e102(repeat,:)=moments(2,:)-skew_theory(repeat,:);
            e103(repeat,:)=moments(3,:)-kurt_theory(repeat,:);
        elseif threshold==0.3
            e111(repeat,:)=moments(1,:)-sigma2_theory(repeat,:);
            e112(repeat,:)=moments(2,:)-skew_theory(repeat,:);
            e113(repeat,:)=moments(3,:)-kurt_theory(repeat,:);
        else
            e121(repeat,:)=moments(1,:)-sigma2_theory(repeat,:);
            e122(repeat,:)=moments(2,:)-skew_theory(repeat,:);
            e123(repeat,:)=moments(3,:)-kurt_theory(repeat,:);
        end
    end
end
save("./results/new_DGP3_garch_snp.mat")


function h=H(k,x)
if k==0
    h=1;
elseif k==1
    h=x;
elseif k==2
    h=(x.^2-1)./sqrt(2);
elseif k==3
    h=x.^3/sqrt(6)-(1/sqrt(6)+sqrt(2/3))*x;
else
    h=x.^4/(2*sqrt(6))-1/2*(1/sqrt(6)+sqrt(2/3)+sqrt(1.5))*x.^2+1/2*sqrt(1.5);
end
end
