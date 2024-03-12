clear all
cnts = 100; % set the times of realization
threshold=0.1;
counts=0;
yyy=zeros(1,1000);
sigma2_theory=zeros(1,1000);
skew_theory=zeros(1,1000);
kurt_theory=zeros(1,1000);
sigma2_hat=zeros(1,1000);
skew_hat=zeros(1,1000);
kurt_hat=zeros(1,1000);
for rep = 1:100
    tic;
counts=counts+1;

% generate simulated data by GARCH(1,1)-snp model
n=1500;
c=0.01;
alpha1=0.01;
alpha2=0.07;
beta=0.85;
% phi01=0.132;
% phi11=0.780;
% phi211=0.003;
% phi212=0.07;
% phi31=0.186;
% phi02=0.036;
% phi12=0.865;
% phi221=0.093;
% phi222=0.110;
% phi32=-0.281;
plist = 0.01:0.01:0.99;
sigma2_t=zeros(1,n+1);
epsilon_t=zeros(1,n+1);
zt=zeros(1,n+1);
etat=zeros(1,n+1);
lambdat=zeros(1,n+1);
skew=zeros(1,n+1);
kurt=zeros(1,n+1);
at=zeros(1,n+1);
bt=zeros(1,n+1);
ct=zeros(1,n+1);
Q0=zeros(99,n+1);
skew=zeros(1,n+1);
kurt=zeros(1,n+1);
for t = 2:1500
counts=counts+1;
sigma2_t(t)=c+alpha1*epsilon_t(t-1)^2+alpha2*epsilon_t(t-1)^2*(epsilon_t(t-1)<0)+beta*sigma2_t(t-1);
lambdat(t)=-0.02-0.13*epsilon_t(t-1)-0.1*epsilon_t(t-1)^2;
etat(t)=5.1+(50-5.1)/(1+exp(-lambdat(t)));
ct(t)=gamma((etat(t)+1)/2)/sqrt(pi*(etat(t)-2))/gamma(etat(t)/2);
at(t)=4*lambdat(t)*ct(t)*(etat(t)-2)/(etat(t)-1);
bt(t)=sqrt(1+3*lambdat(t)^2-at(t)^2);
sy=(4*at(t)*(1+lambdat(t)^2)*(etat(t)-2)/(etat(t)-3));
skew(t)=(sy-3*at(t)*bt(t)^2-at(t)^3)/(bt(t)^3);
ky=(6/(etat(t)-4)+3)*(1+10*lambdat(t)^2+5*lambdat(t)^4);
kurt(t)=(ky-4*at(t)*sy+3*at(t)^4+6*at(t)^2*bt(t)^2)/bt(t)^4;

rng(counts)
randr=rand(1) ;
if randr<(1-lambdat(t))/2
s=tinv(randr/(1-lambdat(t)),etat(t));
yy=s*(1-lambdat(t))/sqrt(etat(t)/(etat(t)-2));
else
s=tinv((randr+lambdat(t))/(1+lambdat(t)),etat(t));
yy=s*(1+lambdat(t))/sqrt(etat(t)/(etat(t)-2));
end
%print(c(t,yy))
zt(t)=(yy-at(t))./bt(t);
epsilon_t(t)=zt(t)*sqrt(sigma2_t(t));
for k = 1:99
  if plist(k)<(1-lambdat(t))/2
    s=tinv(plist(k)/(1-lambdat(t)),etat(t));
    yy=s*(1-lambdat(t))/sqrt(etat(t)/(etat(t)-2));
  else 
    s=tinv((plist(k)+lambdat(t))/(1+lambdat(t)),etat(t));
    yy=s*(1+lambdat(t))/sqrt(etat(t)/(etat(t)-2));
  end
  Q0(k,t)=(yy-at(t))./bt(t).*sqrt(sigma2_t(t));
end
end
Y=epsilon_t(501:n);
Q0=Q0(1:99,501:n);
yyy(rep,1:1000)=Y;
Q_theory{rep}=Q0;
sigma2_theory(rep,:)=sigma2_t(501:n);
skew_theory(rep,:)=skew(501:n);
kurt_theory(rep,:)=kurt(501:n);
row=@(p)[1,norminv(p),(norminv(p)^2-1),(norminv(p)^3-3*norminv(p))];

x_mat=zeros(99,4);
for i = 1:99
  x_mat(i,:)=row(plist(i));
end
sigma_1=zeros(1,1000);
skew_1=zeros(1,1000);
kurt_1=zeros(1,1000);
for t = 1:1000
  w_theta_hat=inv(x_mat'*x_mat)*(x_mat')*Q0(:,t);
  sigma_1(t)=w_theta_hat(2);
  skew_1(t)=w_theta_hat(3)/w_theta_hat(2)*6;
  kurt_1(t)=w_theta_hat(4)/w_theta_hat(2)*24+3;
end
sigma2_hat(rep,:)=sigma_1.^2;
skew_hat(rep,:)=skew_1;
kurt_hat(rep,:)=kurt_1;
   %clear Q_hat
    pp = 0.01:0.01:0.99;
    cnt = 0;
    for ppp = pp
        cnt = cnt + 1;
        for method = 1:4
            if ppp<0.5
                [~, outputs] = new_new_CAViaROptimisation(Y', method, ppp);
                Q_hat{method}(cnt,:) = -outputs.VaR(:,1);

                dqtest_in{rep}(method,cnt) = outputs.DQinSample;
            else
                [~, outputs] = new_new_CAViaROptimisation(-Y', method, 1-ppp);
                Q_hat{method}(cnt,:) = outputs.VaR(:,1);
                dqtest_in{rep}(method,cnt) = outputs.DQinSample;
            end
        end
    end
    estimate_Q_hat{rep,1}=Q_hat{1};
    estimate_Q_hat{rep,2}=Q_hat{2};
    estimate_Q_hat{rep,3}=Q_hat{3};
    estimate_Q_hat{rep,4}=Q_hat{4};
    [moments,num] = new_new_new_estimate_moments(Y',Q_hat,dqtest_in{rep},pp,threshold);
    estimate_moments{rep}=moments;
    toc;
end
Mdl=garch(1,1);
for rep=1:100
    EstMdl = estimate(Mdl,yyy(rep,:)');
    ht=infer(EstMdl,yyy(rep,:)');
    newQ=sqrt(ht).*quantile(yyy(rep,:)'./sqrt(ht),0.01:0.01:0.99);
    moments2=inv(x_mat'*x_mat)*x_mat'*newQ';
    sigma2hat2(rep,:)=moments2(2,:).^2;
    skewhat2(rep,:)=moments2(3,:)./moments2(2,:)*6;
    kurthat2(rep,:)=moments2(4,:)./moments2(2,:)*24+3;
    e21(rep,:)=sigma2hat2(rep,:)-sigma2_theory(rep,:);
    e22(rep,:)=skewhat2(rep,:)-skew_theory(rep,:);
    e23(rep,:)=kurthat2(rep,:)-kurt_theory(rep,:);

end
muhat=readtable("new_dgp4_case5_muhat.csv");
sigma2hat=readtable("new_dgp4_case5_sigma2hat.csv");
skewhat=readtable("new_dgp4_case5_skewhat.csv");
kurthat=readtable("new_dgp4_case5_kurthat.csv");
Mdl = garch(1,1);
for rep=1:100
    EstMdl = estimate(Mdl,yyy(rep,:)');
    ht=infer(EstMdl,yyy(rep,:)');
    newskew(rep,:)=estimate_moments{rep}(2,:).*sqrt(estimate_moments{rep}(1,:))./sqrt(ht');
    newkurt(rep,:)=(estimate_moments{rep}(3,:)-3).*sqrt(estimate_moments{rep}(1,:))./sqrt(ht')+3;
    e31(rep,:)=ht'-sigma2_theory(rep,:);
    e32(rep,:)=newskew(rep,:)-skew_theory(rep,:);
    e33(rep,:)=newkurt(rep,:)-kurt_theory(rep,:);
    e11(rep,:)=sigma2_hat(rep,:)-sigma2_theory(rep,:);
    e12(rep,:)=skew_hat(rep,:)-skew_theory(rep,:);
    e13(rep,:)=kurt_hat(rep,:)-kurt_theory(rep,:);
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
        Y=yyy(repeat,:);
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
save("new_dpg4_gir_skew_studentt.mat")

