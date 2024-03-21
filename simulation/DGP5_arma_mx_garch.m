%clear all

%initialize parameters
cnts = 110; % set the times of realization
threshold=0.1;
for repeat = 1:cnts
    
    try
    tic;
    %% generate simulated data by ARMA-MN-GARCH(1,1) model
    alpha0 = [0.1; 0.3];
    alpha1 = [0.05; 0.1];
    beta = diag([0.85, 0.8]);
    % ARMA parameters
    a0 = 0.5;
    a1 = 0.4;
    b1 = -0.3;
    % MN parameters
    lambda1 = 0.2;
    lambda2 = 1-lambda1;
    lambda = [lambda1; lambda2];

    mu1 = 0.4;
    mu2 = - lambda1 / lambda2 * mu1;
    mu = [mu1; mu2];

    % Generate sigma
    n = 1000;
    ncut = 500;


    plist = 0.01:0.01:0.99;
    starts = (repeat-1)*(n+ncut)+1;
    ends = repeat*(n+ncut);
    seeds = starts:1:ends;
    %plist = [0.9,0.925,0.95,0.975];
    %##########################################################################
    %######################generate ARMA_MN_GARCH############################## 
    %##################return theoratical quantiles and moments################
    %##########################################################################
    output=generate_simulated_data(alpha0,alpha1,beta,a0,a1,b1,lambda,lambda2, mu,n,ncut,plist,seeds);
    Q0 = output.Q;
    Y = output.Y;
    yy(repeat,:)=Y;
    Q_theory{repeat}=Q0;
    mu_0 = output.mu_0;
    sigma_0 = output.sigma_0;
    skew_0 = output.skew_0;
    kurt_0 = output.kurt_0;
    moments_theory{repeat} = [mu_0;sigma_0;skew_0;kurt_0];
%     for i=1:n
%         [mu,h,s,k]=approximate(Q0(:,i));
%         mu_ap(i)=mu;
%         h_ap(i)=h;
%         s_ap(i)=s;
%         k_ap(i)=k;
%     end
    %##########################################################################
    %####################implied conditional moments###########################
    %##########################################################################

    %use theoratical quantiles to estimate four conditional moments 
    row = @(p) [1,norminv(p),(norminv(p)^2 - 1),(norminv(p)^3 - 3 * norminv(p))];
    x_mat = zeros(length(plist),4);
    cnt = 0;
    for p = plist
      cnt = cnt + 1;
      x_mat(cnt,:) = row(p);
    end
    for t=1:length(Y)
        %[w_theta_hat,estCov] = white_weight_regression(Q0(:,t),x_mat);
        w_theta_hat = inv(x_mat.' * x_mat) * x_mat.' * Q0(:,t);
    % theta = inv(x_mat) * Q

        mu_1(t) = w_theta_hat(1);
        sigma_1(t) = w_theta_hat(2);
        skew_1(t) = w_theta_hat(3)/w_theta_hat(2)*6;
        kurt_1(t) = w_theta_hat(4)/w_theta_hat(2)*24+3;
        %[W1,pvalue1]=point_wise_checking(Q0(:,t)-x_mat*w_theta_hat,x_mat);
        %wstat{repeat}(1,t)=W1;
        %wpvalues{repeat}(1,t)=pvalue1;
        %J = [1 0 0 0;0 2*sigma_1(t) 0 0;0 -6*w_theta_hat(3)/(sigma_1(t).^2) 6/sigma_1(t) 0;0 -24*(w_theta_hat(4)+3)/(sigma_1(t).^2) 0 24/sigma_1(t)];
        %cov_matrix = J*estCov*J';
        %se_mu1(t) = sqrt(cov_matrix(1,1));
        %se_sigma1(t) = sqrt(cov_matrix(2,2));
        %se_skew1(t) = sqrt(cov_matrix(3,3));
        %se_kurt1(t) = sqrt(cov_matrix(4,4));
    end
    mu_hat{repeat}(1,:)=mu_1;
    sigma_hat{repeat}(1,:)=sigma_1;
    skew_hat{repeat}(1,:)=skew_1;
    kurt_hat{repeat}(1,:)=kurt_1;
    %se_mu_hat{repeat}(1,:)=se_mu1;
    %se_h_hat{repeat}(1,:)=se_sigma1;
    %se_s_hat{repeat}(1,:)=se_skew1;
    %se_k_hat{repeat}(1,:)=se_kurt1;
    %#########################################
    %###############error#####################
    %#########################################
    disp(' ')
    disp('********* error **********')
    errors{repeat}(1, :) = [sum((mu_0 - mu_1).^2) / n, sum((sigma_0 - sigma_1).^2) / n, sum((sigma_0.^2 - sigma_1.^2).^2) / n, sum((skew_0 - skew_1).^2) / n,...
    sum((kurt_0 - kurt_1).^2) / n];
    relative_errors{repeat}(1, :) = [sum((mu_0 - mu_1).^2) / sum((mean(mu_0)-mu_0).^2), sum((sigma_0 - sigma_1).^2) / sum((mean(sigma_0) - sigma_0).^2), ...
        sum((sigma_0.^2 - sigma_1.^2).^2) / sum((sigma_0.^2 - mean(sigma_0.^2).^2)), sum((skew_0 - skew_1).^2) / sum((skew_0 - mean(skew_0)).^2),...
        sum((kurt_0 - kurt_1).^2) / sum((kurt_0 - mean(kurt_1)).^2)];
%%
    clear Q_hat
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
    estimate_Q_hat{repeat,1}=Q_hat{1};
    estimate_Q_hat{repeat,2}=Q_hat{2};
    estimate_Q_hat{repeat,3}=Q_hat{3};
    estimate_Q_hat{repeat,4}=Q_hat{4};
    %save(sprintf('new_new_%s_quantiles.mat',indexes{i}),"dqtest_in","Q_hat","y",'res')
    [moments,moments_approximation,num] = new_new_new_estimate_moments(Y',Q_hat,dqtest_in{repeat},pp,threshold);
    estimate_moments{repeat}=moments;
    
    catch
        fprintf("no")
    end
end
Mdl=garch(1,1);
for rep=[1:66,68:100]
    EstMdl = estimate(Mdl,yy(rep,:)');
    ht=infer(EstMdl,yy(rep,:)');
    newQ=sqrt(ht).*quantile(yy(rep,:)'./sqrt(ht),0.01:0.01:0.99);
    moments2=inv(x_mat'*x_mat)*x_mat'*newQ';
    sigma2hat2(rep,:)=moments2(2,:).^2;
    skewhat2(rep,:)=moments2(3,:)./moments2(2,:)*6;
    kurthat2(rep,:)=moments2(4,:)./moments2(2,:)*24+3;
    e21(rep,:)=sigma2hat2(rep,:)-moments_theory{rep}(2,:).^2;
    e22(rep,:)=skewhat2(rep,:)-moments_theory{rep}(3,:);
    e23(rep,:)=kurthat2(rep,:)-moments_theory{rep}(4,:);

end
muhat=readtable("dgp5_case5_muhat.csv");
sigma2hat=readtable("dgp5_case5_sigma2hat.csv");
skewhat=readtable("dgp5_case5_skewhat.csv");
kurthat=readtable("dgp5_case5_kurthat.csv");
Mdl = garch(1,1);
for rep=[1:66,68:100]
    EstMdl = estimate(Mdl,yy(rep,:)');
    ht=infer(EstMdl,yy(rep,:)');
    newskew(rep,:)=estimate_moments{rep}(3,:).*estimate_moments{rep}(2,:)./sqrt(ht');
    newkurt(rep,:)=(estimate_moments{rep}(4,:)-3).*estimate_moments{rep}(2,:)./sqrt(ht')+3;
    e31(rep,:)=ht'-moments_theory{rep}(2,:).^2;
    e32(rep,:)=newskew(rep,:)-moments_theory{rep}(3,:);
    e33(rep,:)=newkurt(rep,:)-moments_theory{rep}(4,:);
    e11(rep,:)=sigma_hat{rep}(1,:).^2-moments_theory{rep}(2,:).^2;
    e12(rep,:)=skew_hat{rep}(1,:)-moments_theory{rep}(3,:);
    e13(rep,:)=kurt_hat{rep}(1,:)-moments_theory{rep}(4,:);
    e41(rep,:)=estimate_moments{rep}(2,:).^2-moments_theory{rep}(2,:).^2;
    e42(rep,:)=estimate_moments{rep}(3,:)-moments_theory{rep}(3,:);
    e43(rep,:)=estimate_moments{rep}(4,:)-moments_theory{rep}(4,:);
    e51(rep,:)=sigma2hat{rep,1:end}-moments_theory{rep}(2,:).^2;
    e52(rep,:)=skewhat{rep,1:end}-moments_theory{rep}(3,:);
    e53(rep,:)=kurthat{rep,1:end}-moments_theory{rep}(4,:);
end
count=0;

for threshold=[0,0.1,0.3,0.5]
    count=count+1;
    for repeat=[1:66,68:100]
        Q_hat{1}=estimate_Q_hat{repeat,1};
        Q_hat{2}=estimate_Q_hat{repeat,2};
        Q_hat{3}=estimate_Q_hat{repeat,3};
        Q_hat{4}=estimate_Q_hat{repeat,4};
        Y=yy(repeat,:);
        pp = 0.01:0.01:0.99;
        [moments,num] = new_new_new_estimate_moments(Y',Q_hat,dqtest_in{repeat},pp,threshold);
        estimate_moments_caviar_threshold{count,repeat}=moments;
        if threshold==0
            e91(repeat,:)=moments(1,:)-moments_theory{repeat}(2,:).^2;
            e92(repeat,:)=moments(2,:)-moments_theory{repeat}(3,:);
            e93(repeat,:)=moments(3,:)-moments_theory{repeat}(4,:);
        elseif threshold==0.1
            e101(repeat,:)=moments(1,:)-moments_theory{repeat}(2,:).^2;
            e102(repeat,:)=moments(2,:)-moments_theory{repeat}(3,:);
            e103(repeat,:)=moments(3,:)-moments_theory{repeat}(4,:);
        elseif threshold==0.3
            e111(repeat,:)=moments(1,:)-moments_theory{repeat}(2,:).^2;
            e112(repeat,:)=moments(2,:)-moments_theory{repeat}(3,:);
            e113(repeat,:)=moments(3,:)-moments_theory{repeat}(4,:);
        else
            e121(repeat,:)=moments(1,:)-moments_theory{repeat}(2,:).^2;
            e122(repeat,:)=moments(2,:)-moments_theory{repeat}(3,:);
            e123(repeat,:)=moments(3,:)-moments_theory{repeat}(4,:);
        end
    end
    
end
save('./results/DGP5_arma_mx_garch.mat')
