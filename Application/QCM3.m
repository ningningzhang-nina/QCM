clc
clear all;
indexes = {'AUD_USD','NZD_USD','USD_CAD'};

for i = 1:length(indexes)
    yyyy = readtable(sprintf('./data/%s.xlsx',indexes{i}));
    yyy = yyyy(:,2);
    yy = yyy{1:end,1};
    if i==3
        yy=1./wrev(yy);
    else
        yy = wrev(yy);
    end
    y = (log(yy(2:end))-log(yy(1:end-1)))*100;

    pp = 0.01:0.01:0.99;
    cnt = 0;
    for ppp = pp
        cnt = cnt + 1;
        for method = 1:4
            tic;
            if ppp<0.5
                [~, outputs] = new_new_CAViaROptimisation(y, method, ppp);
                Q_hat{method}(cnt,:) = -outputs.VaR(:,1);
                dqtest_in(method,cnt) = outputs.DQinSample;
            else
                [~, outputs] = new_new_CAViaROptimisation(-y, method, 1-ppp);
                Q_hat{method}(cnt,:) = outputs.VaR(:,1);
                dqtest_in(method,cnt) = outputs.DQinSample;
            end
            toc;
        end
    end
    save(sprintf('final_%s_quantiles.mat',indexes{i}),"dqtest_in","Q_hat","y",'res')
    [moments] = estimate_moments(y,Q_hat,dqtest_in,pp);
    save(sprintf('final_%s_moments.mat',indexes{i}))
    clear dqtest_in pp Q_hat se_hat y res moments DnI_star cvm_pvalues DnI DnC_star DnC outputs
end


%%
clc
clear all;
indexes = {'AUD_USD','NZD_USD','USD_CAD'};
names={'AUD_USD','NZD_USD','CAD_USD'};

for i = 3
    yyyy = readtable(sprintf('./data/%s.xlsx',indexes{i}));
    yyy = yyyy(:,2);
    yy = yyy{1:end,1};
    yy = 1./wrev(yy);
    y = (log(yy(2:end))-log(yy(1:end-1)))*100;
    res = CAD_USD_residual(y');
    [~,DnC,DnC_star,DnI,DnI_star,cvm_pvalues] = CvM_test_tar_CAD_USD(y',100);
    save(sprintf('new_new_CAD_USD_cvm_test.mat',indexes{i}))
%     pp = 0.01:0.01:0.99;
%     cnt = 0;
%     for ppp = pp
%         cnt = cnt + 1;
%         for method = 1:4
%             tic;
%             if ppp<0.5
%                 [~, outputs] = new_new_CAViaROptimisation(y, method, ppp);
%                 Q_hat{method}(cnt,:) = -outputs.VaR(:,1);
%                 dqtest_in(method,cnt) = outputs.DQinSample;
%             else
%                 [~, outputs] = new_new_CAViaROptimisation(-y, method, 1-ppp);
%                 Q_hat{method}(cnt,:) = outputs.VaR(:,1);
%                 dqtest_in(method,cnt) = outputs.DQinSample;
%             end
%             toc;
%         end
%     end
%     save(sprintf('final_%s_quantiles.mat',names{i}),"dqtest_in","Q_hat","y",'res')
%     [moments,pvalues] = estimate_moments(y,Q_hat,dqtest_in,pp);
%     save(sprintf('final_%s_moments.mat',names{i}))
%     clear dqtest_in pp pvalues Q_hat se_hat y res moments DnI_star cvm_pvalues DnI DnC_star DnC outputs
end