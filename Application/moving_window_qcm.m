clc
clear all;
indexes = {'AUD_USD','NZD_USD','USD_CAD'};

for i = 2
    yyyy = readtable(sprintf('./data/%s.xlsx',indexes{i}));
    yyy = yyyy(:,2);
    yy = yyy{1:end,1};
    yy = wrev(yy);
    y = (log(yy(2:end))-log(yy(1:end-1)))*100;
    for j=481:500
        newy=y(1:end-500+j);
    pp = 0.01:0.01:0.99;
    cnt = 0;
    for ppp = pp
        cnt = cnt + 1;
        for method = 1:4
            tic;
            if ppp<0.5
                [~, outputs] = new_new_new_CAViaROptimisation(newy, method, ppp);
                Q_hat{method}(cnt,:) = -outputs.VaR(:,1);
                dqtest_in(method,cnt) = outputs.DQinSample;
            else
                [~, outputs] = new_new_new_CAViaROptimisation(-newy, method, 1-ppp);
                Q_hat{method}(cnt,:) = outputs.VaR(:,1);
                dqtest_in(method,cnt) = outputs.DQinSample;
            end
            toc;
        end
    end
    [moments] = newest_estimate_moments(Q_hat,dqtest_in,pp);
    save(sprintf('./results/final_%s_moments_%d.mat',indexes{i},j))
    clear dqtest_in pp pvalues Q_hat se_hat newy moments DnI_star cvm_pvalues DnI DnC_star DnC outputs
    end
    
end
