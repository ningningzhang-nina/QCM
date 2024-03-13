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
%     y=y';
%     X0=ones(length(y),1);
%     threshold_variable=[0;y(1:end-1)'];
%     X1 = [0;y(1:end-1)'];
%     X2 = [0;0;y(1:end-2)'];
%     X3 = [0;0;0;y(1:end-3)'];
%     X4 = [0;0;0;0;y(1:end-4)'];
%     X5 = [0;0;0;0;0;y(1:end-5)'];
%     X6 = [0;0;0;0;0;0;y(1:end-6)'];
%     X7 = [0;0;0;0;0;0;0;y(1:end-7)'];
%     X8 = [0;0;0;0;0;0;0;0;y(1:end-8)'];
%     X9 = [0;0;0;0;0;0;0;0;0;y(1:end-9)'];
%     XX0 = X0.*(threshold_variable<0);
%     XXX0 = X0.*(threshold_variable>=0);
%     XX1 = X1.*(threshold_variable<0);
%     XXX1 = X1.*(threshold_variable>=0);
%     XX2 = X2.*(threshold_variable<0);
%     XXX2 = X2.*(threshold_variable>=0);
%     XX3 = X3.*(threshold_variable<0);
%     XXX3 = X3.*(threshold_variable>=0);
%     XX4 = X4.*(threshold_variable<0);
%     XXX4 = X4.*(threshold_variable>=0);
%     XX5 = X5.*(threshold_variable<0);
%     XXX5 = X5.*(threshold_variable>=0);
%     XX6 = X6.*(threshold_variable<0);
%     XXX6 = X6.*(threshold_variable>=0);
%     XX7 = X7.*(threshold_variable<0);
%     XXX7 = X7.*(threshold_variable>=0);
%     XX8 = X8.*(threshold_variable<0);
%     XXX8 = X8.*(threshold_variable>=0);
%     XX9 = X9.*(threshold_variable<0);
%     XXX9 = X9.*(threshold_variable>=0);
%     mdl=fitlm([XX0,XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XXX0,XXX1,XXX2,XXX3,XXX4,XXX5,XXX6,XXX7,XXX8,XXX9],y,'Intercept',false);
%     mdl
    nboot=100;
    if i==1
        res = AUD_USD_residual(y');
        [~,DnC,DnC_star,DnI,DnI_star,cvm_pvalues] = CvM_test_tar_AUD_USD(y',nboot);
    elseif i==2
        res = NZD_USD_residual(y');
        [~,DnC,DnC_star,DnI,DnI_star,cvm_pvalues] = CvM_test_tar_NZD_USD(y',nboot);
    else
        res = USD_CAD_residual(y');
        [~,DnC,DnC_star,DnI,DnI_star,cvm_pvalues] = CvM_test_tar_CAD_USD(y',nboot);
    end
    save(sprintf('new_new_%s_cvm_test.mat',indexes{i}))
end
%%

% clc
% clear all;
% indexes = {'AUD_USD','NZD_USD','USD_CAD'};
% 
% for i = 3
%     yyyy = readtable(sprintf('./data/%s.xlsx',indexes{i}));
%     yyy = yyyy(:,2);
%     yy = yyy{1:end,1};
%     yy = 1./wrev(yy);
%     y = (log(yy(2:end))-log(yy(1:end-1)))*100;
% %     y=y';
% %     X0=ones(length(y),1);
% %     threshold_variable=[0;y(1:end-1)'];
% %     X1 = [0;y(1:end-1)'];
% %     X2 = [0;0;y(1:end-2)'];
% %     X3 = [0;0;0;y(1:end-3)'];
% %     X4 = [0;0;0;0;y(1:end-4)'];
% %     X5 = [0;0;0;0;0;y(1:end-5)'];
% %     X6 = [0;0;0;0;0;0;y(1:end-6)'];
% %     X7 = [0;0;0;0;0;0;0;y(1:end-7)'];
% %     X8 = [0;0;0;0;0;0;0;0;y(1:end-8)'];
% %     X9 = [0;0;0;0;0;0;0;0;0;y(1:end-9)'];
% %     XX0 = X0.*(threshold_variable<0);
% %     XXX0 = X0.*(threshold_variable>=0);
% %     XX1 = X1.*(threshold_variable<0);
% %     XXX1 = X1.*(threshold_variable>=0);
% %     XX2 = X2.*(threshold_variable<0);
% %     XXX2 = X2.*(threshold_variable>=0);
% %     XX3 = X3.*(threshold_variable<0);
% %     XXX3 = X3.*(threshold_variable>=0);
% %     XX4 = X4.*(threshold_variable<0);
% %     XXX4 = X4.*(threshold_variable>=0);
% %     XX5 = X5.*(threshold_variable<0);
% %     XXX5 = X5.*(threshold_variable>=0);
% %     XX6 = X6.*(threshold_variable<0);
% %     XXX6 = X6.*(threshold_variable>=0);
% %     XX7 = X7.*(threshold_variable<0);
% %     XXX7 = X7.*(threshold_variable>=0);
% %     XX8 = X8.*(threshold_variable<0);
% %     XXX8 = X8.*(threshold_variable>=0);
% %     XX9 = X9.*(threshold_variable<0);
% %     XXX9 = X9.*(threshold_variable>=0);
% %     mdl=fitlm([XX4,XX5,XXX3,XXX6],y,'Intercept',false);
% %     mdl
%     nboot=100;
%     res = CAD_USD_residual(y');
%     [~,DnC,DnC_star,DnI,DnI_star,cvm_pvalues] = CvM_test_tar_CAD_USD(y',nboot);
% 
%     save(sprintf('new_new_CAD_USD_cvm_test.mat',indexes{i}))
% end