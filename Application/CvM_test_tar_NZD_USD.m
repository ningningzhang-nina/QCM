function [res,DnC,DnC_star,DnI,DnI_star,pvalues] = CvM_test_tar_NZD_USD(y,nboot)
X0=ones(length(y),1);
    threshold_variable=[0;y(1:end-1)'];
    X1 = [0;y(1:end-1)'];
    X2 = [0;0;y(1:end-2)'];
    X3 = [0;0;0;y(1:end-3)'];
    X4 = [0;0;0;0;y(1:end-4)'];
    X5 = [0;0;0;0;0;y(1:end-5)'];
    X6 = [0;0;0;0;0;0;y(1:end-6)'];
    X7 = [0;0;0;0;0;0;0;y(1:end-7)'];
    X8 = [0;0;0;0;0;0;0;0;y(1:end-8)'];
    X9 = [0;0;0;0;0;0;0;0;0;y(1:end-9)'];
    XX0 = X0.*(threshold_variable<0);
    XXX0 = X0.*(threshold_variable>=0);
    XX1 = X1.*(threshold_variable<0);
    XXX1 = X1.*(threshold_variable>=0);
    XX2 = X2.*(threshold_variable<0);
    XXX2 = X2.*(threshold_variable>=0);
    XX3 = X3.*(threshold_variable<0);
    XXX3 = X3.*(threshold_variable>=0);
    XX4 = X4.*(threshold_variable<0);
    XXX4 = X4.*(threshold_variable>=0);
    XX5 = X5.*(threshold_variable<0);
    XXX5 = X5.*(threshold_variable>=0);
    XX6 = X6.*(threshold_variable<0);
    XXX6 = X6.*(threshold_variable>=0);
    XX7 = X7.*(threshold_variable<0);
    XXX7 = X7.*(threshold_variable>=0);
    XX8 = X8.*(threshold_variable<0);
    XXX8 = X8.*(threshold_variable>=0);
    XX9 = X9.*(threshold_variable<0);
    XXX9 = X9.*(threshold_variable>=0);
    mdl=fitlm([XX3,XX4,XX8,XX9,XXX5],y,'Intercept',false);

x0=mdl.Coefficients.Estimate;
[x0,mdl.Coefficients.pValue]
res=y'-[XX3,XX4,XX8,XX9,XXX5]*x0;
tic;
[DnC,DnI]=EscancianoTest(y,res')
toc;
WB = WildBoot(res,nboot); 
DnC_star = zeros(1,nboot);
DnI_star = zeros(1,nboot);
for i=1:nboot
    tic;
    res_star = WB(:,i)';
    f_hat=[XX3,XX4,XX8,XX9,XXX5]*x0;
    y_star = f_hat'+res_star;
    mdl1=fitlm([XX3,XX4,XX8,XX9,XXX5],y_star,'Intercept',false);
    phat=mdl1.Coefficients.Estimate;
    [phat,mdl1.Coefficients.pValue]
    res_star_star=y_star'-[XX3,XX4,XX8,XX9,XXX5]*phat;
    [DnC_star(i),DnI_star(i)]=EscancianoTest(y,res_star_star');
    DnC_star(1:i)
    DnI_star(1:i)
    toc;
end
pvalues=[sum(DnC<DnC_star)/nboot,sum(DnI<DnI_star)/nboot];
end