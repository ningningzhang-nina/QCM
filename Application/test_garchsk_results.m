clc
clear all
results = readtable('garchsk_AUD_USD_results.csv');
%%
y = results.y;
sigma0 = sqrt(results.h);
%size((sigma0').^2)
%size((y).^2)
[p2,~] = waldtest_sd((sigma0).^2,(y).^2);
% test skewness and kurtosis
skew0 = results.sk;
kurt0 = results.ku;
[p3,~,p4,~] = waldtest_skewness_kurtosis((skew0),(y).^3./((sigma0).^3),...
                 ((kurt0)),(y).^4./((sigma0).^4));

% p values of four wald test
pvalues = [p2,p3,p4]
%%
clc
clear all
results = readtable('garchsk_NZD_USD_results.csv');
%%
y = results.y;
sigma0 = sqrt(results.h);
%size((sigma0').^2)
%size((y).^2)
[p2,~] = waldtest_sd((sigma0).^2,(y).^2);
% test skewness and kurtosis
skew0 = results.sk;
kurt0 = results.ku;
[p3,~,p4,~] = waldtest_skewness_kurtosis((skew0),(y).^3./((sigma0).^3),...
                 ((kurt0)),(y).^4./((sigma0).^4));

% p values of four wald test
pvalues = [p2,p3,p4]

%%
clc
clear all
results = readtable('garchsk_CAD_USD_results.csv');
%%
y = results.y;
sigma0 = sqrt(results.h);
%size((sigma0').^2)
%size((y).^2)
[p2,~] = waldtest_sd((sigma0).^2,(y).^2);
% test skewness and kurtosis
skew0 = results.sk;
kurt0 = results.ku;
[p3,~,p4,~] = waldtest_skewness_kurtosis((skew0),(y).^3./((sigma0).^3),...
                 ((kurt0)),(y).^4./((sigma0).^4));

% p values of four wald test
pvalues = [p2,p3,p4]