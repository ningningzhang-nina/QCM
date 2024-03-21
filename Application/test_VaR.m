clc
clear all;
indexes={'AUD_USD','NZD_USD','USD_CAD'};
%indexes={'Copper','USD_CAD','NZD_USD','AUD_USD'};
%indexes = {'EUR_USD','GBP_USD','USD_JPY','USD_CHF','AUD_USD','EUR_GBP','USD_CAD','NZD_USD','EUR_JPY','GBP_JPY'};
%dqtest_in = zeros(4,99);
%hill=zeros(200,10);
for i =1:length(indexes)
    yyyy = readtable(sprintf('./data/%s.xlsx',indexes{i}));
    yyy = yyyy(:,2);
    yy = yyy{1:end,1};
    yy = wrev(yy);

    y = (log(yy(2:end))-log(yy(1:end-1)))*100;
    s=201;
    for j=s:500
        load(sprintf('./results/final_%s_moments_%d.mat',indexes{i},j))
        newy=y(1:end-500+j)';
        X0=ones(length(newy),1);
        threshold_variable=[0;newy(1:end-1)'];
        X1 = [0;newy(1:end-1)'];
        X2 = [0;0;newy(1:end-2)'];
        X3 = [0;0;0;newy(1:end-3)'];
        X4 = [0;0;0;0;newy(1:end-4)'];
        X5 = [0;0;0;0;0;newy(1:end-5)'];
        X6 = [0;0;0;0;0;0;newy(1:end-6)'];
        X7 = [0;0;0;0;0;0;0;newy(1:end-7)'];
        X8 = [0;0;0;0;0;0;0;0;newy(1:end-8)'];
        X9 = [0;0;0;0;0;0;0;0;0;newy(1:end-9)'];
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
        allvariables=[XX0,XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XXX0,XXX1,XXX2,XXX3,XXX4,XXX5,XXX6,XXX7,XXX8,XXX9];
        mdl=fitlm([XX0,XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XXX0,XXX1,XXX2,XXX3,XXX4,XXX5,XXX6,XXX7,XXX8,XXX9],newy,'Intercept',false);
        mdl
        newy=newy';
        newmdl=fitlm(allvariables(:,mdl.Coefficients.pValue<0.10),newy,'Intercept',false)
        x0=newmdl.Coefficients.Estimate;
        [x0,newmdl.Coefficients.pValue]
        var=[newy(end)<0,newy(end)>=0,newy(end-1)*(newy(end)<0),newy(end-1)*(newy(end)>=0),...
            newy(end-2)*(newy(end)<0),newy(end-2)*(newy(end)>=0),...
            newy(end-3)*(newy(end)<0),newy(end-3)*(newy(end)>=0),...
            newy(end-4)*(newy(end)<0),newy(end-4)*(newy(end)>=0),...
            newy(end-5)*(newy(end)<0),newy(end-5)*(newy(end)>=0),...
            newy(end-6)*(newy(end)<0),newy(end-6)*(newy(end)>=0),...
            newy(end-7)*(newy(end)<0),newy(end-7)*(newy(end)>=0),...
            newy(end-8)*(newy(end)<0),newy(end-8)*(newy(end)>=0),...
            newy(end-9)*(newy(end)<0),newy(end-9)*(newy(end)>=0)];
        mu=var(1,mdl.Coefficients.pValue<0.10)*x0;
        
        cnt=0;
        for p=[0.05,0.04,0.03,0.025,0.02,0.01]
            cnt=cnt+1;
            VaR3(cnt,j-s+1) = -mu+moments(2,:).*(-norminv(p)-1/6*(norminv(p)^2-1)*moments(3,:)-1/24*(norminv(p)^3-3*norminv(p))*(moments(4,:)-3)+1/36*(2*norminv(p)^3-5*norminv(p))*moments(3,:).^2);
        end
        VaR1(1:6,j-s+1) = FHS(-y(1:end-500+j),[0.95,0.96,0.97,0.975,0.98,0.99]);
        VaR2(1:6,j-s+1) = EVT(-y(1:end-500+j),[0.95,0.96,0.97,0.975,0.98,0.99]);
    end
    cnt=0;
    for p=[0.95,0.96,0.97,0.975,0.98,0.99]
        cnt=cnt+1;
        vbt1 = varbacktest(y(end-(500-s):end),VaR1(cnt,:)','VaRLevel',p);
        Testresults11=tl(vbt1);
        Testresults12=bin(vbt1);
        Testresults13=pof(vbt1);
        Testresults14=tuff(vbt1);
        Testresults15=cc(vbt1);
        Testresults16=cci(vbt1);
        Testresults17=tbf(vbt1);
        Testresults18=tbfi(vbt1);
        vbt2 = varbacktest(y(end-(500-s):end),VaR2(cnt,:)','VaRLevel',p);
        Testresults21=tl(vbt2);
        Testresults22=bin(vbt2);
        Testresults23=pof(vbt2);
        Testresults24=tuff(vbt2);
        Testresults25=cc(vbt2);
        Testresults26=cci(vbt2);
        Testresults27=tbf(vbt2);
        Testresults28=tbfi(vbt2);
        vbt3 = varbacktest(y(end-(500-s):end),VaR3(cnt,:)','VaRLevel',p);
        Testresults31=tl(vbt3);
        Testresults32=bin(vbt3);
        Testresults33=pof(vbt3);
        Testresults34=tuff(vbt3);
        Testresults35=cc(vbt3);
        Testresults36=cci(vbt3);
        Testresults37=tbf(vbt3);
        Testresults38=tbfi(vbt3);
        pvalues_FHS_tl(cnt,1)=Testresults11.TypeI;
        pvalues_mVaR_tl(cnt,1)=Testresults21.TypeI;
        pvalues_EVT_tl(cnt,1)=Testresults31.TypeI;
        pvalues_FHS_bin(cnt,1)=Testresults12.PValueBin;
        pvalues_mVaR_bin(cnt,1)=Testresults22.PValueBin;
        pvalues_EVT_bin(cnt,1)=Testresults32.PValueBin;
        pvalues_FHS_pof(cnt,1)=Testresults13.PValuePOF;
        pvalues_mVaR_pof(cnt,1)=Testresults23.PValuePOF;
        pvalues_EVT_pof(cnt,1)=Testresults33.PValuePOF;
        pvalues_FHS_tuff(cnt,1)=Testresults14.PValueTUFF;
        pvalues_mVaR_tuff(cnt,1)=Testresults24.PValueTUFF;
        pvalues_EVT_tuff(cnt,1)=Testresults34.PValueTUFF;
        pvalues_FHS_cc(cnt,1)=Testresults15.PValueCC;
        pvalues_mVaR_cc(cnt,1)=Testresults25.PValueCC;
        pvalues_EVT_cc(cnt,1)=Testresults35.PValueCC;
        pvalues_FHS_cci(cnt,1)=Testresults16.PValueCCI;
        pvalues_mVaR_cci(cnt,1)=Testresults26.PValueCCI;
        pvalues_EVT_cci(cnt,1)=Testresults36.PValueCCI;
        pvalues_FHS_tbf(cnt,1)=Testresults17.PValueTBF;
        pvalues_mVaR_tbf(cnt,1)=Testresults27.PValueTBF;
        pvalues_EVT_tbf(cnt,1)=Testresults37.PValueTBF;
        pvalues_FHS_tbfi(cnt,1)=Testresults18.PValueTBFI;
        pvalues_mVaR_tbfi(cnt,1)=Testresults28.PValueTBFI;
        pvalues_EVT_tbfi(cnt,1)=Testresults38.PValueTBFI;
    end
    alpha=[0.95;0.96;0.97;0.975;0.98;0.99];
    %T1=table(alpha,pvalues_FHS_tl,pvalues_EVT_tl,pvalues_mVaR_tl)
    %T2=table(alpha,pvalues_FHS_bin,pvalues_EVT_bin,pvalues_mVaR_bin)
    %T3=table(alpha,pvalues_FHS_pof,pvalues_EVT_pof,pvalues_mVaR_pof)
    %T4=table(alpha,pvalues_FHS_tuff,pvalues_EVT_tuff,pvalues_mVaR_tuff)
    T5=table(alpha,pvalues_FHS_cc,pvalues_EVT_cc,pvalues_mVaR_cc)
    %T6=table(alpha,pvalues_FHS_cci,pvalues_EVT_cci,pvalues_mVaR_cci)
    %T7=table(alpha,pvalues_FHS_tbf,pvalues_EVT_tbf,pvalues_mVaR_tbf)
    %T8=table(alpha,pvalues_FHS_tbfi,pvalues_EVT_tbfi,pvalues_mVaR_tbfi)
end
% data1
% T5 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cc    pvalues_EVT_cc    pvalues_mVaR_cc
%     _____    ______________    ______________    _______________
% 
%      0.95       0.038325           0.7398             0.3547    
%      0.96      0.0051965          0.60539            0.62823    
%      0.97       0.019289           0.5302            0.75627    
%     0.975        0.10089          0.65436            0.74987    
%      0.98        0.38911          0.77991            0.83935    
%      0.99         0.8133          0.52418           0.049041   
% data2
% T5 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cc    pvalues_EVT_cc    pvalues_mVaR_cc
%     _____    ______________    ______________    _______________
% 
%      0.95       0.29296            0.58216           0.58216    
%      0.96       0.66756            0.15909           0.77716    
%      0.97       0.58698           0.099902           0.52702    
%     0.975         0.456            0.34027           0.83085    
%      0.98       0.31755            0.88437           0.88437    
%      0.99       0.97005             0.8133           0.81528  
% data3
% T5 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cc    pvalues_EVT_cc    pvalues_mVaR_cc
%     _____    ______________    ______________    _______________
% 
%      0.95       0.93958           0.47785            0.50089    
%      0.96       0.58872           0.33244            0.55601    
%      0.97       0.66002           0.25854            0.69332    
%     0.975       0.74987           0.49734            0.74987    
%      0.98       0.88437           0.77991            0.83935    
%      0.99        0.8133           0.52418            0.81528 
% %%%%%%%%%%%%%%%%%%%%% data1
% T1 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tl    pvalues_EVT_tl    pvalues_mVaR_tl
%     _____    ______________    ______________    _______________
% 
%      0.95       0.52865           0.29612            0.96571   % 
%      0.96       0.20838           0.62204             0.7842   % 
%      0.97        0.3344           0.53591            0.73618   % 
%     0.975       0.48174           0.27311            0.80195   % 
%      0.98       0.66946           0.54334             0.8724   % 
%      0.99       0.73638           0.56039                  1   % 
% 
% 
% T2 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_bin    pvalues_EVT_bin    pvalues_mVaR_bin
%     _____    _______________    _______________    ________________
% 
%      0.95              1            0.53817             0.10068    
%      0.96        0.36131            0.81948             0.49356    
%      0.97        0.60005                  1             0.60005    
%     0.975        0.88611            0.47392             0.47392    
%      0.98        0.74939                  1              0.3379    
%      0.99         0.6531                  1            0.024619    
% 
% 
% T3 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_pof    pvalues_EVT_pof    pvalues_mVaR_pof
%     _____    _______________    _______________    ________________
% 
%      0.95              1            0.54553            0.082169    
%      0.96        0.37566            0.81803             0.48256    
%      0.97        0.60752                  1             0.59175    
%     0.975        0.88684             0.4873             0.45835    
%      0.98        0.74527                  1             0.31136    
%      0.99        0.64143                  1           0.0015232    
% 
% 
% T4 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tuff    pvalues_EVT_tuff    pvalues_mVaR_tuff
%     _____    ________________    ________________    _________________
% 
%      0.95         0.52041            0.52041              0.52041     %
%      0.96        0.034754            0.41274              0.41274     %
%      0.97       0.0074612            0.30328              0.30328     %
%     0.975        0.020682            0.24848              0.24848     %
%      0.98        0.032805            0.19412              0.31285     %
%      0.99         0.12039            0.86421             0.027977     
% 
% 
% T5 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cc    pvalues_EVT_cc    pvalues_mVaR_cc
%     _____    ______________    ______________    _______________
% 
%      0.95       0.26736           0.73004             0.12118   
%      0.96       0.20077           0.45898             0.42902   
%      0.97       0.48113           0.62817             0.61166   
%     0.975       0.69914           0.49352             0.61911   
%      0.98       0.80405           0.81504             0.54225   
%      0.99        0.8687           0.95065           0.0065705   
% 
% 
% T6 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cci    pvalues_EVT_cci    pvalues_mVaR_cci
%     _____    _______________    _______________    ________________
% 
%      0.95        0.10431            0.60744            0.27344     
%      0.96        0.11931            0.21997            0.27344     %
%      0.97        0.27344            0.33489            0.40428     %
%     0.975        0.40428            0.33489            0.52246     %
%      0.98        0.56529            0.52246            0.65537     %
%      0.99         0.7993            0.75037                  1     %
% 
% 
% T7 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tbf    pvalues_EVT_tbf    pvalues_mVaR_tbf
%     _____    _______________    _______________    ________________
% 
%      0.95        0.28259             0.5482              0.70048   %
%      0.96        0.15016            0.83271              0.93039   %
%      0.97        0.20386            0.98217              0.98218   %
%     0.975        0.33038            0.95557              0.82521   
%      0.98        0.41865            0.88017              0.80435   
%      0.99        0.41064            0.95447           0.00011458   
% 
% 
% T8 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tbfi    pvalues_EVT_tbfi    pvalues_mVaR_tbfi
%     _____    ________________    ________________    _________________
% 
%      0.95        0.23787             0.51441              0.83441    % 
%      0.96         0.1401             0.78931               0.9213    % 
%      0.97        0.16984             0.97072              0.97622    % 
%     0.975        0.26538             0.94935              0.80515     
%      0.98        0.33874             0.82387              0.83249     %
%      0.99         0.3056             0.90444             0.027977  
%%%%%%%%%%%%%%%%%%%%% data2

% T1 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tl    pvalues_EVT_tl    pvalues_mVaR_tl
%     _____    ______________    ______________    _______________
% 
%      0.95       0.52865            0.36861           0.36861    
%      0.96       0.89983            0.10764           0.53156    
%      0.97       0.93307           0.079797           0.81928    
%     0.975       0.87821            0.48174           0.87821    %
%      0.98       0.78251            0.78251            0.8724    %
%      0.99       0.87661            0.96025           0.96025    %
% 
% 
% T2 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_bin    pvalues_EVT_bin    pvalues_mVaR_bin
%     _____    _______________    _______________    ________________
% 
%      0.95              1            0.68152            0.68152     
%      0.96        0.25383             0.1709                  1     %
%      0.97        0.18992            0.11573            0.43158     %
%     0.975        0.31607            0.88611            0.31607     
%      0.98         0.5229             0.5229             0.3379     
%      0.99        0.36869            0.17753            0.17753     
% 
% 
% T3 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_pof    pvalues_EVT_pof    pvalues_mVaR_pof
%     _____    _______________    _______________    ________________
% 
%      0.95              1             0.6852             0.6852     
%      0.96        0.23317            0.18992                  1     %
%      0.97        0.16344            0.13744            0.41548     %
%     0.975        0.29165            0.88684            0.29165     
%      0.98        0.50816            0.50816            0.31136     
%      0.99        0.33148            0.12504            0.12504     
% 
% 
% T4 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tuff    pvalues_EVT_tuff    pvalues_mVaR_tuff
%     _____    ________________    ________________    _________________
% 
%      0.95          0.52041           0.52041              0.52041    % 
%      0.96       0.00030056           0.41274              0.41274    % 
%      0.97        0.0024969           0.30328              0.47095    % 
%     0.975        0.0083955            0.6193               0.6193     %
%      0.98         0.027709           0.80829              0.80829     %
%      0.99          0.12039           0.12039              0.12039     %
% 
% 
% T5 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cc    pvalues_EVT_cc    pvalues_mVaR_cc
%     _____    ______________    ______________    _______________
% 
%      0.95       0.97034           0.84154            0.84154    
%      0.96       0.37552           0.40033            0.97546    %
%      0.97       0.15777           0.32916            0.40351    %
%     0.975       0.19805            0.6268            0.48611    
%      0.98       0.22268            0.7052            0.54225    
%      0.99       0.61281           0.30589            0.30589    
% 
% 
% T6 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cci    pvalues_EVT_cci    pvalues_mVaR_cci
%     _____    _______________    _______________    ________________
% 
%      0.95        0.80616            0.67075            0.67075     
%      0.96         0.4635            0.73711             0.8236    % 
%      0.97        0.18573            0.89924            0.28312     
%     0.975        0.14477            0.33905            0.56529     %
%      0.98        0.10917            0.60964            0.65537     %
%      0.99        0.84892            0.89904            0.89904     %
% 
% 
% T7 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tbf    pvalues_EVT_tbf    pvalues_mVaR_tbf
%     _____    _______________    _______________    ________________
% 
%      0.95        0.24664            0.73967            0.73967     %
%      0.96       0.033274             0.5583            0.65038     %
%      0.97       0.013711            0.54127            0.15894     
%     0.975       0.017406            0.34185             0.3641     %
%      0.98       0.013684            0.74135            0.46673     
%      0.99        0.21725            0.18975            0.18975     
% 
% 
% T8 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tbfi    pvalues_EVT_tbfi    pvalues_mVaR_tbfi
%     _____    ________________    ________________    _________________
% 
%      0.95         0.20542            0.70089              0.70089     %
%      0.96        0.034092            0.60053              0.58884     
%      0.97        0.015973            0.61957              0.13867     
%     0.975        0.015246            0.27567              0.36658     %
%      0.98       0.0091302            0.69774              0.46714     
%      0.99         0.18517            0.29923              0.29923     %
%%%%%%%%%%%%%%%%%% data3
% T1 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tl    pvalues_EVT_tl    pvalues_mVaR_tl
%     _____    ______________    ______________    _______________
% 
%      0.95       0.82115           0.61007            0.82115    %
%      0.96       0.89983            0.4409            0.53156    
%      0.97       0.98323           0.17728            0.73618    
%     0.975       0.96703           0.48174            0.87821    
%      0.98        0.8724           0.66946             0.8724    %
%      0.99       0.73638           0.56039            0.87661    %
% 
% 
% T2 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_bin    pvalues_EVT_bin    pvalues_mVaR_bin
%     _____    _______________    _______________    ________________
% 
%      0.95        0.41177            0.83742            0.41177     
%      0.96        0.25383            0.81948                  1     %
%      0.97       0.066487            0.29434            0.60005     %
%     0.975        0.11515            0.88611            0.31607     
%      0.98         0.3379            0.74939             0.3379     
%      0.99         0.6531                  1            0.36869     
% 
% 
% T3 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_pof    pvalues_EVT_pof    pvalues_mVaR_pof
%     _____    _______________    _______________    ________________
% 
%      0.95         0.3992            0.83639             0.3992     
%      0.96        0.23317            0.82087                  1     %
%      0.97        0.04436             0.3135            0.59175     %
%     0.975       0.086178            0.88684            0.29165     
%      0.98        0.31136            0.74527            0.31136     
%      0.99        0.64143                  1            0.33148     
% 
% 
% T4 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tuff    pvalues_EVT_tuff    pvalues_mVaR_tuff
%     _____    ________________    ________________    _________________
% 
%      0.95        0.92091             0.92091              0.92091     %
%      0.96        0.82012             0.89825              0.89825     %
%      0.97        0.42395             0.69242              0.94194     %
%     0.975        0.56843              0.8042               0.8042     %
%      0.98        0.75617             0.65404              0.65404     
%      0.99        0.10596             0.70571              0.70571     %
% 
% 
% T5 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cc    pvalues_EVT_cc    pvalues_mVaR_cc
%     _____    ______________    ______________    _______________
% 
%      0.95       0.69532           0.96735            0.27852    
%      0.96       0.30861           0.38731            0.43374    %
%      0.97       0.11628            0.2836            0.61166    %
%     0.975       0.20768           0.69914            0.48611    
%      0.98       0.54225           0.80405            0.54225    
%      0.99        0.8687           0.95065            0.61281    
% 
% 
% T6 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_cci    pvalues_EVT_cci    pvalues_mVaR_cci
%     _____    _______________    _______________    ________________
% 
%      0.95        0.89924            0.87753            0.17427     
%      0.96        0.33489            0.17427            0.19617     
%      0.97        0.60964            0.21997            0.40428     
%     0.975        0.65537            0.40428            0.56529     
%      0.98        0.65537            0.56529            0.65537     %
%      0.99         0.7993            0.75037            0.84892     %
% 
% 
% T7 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tbf    pvalues_EVT_tbf    pvalues_mVaR_tbf
%     _____    _______________    _______________    ________________
% 
%      0.95        0.28209             0.7481            0.86529     %
%      0.96        0.56462            0.91131            0.88471     
%      0.97         0.0419            0.66186            0.65109     
%     0.975       0.020018            0.56772            0.64405     %
%      0.98       0.057151            0.76178             0.3654     
%      0.99       0.086294            0.87461             0.4558     
% 
% 
% T8 =
% 
%   6×4 table
% 
%     alpha    pvalues_FHS_tbfi    pvalues_EVT_tbfi    pvalues_mVaR_tbfi
%     _____    ________________    ________________    _________________
% 
%      0.95         0.26455            0.70057              0.86007     %
%      0.96         0.59946            0.88364              0.84845     
%      0.97        0.098095            0.66746              0.59796     
%     0.975        0.033265            0.49026              0.66497     %
%      0.98        0.049789            0.68868              0.35912     
%      0.99        0.051477            0.78481              0.43947 
