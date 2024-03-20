clear all
clc
indexes = {'AUD_USD','NZD_USD','USD_CAD'};
results = zeros(5,3);
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
    n = length(y);
    m = mean(y);
    vars = var(y);
    s = skewness(y);
    k = kurtosis(y);

    results(:,i) = [n;m;vars;s;k];
end
%%
clear all
indexe = {'AUD_USD','NZD_USD','CAD_USD'};

for ii = 1:3
    load(sprintf('final_%s_moments.mat',indexe{ii}))
    n = length(y);
    m = mean(y);
    vars = var(y);
    s = skewness(y);
    k = kurtosis(y);
    
    results(:,ii) = [n;m;vars;s;k];
    h = moments(1,:).^2;
    s = moments(2,:);
    k = moments(3,:);
    [~,p1,jbstat1,~] = lbqtest(h);
    [~,p2,jbstat2,~] = lbqtest(s);
    [~,p3,jbstat3,~] = lbqtest(k);
    r1=iqr(h);
    r2=iqr(s);
    r3=iqr(k);
    h_results(1:4,ii) = [mean(h),max(h),min(h),p1];
    s_results(1:4,ii) = [mean(s),max(s),min(s),p2];
    k_results(1:4,ii) = [mean(k),max(k),min(k),p3];
end



