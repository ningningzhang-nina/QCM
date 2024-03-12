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
    y1=sort(abs(y),"descend");
    
    cnt=0;
    for k=10:100
        cnt=cnt+1;
        hill(cnt,i)=1/(1/k*sum(log(y1(1:k))-log(y1(k+1))));

    end
end
% plot hill index
plot(10:100,hill(1:end,1:3))
legend({'AUD-USD','NZD-USD','USD-CAD'})

%% 
clear all

yyyy = readtable('./data/USD_CAD.xlsx');
yyy = yyyy(:,2);
yy = yyy{1:end,1};
yy = 1./wrev(yy);

y = (log(yy(2:end))-log(yy(1:end-1)))*100;
y1=sort(abs(y),"descend");

cnt=0;
for k=10:100
    cnt=cnt+1;
    hill(cnt)=1/(1/k*sum(log(y1(1:k))-log(y1(k+1))));

end
