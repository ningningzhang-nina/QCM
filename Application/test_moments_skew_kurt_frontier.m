clear all
index_names = {'AUD_USD','NZD_USD','CAD_USD'};
for num =1:3
    load(sprintf("final_%s_moments.mat",index_names{num}))
    sigma0 = moments(1,:);
    skew0 = moments(2,:);
    kurt0 = moments(3,:);
    epsilon = res;
    results(num,:)=kurt0-skew0.^2-1;
end