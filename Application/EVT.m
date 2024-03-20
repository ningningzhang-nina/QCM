function yquantile_hat = EVT(returns,p)
Mdl=garch(1,1);
EstMdl=estimate(Mdl,returns);
Mdl=garch(constant=EstMdl.Constant,GARCH=EstMdl.GARCH{1},ARCH=EstMdl.ARCH{1});
V=infer(Mdl,returns);
zt=returns./sqrt(V);
VF = forecast(EstMdl,1,returns);

newzt=sort(zt,'descend');
k=100;
xi_hat=1/k*sum(log(newzt(1:k)))-log(newzt(k+1));
yquantile_hat=zeros(3,1);
cnt=0;
for q=p
    cnt=cnt+1;
    zquantile_hat=newzt(k+1)*((1-q)/k*length(zt))^(-xi_hat);
    yquantile_hat(cnt,:)=sqrt(VF).*zquantile_hat;
end
% returns=-returns;
% Mdl=garch(1,1);
% EstMdl=estimate(Mdl,returns);
% Mdl=garch(constant=EstMdl.Constant,GARCH=EstMdl.GARCH{1},ARCH=EstMdl.ARCH{1});
% V=infer(Mdl,returns);
% zt=returns./sqrt(V);
% 
% newzt=sort(zt,'descend');
% k=100;
% xi_hat=1/k*sum(log(newzt(1:k)))-log(newzt(k+1));
% cnt=50;
% for q=0.49:-0.01:0.01
%     cnt=cnt+1;
%     zquantile_hat=newzt(k+1)*((q)/k*length(zt))^(-xi_hat);
%     yquantile_hat(cnt,:)=-sqrt(V).*zquantile_hat;
% end
end
