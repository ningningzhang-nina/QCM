function [DnC,DnI]=EscancianoTest(y,res)
[~,n]=size(y);
sigma_e=sum(res.^2)/n;
DnC=0;
DnI=0;
for j=1:n-1
    
    fenmu=sigma_e*j^2*pi^2*(n-j+1);
    res1=res(j+1:n);
    a=(repmat(res1,length(res1),1))'.*(repmat(res1,length(res1),1));
    y1=y(1:n-j);
    b=-((repmat(y1,length(y1),1))'-(repmat(y1,length(y1),1))).^2/2;
    b=exp(b);
    DnC=DnC+sum(sum(a.*b))/fenmu;
    b1=repmat(y,length(y1),1)';
    a1=repmat(y1,n,1);
   ind=(a1<=b1);
   a0=ind*res1';
   DnI=DnI+sum(a0.^2)/fenmu/n;
   
end
