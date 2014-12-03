function [ wt, model, step_size, finalerror] = train_gd( train_d,mu,sigma,train_t )
model=20;

r=0:(0.3/(46*model)):0.3;
s=0:(0.5/(46*model)):0.5;
[n,m]=size(train_d);
w1=0;

sigma=sigma.^2;
r=r(2:length(r));
s=s(2:length(s));

step_size=1;

train_d=repmat(train_d,1,model);
[n,m]=size(train_d);
%replicating the mean to model complexity times
mu=repmat(mu,1,model);
%mu=repmat(mu,model,1);
sigma=repmat(sigma,1,model);

[n1,m1]=size(mu);
n1;
m1;
I=eye(m,m);

sigma=sigma+s;
mu=mu+r;

mu=repmat(mu,n,1);
sigma=repmat(sigma,n,1);

train_m=zeros(n,m);
expo=zeros(n,m);
phi=zeros(n,m);

for i=1:n:n
    for j=1:m
        expo(i,j)=(train_d(i,j)-mu(i,j)).^2;
        phi(i,j)=exp(-1*(expo(i,j)/(2*sigma(i,j))));
    end
end

wt=0.01.*ones((46*model),1);
step = 1;
loop=0;
e=zeros(1,100);
loop=true;
while loop
    %    i = 200;
    
    wtnext = wt - step_size* (phi'*(train_t-(phi*wt)));
    
    %wterms
    train_ts=phi*wt;
    err_sum=(sum((train_t-train_ts).^2))/2;
    err_sum=(2*err_sum/n);
    wterror=sqrt(err_sum);
    
    %wtnexterms
    train_ts=phi*wtnext;
    err_sum=(sum((train_t-train_ts).^2))/2;
    err_sum=(2*err_sum/n);
    wtnexterror=sqrt(err_sum);
    
    
            if  (wterror>wtnexterror)
                wt = wtnext;
                %i=i+1;
                if wterror-wtnexterror<0.01
                    loop=false;
                    break;
                end
            else
                break;
                %i=1;
            end
    
end

finalerror=wterror;

%w1=pinv((phi'*phi+lambda*I))*phi'*train_t;
%save W_cfs.mat
%train_mean=train_d-mu;

save W_gd.mat wt;
save mu_gd.mat mu;
save s_gd.mat sigma;


end











