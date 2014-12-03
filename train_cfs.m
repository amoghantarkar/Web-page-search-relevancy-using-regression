function [ w,model, lambda, error] = train_cfs( train_d,mu,sigma,train_t )

    model=10;

    r=0:(0.3/(46*model)):0.3;
    s=0:(0.5/(46*model)):0.5;
    [n,m]=size(train_d);
    n;
    m;
    
    w=0;
    lambda=20;
    sigma=sigma.^2;
    r=r(2:length(r));
    s=s(2:length(s));
    
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
   
    for i=1:n
        for j=1:m
            expo(i,j)=(train_d(i,j)-mu(i,j)).^2;
            phi(i,j)=exp(-1*(expo(i,j)/(2*sigma(i,j))));
        end
    end
    
    w=pinv((phi'*phi+lambda*I))*phi'*train_t;
    %save W_cfs.mat
    %train_mean=train_d-mu;
    
    tar_main=phi*w;
    %finding the error.
    error=(tar_main-train_t).^2;
    err_sum=(sum(error))/2;
    error=sqrt(2*err_sum/n);
    i=i+1;
    

M_cfs=model;
save W_cfs.mat w;
save mu_cfs.mat mu;
save s_cfs.mat sigma;

end











