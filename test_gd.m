function [ model, step_size,finalerror2 ] = test_gd(validate,validate_t,test,test_t,mu,sigma,wt)

mu_test=mu;
sigma_test=sigma;

%                          VALIDATION PHASE
model=20;

r=0:(0.3/(46*model)):0.3;
s=0:(0.5/(46*model)):0.5;
[n,m]=size(validate);

sigma=sigma.^2;
r=r(2:length(r));
s=s(2:length(s));

step_size=1;

validate=repmat(validate,1,model);
[n,m]=size(validate);
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

expo=zeros(n,m);
phi=zeros(n,m);

for i=1:n:n
    for j=1:m
        expo(i,j)=(validate(i,j)-mu(i,j)).^2;
        phi(i,j)=exp(-1*(expo(i,j)/(2*sigma(i,j))));
    end
end
step_size=1;
e=zeros(1,100);
loop=true;
while loop
    %    i = 200;
    
    wtnext = wt + step_size.* (phi'*(validate_t-(phi*wt)));
    
    %wterms
    validate_ts=phi*wt;
    err_sum=(sum((validate_t-validate_ts).^2))/2;
    err_sum=(2*err_sum/n);
    wterror=sqrt(err_sum);
    
    %wtnexterms
    validate_ts=phi*wtnext;
    err_sum=(sum((validate_t-validate_ts).^2))/2;
    err_sum=(2*err_sum/n);
    wtnexterror=sqrt(err_sum);
        
            if  (wterror>wtnexterror)
                %wt = wtnext;
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

finalerror1=wterror;

%                             TESTING

model=20;

r=0:(0.3/(46*model)):0.3;
s=0:(0.5/(46*model)):0.5;
[n,m]=size(test);

sigma_test=sigma_test.^2;
r=r(2:length(r));
s=s(2:length(s));

step_size=1;

test=repmat(test,1,model);
[n,m]=size(test);
%replicating the mean to model complexity times
mu_test=repmat(mu_test,1,model);
%mu=repmat(mu,model,1);
sigma_test=repmat(sigma_test,1,model);

[n1,m1]=size(mu_test);
n1;
m1;
I=eye(m,m);

sigma_test=sigma_test+s;
mu_test=mu_test+r;

mu_test=repmat(mu_test,n,1);
sigma_test=repmat(sigma_test,n,1);

expo=zeros(n,m);
phi=zeros(n,m);

for i=1:n:n
    for j=1:m
        expo(i,j)=(test(i,j)-mu_test(i,j)).^2;
        phi(i,j)=exp(-1*(expo(i,j)/(2*sigma_test(i,j))));
    end
end

step_size = 1;
loop=true;
while loop
    %    i = 200;
    
    wtnext = wt + step_size* (phi'*(test_t-(phi*wt)));
    
    %wterms
    test_ts=phi*wt;
    err_sum=(sum((test_t-test_ts).^2))/2;
    err_sum=(2*err_sum/n);
    wterror=sqrt(err_sum);
    
    %wtnexterms
    test_ts=phi*wtnext;
    err_sum=(sum((test_t-test_ts).^2))/2;
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

finalerror2=wterror;


%  
% %sprintf('My ubit name is %s\n',yourubitname);
% %sprintf('My student number is %d\n',yournumber);
% %sprintf('the model complexity M cfs is %d\n',M_cfs);
% %sprintf('the model complexity M gd is %d\n',M_gd);
% %sprintf('the regularization parameters lambda cfs is %4.2f\n',lambda_cfs);
% %sprintf('the regularization parameters lambda gd is %4.2f\n',lambda_gd);
% %sprintf('the root mean square error for the closed form solution is %4.2f\n',rms_cfs);
% %sprintf('the root mean square error for the gradient descent method is %4.2fnn',rms_gd);
% 
% sprintf('the model complexity M for the linear regression model is %d', model)
% sprintf('the regularization parameters lambda for the linear regression model is %f', lambda)
% sprintf('the root mean square error for the linear regression model  erms_test is %f', erms_t)
% sprintf('the root mean square error for the linear regression model is erms_validation%f', erms_va)
%  
end
