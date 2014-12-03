function [ ] = test_cfs(validate,validation_t,test,test_t,mu,sigma,w)


%                          VALIDATION PHASE
model=10;

r=0:(0.3/(46*model)):0.3;
sigma=sigma.^2;
lambda=20;
s=0:(0.5/(46*model)):0.5;

[n,m]=size(validate);
n;
m;
r=r(2:length(r));

validate=repmat(validate,1,model);
[n,m]=size(validate);
%replicating the mean to model complexity times
mu_v=repmat(mu,1,model);
%mu_v=repmat(mu,model,1);
sigma_v=repmat(sigma,1,model);
I=eye(m,m);
[n1,m1]=size(mu_v);
n1;
m1;
I=eye(m,m);
s=s(2:length(s));

%creating the random number to be added to the mean
sigma_v=sigma_v+s;
mu_v=mu_v+r;

mu_v=repmat(mu_v,n,1);
sigma_v=repmat(sigma_v,n,1);

validation_m=zeros(n,m);      

expo_v=zeros(n,m);
phi_v=zeros(n,m);

for i=1:n
   for j=1:m  
       expo_v(i,j)=(validate(i,j)-mu_v(i,j)).^2;
       phi_v(i,j)=exp(-1*(expo_v(i,j)/(2*sigma_v(i,j))));
       %phi_v(i,j)=exp(-1*expo_v(i,j));
   end
end

target_v=phi_v*w;
esum=(sum((target_v-validation_t).^2))/2;
erms_va=sqrt(2*esum/n);

%                             TESTING

r=0:(0.3/(46*model)):0.3;
mu_t=repmat(mu,1,model);
%mu_v=repmat(mu,model,1);
[n,m]=size(test);
n;
m;
s=0:(0.5/(46*model)):0.5;
sigma_t=repmat(sigma,1,model);
test=repmat(test,1,model);
[n,m]=size(test);
%replicating the mean to model complexity times
r=r(2:length(r));
[n1,m1]=size(mu_t);
n1;
m1;
I=eye(m,m);
s=s(2:length(s));

%creating the random number to be added to the mean


sigma_t=sigma_t+s;
mu_t=mu_t+r;

mu_t=repmat(mu_t,n,1);
sigma_t=repmat(sigma_t,n,1);

test_m=zeros(n,m);
expo_t=zeros(n,m);
phi_t=zeros(n,m);

for i=1:n
   for j=1:m    

       %test_m(i,j)=test(i,j)-mu_t(i,j);
       expo_t(i,j)=(test(i,j)-mu_t(i,j)).^2;
       phi_t(i,j)=exp(-1*(expo_t(i,j)/(2*sigma_t(i,j))));

       %expo_t(i,j)=(((test(i,j)-mu_t(i,j)).^2)/(2*sigma_t(i,j)));
       %phi_t(i,j)=exp(expo_t(i,j)*-1);
   end
end
target_t=phi_t*w;

err_sum=(sum((target_t-test_t).^2))/2;
err_sum=(2*err_sum/n);
erms_t=sqrt(err_sum);
 


 
%sprintf('My ubit name is %s\n',yourubitname);
%sprintf('My student number is %d\n',yournumber);
%sprintf('the model complexity M cfs is %d\n',M_cfs);
%sprintf('the model complexity M gd is %d\n',M_gd);
%sprintf('the regularization parameters lambda cfs is %4.2f\n',lambda_cfs);
%sprintf('the regularization parameters lambda gd is %4.2f\n',lambda_gd);
%sprintf('the root mean square error for the closed form solution is %4.2f\n',rms_cfs);
%sprintf('the root mean square error for the gradient descent method is %4.2fnn',rms_gd);

% sprintf('the model complexity M for the linear regression model is %d', model)
% sprintf('the regularization parameters lambda for the linear regression model is %f', lambda)
% sprintf('the root mean square error for the linear regression model  erms_test is %f', erms_t)
% sprintf('the root mean square error for the linear regression model is erms_validation%f', erms_va)
 
end
