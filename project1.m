clc 
clear all 
close all 
yourubitname='amoghsun';
yournumber=50134359;
load project1_data.mat
[w1,M_cfs,lambda_cfs, rms_cfs]=train_cfs(train_d,mu,sigma,train_t); 
test_cfs(validate,validation_t,test,test_t,mu,sigma,w1); 
[w2,M_gd,lambda_gd, rms_gd]=train_gd(train_d,mu,sigma,train_t); 
[M_gd,lambda_gd, rms_gd]=test_gd(test,test_t,validate,validation_t,mu,sigma,w2); 

fprintf('My ubit name is %s\n',yourubitname);
fprintf('My student number is %d\n',yournumber);
fprintf('the model complexity M_cfs is %d\n',M_cfs);
fprintf('the model complexity M_gd is %d\n',M_gd);
fprintf('the regularization parameters lambda_cfs is %4.2f\n',lambda_cfs);
fprintf('the regularization parameters lambda_gd is %4.2f\n',lambda_gd);
fprintf('the root mean square error for the closed form solution is %4.2f\n',rms_cfs);
fprintf('the root mean square error for the gradient descent method is %4.2f',rms_gd);

save W_cfs.mat w1;
save mu_cfs.mat mu;
save s_cfs.mat sigma;


save W_gd.mat w2;
save mu_gd.mat mu;
save s_gd.mat sigma;