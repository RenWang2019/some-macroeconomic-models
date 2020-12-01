
clear all; close all

disp('--- Solve G.E. Model ------');
% Parameters: {beta,sigma,A,alpha,delta} - exog. var
beta    = .99^30;
sigma   = 2; 
A       = 1.0;
alpha   = .33;
delta   = 1-(1-.05)^30;

% endo. var: {c1,c1,s,w,r,K,L,Y}
X0    = [.1,.1,.1,.1,1,.1,1,.5];        
%
% Call the function to solve for the ge model
[c_1,c_2,s,w,r,K,L,Y] = sol_GEmodel_sys_f(beta,sigma,A,alpha,delta,X0);


disp('---- Results---------------- ');
disp(['c1    =' num2str(c_1)]);
disp(['c2    =' num2str(c_2)]);
disp(['s    =' num2str(s)]);
disp(['Y    =' num2str(Y)]);
disp(['K    =' num2str(K)]);
disp(['L    =' num2str(L)]);
disp(['w    =' num2str(w)]);
disp(['R    =' num2str((1+r)^(1/30))]);

disp('-------------------------- ');