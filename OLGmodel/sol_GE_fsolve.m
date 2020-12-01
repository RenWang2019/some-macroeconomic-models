clear all; close all

disp('--- New run ------');
% Para
% exog. var: {?,?,A,?,?}
beta    = .98^30;
sigma   = 2; 
A       = 1;
alpha   = .33;
delta   = 1-(1-.05)^30;
% Initials
% endo. var: {c1,c2,s,w,r,K,N,Y}
X0    = [.1,.1,.1,.1,1,1,1,1];        
%
% Call the f to solve for a ge model
[c_1,c_2,s,w,r,K,N,Y] = sol_GEsys_f(beta,sigma,A,alpha,delta,X0);
%
disp('---- Results---------------- ');
disp(['Y    =' num2str(Y)]);
disp(['K    =' num2str(K)]);
disp(['N    =' num2str(N)]);
disp(['R    =' num2str((1+r)^(1/30))]);
disp(['w    =' num2str(w)]);
disp('-------------------------- ');