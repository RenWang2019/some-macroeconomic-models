%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solving a g.e. model  w/ Gauss-Seidel algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all
tic
%disp('--- New run ------');
% Para
beta    = .98^30;
sigma   = 2; 
A       = 1;
alpha   = .33;
delta   = 1-(1-.05)^30;

% Initials
Kold    = 1;
Nold    = 1; 
w       = (1-alpha)*A*Kold^alpha*Nold^(-alpha);
R       = 1+alpha*A*Kold^(alpha-1)*Nold^(1-alpha)-delta;
%  
error   = 100;
errorv  = 100;
iter    = 0;
itermax = 50;
tol     = .001;
update  = .5;

while (iter<itermax)&(error>tol)
     %
     % Household
     lambda_sig=w/(1+(1/R)*(1/(R*beta))^(-1/sigma));
     %
     c_1    = lambda_sig;
     s      = w - c_1; 
     c_2    = lambda_sig*(1/(R*beta))^(-1/sigma);
     % Labor and capital markets
     N      = 1; 
     Knew   = s;
%      K = Knew;
     K      = update*Kold + (1-update)*Knew;   % convex update
     % Factor prices 
     w      = (1-alpha)*A*K^alpha*N^(-alpha);
     q      = alpha*A*K^(alpha-1)*N^(1-alpha);
     r      = q - delta;
     R      = 1 + r;
     % Output
     Y      = A*K^(alpha)*N^(1-alpha);	
     % update
     error  = 100*abs(K - Kold)/Kold;   % error in percentage
     errorv = [errorv error];
     Kold   = K;
     iter   = iter+1    
     
end
             
          
disp('---- Results---------------- ');
disp(['Y    =' num2str(Y)]);
disp(['K    =' num2str(K)]);
disp(['N    =' num2str(N)]);
disp(['R    =' num2str(R^(1/30))]);
disp(['w    =' num2str(w)]);
disp(['error=' num2str(error)]);
disp('-------------------------- ');


toc