%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving a g.e. model  w/ fsolve 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
function [c_1,c_2,s,w,r,K,N,Y] = sol_GEsys_f(beta,sigma,A,alpha,delta,X0);
% Exog. variables: {?,?,A,?,?}
% Endo. variables: {c?,c?,s,w,r,K,L,Y}
options = optimset('Display', 'off'); % Turn off Display
fsolve(@cFOCs_f,X0,options);                    % call fsolve to solve the func csys_f
   function F = cFOCs_f(X)              % define the function F with X 
          % define variables
          %{c?,c?,s,w,r,K,L,Y}
          c_1       = X(1);               % consumption 1
          c_2       = X(2);               % consumption 2
          s         = X(3);               % savings
          w         = X(4);               % wage  
          r         = X(5);               % interest rate
          K         = X(6);               % Capital
          N         = X(7);               % Labor  
          Y         = X(8);               % Output
          %
          % system of equations  
          F(1)      = beta*c_2^(-sigma) - c_1^(-sigma)/(1+r);
          F(2)      = c_1 + c_2/(1+r)- w;
          F(3)      = s + c_1 - w;       
          F(4)      = w - (1-alpha)*A*K^alpha*N^(-alpha);
          F(5)      = r - alpha*A*K^(alpha-1)*N^(1-alpha) + delta;   
          F(6)      = N - 1;
          F(7)      = K - s;
          F(8)      = Y - A*K^alpha*N^(1-alpha);
   end
end


