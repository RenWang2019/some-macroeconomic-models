%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the Aiyagari-type model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear
close all; clear all;

tic
diary model2.out; 
disp(' --------- Aiyagari Type Model ------------------');
disp('');

% t1 = clock; 

%--------------------------------------------------%
% [1.] Parameter values
%--------------------------------------------------%
%
sigma  = 2.0;            % risk aversion              
beta   = 0.9;            % subjective discount factor 
prob   = [.7 .3; .2 .8];  ; % prob(i,j) = probability (s(t+1)=sj | s(t) = si)
delta  = 0.975;            % 1 - depreciation
A      = 1.00;            % production technology
alpha  = 0.36;            % capital's share of income
theta  = 0.05;            % non-rental income if unemployed is theta*wage
Kstart = 3.0;            % initial value for aggregate capital stock
g      = 0.20;            % relaxation parameter
gzv    = [.1 .5];             % level of shock 
%   calculate aggregate labor supply
gzn = 2;
gPv =1/gzn*ones(1,gzn);
for i=1:30; gPv=gPv*prob; end
N = gPv*gzv';
%   form capital grid
%   
maxkap = 20;                      % maximum value of capital grid   
inckap = 0.01;                   % size of capital grid increments
nkap   = round(maxkap/inckap+1);  % number of grid points
% %

%--------------------------------------------------%
%  Search  for equi aggregate capital stock
%--------------------------------------------------%
K       = Kstart;
Kold    = Kstart;
%
iter    = 1;
maxiter = 50;
toler   = 10.0^(-3);
errorK  = 100;

disp('------------------------------------------');
disp('    iter     erorrK     meanK      Kold');

while  (errorK>toler) & (iter<=maxiter);
   %
   % -------------------------%
   %  calculate rental and wage rates
   % -------------------------%
   wage = (1-alpha)*A*K^(alpha)*N^(-alpha);
   rent = (alpha)*A*K^(alpha-1)*N^(1-alpha);
   % 
   %  tabulate the utility function such that for zero or negative
   %  consumption utility remains a large negative number so that
   %  such values will never be chosen as utility maximizing      
   %
   util1=-10000*ones(nkap,nkap);  % utility when employed     
   util2=-10000*ones(nkap,nkap);  % utility when unemployed   
   for i=1:nkap;
         kap=(i-1)*inckap;
         for j=1:nkap; 
               kapp = (j-1)*inckap;
               cons1 = gzv(1)*wage + (rent+delta)*kap - kapp; 
               if cons1>0;
                  util1(j,i)=(cons1)^(1-sigma)/(1-sigma);
               end;
               cons2 = gzv(2)*wage + (rent+delta)*kap - kapp;
               if cons2 > 0;
                  util2(j,i)=(cons2)^(1-sigma)/(1-sigma);
               end;
         end;
   end;
   % -------------------------%
   % ----Find decision rules
   % -------------------------%
   %
   %  initialize some variables
   %
   v       = zeros(nkap,2);
   decis   = zeros(nkap,2);
   test    = 10;
   [rs,cs] = size(util1);
   %
   %  iterate on Bellman's equation and get the decision 
   %  rules and the value function at the optimum         
   %
   while test > 10^-8
       for i=1:cs;
           r1(:,i)= util1(:,i)+ beta*(prob(1,1)*v(:,1)+ prob(1,2)*v(:,2));
           r2(:,i)= util2(:,i)+ beta*(prob(2,1)*v(:,1)+ prob(2,2)*v(:,2));
       end;

       [tv1,tdecis1]= max(r1);
       [tv2,tdecis2]= max(r2);
       tdecis       = [tdecis1' tdecis2'];
       tv           = [tv1' tv2'];
       %
       test         = max(any(tdecis-decis));
       v            = tv;
       decis        = tdecis;
   end;
   decis=(decis-1)*inckap;
   %
   % -------------------------------%
   % [.]-Find a stationary distribution
   % -------------------------------%
   %   form transition matrix
   %   trans is the transition matrix from state at t (row)
   %   to the state at t+1 (column) 
   %   The eigenvector associated with the unit eigenvalue
   %   of trans' is  the stationary distribution. 
   % 
   g2=sparse(cs,cs);
   g1=sparse(cs,cs);
   for i=1:cs
       g1(i,tdecis1(i))=1;
       g2(i,tdecis2(i))=1;
   end
   trans=[ prob(1,1)*g1 prob(1,2)*g1; prob(2,1)*g2 prob(2,2)*g2];
   trans=trans';
   probst = (1/(2*nkap))*ones(2*nkap,1);
   test=1;
   while test>10^(-8);
      probst1 = trans*probst;
      test = max(abs(probst1-probst));
      probst = probst1;
   end;
   %
   %   vectorize the decision rule to be conformable with probst
   %   calculate new aggregate capital stock  meanK
   %
   kk=decis(:);
   meanK=probst'*kk;
   %
   %  calculate measure over (k,s) pairs
   %  lambda has same dimensions as decis
   %
   lambda   = zeros(cs,2);
   lambda(:)= probst;
   %
   %   calculate stationary distribution of k
   %
   [v1,d1]      = eig(prob');
   [dmax,imax]  = max(diag(d1));
   probst1      = v1(:,imax);
   ss           = sum(probst1);
   probst1      = probst1/ss;
   probk        = sum(lambda');     %  stationary distribution of `captal' 
   probk        = probk';
   %
   %   form metric and update K
   %
   errorK   = abs((Kold-meanK)/Kold);
   K        = g*meanK + (1-g)*Kold;
   Kold     = K;
   
   disp([  iter errorK meanK Kold ]);
   iter = iter+1;
end;

%--------------------------------------------------%
%   print out results
%--------------------------------------------------%
disp('PARAMETER VALUES');
disp('');
disp('    sigma      beta      delta       A      alpha      theta'); 
disp([ sigma beta delta A alpha theta]);
disp(''); 
disp('EQUILIBRIUM RESULTS ');
disp('');
disp('      K         N        wage      rent');
disp([ Kold N wage rent ]);
disp('-----------------------------------------');

Y       = A*K^(alpha)*N^(1-alpha);
H       = N;
q       = rent;
w       = wage;
gdelta  = 1. - delta;
r       = q - gdelta;
R       = 1 + r;
tauL    = 0.0;
tauK    = 0.0;

fprintf('------- Results -------\n');

fprintf('Y= %12.8f\n', Y);
fprintf('K= %12.8f\n', K);
fprintf('H= %12.8f\n', H);
fprintf('N= %12.8f\n', N);
fprintf('R= %12.8f\n', R);
fprintf('w= %12.8f\n', w);
fprintf('tauL= %12.8f\n', tauL);
fprintf('tauK= %12.8f\n', tauK);
fprintf('-----------------------\n');
fprintf('K/Y= %12.8f\n', K/Y);
fprintf('-----------------------\n');

% break
%
%    simulate life histories of the agent
%
% disp('SIMULATING LIFE HISTORY');
% k = Kold;               % initial level of capital 
% n = 100;                % number of periods to simulate
% s0 = 1;                 % initial state 
% hist = zeros(n-1,2);
% cons = zeros(n-1,1);
% invest = zeros(n-1,1);
% grid = [ (0:inckap:maxkap)' ];  
% [chain,state] = markov(prob,n,s0);
% for i = 1:n-1;
%     hist(i,:) = [ k chain(i) ];
%     I1 = round(k/inckap) ;
%     I2 = round(k/inckap) + 1;
%     if I1 == 0;
%        I1=1;
%        disp('N.B.  I1 = 0');
%     end;
%     if I2 > nkap;
%        I2 = nkap;
%        disp('N.B.  I2 > nkap');
%     end;
%     weight = (grid(I2,1) - k)/inckap; 
%     kprime = weight*(decis(I1,chain(i))) +  (1-weight)*(decis(I2,chain(i)));
%     if chain(i) == 1;
%        cons(i) = wage + (rent + delta)*k - kprime;
%     elseif chain(i) == 2;
%        cons(i) = wage*theta + (rent + delta)*k - kprime;
%     else;
%       disp('something is wrong with chain');
%       chain
%     end;
%     k = kprime;
%     invest(i) = kprime;
% end;
% 
% 
% figure(1);
% subplot(2,2,1);
% plot((1:n-1)',invest,(1:n-1)',cons);
% title('MODEL 2: INVESTMENT AND CONSUMPTION');
% print histmod2
% disp('Covariance matrix');
% disp([cov(cons,invest)]);
% %
% %     calculate income distribution
% %
% income =  [ (rent*grid + wage)  (rent*grid + wage*theta) ]  ; 
% [ pinc, index ] = sort(income(:));
% plambda = lambda(:);
% 
% subplot(2,2,2);
% plot(pinc,plambda(index));
% title('MODEL 2: INCOME DISTRIBUTION');
% xlabel('INCOME LEVEL');
% ylabel('% OF AGENTS');
% print distmod2
% %
% %    calculate capital distribution
% %
% subplot(2,2,3);
% plot(grid,probk);
% title('MODEL 2: CAPITAL DISTRIBUTION');
% xlabel('CAPITAL GRID');
% ylabel('% OF AGENTS');
% print capdmod2
% timer=etime(clock,t1);
% disp([timer]);
% diary off
% % save mod2 grid v lambda probk income hist invest cons

toc



