%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Huggett-type model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;
tic

%%----------------------------------------------------------%%
%%[1.]-Parameters
%%----------------------------------------------------------%%
% Preference
gsigma      = 2.0;                          % intertemporal elasticity of consumption
gbeta       = .98;                           % discount rate  
% Shock and markov Switching 
yv          = [.25 3];       % level of shock 
P           = [.6 .4 ; 
               .3 .7];       % markov Switching      
% Stationary measure of types
Pv =1/2*ones(1,2);
for i=1:30; Pv=Pv*P; end
% 
% Asset space
bmin        = -2; %-2.; %-4;                      % minimum asset holding
bmax        = 2;                           % maximum asset holding 
grid        = round(100*(bmax-bmin));        % number of grid points     
bv          = linspace(bmin,bmax,grid); 	% vector of k grid points
gridstep    = bv(2)-bv(1);      			% gridstep
distance    = bv(end)-bv(1); 

tolerance   = 10^-3;
Lambda      = .5;
itermax     = 100;  

%%----------------------------------------------------------%%
%%[2.] Start a big loop 
%%----------------------------------------------------------%% 
r           = 0.03;    
R           = 1.+r;
Rold        = 1.; 
iter        = 0;
error       = 100; 
errorv      = [];

%%-----------------------------------------------% 
%%[2.1]-Decision Rule of saving 
%%-----------------------------------------------% 
%  initial guess
Vnext       = [(bv'-bmin+.1).^(1-gsigma)/(1-gsigma), (bv'-bmin+.1).^(1-gsigma)/(1-gsigma)];    
Vnow        = Vnext;
EV          = (P*Vnext')';
bbv         = bv';

% errors  
iter1        = 0;
error1       = 100; 
error1v      = 100;
%
while (iter1 < 300 & error1>10^-4)
    % expected value    
    EV           = (P*Vnext')';
    % loop over shocks and assets
    for i = 1: grid     % loop over grid points
        for z = 1:2     % number of shocks 
            %
            income  = yv(z) + R*bv(i);     
            %           
            cv  = income-bbv;
            cv  = (cv>0).*cv + (cv<=0)*10^-10;
            vv  = cv.^(1-gsigma)/(1-gsigma) + gbeta*EV(:,z);
            % find max
            [val pos]  = max(vv);            
            %            
            Bopt(i,z)   = bbv(pos);             % optimal capital
            Copt(i,z)   = cv(pos);              % optimal consumption 
            Vnow(i,z)   = vv(pos);              % value function            
            %            
        end% for z: over shock                        
    end % for i: over grid    
    % error 
    error1       = 100*sum(sum(abs(Vnow-Vnext)))/sum(sum(Vnext));
    iter1        = iter1+1;
    % update
    Vnext        = Vnow;    
      
end % while iter1

figure(1)
subplot(2,1,1);
plot(bv,Vnow);
legend('Low','High');
xlabel('b');
title('Value Function');

subplot(2,1,2);
plot(bv,Bopt);
xlabel('b_t');
ylabel('b_{t+1}');
title('Decision rule');


  
%%-----------------------------------------------% 
%% [2.2] Stationary distribution 
%%-----------------------------------------------% 
 
% initial 
Mnow         = ones(grid,2)/(grid*2);
%
iter2       = 0;
error2      = 10;   

while (iter2 <1000 & error2>10^-10)   
    %
    Mnext   = zeros(grid,2);
    %
    for i = 1: grid
        for z = 1:2
                % a trick to find positions in grid points
                posL=min(floor((Bopt(i,z)-bmin)/distance*grid)+1,grid);
                posL=round(posL);
                % adjust correct position
                if bv(posL)>Bopt(i,z)
                    posL=posL-1;
                end
                posH		= min(posL+1,grid);
                weight      = ((Bopt(i,z) - bv(posL))/gridstep);
                % probabilities in the next period     
                transp =  Mnow(i,z)*P(z,:); 
                % fill in the probabilities in coresponding positions
                for zz=1:2
                    Mnext(posL,zz)= Mnext(posL,zz) + (1-weight)*transp(zz);
                    Mnext(posH,zz)= Mnext(posH,zz) + weight*transp(zz);
                end
        end % end for z   
    end % end for i: grid
    error2      = sum(sum(abs(Mnext-Mnow)));
    Mnow        = Mnext;
    iter2       = iter2+1;
end % while iter2


%%-----------------------------------------------% 
%% [3] Reporting results 
%%-----------------------------------------------% 

[rs,cs] = size(Mnow);
F       = zeros(rs,cs);
F(1,:)  = Mnow(1,:);
for z = 1:2
    for i =2:rs
    F(i,z) = F(i-1,z)+ Mnow(i,z);
    end
end

figure(2);
subplot(2,1,1);
plot(bv,Mnow);
legend('Low','High');
xlabel('b');
title('Stationary Distribution:PDF');

subplot(2,1,2);
plot(bv',sum(F,2));
xlabel('b');
title('Stationary Distribution:CDF');


toc