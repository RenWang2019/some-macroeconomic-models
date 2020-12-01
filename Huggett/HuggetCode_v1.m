%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Huggett-type model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
tic

%%----------------------------------------------------------%%
%%[1.]-Parameters
%%----------------------------------------------------------%%
% Preference
gsigma      = 2;                          % intertemporal elasticity of consumption
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
bmin        = -2; %-4;                      % minimum asset holding
bmax        = 2;                           % maximum asset holding 
grid        = round(100*(bmax-bmin));        % number of grid points     
bv          = linspace(bmin,bmax,grid); 	% vector of k grid points
gridstep    = bv(2)-bv(1);      			% gridstep
distance    = bv(end)-bv(1); 

%%----------------------------------------------------------%%
%%[2.] Start a big loop 
%%----------------------------------------------------------%% 
% intials 
r           = 0.1;    
R           = 1+r;
Rold        = 1.; 
% 
tolerance   = 10^-5;
Lambda      = .5;
itermax     = 100;  
iter        = 0;
error       = 100; 
errorv      = [];

while iter < itermax&error>tolerance    
    if rem(iter,1) == 0
       fprintf('%6.0f %3.4f %3.2f...\n',iter,error,toc );        % tells you how far we are in the loop
    end  

    %%-----------------------------------------------% 
    %%[2.1]-Decision Rule of saving 
    %%-----------------------------------------------% 
    %  initial guess
    Vnext       = [(bv'-bmin+.1).^(1-gsigma)/(1-gsigma), (bv'-bmin+0.1).^(1-gsigma)/(1-gsigma)];    
    Vnow        = Vnext;
    EV          = (P*Vnext')';
    bbv         = bv';
    % errors  
    iter1        = 0;
    error1       = 100; 
    error1v      = 100;
    %
    while (iter1 < 300 & error1>.0001)
        % loop over shocks and assets
        for i = 1: grid     % loop over grid points
            for z = 1:2     % number of shocks 
                %
                income  = yv(z) + R*bv(i);     
                %           
                % start w/ all possible
                %[.1]
    %             for jj=1:grid
    %                cjj      = income -bbv(jj);
    %                cv(jj)   = (cjj>0).*cjj + (cjj<=0)*10^-3;
    %                vv(jj)   = cjj^(1-sigma)/(1-sigma) + beta*EV(jj,z);
    %             end
                %[.2]
                cv  = income-bbv;
                cv  = (cv>0).*cv + (cv<=0)*10^-3;
                vv  = cv.^(1-gsigma)/(1-gsigma) + gbeta*EV(:,z);
                % find max
                [val pos]  = max(vv);            
                %
    %             bt   = bbv(pos);  
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
        % expected value    
        EV           = (P*Vnext')';  
    end % while iter1
    % if iter ==10;  keyboard; end

    %%-----------------------------------------------% 
    %% [2.2] Stationary distribution 
    %%-----------------------------------------------% 

    % initial 
    Mnow         = ones(grid,2)/(grid*2);
    %
    iter2       = 0;
    error2      = 10;   

    while (iter2 <500 & error2>10^-10)   
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
    %% [2.3] Clear bond market 
    %%-----------------------------------------------% 
    surplus     = sum(sum(Mnow.*Bopt));
    % % 
    R           = R - 0.0001*surplus;    
    R           = (1-Lambda)*Rold + Lambda*R;
    r           = R-1;
    %
    error       = 100*abs(R-Rold);
    errorv      = [errorv error];
    Rold        = R;
    iter        = iter+1;

end % while

%%-----------------------------------------------% 
%% [3] Reporting results 
%%-----------------------------------------------% 
% pdf
F       = zeros(length(Mnow(:,1)),2);
F(1,:)  = Mnow(1,:);
for z = 1:2
    for i =2:length(Mnow(:,1));
    F(i,z) = F(i-1,z)+Mnow(i,z);
    end
end
disp('----Results------'); 
disp(['Inteterest rate: r =' num2str((R-1)*100) '%']);
disp(['Excess supply      =' num2str(surplus)]);

% disp('--- Annual net interest rates: percent -----')
%%--[3.] Plot 
subplot(2,2,1);
plot(bv,Vnow);
legend('Low','High');
xlabel('b');
title('Value Function');

subplot(2,2,2);
plot(bv,Bopt);
xlabel('b_t');
ylabel('b_{t+1}');
title('Decision rule');

subplot(2,2,3);
plot(bv',F);
xlabel('b');
title('Stationary Distribution:CDF');

subplot(2,2,4);
plot(bv,Mnow);
legend('Low','High');
xlabel('b');
title('Stationary Distribution:PDF');


% figure(2)
% plot(errorv)
toc