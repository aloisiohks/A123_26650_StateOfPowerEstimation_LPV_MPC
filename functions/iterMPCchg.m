%   iterMPC computes SOP charge at time k
%
%   xk_1: current state
%   mpcData: mpcData configuration
%   Nsim: prediction interval
%   deltaT: sampling time

function [Y,X,U,DelU,Pavg,Pinst,mpcData] = iterMPCchg(xk_1,mpcData,Nsim,deltaT)

    % Load MPC data
    Np = mpcData.Np;
    Nc = mpcData.Nc;
    Tfk = mean(mpcData.Tfk);    
    model = mpcData.model;
    DelU = zeros(1,Nsim+1);
    Ru  = mpcData.Ru;
    Ref = mpcData.const.z_max*ones(Np,1);
    u = mpcData.const.u_min;
    mpcData.uk_1 = u;

if ((xk_1(1) < mpcData.const.z_max) && (xk_1(4) < (mpcData.const.tc_max )) )
   
   % Calculate intial state and output values
   [v,~,~,~] = iterModel(xk_1, u,Tfk, model,deltaT);
   Y(:,1) = v;
   Xf = [xk_1;Tfk; 1; sign(u); u];
   X(:,1) = Xf;
   U(:,1) = u;
   xold = xk_1;

%    u=0;   
   %*************************************************************************
   % MAIN LOOP
   
   for kk = 2:Nsim
       
    mpcData.k = kk;   
    
    % Augmented state vector
    dx = [xold; Tfk; 1; sign(u); u];
    
   % Y = G*x + Phi*U
   %  Compute SOC prediction matrices
   [Phi,G,~] = predMat_LPV(xold,u,mpcData,'soc');
   
      
   F = -G'*(Ref - Phi*X(:,kk-1) );
   if mpcData.adap == 1
       [~,S,~] = svd(G'*G);
       [m,n] = size(S);
       Ru = (norm(F,2)/(2*mpcData.const.du_max*sqrt(Nc)))-(S(m,n)); 
   end
   E = (G'*G + Ru*eye(mpcData.Nc,mpcData.Nc)); 
   
   DU = -E\F;
   
   mpcData.Ru_store(kk) = Ru;

   
   [M,gamma] = constraints(dx,mpcData,'charge');
   
  
    if sum(M*DU - gamma > 0) > 0
        [DU,~,~] = hildreth(E,F,M,gamma,[],150);
    end
      
    du = DU(1);
    DelU(:,kk) = du;
    u = du + U(:,kk-1) ;
    U(:,kk) = u;
  
    [v,x,~,~] = iterModel(xold, u,Tfk, model,deltaT);
    xold = x;
    Y(:,kk) = v;
    X(:,kk) = [x;Tfk; 1; sign(u) ;u];
    mpcData.uk_1 = u;
    mpcData.v = v;
    
    
    % LPV
    ik = zeros(1,Np); OCV = ik;
    ik(1) = u + DU(2);
    [~,xp,OCV(1)] = iterModel(x,ik(1),Tfk,mpcData.model,mpcData.deltaT);
    xf(:,1)=xp;
    
    for i=2:Np
            if i < Nc
                ik(i) = ik(i-1) + DU(i+1);
            else
                ik(i) = ik(i-1);
            end
            
            [~,xp,OCV(i)] = iterModel(xp,ik(i),Tfk,mpcData.model,mpcData.deltaT);
            xf(:,i) = xp;
            
    end
      
    mpcData.xf = xf;
    mpcData.uk_f = ik;
    mpcData.OCV = OCV;
   
   end
   %*************************************************************************
   
    %POWER COMPUTATION
    Pinst = Y(2:Nsim).*U(1,2:Nsim);  % Compute instantaneous power over prediction horizon
    Etot = trapz(Pinst);             % Compute corresponding total energy
    Pavg = Etot/10;                   % Compute max (average) sustained power over interval


else

u = 0;
[v,x,~,~] = iterModel(xk_1, u,Tfk, model,deltaT);

mpcData.Ru_store = repmat(0,1,Nsim);
Y(:,1:Nsim) = repmat(v,1,Nsim);
X(:,1:Nsim) = repmat([x; Tfk; 1; sign(u); u],1,Nsim);
U(1:Nsim) = u*ones(1,Nsim);
Pavg = 0;
Pinst = [];    

end

end

