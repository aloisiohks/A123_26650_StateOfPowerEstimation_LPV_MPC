function [ M, gamma ] = constraints( dx, mpcData,rflag)


    Np = mpcData.Np;
    Nc = mpcData.Nc;
    Sigma = mpcData.Sigma;
    model = mpcData.model; 
    
    [Phi_Tc,G_Tc] = predMat_LPV(dx(1:5),dx(end),mpcData,'temp');
    [Phi_v,G_v] = predMat_LPV(dx(1:5),dx(end),mpcData,'voltage');

    
    if mpcData.k > 2
        OCV = mpcData.OCV;
    else
        OCV =  OCVfromSOCtemp(dx(1),dx(4),model)*ones(1,Np);
    end
    
    %****************************************************************************        
        % CHARGE
    if strcmp(rflag,'charge')
        umax = 0;
        umin = mpcData.const.u_min;
        
        M = [
                -Sigma;         % i  > i_min
                  Sigma;         % i  < i_max
                  G_v;         % v  < v_max
                  G_Tc;        % Tc < Tc_max
                 
        ] ;
        
        gamma = [ 
            
            % u constraints
            -umin*ones(Nc,1) + mpcData.uk_1*ones(Nc,1);
              umax*ones(Nc,1) - mpcData.uk_1*ones(Nc,1);
       
            % v constraints
              mpcData.const.v_max*ones(Np ,1) - Phi_v*dx - OCV' ;
%              -mpcData.const.v_min + G_v*dx  + (OCV -C_v(1)*dx(1))*ones(Np,1);

             % Tc constraints
              mpcData.const.tc_max - Phi_Tc*dx ;
%              -mpcData.const.tc_min + G_Tc*dx 
            ];
        
    end
        
%****************************************************************************        
        % DISCHARGE
 if strcmp(rflag,'discharge')         
        umax = mpcData.const.u_max;
        umin = 0;
         
                M = [
                   -Sigma;         % i  > i_min
                   Sigma;         % i  < i_max
                     -G_v  ;        % v  > v_min
                   G_Tc;        % Tc < Tc_max 
%                  -eye(Nc);       % dig[k] > dig_min 
        ] ;
        
         gamma = [ 
                 -umin*ones(Nc,1) + mpcData.uk_1*ones(Nc,1);     % i  > i_min
              umax*ones(Nc,1) - mpcData.uk_1*ones(Nc,1);     % i  < i_max

                      
             -mpcData.const.v_min + Phi_v*dx  + OCV';   % v  > v_min
             mpcData.const.tc_max - Phi_Tc*dx ;      % Tc < Tc_max
%        -mpcData.const.du_min*ones(Nc,1);      % dig[k] > dig_min 
                           
            ];
        
        
    end
       

end

