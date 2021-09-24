function [Phi,G,augMatrices] = predMat_LPV(xhat,ik,mpcData,option)


% tic
% n = number of rows
% m = number of columms;

Np = mpcData.Np;
Nc = mpcData.Nc;

A_store = zeros(9,9,Np);
C_store = zeros(1,9,Np);


% G = zeros(Np,mA+mB+mBt);
% Phi = zeros(Np,Nc);

if mpcData.k-1 > 1    
    
    for j=1:Np
         
          xp = mpcData.xf(:,j);
          uk = mpcData.uk_f(j);
      
          [linMatrices] = linMat_LPV(xp,uk,mpcData);
          Am = linMatrices.A;
          Bm = linMatrices.B;
          
          if strcmp(option,'soc')
               Cm = linMatrices.C_z ;
               D = linMatrices.D_z ;
          end

          if strcmp(option,'voltage')
              Cm = linMatrices.C_v ;
              D = linMatrices.D_v ;
          end

          if strcmp(option,'temp')
              Cm = linMatrices.C_Tc ;
              D = linMatrices.D_Tc ;
          end
          
          [nA,mA] = size(Am);
          [~,mB] = size(Bm);
          [~,mC] = size(Cm);
          

          A = zeros(nA+mB,mA+mB);
          A(1:nA,1:mA) = Am;
          A(1:nA,nA+1:nA+mB)=Bm;
          A(nA+1:end,mA+1:end)=eye(mB);
          B = zeros(nA+mB,1);
          B(nA+mB,end) = 1;
          C = zeros(1,mA+mB);
          C(1,1:mC) = Cm;
          C(1,mA+1:end) = D;
        
          A_store(:,:,j) = A;
          C_store(:,:,j) = C;
        
    end
    
    % Make Phi 
    for i=1:Np
        if i==1
            Aux = A_store(:,:,i);
            Phi(i,:) = C_store(:,:,i)*Aux;
        else
            Aux = A_store(:,:,i)*Aux;
            Phi(i,:) = C_store(:,:,i)*Aux;
        end
    end  
    
    %Make G 
    
    for i = 1:Nc
        if i ==1
            for j = 1:Np
                if j == 1
                    Aux = B;
                else
                    Aux = A_store(:,:,j)*Aux;
                end
                G(j,i) = C_store(:,:,j)*Aux;
            end
        else
            for j=i:Np
                if j==i
                   Aux = B;
                else
                    Aux = A_store(:,:,j)*Aux;
                end
                G(j,i) = C_store(:,:,j)*Aux;
            end
        end
    end
                  
       
    
%%
       else
  
          xp = xhat;
          uk = ik;
      
          [linMatrices] = linMat_LPV(xp,uk,mpcData);
          Am = linMatrices.A;
          Bm = linMatrices.B;
          
          if strcmp(option,'soc')
               Cm = linMatrices.C_z ;
               D = linMatrices.D_z ;
          end

          if strcmp(option,'voltage')
              Cm = linMatrices.C_v ;
              D = linMatrices.D_v ;
          end

          if strcmp(option,'temp')
              Cm = linMatrices.C_Tc ;
              D = linMatrices.D_Tc ;
          end
          
          [nA,mA] = size(Am);
          [~,mB] = size(Bm);
          [~,mC] = size(Cm);
          

          A = zeros(nA+mB,mA+mB);
          A(1:nA,1:mA) = Am;
          A(1:nA,nA+1:nA+mB)=Bm;
          A(nA+1:end,mA+1:end)=eye(mB);
          B = zeros(nA+mB,1);
          B(nA+mB,end) = 1;
          C = zeros(1,mA+mB); 
          C(1,1:mC) = Cm;
          C(1,mA+1:end) = D;
        
        
            Phi = C*A;
            cy = C*B;  r = cy;
        for j = 2:Np
            Phi = [Phi; C*A^(j)];
            cy = [cy; C*(A^(j-1))*B];
        end
        
        ry = zeros(1,Nc);
        ry(1,1) = r;

        G = toeplitz(cy,ry);
        
     
end





augMatrices.A = A;
augMatrices.B = B;
augMatrices.C = C;
% toc
% disp('predMat')
end
