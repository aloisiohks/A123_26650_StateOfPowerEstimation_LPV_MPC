function [vk,x,OCV,dOCVdT] = iterModel(x, ik,Tfk, model,deltaT)
    % Load constants
    
  T = x(4) ;  
  z0 = x(1);
  iR0 = x(2);
  h0 = x(3);
  
  ik = ik(:); iR0 = iR0(:);
  
  % Get model parameters from model structure
  RCfact = exp(-deltaT./abs(getParamESC('RCParam',T,model)))';
  G = getParamESC('GParam',T,model);
  Q = getParamESC('QParam',T,model);
  M = getParamESC('MParam',T,model);
  M0 = getParamESC('M0Param',T,model);
  RParam = getParamESC('RParam',T,model);
  R0Param = getParamESC('R0Param',T,model);
  etaParam = getParamESC('etaParam',T,model);
  Rc = model.Rc;
  Ru = model.Ru;
  Cc = model.Cc;
  Cs = model.Cs;
  
  if ik<0
      eta = etaParam; 
  else
      eta = 1;
  end
  
  % Simulate the dynamic states of the model
%   if exist('ss','file'), % use control-system-toolbox method, if available
%     sysd= ss(diag(RCfact),1-RCfact,eye(length(RCfact)),0,-1);
%     irk = lsim(sysd,etaik,[],iR0);
%   else

      irk = RCfact*iR0 + (1-RCfact)*eta*ik;
    
%   end
  zk = z0-eta*ik*deltaT/(Q*3600); 
  if any(zk>1.1),
    warning('Current may have wrong sign as SOC > 110%');
  end
  
  % Hysteresis stuff
  hk = h0; 
  fac=exp(-abs(G*eta*ik*deltaT/(3600*Q)));
  hk=fac*hk-(1-fac)*sign(ik);
  
 
    
  % Compute output equation
  OCV = OCVfromSOCtemp(zk,T,model);
  vk = OCV - irk*RParam' - ik*R0Param + M*hk - M0*sign(ik);

    
  dOCVdT = interp1(model.SOCs,model.dOCVdT,z0); 
  
  
  A = [ (1-deltaT/Rc/Cc+(deltaT/Cc)*dOCVdT*ik) (deltaT/Rc/Cc);...
        (deltaT/Rc/Cs) (1 - deltaT/Rc/Cs - deltaT/Ru/Cs)];
  
  B = [(deltaT/Cc)*(OCV-vk) 0; 0  deltaT/(Ru*Cs)];
  
  xT = A*[x(4);x(5)] + B*[(ik) ; Tfk];
  
  x = [zk; irk; hk; xT];
  
end  
    
    
    
    

























% function [v,x,OCV] = iterModel(x, ik,Tfk, model,deltaT)
%     % Load constants
%     
% temp = x(4) ; % Core temperature   Tc[k-1]    
% 
% R = getParamESC('RParam',temp,model);
% RC = getParamESC('RCParam',temp,model);    
% eta = getParamESC('etaParam',temp,model);  
% Q = getParamESC('QParam',temp,model); 
% M = getParamESC('MParam',temp,model); 
% M0 = getParamESC('M0Param',temp,model); 
% R0 = getParamESC('R0Param',temp,model);
% Rc = model.Rc;
% Ru = 1;
% Cc = model.Cc;
% Cs = model.Cs;
% 
% 
% zk = x(1);
% ir = x(2);
% hk = x(3);
% temp = x(4);  % Core temperature   Tc[k]   
% 
% OCV = OCVfromSOCtemp(zk,temp,model);
% v = OCV + M*hk -M0*sign(ik) - R*ir - R0*ik;
% 
% A_RC = exp(-deltaT/RC);
% A_HK = exp(-abs(ik*eta*deltaT/(3600*Q)));
% 
% dOCVdT = entropyvariation(zk);
% 
% A = [1 0 0  0 0 ;...
%      0 A_RC 0 0 0;...
%      0 0 A_HK 0 0;...
%      0 (ik*R*deltaT/Cc)  (-ik*M*deltaT/Cc) (1-deltaT/Rc/Cc+(deltaT/Cc)*ik*dOCVdT)  (deltaT/Rc/Cc);...
%      0 0 0 (deltaT/Rc/Cs)     (1 - deltaT/Rc/Cs - deltaT/Ru/Cs)];
% 
% 
% 
% B = [-deltaT*eta/(3600*Q) 0  0; ...
%     (1-A_RC) 0 0; ...
%     0 (A_HK-1) 0;...
%     (-M*hk + M0*sign(ik) + R*ir + R0*ik)*deltaT/Cc 0 0;...
%     0  0  deltaT/Ru/Cs];
%     
% x = A*x + B*[ik;sign(ik);Tfk];
% 
% 
%     
%     
%     
% end  
%     
%     
%     
%     
