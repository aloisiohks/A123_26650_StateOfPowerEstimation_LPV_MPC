function [linMatrices,xekf,ekfData] = initEKF(v0,SOC0,ir0,hk0,Tc0,Ts0,SigmaX0,SigmaV,SigmaW,model,deltaT)
 
    
% Initial state description
ekfData.zkInd = 1;
ekfData.irInd = 2;
ekfData.hkInd = 3;
ekfData.TcInd = 4;
ekfData.TsInd = 5;

xekf = [SOC0 ir0 hk0  Tc0 Ts0]'; % initial state
ekfData.xhat = xekf; % initial state
% Covariance values
ekfData.SigmaX = SigmaX0;
ekfData.SigmaV = SigmaV;
ekfData.SigmaW = SigmaW;
ekfData.Qbump = 5;

ekfData.SOC = SOC0;
ekfData.prior_SOC = ekfData.SOC;    % Past SOC guess
OCV = OCVfromSOCtemp(SOC0,Tc0,model);
ekfData.OCV_1 = OCV;
ekfData.v_1 = v0;
% previous value of current
ekfData.uk_1 = 0;   % Current applied before first time step
ekfData.signIk = 0;
% store model data structure too
ekfData.model = model;
ekfData.deltaT = deltaT;
linMatrices = linMat_Ts(xekf,model,ekfData.uk_1,deltaT);
ekfData.vkbounds = 3*sqrt(linMatrices.Chat*SigmaX0*linMatrices.Chat' + linMatrices.Dhat*SigmaV*linMatrices.Dhat');
ekfData.soc_bounds = 3*sqrt(SigmaX0(1,1));
ekfData.Tc_bounds = 3*sqrt(SigmaX0(4,4));
ekfData.Ts_bounds = 3*sqrt(SigmaX0(5,5));
 
end