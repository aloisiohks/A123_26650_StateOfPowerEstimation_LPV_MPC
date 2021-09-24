function [spkfData,vk] = iterSPKF(vk,ik,Tk,Tfk,Tsk,deltat,spkfData)
model = spkfData.model;
% Load the cell model parameters
Q = getParamESC('QParam',Tk,model);
G = getParamESC('GParam',Tk,model);
M = getParamESC('MParam',Tk,model);
M0 = getParamESC('M0Param',Tk,model);
RC = exp(-deltat./abs(getParamESC('RCParam',Tk,model)))';
R = getParamESC('RParam',Tk,model)';
R0 = getParamESC('R0Param',Tk,model);
Rc = model.Rc;
Ru = model.Ru;
Cc = model.Cc;
Cs = model.Cs;
eta = getParamESC('etaParam',Tk,model);
if ik<0, ik=ik*eta; end;
% Get data stored in spkfData structure
I = spkfData.priorI;
SigmaX = spkfData.SigmaX;
xhat = spkfData.xhat;
Nx = spkfData.Nx;
Nw = spkfData.Nw;
Nv = spkfData.Nv;
Na = spkfData.Na;
Snoise = spkfData.Snoise;
Wc = spkfData.Wc;
irInd = spkfData.irInd;
hkInd = spkfData.hkInd;
zkInd = spkfData.zkInd;
TcInd = spkfData.TcInd;
TsInd = spkfData.TsInd;
if abs(ik)>Q/100, spkfData.signIk = sign(ik); end;
signIk = spkfData.signIk;
% Step 1a: State estimate time update
% -Create xhatminus augmented SigmaX points
% -Extract xhatminus state SigmaX points
% -Compute weighted average xhatminus(k)
% Step 1a-1: Create augmented SigmaX and xhat

[sigmaXa,p] = chol(SigmaX,'lower');
if p>0,
fprintf('Cholesky error. Recovering...\n');
theAbsDiag = abs(diag(SigmaX));
sigmaXa = diag(max(SQRT(theAbsDiag),(SQRT(spkfData.SigmaW(1,1)))));
end
% real(sigmaXa)
% zeros([Nx Nw+Nv])
% zeros([Nw+Nv Nx])
% Snoise
sigmaXa=[real(sigmaXa) zeros([Nx Nw+Nv]); zeros([Nw+Nv Nx]) Snoise];
xhata = [xhat; zeros([Nw+Nv 1])];
% NOTE: sigmaXa is lower-triangular
% Step 1a-2: Calculate SigmaX points (strange indexing of xhata to
% avoid "repmat" call, which is very inefficient in MATLAB)
Xa = xhata(:,ones([1 2*Na+1])) + ...
spkfData.h*[zeros([Na 1]), sigmaXa, -sigmaXa];
% Step 1a-3: Time update from last iteration until now
% stateEqn(xold,current,xnoise)

Xx = stateEqn(Xa(1:Nx,:),I,Xa(Nx+1:Nx+Nw,:),Tfk);
xhat = Xx*spkfData.Wm;
% Step 1b: Error covariance time update
% -Compute weighted covariance sigmaminus(k)
% (strange indexing of xhat to avoid "repmat" call)
Xs = Xx - xhat(:,ones([1 2*Na+1]));
SigmaX = Xs*diag(Wc)*Xs';

% Step 1c: Output estimate
% -Compute weighted output estimate yhat(k)
I = ik; yk = [vk;Tsk];

Y = outputEqn(Xx,I,Xa(Nx+Nw+1:end,:),Tk,model);
yhat = Y*spkfData.Wm;

% Step 2a: Estimator gain matrix
Ys = Y - yhat(:,ones([1 2*Na+1]));
SigmaXY = Xs*diag(Wc)*Ys';
SigmaY = Ys*diag(Wc)*Ys';
L = SigmaXY/SigmaY;
% Step 2b: State estimate measurement update

r = [yk - yhat]; % residual. Use to check for sensor errors...
% if r^2 > 100*SigmaY, L(:,1)=0.0; end
xhat = xhat + L*r;
xhat(zkInd)=min(1.05,max(-0.05,xhat(zkInd)));
% Step 2c: Error covariance measurement update
SigmaX = SigmaX - L*SigmaY*L';
[~,S,V] = svd(SigmaX);
HH = V*S*V';
SigmaX = (SigmaX + SigmaX' + HH + HH')/4; % Help maintain robustness
% Q-bump code
if r(1)^2>4*SigmaY(1), % bad voltage estimate by 2-SigmaX, bump Q
fprintf('Bumping sigmax\n');
SigmaX(zkInd,zkInd) = SigmaX(zkInd,zkInd)*spkfData.Qbump;
end
% Save data in spkfData structure for next time...
spkfData.priorI = ik;
spkfData.SigmaX = SigmaX;
spkfData.xhat = xhat;

spkfData.zkBounds = 3*sqrt(SigmaX(zkInd,zkInd));
spkfData.hkBounds = 3*sqrt(SigmaX(hkInd,hkInd));
spkfData.irckBounds = 3*sqrt(SigmaX(irInd,irInd));
spkfData.TckBounds = 3*sqrt(SigmaX(TcInd,TcInd));
spkfData.TskBounds = 3*sqrt(SigmaX(TsInd,TsInd));



vk = yhat(1);
Tck = xhat(4);
Tsk = xhat(5);


% Calculate new states for all of the old state vectors in xold.
function xnew = stateEqn(xold,current,xnoise,Tfk)
current = current + xnoise(1,:); % noise adds to current
Tfk = Tfk + xnoise(2,:);
xnew = 0*xold;

xnew(irInd,:) = RC*xold(irInd,:) + (1-RC)*current;
Ah = exp(-abs(current*G*deltat/(3600*Q))); % hysteresis factor
xnew(hkInd,:) = Ah.*xold(hkInd,:) - (1-Ah).*sign(current);
xnew(zkInd,:) = xold(zkInd,:) - current/3600/Q;

xnew(hkInd,:) = min(1,max(-1,xnew(hkInd,:)));
xnew(zkInd,:) = min(1.05,max(-0.05,xnew(zkInd,:)));

xnew(TcInd,:) = (1-deltat/Rc/Cc)*xold(TcInd,:) + (deltat/Rc/Cc)*xold(TsInd,:)...
    + (current*R*deltat/Cc).*xold(irInd,:) -(current*M*deltat/Cc).*xold(hkInd,:) + (current.^2)*R0*deltat/Cc;
xnew(TsInd,:) = (deltat/Rc/Cs)*xold(TcInd,:) +(1 - deltat/Rc/Cs - deltat/Ru/Cs)*xold(TsInd,:) + Tfk*deltat/Ru/Cs;


end
% Calculate cell output voltage for all of state vectors in xhat
function [yhat] = outputEqn(xhat,current,ynoise,T,model)
yhat = OCVfromSOCtemp(xhat(zkInd,:),T,model);
yhat = yhat + M*xhat(hkInd,:) + M0*signIk;
yhat = yhat - R*xhat(irInd,:) - R0*current + ynoise(1,:);
Tshat = xhat(TsInd,:)+ynoise(2,:);
yhat = [yhat;Tshat];
end
% "Safe" square root
function X = SQRT(x)
X = sqrt(max(0,x));
end
end