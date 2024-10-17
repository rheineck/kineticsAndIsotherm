%% Pseudo Second Order Model

function [qePSO,k2,r2PSO,PSO] = pseudoSecondOrder(txqt,x0)

t = txqt(:,1);
qt = txqt(:,2);

PSO = @(x,t) ((x(1)^2)*x(2).*t)./((x(1)*x(2).*t)+1);
% x0(1) = input('Initial Values: ');
% x0(2) = input('Initial Values: ');
% % x0 = [1,1];
par = lsqcurvefit(PSO,x0,t,qt);

qePSO = par(1);
k2 = par(2);

r2PSO = determinationCoefficient(t,qt,PSO,par);

end