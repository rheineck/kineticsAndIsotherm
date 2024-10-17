%% Pseudo First Order Model

function [qePFO,k1,r2PFO,PFO] = pseudoFirstOrder(txqt,x0)

t = txqt(:,1);
qt = txqt(:,2);

PFO = @(x,t) x(1).*(1-exp(-x(2).*t));
% x0(1) = input('Initial Values: ');
% x0(2) = input('Initial Values: ');
% % x0 = [1,1];
par = lsqcurvefit(PFO,x0,t,qt);

qePFO = par(1);
k1 = par(2);

r2PFO = determinationCoefficient(t,qt,PFO,par);

end