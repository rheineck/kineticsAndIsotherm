% Sips Model
function [sips,r2S,sdS] = sipsModel(cexqe,x0)

[~,j] = size(cexqe);
if j ~= 2
    error('Just two columns are supported!');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

S = @(x,ce) (x(1)*x(2).*ce^(1/x(3)))./(1+x(2).*ce^(1/x(3)));
par = lsqcurvefit(S,x0,ce,qe);

qm = par(1);
kS = par(2);
nS = par(3);

sips = par;

r2S = determinationCoefficient(ce,qe,S,par);
sdS = standardDeviation(ce,qe,S,par);
end