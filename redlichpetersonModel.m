% Redlich-Peterson Isotherm Model
function [kRP,aRP,beta,r2RP,sdRP]=redlichpetersonModel(cexqe,x0)

[~,j] = size(cexqe);
if j ~= 2
    error('Just two columns are supported!');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

RP = @(x,ce) (x(1).*ce)/(1+x(2).*ce.^(x(3)));
par = lsqcurvefit(RP,x0,ce,qe);

kRP = par(1);
aRP = par(2);
beta = par(3);

r2RP = determinationCoefficient(ce,qe,RP,par);
sdRP = standardDeviation(ce,qe,RP,par);
end