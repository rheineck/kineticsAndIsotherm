function [qmax,kL,r2L]=langmuirModel(cexqe)

[~,j]=size(cexqe);
if j~=2
    error('Just two columns are supported');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

L = @(x,ce) (x(1)*x(2).*ce)./(1+x(2).*ce);
x0(1) = input('Initial Values: ');
x0(2) = input('Initial Values: ');
% x0 = [1 1];
par = lsqcurvefit(L,x0,ce,qe);

qmax = par(1);
kL = par(2);

r2L = determinationCoefficient(ce,qe,L,par);

end