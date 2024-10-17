function [kF,n,r2F]=freundlichModel(cexqe)

[~,j]=size(cexqe);
if j~=2
    error('Just two columns are supported!');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

F = @(x,ce) (x(1).*ce.^(1/x(2)));
x0(1) = input('Initial Values: ');
x0(2) = input('Initial Values: ');
% x0 = [1 1];
par = lsqcurvefit(F,x0,ce,qe);

kF = par(1);
n = par(2);

r2F = determinationCoefficient(ce,qe,F,par);

end