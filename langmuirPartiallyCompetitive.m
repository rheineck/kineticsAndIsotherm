function [cexqeAPC,cexqeBPC,erLPC] = langmuirPartiallyCompetitive(langmuirSingle,cexqeA,cexqeB)

qmaxA = langmuirSingle(1,1);
qmaxB = langmuirSingle(1,2);
kLA = langmuirSingle(2,1);
kLB = langmuirSingle(2,2);

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

n = length(ceA);
m = length(ceB);
if m == n
    qeAPC = zeros(n,1);
    qeBPC = zeros(n,1);
    for i = 1:n
        qeAPC(i) = (((qmaxA-qmaxB)*kLA*ceA(i))/1+kLA*ceA(i))+((qmaxA*kLA*ceA(i))/(1+kLA*ceA(i)+kLB*ceB(i)));
        qeBPC(i) = (qmaxB*kLB*ceB(i))/(1+kLA*ceA(i)+kLB*ceB(i));
    end
else
    error('Dimensions dont match!');
end

erA = avgerror(qeAPC,qeA);
erB = avgerror(qeBPC,qeB);
% [~,~,r2LA] = linRegression(ceA,qeAPC);
% [~,~,r2LB] = linRegression(ceB,qeBPC);
% [r2LA] = determinationCoefficient(ceA,qeAPC);
% [r2LB] = determinationCoefficient(ceB,qeBPC);

cexqeAPC = [ceA qeAPC];
cexqeBPC = [ceB qeBPC];
erLPC = [erA; erB];

end