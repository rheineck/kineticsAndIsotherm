function [cexqeAM,cexqeBM,erLM] = langmuirModified(langmuirSingle,cexqeA,cexqeB)

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
if n == m
    qeAM = zeros(n,1);
    qeBM = zeros(n,1);
    for i = 1:n
        qeAM(i) = (qmaxA*kLA*ceA(i))/(1+kLA*ceA(i)+kLB*ceB(i));
        qeBM(i) = (qmaxB*kLB*ceB(i))/(1+kLA*ceA(i)+kLB*ceB(i));
    end
else
    error('Dimensions dont match!');
end

erA = avgerror(qeAM,qeA);
erB = avgerror(qeBM,qeB);

% [~,~,r2LA] = linRegression(ceA,qeAM);
% [~,~,r2LB] = linRegression(ceB,qeBM);
% [r2LA] = determinationCoefficient(ceA,qeAM);
% [r2LB] = determinationCoefficient(ceB,qeBM);

cexqeAM = [ceA qeAM];
cexqeBM = [ceB qeBM];
erLM = [erA; erB];

end