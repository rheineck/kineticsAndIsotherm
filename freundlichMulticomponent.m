function [cexqeAF,cexqeBF,erLF] = freundlichMulticomponent(freundlichSingle,cexqeA,cexqeB)

kFA = freundlichSingle(1,1);
kFB = freundlichSingle(1,2);
nA = freundlichSingle(2,1);
nB = freundlichSingle(2,2);

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

[a,~,~] = linRegression(ceA,ceB);

a12 = a;
a21 = 1/a12;

n = length(ceA);
m = length(ceB);
if n == m
    qeAF = zeros(n,1);
    qeBF = zeros(n,1);
    for i = 1:n
        qeAF(i) = kFA*ceA(i)*(ceA(i)+a12*ceB(i))^((1/nA)-1);
        qeBF(i) = kFB*ceB(i)*(ceA(i)*a21+ceB(i))^((1/nB)-1);
    end
else
    error('Dimensions dont match!');
end

erA = avgerror(qeAF,qeA);
erB = avgerror(qeBF,qeB);
% [~,~,r2LA] = linRegression(ceA,qeAF);
% [~,~,r2LB] = linRegression(ceB,qeBF);
% [r2LA] = determinationCoefficient(ceA,qeAF);
% [r2LB] = determinationCoefficient(ceB,qeBF);

cexqeAF = [ceA qeAF];
cexqeBF = [ceB qeBF];
erLF = [erA; erB];

end