function [cexqeAIF,cexqeBIF,erLIF] = iastFreunclich(freundlichSingle,cexqeA,cexqeB)

kFA = freundlichSingle(1,1);
kFB = freundlichSingle(1,2);
nA = freundlichSingle(2,1);
nB = freundlichSingle(2,2);

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

n = length(qeA);
m = length(qeB);
if n == m
    ceAIF = zeros(n,1);
    ceBIF = zeros(n,1);
    for i = 1:n
        ceAIF(i) = (qeA(i)/(qeA(i)+qeB(i)))*((nA*qeA(i)+nB*qeB(i))/(nA*kFA))^(nA);
        ceBIF(i) = (qeB(i)/(qeA(i)+qeB(i)))*((nA*qeA(i)+nB*qeB(i))/(nB*kFB))^(nB);
    end
else
    error('Dimension dont match!');
end

erA = avgerror(ceAIF,ceA);
erB = avgerror(ceBIF,ceB);
% [~,~,r2LA] = linRegression(ceAIF,qeA);
% [~,~,r2LB] = linRegression(ceBIF,qeB);
% [r2LA] = determinationCoefficient(ceAIF,qeA);
% [r2LB] = determinationCoefficient(ceBIF,qeB);

cexqeAIF = [ceAIF qeA];
cexqeBIF = [ceBIF qeB];
erLIF = [erA; erB];

end