% Extended Sips Model
function [cexqeAES,cexqeBES,erS,sdSM] = extendedSipsModel(cexqeA,cexqeB,sips)

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

n = length(ceA);
m = length(ceB);

qeAES = zeros(n,1);
qeBES = zeros(n,1);

if n == m
    ESA = @(sips,ceA,ceB) (sips(1,1)*sips(2,1).*ceA.^sips(3,1))/(1+(sips(2,1).*ceA.^sips(3,1))+(sips(2,2).*ceB.^sips(3,2)));
    ESB = @(sips,ceA,ceB) (sips(1,2)*sips(2,2).*ceB.^sips(3,2))/(1+(sips(2,1).*ceA.^sips(3,1))+(sips(2,2).*ceB.^sips(3,2)));

    for i = 1:n
        qeAES(i) = ESA(sips,ceA,ceB);
        qeBES(i) = ESB(sips,ceA,ceB);
    end
    
    erSA = avgerror(qeAES,qeA);
    erSB = avgerror(qeBED,qeB);
    sdSMA = standardDeviation(ceA,qeA,ESA,sips);
    sdSMB = standardDeviation(ceB,qeB,ESB,sips);
    
    erS = [erSA; erSB];
    sdSM = [sdSMA; sdSMB];
    cexqeAES = [ceA qeAES];
    cexqeBES = [ceB qeBES];
    
else
    error('Dimensions dont  match!');
end
end