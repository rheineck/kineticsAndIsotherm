% Non modified Competitive Redlich Peterson Isotherm
function [cexqeANRP,cexqeBNRP,erNRP,sdNRP] = nModifiedRP(cexqeA,cexqeB,redlichPeterson)

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

n = length(ceA);
m = length(ceB);

qeANRP = zeros(n,1);
qeBNRP = zeros(n,1);

if n == m
    NRPA = @(redlichPeterson,ceA,ceB) (redlichPeterson(1,1).*ceA)/(1+(redlichPeterson(2,1).*ceA.^redlichPeterson(3,1))+(redlichPeterson(2,2).*ceB.^redlichPeterson(3,2)));
    NRPB = @(redlichPeterson,ceA,ceB) (redlichPeterson(1,2).*ceB)/(1+(redlichPeterson(2,1).*ceA.^redlichPeterson(3,1))+(redlichPeterson(2,2).*ceB.^redlichPeterson(3,2)));
    
    for i = 1:n
        qeANRP = NRPA(redlichPeterson,ceA,ceB);
        qeBNRP = NRPB(redlichPeterson,ceA,ceB);
    end
    
    erANRP = avgerror(qeA,qeANRP);
    erBNRP = avgerror(qeB,qeBNRP);
    sdANRP = standardDeviation(ceA,qeA,NRPA,redlichPeterson);
    sdBNRP = standardDeviation(ceB,qeB,NRPB,redlichPeterson);
    
    erNRP = [erANRP; erBNRP];
    sdNRP = [sdANRP; sdBNRP];
    cexqeANRP = [ceA qeANRP];
    cexqeBNRP = [ceB qeBNRP];
else
    error('Dimensions dont match!');
end    
end