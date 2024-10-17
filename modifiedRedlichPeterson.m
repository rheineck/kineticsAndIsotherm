% Modified Competitive Redlich Peterson Isotherm
function [cexqeAMRP,cexqeBMRP,erMRP,sdMRP] = modifiedRedlichPeterson(cexqeA,cexqeB,redlichPeterson)

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

n = length(ceA);
m = length(ceB);

if n == m
    MRPA = @(ceA,ceB,redlichPeterson) (redlichPeterson(1,1).*(ceA/redlichPeterson(4,1)))./(1+(redlichPeterson(2,1).*(ceA/redlichPeterson(4,1)).^redlichPeterson(3,1))+(redlichPeterson(2,2).*(ceB/redlichPeterson(4,2)).^redlichPeterson(3,2)));
    MRPB = @(ceA,ceB,redlichPeterson) (redlichPeterson(2,1).*(ceB/redlichPeterson(4,2)))./(1+(redlichPeterson(2,1).*(ceA/redlichPeterson(4,1)).^redlichPeterson(3,1))+(redlichPeterson(2,2).*(ceB/redlichPeterson(4,2)).^redlichPeterson(3,2)));
    
    parA = lsqcurvefit(MRPA,x0,ceA,qeA);
    parB = lsqcurvefit(MRPB,x0,ceB,qeB);
    
    redlichPeterson(4,1) = parA(1);
    redlichPeterson(4,2) = parB(2);
    
    for i = 1:n
        qeAMRP = MRPA(ceA,ceB,redlichPeterson);
        qeBMRP = MRPB(ceA,ceB,redlichPeterson);
    end
    
    erAMRP = avgerror(qeA,qeAMRP);
    erBMRP = avgerror(qeB,qeBMRP);
    sdAMRP = standardDeviation(ceA,qeA,MRPA,redlichPeterson);
    sdBMRP = standardDeviation(ceB,qeB,NRPB,redlichPeterson);
    
    cexqeAMRP = [ceA qeAMRP];
    cexqeBMRP = [ceB qeBMRP];
    erMRP = [erAMRP; erBMRP];
    sdMRP = [sdAMRP; sdBMRP];
else
    error('Dimensions dont match!');
end

end