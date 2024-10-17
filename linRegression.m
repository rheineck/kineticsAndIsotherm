function [a,b,r2] = linRegression(x,y)

sumx = 0;
sumy = 0;
sumxy = 0;
sumx2 = 0;
sumy2 = 0;
r = 0;

nc = length(x);
na = length(y);
if na == nc
    n = nc;
    for i = 1:n
        sumx = sumx+x(i);
    end
    for i = 1:n
        sumy = sumy+y(i);
    end
    for i = 1:n
        sumxy = sumxy+(x(i)*y(i));
    end
    for i = 1:n
        sumx2 = sumx2+(x(i))^2;
    end
    for i = 1:n
        sumy2 = sumy2+(y(i)^2);
    end
    a=((sumx2*sumy)-(sumxy*sumx))/((n*sumx2)-(sumx)^2);
    b=((n*sumxy)-(sumx*sumy))/((n*sumx2)-(sumx^2));
    r=(n*sumxy-(sumx*sumy))/(sqrt(n*sumx2-((sumx)^2))*sqrt(n*sumy2-((sumy)^2)));
    r2=r^2;
else
    error('Dimensions dont match!');
end
end