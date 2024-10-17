function [r2] = determinationCoefficient(x,y,fun,par)

sumy2 = 0;
sumd = 0;
sumd2 = 0;

nc = length(x);
na = length(y);

ycalc = fun(par,x);
sumycalc = sum(ycalc);

if na == nc
    n = nc;
    for i = 1:n
        sumy2 = sumy2 + (ycalc(i)^2);
    end
    for i = 1:n
        d(i) = y(i) - ycalc(i);
        sumd = sumd + d(i);
    end
    for i = 1:n
        sumd2 = sumd2 + (d(i))^2;
    end
else
    error('Dimensions dont match!');
end

n = length(x);
r2 = 1 - (sumd2) / (sumy2 - ((sumycalc)^2/n));

end