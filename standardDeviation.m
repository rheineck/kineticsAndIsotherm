% Standard Deviation
function sd = standardDeviation(x,y,E,par)

sum = 0;
n = length(y);
ycalc = zeros(n,1);

for n = 1:n
    ycalc(n) = E(par,x);
end

for n = 1:n
    sum = sum + ((y-ycalc)/y)^2;
end

sd = sqrt(sum/(n-1));

end