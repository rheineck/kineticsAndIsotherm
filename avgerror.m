function er = avgerror(xmod,xexp)

sum = 0;
n = length(xmod);
m = length(xexp);

if n == m
    for i=1:n
        sum = sum + (abs(xmod(i)-xexp(i))/xexp(i));
    end
    n = length(xmod);
    er = 100*(1/n)*sum;
else
    error('Dimensions dont match!');
end
end