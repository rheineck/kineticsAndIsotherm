function cexqe = adsorptionCapacity(par)

ci = par(:,1);
cf = par(:,2);
m = par(:,3);
V = xlsread('adsorptionData.xlsx','dadosIsoterma','D3');
n = length(ci);
qe = zeros(n,1);

for n = 1:n
    qe(n) = ((ci(n) - cf(n))*V)/m(n);
end

ce = cf;
cexqe = [ce qe];

end