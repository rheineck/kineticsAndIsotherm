function [qmax,kL,ceqe,ce,r2L]=modeloLangmuir(cexqe)

% a = 1/qmax
% b = 1/kL*qmax
% x = Ce
% y = Ce/qe (variável ceqe)

[i,~]=size(cexqe);
% if j~=2
%     error('Apenas duas colunas são suportadas');
% end
for i=1:i
    ce(i)=cexqe(i,1);
    qe(i)=cexqe(i,2);
end
n=length(ce);
for i=1:n
    y(i)=ce(i)/qe(i);
    x(i)=ce(i);
end
[a,b,r2]=linRegression(y,x);
qmax=1/a;
kL=1/(b*qmax);
ceqe = y';
ce = x';
r2L=r2;
%disp(r2L);
clear x y
end