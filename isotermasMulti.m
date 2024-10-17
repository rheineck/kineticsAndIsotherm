% Isotermas de Adsorção Multicomponentes
% Raphael Gilioli Heineck

%% Desenvolvimento
% 26/01/2021 - Primeira Versão
% 27/01/2021 - Código pronto, sem testes
%            - Ajustar carregaParametros para coincidir com a tabela
%               do excel 'adsorptionData.xlsx'
% 28/01/2021 - Código Pronto para ser utilizado, devidamente testado


%% Função Principal - Modelos de Adsorção para Multicomponentes

function isotermasMulti()

diary resultadosIsotermasMulti.txt;

close;
close all;
clear;
clearvars -global;
clc;
tic;

disp('Cálculo de Isotermas Multicomponentes- v1.0');
disp('Autor: Raphael Gilioli Heineck, Eng.º Químico');

global par V cexqeA cexqeB qmaxLA kLA ceqeA ceA r2LA qmaxLB ...
    kLB ceqeB ceB r2LB kFA nFA lnqeA lnceA r2FA kFB nFB lnqeB lnceB r2FB ...
    qmaxSA kSA qeSA ceSA r2SA qmaxSB kSB qeSB ceSB r2SB qeA qeB ...
    ajuste maior langmuir freundlich sips plotL plotF plotS plotM ceLA ceLB ...
    qmaxL kL r2L kF nF r2F qmaxS kS r2S qA qB qepA qepB resultadoA resultadoB

% Carrega Parametros
    disp('Etapa 01 - Carregamento de Parametros...');
    [par,V]=carregaParametros();
    disp('Etapa Concluída');
    
% Cálculos Preliminares
    disp('Etapa 02 - Cálculos Preliminares...');
    [cexqeA,cexqeB]=calculosPreliminares(par,V);
    disp('Etapa Concluída');
    
% Ajuste pelos Modelos Matemáticos
    disp('Etapa 03 - Ajustando pelos Modelos Matemático...');
    disp('Etapa 03.1 - Modelo de Langmuir...');
    [qmaxLA,kLA,ceqeA,ceLA,r2LA]=modeloLangmuir(cexqeA);
    [qmaxLB,kLB,ceqeB,ceLB,r2LB]=modeloLangmuir(cexqeB);
    qmaxL=[qmaxLA qmaxLB];
    kL=[kLA kLB];
    r2L=[r2LA r2LB];
    disp('Etapa 03.2 - Modelo de Freundlich...');
    [kFA,nFA,lnqeA,lnceA,r2FA]=modeloFreundlich(cexqeA);
    [kFB,nFB,lnqeB,lnceB,r2FB]=modeloFreundlich(cexqeB);
    kF=[kFA kFB];
    nF=[nFA nFB];
    r2F=[r2FA r2FB];
    disp('Etapa 03.3 - Modelo de Sips...');
    [qmaxSA,kSA,qeSA,ceSA,r2SA]=modeloSips(cexqeA);
    [qmaxSB,kSB,qeSB,ceSB,r2SB]=modeloSips(cexqeB);
    qmaxS=[qmaxSA qmaxSB];
    kS=[kSA kSB];
    r2S=[r2SA r2SB];
    disp('Etapa 03.4 - Modelo Competitivo de Langmuir...');
    [qeA,qeB,ceA,ceB,qA,qB]=modeloLangmuirCompetitivo(qmaxLA,qmaxLB,kLA,kLB,ceLA,ceLB,qepA,qepB);
    disp('Etapa Concluida');
    
% Escolha do melhor ajuste
    disp('Etapa 04 - Escolhendo Melhor Ajuste...');
    [ajuste,maior,langmuir,freundlich,sips]=melhorAjustev2(qmaxL,kL,r2L,kF,nF,r2F,qmaxS,kS,r2S);
    disp('Etapa Concluída');
    
% Interpretação Modelo COmpetitivo de Langmuir
    disp('Etapa 05 - Interpretando Modelo Binário Langmuir...');
    [qA,qB,resultadoA,resultadoB]=interpretacaoLangmuir(qA,qB);
    disp('Etapa Concluída');
    
% Plotagem dos gráficos
    disp('Etapa 05 - Plotando Gráficos...');
    [plotL,plotF,plotS,plotM]=plotaGraficos;
    disp('Etapa Concluída');
    
%     close all
    
% Armazenamento de dados
    disp('Etapa 06 - Armazenando dados...');
    armazenaDados(ajuste,maior,langmuir,freundlich,qA,qB,resultadoA,resultadoB,sips,plotL,plotF,plotS,plotM);
    disp('Etapa Concluída');
    
    toc
    diary off;
        
end

%% Armazenamento de Dados

function armazenaDados(ajuste,maior,langmuir,freundlich,qA,qB,resultadoA,resultadoB,sips,plotL,plotF,plotS,plotM)

xlswrite('adsorptionData.xlsx',maior(1),'resultadoIsotermasMulti','N3');
xlswrite('adsorptionData.xlsx',maior(2),'resultadoIsotermasMulti','Q3');
xlswrite('adsorptionData.xlsx',ajuste(1),'resultadoIsotermasMulti','M2');
xlswrite('adsorptionData.xlsx',ajuste(2),'resultadoIsotermasMulti','P2');

xlswrite('adsorptionData.xlsx',qA','resultadoIsotermasMulti','Q6');
xlswrite('adsorptionData.xlsx',qB','resultadoIsotermasMulti','R6');
xlswrite('adsorptionData.xlsx',resultadoA,'resultadoIsotermasMulti','S4');
xlswrite('adsorptionData.xlsx',resultadoB,'resultadoIsotermasMulti','S5');

[n,m]=size(langmuir);
for n=1:n
    if n==1
        qmaxL(n)=langmuir(n,m-1);
        kL(n)=langmuir(n+1,m-1);
        r2L(n)=langmuir(n+2,m-1);
    else
        if n==2
            qmaxL(n)=langmuir(n-1,m);
            kL(n)=langmuir(n,m);
            r2L(n)=langmuir(n+1,m);
        end
    end
end
clear n m
[n,m]=size(freundlich);
for n=1:n
    if n==1
        kF(n)=freundlich(n,m-1);
        nF(n)=freundlich(n+1,m-1);
        r2F(n)=freundlich(n+2,m-1);
    else
        if n==2
            kF(n)=freundlich(n-1,m);
            nF(n)=freundlich(n,m);
            r2F(n)=freundlich(n+1,m);
        end
    end
end
clear n m
[n,m]=size(sips);
for n=1:n
    if n==1
        qmaxS(n)=sips(n,m-1);
        kS(n)=sips(n+1,m-1);
        r2S(n)=sips(n+2,m-1);
    else
        if n==2
            qmaxS(n)=sips(n-1,m);
            kS(n)=sips(n,m);
            r2S(n)=sips(n+1,m);
        end
    end
end

% Escreve os resultados
% Langmuir
xlswrite('adsorptionData.xlsx',qmaxL(1),'resultadoIsotermasMulti','B4');
xlswrite('adsorptionData.xlsx',qmaxL(2),'resultadoIsotermasMulti','C4');
xlswrite('adsorptionData.xlsx',kL(1),'resultadoIsotermasMulti','B5');
xlswrite('adsorptionData.xlsx',kL(2),'resultadoIsotermasMulti','C5');
xlswrite('adsorptionData.xlsx',r2L(1),'resultadoIsotermasMulti','B6');
xlswrite('adsorptionData.xlsx',r2L(2),'resultadoIsotermasMulti','C6');

% Freundlich
xlswrite('adsorptionData.xlsx',kF(1),'resultadoIsotermasMulti','F4');
xlswrite('adsorptionData.xlsx',kF(2),'resultadoIsotermasMulti','G4');
xlswrite('adsorptionData.xlsx',nF(1),'resultadoIsotermasMulti','F5');
xlswrite('adsorptionData.xlsx',nF(2),'resultadoIsotermasMulti','G5');
xlswrite('adsorptionData.xlsx',r2F(1),'resultadoIsotermasMulti','F6');
xlswrite('adsorptionData.xlsx',r2F(2),'resultadoIsotermasMulti','G6');

% sips
xlswrite('adsorptionData.xlsx',qmaxS(1),'resultadoIsotermasMulti','J4');
xlswrite('adsorptionData.xlsx',qmaxS(2),'resultadoIsotermasMulti','K4');
xlswrite('adsorptionData.xlsx',kS(1),'resultadoIsotermasMulti','J5');
xlswrite('adsorptionData.xlsx',kS(2),'resultadoIsotermasMulti','K5');
xlswrite('adsorptionData.xlsx',r2S(1),'resultadoIsotermasMulti','J6');
xlswrite('adsorptionData.xlsx',r2S(2),'resultadoIsotermasMulti','K6');

% Armazena Dados para posterior plotagem
xlswrite('adsorptionData.xlsx',plotL,'resultadoIsotermasMulti','A10');
xlswrite('adsorptionData.xlsx',plotF,'resultadoIsotermasMulti','E10');
xlswrite('adsorptionData.xlsx',plotS,'resultadoIsotermasMulti','I10');
xlswrite('adsorptionData.xlsx',plotM,'resultadoIsotermasMulti','M6');
end

%% Plotagem de Gráficos

function [plotL,plotF,plotS,plotM]=plotaGraficos

global ceqeA ceLA lnqeA lnceA qeSA ceSA ceqeB ceLB lnqeB lnceB qeSB ceSB ...
    ceA ceB qeA qeB

qeA=qeA';
qeB=qeB';

plotL=[ceLA,ceqeA,ceLB,ceqeB];
plotF=[lnceA,lnqeA,lnceB,lnqeB];
plotS=[ceSA,qeSA,ceSB,qeSB];
plotM=[ceA,qeA,ceB,qeB];

[~,j]=size(plotL);
[~,n]=size(plotF);
[~,m]=size(plotS);
[~,o]=size(plotM);
if j==n && n==m && n==o
    for j=1:j
        if j==3
            break
        end
        if j==2
            j=j+1;
        end
        if j==4
            break
        end
        figure();
        if j==1
            plot(ceLA,ceqeA);
            title('Modelo de Langmuir - A');
        else
            if j==3
                plot(ceLB,ceqeB);
                title('Modelo de Langmuir - B');
            end
        end
        ylabel('Ce/Qe (mg L^-1 mg^-1 g^-1');
        xlabel('Ce (mg L^-1)');
        grid;

        figure();
        if j==1
            plot(lnceA,lnqeA);
            title('Modelo de Freundlich - A');
        else
            if j==3
                plot(lnceB,lnqeB);
                title('Modelo de Freundlich - B');
            end
        end
        xlabel('ln Ce (mg L^-1)');
        ylabel('ln qe (mg g^-1)');
        grid;

        figure();
        if j==1
            plot(ceSA,qeSA);
            title('Modelo de Sips - A');
        else
            if j==3
                plot(ceSB,qeSB)
                title('Modelo de Sips - B');
            end
        end
        xlabel('(1/Ce)^(1/n) (L mg^-1)');
        ylabel('1/qt (g mg^-1)');
        grid;
        
        figure();
        if j==1
            plot(ceA,qeA);
            title('Modelo de Langmuir Multicomponente - A');
        else
            if j==3
                plot(ceB,qeB);
                title('Modelo de Langmuir Multicomponente - B');
            end
        end
        ylabel('qe (mg g^-1)');
        xlabel('Ce (mg L^-1)');
        grid;
    end
end

clear ceqeA ceA lnqeA lnceA qeSA ceSA ceqeB ceB lnqeB lnceB qeSB ceSB
end

%% Interpretação Modelo Competitivo de Langmuir

function [qA,qB,resultadoA,resultadoB]=interpretacaoLangmuir(qA,qB)

if qA > 1
    disp('O processo de adsorção é aumentado pela');
    disp('presença de outro composto');
    disp(['qA = ',num2str(qA)]);
    resultadoA={'Aumenta'};
else
    if qA == 0
        disp('O processo de adsorção não sofre');
        disp('interferência pelo outro composto');
        disp(['qA = ',num2str(qA)]);
        resultadoA={'Não Interfere'};
    else
        if qA < 1
            disp('O processo de adsorção é diminuído pela');
            disp('presença de outro composto');
            disp(['qA = ',num2str(qA)]);
            resultadoA={'Diminui'};
        end
    end
end

if qB > 1
    disp('O processo de adsorção é aumentado pela');
    disp('presença de outro composto');
    disp(['qB = ',num2str(qB)]);
    resultadoB={'Aumenta'};
else
    if qB == 0
        disp('O processo de adsorção não sofre');
        disp('interferência pelo outro composto');
        disp(['qB = ',num2str(qB)]);
        resultadoB={'Não Interfere'};
    else
        if qB < 1
            disp('O processo de adsorção é diminuído pela');
            disp('presença de outro composto');
            disp(['qB = ',num2str(qB)]);
            resultadoB={'Diminui'};
        end
    end
end

end

%% Escolha do Melhor Ajuste v2

function [ajuste,maior,langmuir,freundlich,sips]=melhorAjustev2(qmaxL,kL,r2L,kF,nF,r2F,qmaxS,kS,r2S)

% Assumindo que todas estão com as mesmas dimensões
langmuir=[qmaxL; kL; r2L];
% disp(langmuir);
% disp(kF);
% disp(nF);
% disp(r2F);
freundlich=[kF; nF; r2F];
sips=[qmaxS; kS; r2S];
% disp(freundlich);
% disp(sips);
ajuste=0;
maior=0;
n=length(r2L);
m=length(r2F);
o=length(r2S);
if n~=m && m~=o
    error('Dimensões não condizem');
else
    for i=1:n
        if i==1
            if r2L(i) > r2F(i) && r2L(i) > r2S(i)
                maiorA=r2L(i);
                ajusteA={'Modelo de Langmuir'};
            else
                if r2F(i) > r2L(i) && r2F(i) > r2S(i)
                    maiorA=r2F(i);
                    ajusteA={'Modelo de Freundlich'};
                else
                    if r2S(i) > r2L(i) && r2S(i) > r2F(i)
                        maiorA=r2S(i);
                        ajusteA={'Modelo de Sips'};
                    end
                end
            end
        else
            if i==2
                if r2L(i) > r2F(i) && r2L(i) > r2S(i)
                    maiorB=r2L(i);
                    ajusteB={'Modelo de Langmuir'};
                else
                    if r2F(i) > r2L(i) && r2F(i) > r2S(i)
                        maiorB=r2F(i);
                        ajusteB={'Modelo de Freundlich'};
                    else
                        if r2S(i) > r2L(i) && r2S(i) > r2F(i)
                            maiorB=r2S(i);
                            ajusteB={'Modelo de Sips'};
                        end
                    end
                end
            end
        end
    end
end
maior=[maiorA maiorB];
ajuste=[ajusteA ajusteB];
% disp(maior);
% disp(ajuste);
for i=1:n
    if maior(i) == r2L(i)
        disp('Modelo de Isoterma de Langmuir');
        disp('se ajusta melhor aos dados estudados');
        disp(['r^2 = ',num2str(r2L(i))]);
        disp(['kL = ',num2str(kL(i))]);
        disp(['qmax = ',num2str(qmaxL(i))]);
    else
        if maior(i) == r2F(i)
            disp('Modelo de Isoterma de Freundlich');
            disp('se ajusta melhor aos dados estudados');
            disp(['r^2 = ',num2str(r2F(i))]);
            disp(['kF = ',num2str(kF(i))]);
            disp(['nF = ',num2str(nF(i))]);
        else
            if maior(i) == r2S(i)
                disp('Modelo de Langmuir-Freundlich (SIPS)');
                disp('se ajusta melhor ao modelo estudado');
                disp(['r^2 = ',num2str(r2S(i))]);
                disp(['kS = ',num2str(kS(i))]);
                disp(['qmaxS = ',num2str(qmaxS(i))]);
            end
        end
    end
end
end

%% Modelo Competitivo de Langmuir

function [qeA,qeB,ceA,ceB,qA,qB]=modeloLangmuirCompetitivo(qmaxLA,qmaxLB,kLA,kLB,ceLA,ceLB,qepA,qepB)

n=length(ceLA);
m=length(ceLB);
if n==m
    for i=1:n
        qeA(i)=(qmaxLA*kLA*ceLA(i))/(1+kLA*ceLA(i)+kLB*ceLB(i));
        qeB(i)=(qmaxLB*kLB*ceLB(i))/(1+kLA*ceLA(i)+kLB*ceLB(i));
    end
else
    error('Dimensões não condizem');
end
clear n m
ceA=ceLA;
ceB=ceLB;

n=length(qepA);
m=length(qepB);
if n==m
    for i=1:n
        qA(i)=qeA(i)/qepA(i);
        qB(i)=qeB(i)/qepB(i);
    end
end

end

%% Modelo de Sips (Langmuir-Freundlich)

function [qmaxS,kS,qeS,ceS,r2S]=modeloSips(cexqe)

% y = 1/qe
% x = (1/ce)^(1/nS)
% a = 1/qmax*kS
% b = 1/qmax

[i,j]=size(cexqe);
if j~=2
    error('Apenas duas colunas são suportadas');
end
for i=1:i
    ce(i)=cexqe(i,1);
    qe(i)=cexqe(i,2);
end
n=length(ce);
for i=1:n
    y(i)=(1/qe(i));
    x(i)=(1/ce(i));
end
[a,b,r2]=linearizacao(y,x);
qmaxS=1/b;
kS=1/(qmaxS*a);
qeS=y';
ceS=x';
r2S=r2;
%disp(r2S);
clear x y
end

%% Modelo de Freundlich

function [kF,nF,lnqe,lnce,r2F]=modeloFreundlich(cexqe)

% y = lnqe
% x = lnce
% a = 1/N
% b = lnkF

[i,j]=size(cexqe);
if j~=2
    error('Apenas duas colunas são suportadas');
end
for i=1:i
    ce(i)=cexqe(i,1);
    qe(i)=cexqe(i,2);
end
n=length(ce);
for i=1:n
    y(i)=log(qe(i));
    x(i)=log(ce(i));
end
[a,b,r2]=linearizacao(y,x);
kF=exp(b);
nF=1/a;
r2F=r2;
lnqe=y';
lnce=x';
%disp(r2F);
clear x y
end

%% Modelo de Langmuir

function [qmax,kL,ceqe,ce,r2L]=modeloLangmuir(cexqe)

% a = 1/qmax
% b = 1/kL*qmax
% x = Ce
% y = Ce/qe (variável ceqe)

[i,j]=size(cexqe);
if j~=2
    error('Apenas duas colunas são suportadas');
end
for i=1:i
    ce(i)=cexqe(i,1);
    qe(i)=cexqe(i,2);
end
n=length(ce);
for i=1:n
    y(i)=ce(i)/qe(i);
    x(i)=ce(i);
end
[a,b,r2]=linearizacao(y,x);
qmax=1/a;
kL=1/(b*qmax);
ceqe = y';
ce = x';
r2L=r2;
%disp(r2L);
clear x y
end

%% Linearização

function [a,b,r2]=linearizacao(y,x)

global somax somay somaxy somax2 somay2 r

somax=0;
somay=0;
somaxy=0;
somax2=0;
somay2=0;
r=0;

nc=length(x);
na=length(y);
if na==nc
    n=nc;
    for i=1:n
        somax=somax+x(i);
    end
    for i=1:n
        somay=somay+y(i);
    end
    for i=1:n
        somaxy=somaxy+(x(i)*y(i));
    end
    for i=1:n
        somax2=somax2+(x(i))^2;
    end
    for i=1:n
        somay2=somay2+(y(i)^2);
    end
    a=((somax2*somay)-(somaxy*somax))/((n*somax2)-(somax)^2);
    b=((n*somaxy)-(somax*somay))/((n*somax2)-(somax^2));
    r=(n*somaxy-(somax*somay))/(sqrt(n*somax2-((somax)^2))*sqrt(n*somay2-((somay)^2)));
    r2=r^2;
else
    error('Dimensões não condizentes!');
end
end

%% Calculos Preliminares

function [cexqeA,cexqeB]=calculosPreliminares(par,V)

global mass qepA qepB ciA ciB cfA cfB
% j = 1 -> Cf(PCT)
% j = 2 -> Cf(DCF)
% j = 3 -> mass

% qepA e qepB são parametros preliminares, necessários no cálculos do
% modelo competitivo de langmuir
% 
% ciA=ci(1);
% ciB=ci(2);
% clear ci

[i,j]=size(par);
if j~=5
    error('Cinco colunas são permitidas!');
end
for i = 1:i
    ciA(i)=par(i,1);
    cfA(i)=par(i,2);
    ciB(i)=par(i,3);
    cfB(i)=par(i,4);
    mass(i)=par(i,5);
end
V=V*(10^-3);
n=length(ciA);
m=length(ciB);
if n==m
    for i=1:n
        qepA(i)=((ciA(i)-cfA(i))*V)/mass(i);
        qepB(i)=((ciB(i)-cfB(i))*V)/mass(i);
    end
else
    error('Dimensões não condizem');
end
cexqeA=[cfA' qepA'];
cexqeB=[cfB' qepB'];
% 
% disp(cexqeA);
% disp(cexqeB);
% 
% [i,j]=size(cexqeA);
% [m,n]=size(cexqeB);
% disp(i);
% disp(j);
% disp(m);
% disp(n);
end

%% Carregamento dos Parametros

function [par,V]=carregaParametros()

global conc qe ce

conc=0;
qe=0;
ce=0;
% ciA=xlsread('adsorptionData.xlsx','dadosIsoterma','N2');
% ciB=xlsread('adsorptionData.xlsx','dadosIsoterma','N3');
V=xlsread('adsorptionData.xlsx','dadosIsoterma','C3');
par=xlsread('adsorptionData.xlsx','dadosIsoterma','M19:Q22');
% ci=[ciA; ciB];

end