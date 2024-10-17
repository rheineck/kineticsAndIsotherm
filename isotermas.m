% Isotermas de Adsorção
% Raphael Gilioli Heineck

%% Função Principal - Modelos de Isoterma Linear

function isotermas()

diary resultadoIsotermas.txt;
% Função Principal

close;
close all;
clear;
clearvars -global;
clc;
tic;

disp('Cálculo de Isotermas - versão 1.0');

global cexqe qmaxL kL ceqe ce r2L kF nF lnqe lnce r2F qmaxS kS qeS ceS r2S
global langmuir freundlich sips plotL plotF plotS par V ajuste maior

% Carrega Parametros de uma planilha de Excel
    disp('Etapa 01 - Carregando Parametros...');
    [par,V]=carregaParametros();
    disp('Etapa Concluída');

% Efetua cálculos preliminares
    disp('Etapa 02 - Efetuando Cálculos Preliminares...');
    [cexqe]=calculosPreliminares(par,V);
    disp(cexqe);
    disp('Etapa Concluída');
    
% Ajusta pelo Modelo de Langmuir
    disp('Etapa 03 - Ajustando pelos Modelos Matemáticos...');
    disp('Etapa 03.1 - Modelo de Langmuir...');
    [qmaxL,kL,ceqe,ce,r2L]=modeloLangmuir(cexqe);
    disp('Etapa 03.2 - Modelo de Freundlich...');
    [kF,nF,lnqe,lnce,r2F]=modeloFreundlich(cexqe);
    disp('Etapa 03.3 - Modelo de Sips...');
    [qmaxS,kS,qeS,ceS,r2S]=modeloSips(cexqe);
    disp('Etapa Concluida');

% Esxolha do melhor ajuste
    disp('Etapa 04 - Escolha melhor ajuste...');
    [ajuste,maior,langmuir,freundlich,sips]=melhorAjustev2(qmaxL,kL,r2L,kF,nF,r2F,qmaxS,kS,r2S);
    disp('Etapa Concluida');
    
% Plota gráficos
    disp('Etapa 05 - Plotagem gráficos...');
    [plotL,plotF,plotS]=plotaGraficos;
    disp('Etapa Concluida');
    
% Armazenando dados
    disp('Etapa 06 - Armazenando Parâmetros...');
    armazenaDados(ajuste,maior,langmuir,freundlich,sips,plotL,plotF,plotS);
    disp('Execução finalizada com sucesso!');
    
    toc;
    diary off
    
end

%% Armazenamento de Dados

function armazenaDados(ajuste,maior,langmuir,freundlich,sips,plotL,plotF,plotS)

xlswrite('adsorptionData.xlsx',maior,'resultadoIsotermas','K3');
xlswrite('adsorptionData.xlsx',ajuste,'resultadoIsotermas','J2');

n=length(langmuir);
for i=1:n
    if i==1
        qmaxL=langmuir(i);
    else
        if i==2
            kL=langmuir(i);
        else
            if i==3
                r2L=langmuir(i);
            end
        end
    end
end
xlswrite('adsorptionData.xlsx',qmaxL,'resultadoIsotermas','B3');
xlswrite('adsorptionData.xlsx',kL,'resultadoIsotermas','B4');
xlswrite('adsorptionData.xlsx',r2L,'resultadoIsotermas','B5');
clear n
n=length(freundlich);
for i=1:n
    if i==1
        kF=freundlich(i);
    else
        if i==2
            nF=freundlich(i);
        else
            if i==3
                r2F=freundlich(i);
            end
        end
    end
end
xlswrite('adsorptionData.xlsx',kF,'resultadoIsotermas','E4');
xlswrite('adsorptionData.xlsx',nF,'resultadoIsotermas','E3');
xlswrite('adsorptionData.xlsx',r2F,'resultadoIsotermas','E5');
clear n
n=length(sips);
for i=1:n
    if i==1
        qmaxS=sips(i);
    else
        if i==2
            kS=sips(i);
        else
            if i==3
                r2S=sips(i);
            end
        end
    end
end
xlswrite('adsorptionData.xlsx',qmaxS,'resultadoIsotermas','H3');
xlswrite('adsorptionData.xlsx',kS,'resultadoIsotermas','H4');
xlswrite('adsorptionData.xlsx',r2S,'resultadoIsotermas','H5');

% Armazena Dados para posterior plotagem
xlswrite('adsorptionData.xlsx',plotL,'resultadoIsotermas','A7');
xlswrite('adsorptionData.xlsx',plotF,'resultadoIsotermas','D7');
xlswrite('adsorptionData.xlsx',plotS,'resultadoIsotermas','G7');
end
%% Plotagem de Gráficos

function [plotL,plotF,plotS]=plotaGraficos

global ceqe ce lnqe lnce qeS ceS

plotL=[ce,ceqe];
plotF=[lnce,lnqe];
plotS=[ceS,qeS];

figure();
plot(ce,ceqe);
title('Modelo de Langmuir');
ylabel('Ce/Qe (mg L^-1 mg^-1 g^-1');
xlabel('Ce (mg L^-1)');
grid;

figure();
plot(lnqe,lnce);
title('Modelo de Freundlich');
xlabel('ln Ce (mg L^-1)');
ylabel('ln qe (mg g^-1)');
grid;

figure();
plot(qeS,ceS);
title('Modelo de Sips');
xlabel('(1/Ce)^(1/n) (L mg^-1)');
ylabel('1/qt (g mg^-1)');
grid;

clear ce ceqe lnqe lnce qeS ceS
end

%% Escolha do Melhor Ajuste v2

function [ajuste,maior,langmuir,freundlich,sips]=melhorAjustev2(qmaxL,kL,r2L,kF,nF,r2F,qmaxS,kS,r2S)

langmuir=[qmaxL; kL; r2L];
freundlich=[kF; nF; r2F];
sips=[qmaxS; kS; r2S];
ajuste=0;
maior=0;
if r2L > r2F && r2L > r2S
    maior=r2L;
    ajuste={'Modelo de Langmuir'};
else
    if r2F > r2L && r2F > r2S
        maior=r2F;
        ajuste={'Modelo de Freundlich'};
    else
        if r2S > r2L && r2S > r2F
            maior=r2S;
            ajuste={'Modelo de Sips'};
        end
    end
end

if maior == r2L
    disp('Modelo de Isoterma de Langmuir');
    disp('se ajusta melhor aos dados estudados');
    disp(['r^2 = ',num2str(r2L)]);
    disp(['kL = ',num2str(kL)]);
    disp(['qmax = ',num2str(qmaxL)]);
else
    if maior == r2F
        disp('Modelo de Isoterma de Freundlich');
        disp('se ajusta melhor aos dados estudados');
        disp(['r^2 = ',num2str(r2F)]);
        disp(['kF = ',num2str(kF)]);
        disp(['nF = ',num2str(nF)]);
    else
        if maior == r2S
            disp('Modelo de Langmuir-Freundlich (SIPS)');
            disp('se ajusta melhor ao modelo estudado');
            disp(['r^2 = ',num2str(r2S)]);
            disp(['kS = ',num2str(kS)]);
            disp(['qmaxS = ',num2str(qmaxS)]);
        end
    end
end
% p={ajuste; maior};
end

%% Escolha do Melhor Ajuste v1
% function [langmuir,freundlich,sips]=melhorAjuste(qmaxL,kL,r2L,kF,nF,r2F,qmaxS,kS,r2S)
% 
% langmuir=[qmaxL; kL; r2L];
% freundlich=[kF; nF; r2F];
% sips=[qmaxS; kS; r2S];
% 
% if r2L>=0.8000 && r2L<=1
%     disp('Modelo de Isoterma de Langmuir');
%     disp('se ajusta melhor aos dados estudados!');
%     disp(['R^2 = ',num2str(r2L)]);
%     disp(['kL = ',num2str(kL)]);
%     disp(['qmax = ',num2str(qmaxL)]);
%     %p1 = kL;
%     %p2 = qmaxL;
%     %p3 = r2L;
%     %p4 = "Modelo de Langmuir";
%     %p=[p1; p2; p3; p4];
% else
%     if r2F>=0.8000 && r2F<=1
%         disp('Modelo de Isoterma de Freundlich');
%         disp('se ajusta melhor aos dados estudados!');
%         disp(['R^2 = ',num2str(r2F)]);
%         disp(['kF = ',num2str(kF)]);
%         disp(['nF = ',num2str(nF)]);
%         %p1 = kF;
%         %p2 = nF;
%         %p3 = r2F;
%         %p4 = "Modelo de Freundlich";
%         %p=[p1; p2; p3; p4];
%     else
%         if r2S>=0.8000 && r2S<=1
%             disp('Modelo de Isoterma de Sips');
%             disp('se ajusta melhor aos dados estudados!');
%             disp(['R^2 = ',num2str(r2S)]);
%             disp(['kS = ',num2str(kS)]);
%             disp(['qmax = ',num2str(qmaxS)]);
%             %p1 = kS;
%             %p2 = qmaxS;
%             %p3 = r2S;
%             %p4 = "Modelo de Sips";
%             %p=[p1; p2; p3; p4];
%         end
%     end
% end
% clear qmaxL kL r2L kF nF r2F qmaxS kS r2S
% end

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

%% Cálculos Preliminares

function [cexqe]=calculosPreliminares(par,V)

global ci cf ce qe mass

% j = 1 -> Ci
% j = 2 -> Cf
% j = 3 -> mass

[i,j]=size(par);
if j~=3
    error('Apenas três colunas de dados são permitidas!');
end
for i=1:i
    ci(i)=par(i,1);
    cf(i)=par(i,2);
    mass(i)=par(i,3);
end
disp(ci);
disp(cf);
disp(mass);
n=length(ci);
for i=1:n
    qe(i)=((ci(i)-cf(i))*V)/mass(i);
    disp(qe);
end
ce=cf';
qe=qe';
cexqe=[ce qe];
end

%% Carregamento de Parametros

function [par,V]=carregaParametros()

global conc qe ce 

conc=0;
qe=0;
ce=0;
V=xlsread('adsorptionData.xlsx','dadosIsoterma','C3');
par=xlsread('adsorptionData.xlsx','dadosIsoterma','M6:O9');
% par=xlsread('adsorptionData.xlsx','dadosIsoterma','P6:R9');

end