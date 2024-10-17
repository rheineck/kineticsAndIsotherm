% Cin�tica de Adsor��o
% Raphael Gilioli Heineck

%% Fun��o Principal - Modelos Cin�ticos Lineares

function cineticaLinear()

diary resultadoCineticas.txt;
% Fun��o Principal
tic;
close;
close all;
clear;
clearvars -global;
clc;

disp('C�lculo de Cin�ticas de Adsor��o - Ajuste Linear - v2.0a');

global txqt lnq ln k r2p t tqt ks r2s C kd r2d qt rt qe plotP plotS plotD
global ajuste maior pseudoPrimeira pseudoSegunda difusaoIntraparticula

% Carregar Parametros de uma planilha do excel
    disp('Etapa 01 - Carregando Par�metros...');
    carregarParametros;
    disp('Etapa Conclu�da!');
    
% Calculos Preliminares a adsor��o
    disp('Etapa 02 - Executando C�lculos Preliminares...');
    [txqt]=calculosPreliminares;
    disp('Etapa Conclu�da!');
    
% Calculos utilizando modelos Cin�ticos
    disp('Etapa 03 - Calculando Cin�ticas...');
    disp('Etapa 03.1 - Pseudo Primeira Ordem...');
    [lnq,ln,k,r2p,t]=pseudoPrimeiraOrdem(txqt);
    disp('Etapa 03.2 - Pseudo Segunda Ordem...');
    [qe,tqt,ks,r2s]=pseudoSegundaOrdem(txqt);
    disp('Etapa 03.3 - Difus�o Intrapart�cula...');
    [C,kd,r2d,qt,rt]=difusaoIntraparticul(txqt);
    disp('Etapa Conclu�da!');
    
% Descobrir qual modelo se ajusta melhor
    disp('Etapa 04 - Escolhendo melhor ajuste...');
    [ajuste,maior,pseudoPrimeira,pseudoSegunda,difusaoIntraparticula]=melhorAjustev2(lnq,k,r2p,qe,ks,r2s,C,kd,r2d);
%   [p1,p2,p3]=melhor_ajuste(k,r2p,ks,r2s,kd,r2d);
    disp('Etapa Conclu�da!');

% Plota Gr�ficos
    disp('Etapa 05 - Plotando Gr�ficos...');
    [plotP,plotS,plotD]=plotaGraficos;
    disp('Etapa Conclu�da!'); 

% Encaminhamento de Parametros
    disp('Etapa 06 - Armazenando Par�metros...');
    armazenaDados(ajuste,maior,pseudoPrimeira,pseudoSegunda,difusaoIntraparticula,plotP,plotS,plotD)
%   [pC]=armazenaParametros(p1,p2,p3);
    disp('Execu��o Finalizada com sucesso!');
    
% Volta para programa de Sele��o
    %selecao;
    
    toc;
    diary off
end

%% Armazenamento de Dados v2

function armazenaDados(ajuste,maior,pseudoPrimeira,pseudoSegunda,difusaoIntraparticula,plotP,plotS,plotD)

xlswrite('adsorptionData.xlsx',maior,'resultadoCineticas','K3');
xlswrite('adsorptionData.xlsx',ajuste,'resultadoCineticas','J2');

n=length(pseudoPrimeira);
for i=1:n
    if i==1
        lnq=pseudoPrimeira(i);
    else
        if i==2
            k=pseudoPrimeira(i);
        else
            if i==3
                r2p=pseudoPrimeira(i);
            end
        end
    end
end
xlswrite('adsorptionData.xlsx',lnq,'resultadoCineticas','B3');
xlswrite('adsorptionData.xlsx',k,'resultadoCineticas','B4');
xlswrite('adsorptionData.xlsx',r2p,'resultadoCineticas','B5');
clear n
n=length(pseudoSegunda);
for i=1:n
    if i==1
        qe=pseudoSegunda(i);
    else
        if i==2
            ks=pseudoSegunda(i);
        else
            if i==3
                r2s=pseudoSegunda(i);
            end
        end
    end
end
xlswrite('adsorptionData.xlsx',qe,'resultadoCineticas','E4');
xlswrite('adsorptionData.xlsx',ks,'resultadoCineticas','E3');
xlswrite('adsorptionData.xlsx',r2s,'resultadoCineticas','E5');
clear n
n=length(difusaoIntraparticula);
for i=1:n
    if i==1
        C=difusaoIntraparticula(i);
    else
        if i==2
            kd=difusaoIntraparticula(i);
        else
            if i==3
                r2d=difusaoIntraparticula(i);
            end
        end
    end
end
xlswrite('adsorptionData.xlsx',C,'resultadoCineticas','H3');
xlswrite('adsorptionData.xlsx',kd,'resultadoCineticas','H4');
xlswrite('adsorptionData.xlsx',r2d,'resultadoCineticas','H5');

% Armazena Dados para posterior plotagem
xlswrite('adsorptionData.xlsx',plotP,'resultadoCineticas','A7');
xlswrite('adsorptionData.xlsx',plotS,'resultadoCineticas','D7');
xlswrite('adsorptionData.xlsx',plotD,'resultadoCineticas','G7');
end


%% Armazenamento de Dados v1

% function [pC] = armazenaParametros(p1,p2,p3)
% 
% % p1 = k
% % p2 = r^2
% % p3 = modelo
% 
% pC = [p1; p2; p3];
% 
% end

%% Plotagem de Gr�ficos

function [plotP,plotS,plotD]=plotaGraficos

global t ln tqt qt rt

plotP=[t,ln];
plotS=[t,tqt];
plotD=[rt,qt];

figure();
plot(t,ln);
title('Pseudo Primeira Ordem');
xlabel('Tempo (min)');
ylabel('ln(qe-qt)');
grid;

figure();
plot(t,tqt);
title('Pseudo Segunda Ordem');
xlabel('Tempo (min)');
ylabel('t/qt');
grid;

figure();
plot(rt,qt);
title('Difus�o Intrapart�cula');
xlabel('Tempo (raiz(min))');
ylabel('qt');
grid;

end

%% Modelo de Pseudo Primeira Ordem

function [lnq,ln,k,r2p,t]=pseudoPrimeiraOrdem(txqt)

clear qt;
clear qe;

% Pseudo Primeira Ordem
% y=ln => ln(qe-qt)
% a=lnq => ln(qe)
% b=k => k
% x=t => t

[i,j]=size(txqt);
if j~=2
    error('Apenas duas colunas de dados s�o permitidas!');
end
for i=1:i
    t(i)=txqt(i,1);
    qt(i)=txqt(i,2);
end
n=length(t);
qe=qt(n);
%disp(t);
for i=1:n
    ln(i)=log(qe-qt(i));
end
n=length(ln);
for i=n
    ln(i)=qe;
end
%disp(ln);
clear qt;
qt=ln';
[a,b,r2]=linearizacao(qt,t);
t=t';
r2p=r2;
ln=qt;
k=b;
lnq=a;
clear qe;

end

%% Modelo de Pseudo Segunda Ordem

function [qe,tqt,ks,r2s]=pseudoSegundaOrdem(txqt)

global t qt

clear t;
clear qe;
clear qt;

% Pseudo Segunda Ordem
% y=tqt => t/qt
% a=ks => 1/(ks.qe^2)
% b=qe => 1/qe
% x=t => t

[i,j]=size(txqt);
if j~=2
    error('Apenas duas colunas de dados s�o permitidas!');
end
for i=1:i
    t(i)=txqt(i,1);
    qt(i)=txqt(i,2);
end
n=length(t);
qe=qt(n);
for i=1:n
    tqt(i)=t(i)/qt(i);
end
clear qt;
qt=tqt';
[a,b,r2]=linearizacao(qt,t);
t=t';
r2s=r2;
tqt=qt;
qe=1/b;
ks=1/(a*qe^2);
end

%% Modelo de Difus�o Intrapart�cula

function [C,kd,r2d,qt,rt]=difusaoIntraparticul(txqt)

global t qe

clear t;
clear qe;

% Difus�o Intrapart�cula
% y=qt => qt
% a=kd => k
% x=rt => sqrt(t)
% b=C => C

[i,j]=size(txqt);
if j~=2
    error('Apenas duas colunas de dados s�o permitidas!');
end
for i=1:i
    t(i)=txqt(i,1);
    qt(i)=txqt(i,2);
end
n=length(t);
qe=qt(n);
for i=1:n
   rt(i)=(sqrt(t(i)));
end
rt=rt';
t=rt;
[a,b,r2]=linearizacao(qt,t);
t=t';
qt=qt';
r2d=r2;
C=b;
kd=a;
end

%% Escolha do melhor ajuste v2

function [ajuste,maior,pseudoPrimeira,pseudoSegunda,difusaoIntraparticula]=melhorAjustev2(lnq,k,r2p,qe,ks,r2s,C,kd,r2d)

pseudoPrimeira=[lnq; k; r2p];
pseudoSegunda=[qe; ks; r2s];
difusaoIntraparticula=[C; kd; r2d];
% disp(pseudoPrimeira);
% disp(pseudoSegunda);
% disp(difusaoIntraparticula);
maior=0;
ajuste=0;
if r2p > r2s && r2p > r2d
    maior=r2p;
    ajuste={'Pseudo Primeira Ordem'};
else
    if r2s > r2p && r2s > r2d
        maior=r2s;
        ajuste={'Pseudo Segunda Ordem'};
    else
        if r2d > r2p && r2d > r2s
            maior=r2d;
            ajuste={'Difus�o Intrapart�cula'};
        end
    end
end

if maior == r2p
    disp('Modelo cin�tico de Pseudo Primeira Ordem');
    disp('se ajusta melhor aos dados estudados');
    disp(['r^2 = ',num2str(r2p)]);
    disp(['ln(qe) = ',num2str(lnq)]);
    disp(['k = ',num2str(k)]);
else
    if maior == r2s
        disp('Modelo cin�tico de Pseudo Segunda Ordem');
        disp('se ajusta melhor aos dados estudados');
        disp(['r^2 = ',num2str(r2s)]);
        disp(['k = ',num2str(ks)]);
        disp(['qe = ',num2str(qe)]);
    else
        if maior == r2d
            disp('Modelo cin�tico de Difus�o Intraparticula');
            disp('se ajusta melhor ao modelo estudado');
            disp(['r^2 = ',num2str(r2d)]);
            disp(['k = ',num2str(kd)]);
            disp(['C = ',num2str(C)]);
        end
    end
end
end

%% Escolha do Melhor Ajuste v1

% function [p1,p2,p3]=melhor_ajuste(k,r2p,ks,r2s,kd,r2d)
% 
% if r2p>=0.9000 && r2p<=1
%     disp('Modelo cin�tico Pseudo Primeira Ordem');
%     disp('se ajusta melhor ao modelo estudado!');
%     disp(['R^2 = ',num2str(r2p)]);
%     disp(['k = ',num2str(k)]);
%     p1 = k;
%     p2 = r2p;
%     p3 = "Pseudo Primeira Ordem";
% else
%     if r2s>=0.9000 && r2s<=1
%         disp('Modelo cin�tico Pseudo Segunda Ordem');
%         disp('se ajusta melhor ao modelo estudado!');
%         disp(['R^2 = ',num2str(r2s)]);
%         disp(['k = ',num2str(ks)]);
%         p1 = ks;
%         p2 = r2s;
%         p3 = "Pseudo Segunda Ordem";
%     else
%         if r2d>=0.9000 && r2d<=1
%             disp('Modelo cin�tico Difus�o Intrapart�cula');
%             disp('se ajusta melhor ao modelo estudado!');
%             disp(['R^2 = ',num2str(r2d)]);
%             disp(['k = ',num2str(kd)]);
%             p1 = kd;
%             p2 = r2d;
%             p3 = "Difus�o Intrapart�cula";
%         end
%     end
% end
% end
%% Lineariza��o

function [a,b,r2]=linearizacao(qt,t)

global somax somay somaxy somax2 somay2 r

somax=0;
somay=0;
somaxy=0;
somax2=0;
somay2=0;
r=0;

nc=length(t);
na=length(qt);
if na==nc
    n=nc;
    for i=1:n
        somax=somax+t(i);
    end
    for i=1:n
        somay=somay+qt(i);
    end
    for i=1:n
        somaxy=somaxy+(t(i)*qt(i));
    end
    for i=1:n
        somax2=somax2+(t(i))^2;
    end
    for i=1:n
        somay2=somay2+(qt(i)^2);
    end
    a=((somax2*somay)-(somaxy*somax))/((n*somax2)-(somax)^2);
    b=((n*somaxy)-(somax*somay))/((n*somax2)-(somax^2));
    r=(n*somaxy-(somax*somay))/(sqrt(n*somax2-((somax)^2))*sqrt(n*somay2-((somay)^2)));
    r2=r^2;
else
    error('Dimens�es n�o condizentes!');
end
end

%% C�lculos Preliminares

function [txqt]=calculosPreliminares()

global tempo conc par concmae V qt mas

[i,j]=size(par);
if j~=3
    error('Apenas duas colunas de dados s�o permitidas!');
end
concmae=par(1,2);
for i=2:i
    tempo(i-1)=par(i,1);
    conc(i-1)=par(i,2);
    mas(i-1)=par(i,3);
end
V=V*(10^-3);
n=length(conc);
for i=1:n
    qt(i)=((concmae-conc(i))*V)/mas(i);
end
tempo=tempo';
qt=qt';
txqt=[tempo qt];
end

%% Carrega Parametros

function carregarParametros()

global par tempo conc qe V concmae

tempo=0;
conc=0;
qe=0;
concmae=xlsread('adsorptionData.xlsx','dadosCinetica','N2');
V=xlsread('adsorptionData.xlsx','dadosCinetica','C3');
par=xlsread('adsorptionData.xlsx','dadosCinetica','M6:O16');
% par=xlsread('adsorptionData.xlsx','dadosCinetica','P6:R16');

end