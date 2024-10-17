% Kinetics and Isotherms non-linear
% Raphael Gilioli Heineck

%% Development
% Version 1.0
% 04/02/2021 - Begin
% 08/02/2021 - All models are ready, best adjust choice (multi) to be
%              implemented. Plots and save data also.
% 09/02/2021 - Better adjust choice algorithm for single and
%              multicomponent systems are ready. Graophic plot is ready too.
% 10/02/2021 - Code is ready, but not tested.
% 18/02/2021 - Single component ready to use.
%            - Multi component ready to use.
%            - Final version.
%            - Validation through Software Minitab 19 - OK
%            - Code aproved

% Version 2.0
% 03/03/2021 - Changes in function graphPlotSingle and graphPlotMulti.
%              Algorithm has to plot the adjusted data as well.
% 08/03/2021 - Graph plot was changed. best adjust no longer available
%            - Best adjust available, new equation
% 30/03/2021 - Average error function implemented, to use in multicomponent
%              models.

%% Main Function

function adsorption()

diary adsorption.txt;

close;
close all;
clear;
clearvars -global;
clc;
tic;

global options PFO PSO

disp('Kinetics and Isotherms - v2.0');
disp('Author: Raphael Gilioli Heineck, Chemical Engineer');
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxIterations',1000,'FunctionTolerance',1e-99);
% disp(options);
p = input('Your data is single-component(1) or multi-component(2)? ');
% p = 2;
if p == 1
    disp('Single component choosen');
    disp('Step 01 - Loading Parameters...');
    [parK,~,parIS,V,concA,~] = loadParameters();
    disp('Loading Complete!');
    disp(['concA = ',num2str(concA)]);
    
    disp('Step 02 - Preliminar Calculation...');
    [txqt,cexqe] = preCalculationSingle(parK,parIS,V,concA);
    disp('Step Complete!');
    
    disp('Step 03 - Performing adjustments...');
    disp('Kinetic Adjusts...');
    disp('  Pseudo First Order...');
    [qePFO,k1,r2PFO,PFO] = pseudoFirstOrder(txqt);
    pseudoFirst = [qePFO; k1; r2PFO];
    disp('  Pseudo Second Order...');
    [qePSO,k2,r2PSO,PSO] = pseudoSecondOrder(txqt);
    pseudoSecond = [qePSO; k2; r2PSO];
    disp('Isotherms Adjusts...');
    disp('  Langmuir Model...');
    [qmax,kL,r2L]=langmuirModel(cexqe);
    langmuir = [qmax; kL; r2L];
    disp('  Freundlich Model...');
    [kF,n,r2F]=freundlichModel(cexqe);
    freundlich = [kF; n; r2F];
    disp('Step Complete');
    
    disp('Step 04 - Choose best fit models...');
    [biggerK,biggerI,adjustK,adjustI]=bestAdjustS(pseudoFirst,pseudoSecond,langmuir,freundlich);
    disp('Step Complete');
    
    disp('Step 05 - Graphics Plot...');
    [plotK,plotI] = graphPlotSingle(txqt,cexqe);
    disp('Plot Conclude');
    
    disp('Step 06 - Saving Data...');
    saveDataSingle(plotK,plotI,biggerK,adjustK,biggerI,adjustI,pseudoFirst,pseudoSecond,langmuir,freundlich);
    disp('Step Complete!');
    
else
    if p == 2
        disp('Multi-component choosen');
        disp('Step 01 - Loading Parameters...');
        [parK,parIM,parIS,V,concA,concB] = loadParameters();
        disp('Loading Complete!');
        
        disp('Step 02 - Preliminar Calculation...');
        [txqtA,txqtB,cexqeA,cexqeB,cexqeAS,cexqeBS,langmuirSingle,freundlichSingle] = preCalculationMulti(parK,parIS,parIM,V,concA,concB);
        plotIS = [cexqeAS cexqeBS];
        plotFS = plotIS;
        disp('Step Complete!');
        
        disp('Step 03 - Performing adjustments...');
        disp('Kinetics Adjusts...');
        disp(' Pseudo First Order...');
        x0(1) = xlsread('adsorptionData.xlsx','resultadoCineticasMulti','N7');
        x0(2) = xlsread('adsorptionData.xlsx','resultadoCineticasMulti','O7');
        [qePFOA,k1A,r2PFOA] = pseudoFirstOrder(txqtA,x0);
        clear x0;
        x0(1) = xlsread('adsorptionData.xlsx','resultadoCineticasMulti','N8');
        x0(2) = xlsread('adsorptionData.xlsx','resultadoCineticasMulti','O8');
        [qePFOB,k1B,r2PFOB,PFO] = pseudoFirstOrder(txqtB,x0);
        clear x0;
        pseudoFirst = [qePFOA qePFOB; k1A k1B; r2PFOA r2PFOB];
        disp(' Pseudo Second Order...');
        x0(1) = xlsread('adsorptionData.xlsx','resultadoCineticasMulti','N9');
        x0(2) = xlsread('adsorptionData.xlsx','resultadoCineticasMulti','O9');
        [qePSOA,k2A,r2PSOA] = pseudoSecondOrder(txqtA,x0);
        clear x0;
        x0(1) = xlsread('adsorptionData.xlsx','resultadoCineticasMulti','N10');
        x0(2) = xlsread('adsorptionData.xlsx','resultadoCineticasMulti','O10');
        [qePSOB,k2B,r2PSOB,PSO] = pseudoSecondOrder(txqtB,x0);
        clear x0;
        pseudoSecond = [qePSOA qePSOB; k2A k2B; r2PSOA r2PSOB];
        disp('Isotherms Adjusts...');
        disp('  Langmuir Modified Model...');
        [cexqeAM,cexqeBM,erLM] = langmuirModified(langmuirSingle,cexqeA,cexqeB);
        disp('  Langmuir Partially Competitive Multicomponent Model...');
        [cexqeAPC,cexqeBPC,erLPC] = langmuirPartiallyCompetitive(langmuirSingle,cexqeA,cexqeB);
        disp('  Freundlich Multicomponent Model...');
        [cexqeAF,cexqeBF,erLF,a] = freundlichMulticomponent(freundlichSingle,cexqeA,cexqeB);
        disp('  IAST-Freundlich Model...');
        [cexqeAIF,cexqeBIF,erLIF] = iastFreunclich(freundlichSingle,cexqeA,cexqeB);
        disp('Step Complete!');
        
        disp('Step 04 - Choose best fit models...');
        [biggerK,biggerI,adjustK,adjustI] = bestAdjustM(pseudoFirst,pseudoSecond,erLM,erLPC,erLF,erLIF);
        disp('Step Complete!');
        
        disp('Step 05 - Graphics Plot...');
        [plotKA,plotKB,plotLM,plotLPC,plotFM,plotIF] = graphPlotMulti(txqtA,txqtB,cexqeAM,cexqeBM,cexqeAPC,cexqeBPC,cexqeAF,cexqeBF,cexqeAIF,cexqeBIF,...
    pseudoFirst,pseudoSecond,PFO,PSO,langmuirSingle,freundlichSingle,a);
        disp('Plot Conclude!');
        
        disp('Step 06 - Saving Data...');
        saveDataMulti(plotIS,plotFS,plotKA,plotKB,plotLM,plotLPC,plotFM,plotIF,adjustK,biggerK,...
    pseudoFirst,pseudoSecond,langmuirSingle,freundlichSingle,erLM,erLPC,erLF,erLIF,biggerI,adjustI);
        disp('Step Conclude!');

    end
end

time = toc;
disp(['Execution Time = ',num2str(time),' seconds.']);
diary off;

end

%% Saving Data Multi - v2.0

function saveDataMulti(plotIS,plotFS,plotKA,plotKB,plotLM,plotLPC,plotFM,plotIF,adjustK,biggerK,...
    pseudoFirst,pseudoSecond,langmuirSingle,freundlichSingle,erLM,erLPC,erLF,erLIF,biggerI,adjustI)

% Kinetics
xlswrite('adsorptionData.xlsx',pseudoFirst(1,1),'resultadoCineticasMulti','B4');
xlswrite('adsorptionData.xlsx',pseudoFirst(2,1),'resultadoCineticasMulti','B5');
xlswrite('adsorptionData.xlsx',pseudoFirst(3,1),'resultadoCineticasMulti','B6');
xlswrite('adsorptionData.xlsx',pseudoFirst(1,2),'resultadoCineticasMulti','H4');
xlswrite('adsorptionData.xlsx',pseudoFirst(2,2),'resultadoCineticasMulti','H5');
xlswrite('adsorptionData.xlsx',pseudoFirst(3,2),'resultadoCineticasMulti','H6');
xlswrite('adsorptionData.xlsx',pseudoSecond(1,1),'resultadoCineticasMulti','E4');
xlswrite('adsorptionData.xlsx',pseudoSecond(2,1),'resultadoCineticasMulti','E5');
xlswrite('adsorptionData.xlsx',pseudoSecond(3,1),'resultadoCineticasMulti','E6');
xlswrite('adsorptionData.xlsx',pseudoSecond(1,2),'resultadoCineticasMulti','K4');
xlswrite('adsorptionData.xlsx',pseudoSecond(2,2),'resultadoCineticasMulti','K5');
xlswrite('adsorptionData.xlsx',pseudoSecond(3,2),'resultadoCineticasMulti','K6');
xlswrite('adsorptionData.xlsx',adjustK(1),'resultadoCineticasMulti','M4');
xlswrite('adsorptionData.xlsx',adjustK(2),'resultadoCineticasMulti','P4');
xlswrite('adsorptionData.xlsx',biggerK(1),'resultadoCineticasMulti','N5');
xlswrite('adsorptionData.xlsx',biggerK(2),'resultadoCineticasMulti','Q5');
xlswrite('adsorptionData.xlsx',plotKA,'resultadoCineticasMulti','B8:C18');
xlswrite('adsorptionData.xlsx',plotKB,'resultadoCineticasMulti','H8:I18');
    
% Isotherms
xlswrite('adsorptionData.xlsx',langmuirSingle,'resultadoIsotermasMulti','B4:C6');
xlswrite('adsorptionData.xlsx',freundlichSingle,'resultadoIsotermasMulti','F4:G6');
xlswrite('adsorptionData.xlsx',plotIS,'resultadoIsotermasMulti','A10:D15');
xlswrite('adsorptionData.xlsx',plotFS,'resultadoIsotermasMulti','E10:H15');
xlswrite('adsorptionData.xlsx',erLM,'resultadoIsotermasMulti','J3:J4');
xlswrite('adsorptionData.xlsx',plotLM,'resultadoIsotermasMulti','I6:L11');
xlswrite('adsorptionData.xlsx',erLPC,'resultadoIsotermasMulti','N3:N4');
xlswrite('adsorptionData.xlsx',plotLPC,'resultadoIsotermasMulti','M6:P11');
xlswrite('adsorptionData.xlsx',erLF,'resultadoIsotermasMulti','R3:R4');
xlswrite('adsorptionData.xlsx',plotFM,'resultadoIsotermasMulti','Q6:T11');
xlswrite('adsorptionData.xlsx',erLIF,'resultadoIsotermasMulti','V3:V4');
xlswrite('adsorptionData.xlsx',plotIF,'resultadoIsotermasMulti','U6:X11');
xlswrite('adsorptionData.xlsx',biggerI(1),'resultadoIsotermasMulti','J15');
xlswrite('adsorptionData.xlsx',biggerI(2),'resultadoIsotermasMulti','M15');
xlswrite('adsorptionData.xlsx',adjustI(1),'resultadoIsotermasMulti','I14');
xlswrite('adsorptionData.xlsx',adjustI(2),'resultadoIsotermasMulti','L14');

end

%% Saving Data Single - v2.0

function saveDataSingle(plotK,plotI,biggerK,adjustK,biggerI,adjustI,pseudoFirst,pseudoSecond,langmuir,freundlich)

% Kinetics
qePFO = pseudoFirst(1);
qePSO = pseudoSecond(1);
k1 = pseudoFirst(2);
k2 = pseudoSecond(2);
r2PFO = pseudoFirst(3);
r2PSO = pseudoSecond(3);
xlswrite('adsorptionData.xlsx',qePFO,'resultadoCineticas','B3');
xlswrite('adsorptionData.xlsx',k1,'resultadoCineticas','B4');
xlswrite('adsorptionData.xlsx',r2PFO,'resultadoCineticas','B5');
xlswrite('adsorptionData.xlsx',qePSO,'resultadoCineticas','E3');
xlswrite('adsorptionData.xlsx',k2,'resultadoCineticas','E4');
xlswrite('adsorptionData.xlsx',r2PSO,'resultadoCineticas','E5');
xlswrite('adsorptionData.xlsx',biggerK,'resultadoCineticas','H3');
xlswrite('adsorptionData.xlsx',adjustK,'resultadoCineticas','G2');
xlswrite('adsorptionData.xlsx',plotK,'resultadoCineticas','B7:C17');

% Isotherms
qmax = langmuir(1);
kL = langmuir(2);
r2L = langmuir(3);
kF = freundlich(1);
n = freundlich(2);
r2F = freundlich(3);
xlswrite('adsorptionData.xlsx',qmax,'resultadoIsotermas','B4');
xlswrite('adsorptionData.xlsx',kL,'resultadoIsotermas','B5');
xlswrite('adsorptionData.xlsx',r2L,'resultadoIsotermas','B6');
xlswrite('adsorptionData.xlsx',n,'resultadoIsotermas','E4');
xlswrite('adsorptionData.xlsx',kF,'resultadoIsotermas','E5');
xlswrite('adsorptionData.xlsx',r2F,'resultadoIsotermas','E6');
xlswrite('adsorptionData.xlsx',plotI,'resultadoIsotermas','A10:B13');
xlswrite('adsorptionData.xlsx',adjustI,'resultadoIsotermas','G3');
xlswrite('adsorptionData.xlsx',biggerI,'resultadoIsotermas','H4');
    
end

%% Saving Data - v1.0

function saveData(p,plotK,plotI,plotIS,plotKA,plotKB,plotLM,plotLPC,plotFM,plotIF,biggerK,biggerI,...
    adjustK,adjustI,pseudoFirst,pseudoSecond,langmuir,freundlich,langmuirSingle,freundlichSingle,...
    r2LM,r2LPC,r2LF,r2LIF)

if p == 1
    % Kinetics
    qePFO = pseudoFirst(1);
    qePSO = pseudoSecond(1);
    k1 = pseudoFirst(2);
    k2 = pseudoSecond(2);
    r2PFO = pseudoFirst(3);
    r2PSO = pseudoSecond(3);
    xlswrite('adsorptionData.xlsx',qePFO,'resultadoCineticas','B3');
    xlswrite('adsorptionData.xlsx',k1,'resultadoCineticas','B4');
    xlswrite('adsorptionData.xlsx',r2PFO,'resultadoCineticas','B5');
    xlswrite('adsorptionData.xlsx',qePSO,'resultadoCineticas','E3');
    xlswrite('adsorptionData.xlsx',k2,'resultadoCineticas','E4');
    xlswrite('adsorptionData.xlsx',r2PSO,'resultadoCineticas','E5');
    xlswrite('adsorptionData.xlsx',biggerK,'resultadoCineticas','H3');
    xlswrite('adsorptionData.xlsx',adjustK,'resultadoCineticas','G2');
    xlswrite('adsorptionData.xlsx',plotK,'resultadoCineticas','B7:C17');
    
    % Isotherms
    qmax = langmuir(1);
    kL = langmuir(2);
    r2L = langmuir(3);
    kF = freundlich(1);
    n = freundlich(2);
    r2F = freundlich(3);
    xlswrite('adsorptionData.xlsx',qmax,'resultadoIsotermas','B4');
    xlswrite('adsorptionData.xlsx',kL,'resultadoIsotermas','B5');
    xlswrite('adsorptionData.xlsx',r2L,'resultadoIsotermas','B6');
    xlswrite('adsorptionData.xlsx',n,'resultadoIsotermas','E4');
    xlswrite('adsorptionData.xlsx',kF,'resultadoIsotermas','E5');
    xlswrite('adsorptionData.xlsx',r2F,'resultadoIsotermas','E6');
    xlswrite('adsorptionData.xlsx',plotI,'resultadoIsotermas','A10:B13');
    xlswrite('adsorptionData.xlsx',adjustI,'resultadoIsotermas','G3');
    xlswrite('adsorptionData.xlsx',biggerI,'resultadoIsotermas','H4');
    
else
    if p == 2
    % Kinetics
    xlswrite('adsorptionData.xlsx',pseudoFirst(1,1),'resultadoCineticasMulti','B4');
    xlswrite('adsortionData.xlsx',pseudoFirst(2,1),'resultadoCineticasMulti','B5');
    xlswrite('adsorptionData.xlsx',pseudoFirst(3,1),'resultadoCineticasMulti','B6');
    xlswrite('adsorptionData.xlsx',pseudoFirst(1,2),'resultadoCineticasMulti','H4');
    xlswrite('adsorptionData.xlsx',pseudoFirst(2,2),'resultadoCineticasMulti','H5');
    xlswrite('adsorptionData.xlsx',pseudoFirst(3,2),'resultadoCineticasMulti','H6');
    xlswrite('adsorptionData.xlsx',pseudoSecond(1,1),'resultadoCineticasMulti','E4');
    xlswrite('adsorptionData.xlsx',pseudoSecond(2,1),'resultadoCineticasMulti','E5');
    xlswrite('adsorptionData.xlsx',pseudoSecond(3,1),'resultadoCineticasMulti','E6');
    xlswrite('adsorptionData.xlsx',pseudoSecond(1,2),'resultadoCineticasMulti','K4');
    xlswrite('adsorptionData.xlsx',pseudoSecond(2,2),'resultadoCineticasMulti','K5');
    xlswrite('adsortionData.xlsx',pseudoSecond(3,2),'resultadoCineticasMulti','K6');
    xlswrite('adsorptionData.xlsx',adjustK(1),'resultadoCineticasMulti','M4');
    xlswrite('adsorptionData.xlsx',adjustK(2),'resultadoCineticasMulti','P4');
    xlswrite('adsorptionData.xlsx',biggerK(1),'resultadoCineticasMulti','N5');
    xlswrite('adsorptionData.xlsx',biggerK(2),'resultadoCineticasMulti','Q5');
    xlswrite('adsorptionData.xlsx',plotKA,'resultadoCineticasMulti','B8:C18');
    xlswrite('adsorptionData.xlsx',plotKB,'resultadoCineticasMulti','H8:I18');
    
    % Isotherms
    xlswrite('adsorptionData.xlsx',langmuirSingle,'resultadoIsotermasMulti','B4:C6');
    xlswrite('adsorptionData.xlsx',freundlichSingle,'resultadoIsotermasMulti','F4:G6');
    xlswrite('adsorptionData.xlsx',plotIS,'resultadoIsotermasMulti','A10:D13');
    xlswrite('adsorptionData.xlsx',r2LM,'resultadoIsotermasMulti','J3:J4');
    xlswrite('adsorptionData.xlsx',plotLM,'resultadoIsotermasMulti','I6:L9');
    xlswrite('adsorptionData.xlsx',r2LPC,'resultadoIsotermasMulti','N3:N4');
    xlswrite('adsorptionData.xlsx',plotLPC,'resultadoIsotermasMulti','M6:P9');
    xlswrite('adsorptionData.xlsx',r2LF,'resultadoIsotermasMulti','R3:R4');
    xlswrite('adsorptionData.xlsx',plotFM,'resultadoIsotermasMulti','Q6:T9');
    xlswrite('adsorptionData.xlsx',r2LIF,'resultadoIsotermasMulti','V3:V4');
    xlswrite('adsorptionData.xlsx',plotIF,'resultadoIsotermasMulti','V6:X9');
    xlswrite('adsorptionData.xlsx',biggerI(1),'resultadoIsotermasMulti','J12');
    xlswrite('adsorptionData.xlsx',biggerI(2),'resultadoIsotermasMulti','M12');
    xlswrite('adsorptionData.xlsx',adjustI(1),'resultadoIsotermasMulti','I11');
    xlswrite('adsorptionData.xlsx',adjustI(2),'resultadoIsotermasMulti','L11');
    end    
end
end

%% Graphical plot - Multi - v2.0

function [plotKA,plotKB,plotLM,plotLPC,plotFM,plotIF] = graphPlotMulti(txqtA,txqtB,cexqeAM,cexqeBM,cexqeAPC,cexqeBPC,cexqeAF,cexqeBF,cexqeAIF,cexqeBIF,...
    pseudoFirst,pseudoSecond,PFO,PSO,langmuirSingle,freundlichSingle,a)

pfoA = pseudoFirst(:,1);
pfoB = pseudoFirst(:,2);
psoA = pseudoSecond(:,1);
psoB = pseudoSecond(:,2);
plotKA = txqtA;
plotKB = txqtB;
plotLM = [cexqeAM cexqeBM];
plotLPC = [cexqeAPC cexqeBPC];
plotFM = [cexqeAF cexqeBF];
plotIF = [cexqeAIF cexqeBIF];
tA = txqtA(:,1);
qtA = txqtA(:,2);
tB = txqtB(:,1);
qtB = txqtB(:,2);
ceAM = cexqeAM(:,1);
qeAM = cexqeAM(:,2);
ceBM = cexqeBM(:,1);
qeBM = cexqeBM(:,2);
ceAPC = cexqeAPC(:,1);
qeAPC = cexqeAPC(:,2);
ceBPC = cexqeBPC(:,1);
qeBPC = cexqeBPC(:,2);
ceAF = cexqeAF(:,1);
qeAF = cexqeAF(:,2);
ceBF = cexqeBF(:,1);
qeBF = cexqeBF(:,2);
ceAIF = cexqeAIF(:,1);
qeAIF = cexqeAIF(:,2);
ceBIF = cexqeBIF(:,1);
qeBIF = cexqeBIF(:,2);

qmaxA = langmuirSingle(1,1);
qmaxB = langmuirSingle(1,2);
kLA = langmuirSingle(2,1);
kLB = langmuirSingle(2,2);

kFA = freundlichSingle(1,1);
kFB = freundlichSingle(1,2);
nA = freundlichSingle(2,1);
nB = freundlichSingle(2,2);

a12 = a(1,1);
a21 = a(1,2);

time = linspace(tA(1),tA(end));

figure();
plot(tA,qtA,'ko',time,PFO(pfoA,time),'b-');
title('Pseudo First Order - Component A');
legend('Data','Ajusted Data');
xlabel('Time (min)');
ylabel('Adsorption Capacity (mg g^-1)');
grid;

figure();
plot(tB,qtB,'ko',time,PFO(pfoB,time),'b-');
title('Pseudo First Order - Component B');
legend('Data','Ajusted Data');
xlabel('Time (min)');
ylabel('Adsorption Capacity (mg g^-1)');
grid;

figure();
plot(tA,qtA,'ko',time,PSO(psoA,time),'b-');
title('Pseudo Second Order - Component A');
legend('Data','Ajusted Data');
xlabel('Time (min)');
ylabel('Adsorption Capacity (mg g^-1)');
grid;

figure();
plot(tB,qtB,'ko',time,PSO(psoB,time),'b-');
title('Pseudo Second Order - Component B');
legend('Data','Ajusted Data');
xlabel('Time (min)');
ylabel('Adsorption Capacity (mg g^-1)');
grid;

% figure();
% scatter(tA,qtA);
% title('Kinetics - Component A');
% xlabel('Time (min)');
% ylabel('Adsorption Capacity (mg g^-1)');
% grid;
%         
% figure();
% scatter(tB,qtB);
% title('Kinetics - Component B');
% xlabel('Time (min)');
% ylabel('Adsorption Capacity (mg g^-1)');
% grid;
%

LMA = @(ceAM,ceBM) (qmaxA*kLA.*ceAM)./(1+kLA.*ceAM+kLB.*ceBM);
LMB = @(ceAM,ceBM) (qmaxB*kLB.*ceBM)./(1+kLA.*ceAM+kLB.*ceBM);

ceA = linspace(ceAM(1),ceAM(end));
ceB = linspace(ceBM(1),ceBM(end));

% qe = linspace(qeAM(1),qeAM(end));
figure();
plot(ceAM,qeAM,'ko',ceA,LMA(ceA,ceB),'b-');
% scatter(ceAM,qeAM);
title('Langmuir Modified Isotherm - Component A');
xlabel('Equilibrium Concentration (mg L^-1');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
legend('Data','Ajusted Data');
grid;

% ce = linspace(ceBM(1),ceBM(end));
       
figure();
plot(ceBM,qeBM,'ko',ceB,LMB(ceA,ceB),'b-');
% scatter(ceBM,qeBM);
title('Langmuir Modified Isotherm - Component B');
xlabel('Equilibrium Concentration (mg L^-1');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
legend('Data','Ajusted Data');
grid;

clear ceA ceB

LPCA = @(ceAPC,ceBPC) (((qmaxA-qmaxB)*kLA.*ceAPC)./(1+kLA.*ceAPC))+((qmaxA*kLA.*ceAPC)./(1+kLA.*ceAPC+kLB.*ceBPC));
LPCB = @(ceAPC,ceBPC) (qmaxB*kLB.*ceBPC)./(1+kLA.*ceAPC+kLB.*ceBPC);

ceA = linspace(ceAPC(1),ceAPC(end));
ceB = linspace(ceBPC(1),ceBPC(end));

figure();
plot(ceAPC,qeAPC,'ko',ceA,LPCA(ceA,ceB),'b-');
% scatter(ceAPC,qeAPC);
title('Langmuir Partially Competitive Multicomponent Model - Component A');
xlabel('Equilibrium Concentration (mg L^-1');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
legend('Data','Ajusted Data');
grid;

% ce = linspace(ceBPC(1),ceBPC(end));
        
figure();
plot(ceBPC,qeBPC,'ko',ceB,LPCB(ceA,ceB),'b-');
% scatter(ceBPC,qeBPC);
title('Langmuir Partially Competitive Multicomponent Model - Component B');
xlabel('Equilibrium Concentration (mg L^-1');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
legend('Data','Ajusted Data');
grid;

ceA = linspace(ceAF(1),ceAF(end));
ceB = linspace(ceBF(1),ceBF(end));

mA = ((1/nA)-1);
mB = ((1/nB)-1);

AFM = @(ceAF,ceBF) kFA.*ceAF.*(ceAF+a12.*ceBF).^mA;
BFM = @(ceAF,ceBF) kFB.*ceBF.*(a21.*ceAF+ceBF).^mB;

figure();
plot(ceAF,qeAF,'ko',ceA,AFM(ceA,ceB),'b-');
% scatter(ceAF,qeAF);
title('Freundlich Multicomponent Isotherm - Component A');
xlabel('Equilibrium Concentration (mg L^-1');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
legend('Data','Ajusted Data');
grid;

% ce = linspace(ceBF(1),ceBF(end));
        
figure();
plot(ceBF,qeBF,'ko',ceB,BFM(ceA,ceB),'b-');
% scatter(ceBF,qeBF);
title('Freunlich Multicomponent Isotherm - Component B');
xlabel('Equilibrium Concentration (mg L^-1)');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
legend('Data','Ajusted Data');
grid;

qeIFA = linspace(qeAIF(1),qeAIF(end));
qeIFB = linspace(qeBIF(1),qeBIF(end));

AIF = @(qeIFA,qeIFB) (qeIFA./(qeIFA+qeIFB)).*((nA.*qeIFA+nB.*qeIFB)/(nA*kFA)).^nA;
BIF = @(qeIFA,qeIFB) (qeIFB./(qeIFA+qeIFB)).*((nA.*qeIFA+nB.*qeIFB)/(nB*kFB)).^nB;

% ceIFA = AIF(qeIFA,qeIFB);
% ceIFB = BIF(qeIFA,qeIFB);

figure();
plot(qeAIF,ceAIF,'ko',qeIFA,AIF(qeIFA,qeIFB),'b-');
% scatter(ceAIF,qeAIF);
title('IAST-Freundlich Isotherm - Component A');
xlabel('Equilibrium Concentration (mg L^-1)');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
legend('Data','Ajusted Data');
grid;
        
figure();
plot(qeBIF,ceBIF,'ko',qeIFB,BIF(qeIFA,qeIFB),'b-');
% scatter(ceBIF,qeBIF);
title('IAST-Freundlich Isotherm - Component B');
xlabel('Equilibrium Concentration (mg L^(-1))');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
legend('Data','Ajusted Data');
grid;

end

%% Graphical plot - Single - v2.0

function [plotK,plotI] = graphPlotSingle(txqt,cexqe)

plotK = txqt;
plotI = cexqe;
t = txqt(:,1);
qt = txqt(:,2);
ce = cexqe(:,1);
qe = cexqe(:,2);
 
figure();
scatter(t,qt);
title('Kinetics');
xlabel('Time (min)');
ylabel('Adsorption Capacity (mg g^-1)');
grid;
 
figure();
scatter(ce,qe);
title('Isotherms');
xlabel('Equilibrium Concentration (mg L^-1');
ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
grid;

end

%% Graphical plot - v1.0

function [plotK,plotI,plotKA,plotKB,plotLM,plotLPC,plotFM,plotIF] = graphPlot(p,txqt,cexqe,txqtA,txqtB,cexqeAM,cexqeBM,cexqeAPC,cexqeBPC,cexqeAF,cexqeBF,cexqeAIF,cexqeBIF)

if p == 1
    plotK = txqt;
    plotI = cexqe;
    t = txqt(:,1);
    qt = txqt(:,2);
    ce = cexqe(:,1);
    qe = cexqe(:,2);
    
    figure();
    scatter(t,qt);
    title('Kinetics');
    xlabel('Time (min)');
    ylabel('Adsorption Capacity (mg g^-1)');
    grid;
    
    figure();
    scatter(ce,qe);
    title('Isotherms');
    xlabel('Equilibrium Concentration (mg L^-1');
    ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
    grid;
    
else
    if p == 2
        plotKA = txqtA;
        plotKB = txqtB;
        plotLM = [cexqeAM cexqeBM];
        plotLPC = [cexqeAPC cexqeBPC];
        plotFM = [cexqeAF cexqeBF];
        plotIF = [cexqeAIF cexqeBIF];
        tA = txqtA(:,1);
        qtA = txqtA(:,2);
        tB = txqtB(:,1);
        qtB = txqtB(:,2);
        ceAM = cexqeAM(:,1);
        qeAM = cexqeAM(:,2);
        ceBM = cexqeBM(:,1);
        qeBM = cexqeBM(:,2);
        ceAPC = cexqeAPC(:,1);
        qeAPC = cexqeAPC(:,2);
        ceBPC = cexqeBPC(:,1);
        qeBPC = cexqeBPC(:,2);
        ceAF = cexqeAF(:,1);
        qeAF = cexqeAF(:,2);
        ceBF = cexqeBF(:,1);
        qeBF = cexqeBF(:,2);
        ceAIF = cexqeAIF(:,1);
        qeAIF = cexqeAIF(:,2);
        ceBIF = cexqeBIF(:,1);
        qeBIF = cexqeBIF(:,2);
        
        figure();
        scatter(tA,qtA);
        title('Kinetics - Component A');
        xlabel('Time (min)');
        ylabel('Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(tB,qtB);
        title('Kinetics - Component B');
        xlabel('Time (min)');
        ylabel('Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(ceAM,qeAM);
        title('Langmuir Modified Isotherm - Component A');
        xlabel('Equilibrium Concentration (mg L^-1');
        ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(ceBM,qeBM);
        title('Langmuir Modified Isotherm - Component B');
        xlabel('Equilibrium Concentration (mg L^-1');
        ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(ceAPC,qeAPC);
        title('Langmuir Partially Competitive Multicomponent Model - Component A');
        xlabel('Equilibrium Concentration (mg L^-1');
        ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(ceBPC,qeBPC);
        title('Langmuir Partially Competitive Multicomponent Model - Component B');
        xlabel('Equilibrium Concentration (mg L^-1');
        ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(ceAF,qeAF);
        title('Freundlich Multicomponent Isotherm - Component A');
        xlabel('Equilibrium Concentration (mg L^-1');
        ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(ceBF,qeBF);
        title('Freunlich Multicomponent Isotherm - Component B');
        xlabel('Equilibrium Concentration (mg L^-1');
        ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(ceAIF,qeAIF);
        title('IAST-Freundlich Isotherm - Component A');
        xlabel('Equilibrium Concentration (mg L^-1');
        ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
        grid;
        
        figure();
        scatter(ceBIF,qeBIF);
        title('IAST-Freundlich Isotherm - Component B');
        xlabel('Equilibrium Concentration (mg L^-1');
        ylabel('Equilibrium Adsorption Capacity (mg g^-1)');
        grid;
    end
end    

end

%% Choose the best adjust (multi);

function [biggerK,biggerI,adjustK,adjustI] = bestAdjustM(pseudoFirst,pseudoSecond,erLM,erLPC,erLF,erLIF)

% pseudoFirst = [qePFOA qePFOB; k1A k1B; r2PFOA r2PFOB];
% pseudoSecond = [qePSOA qePSOB; k2A k2B; r2PSOA r2PSOB];

biggerIA = 0;
biggerIB = 0;
biggerKA = 0;
biggerKB = 0;
adjustKA = 0;
adjustKB = 0;
adjustIA = 0;
adjustIB = 0;

[n,m] = size(pseudoFirst);
[o,p] = size(pseudoSecond);

if n == o && m == p
    for i = 1:n
        if i == 1
            qePFOA = pseudoFirst(i,i);
            qePFOB = pseudoFirst(i,i+1);
            qePSOA = pseudoSecond(i,i);
            qePSOB = pseudoSecond(i,i+1);
        else
            if i == 2
                k1A = pseudoFirst(i,i-1);
                k1B = pseudoFirst(i,i);
                k2A = pseudoSecond(i,i-1);
                k2B = pseudoSecond(i,i);
            else
                if i == 3
                    r2PFOA = pseudoFirst(i,i-2);
                    r2PFOB = pseudoFirst(i,i-1);
                    r2PSOA = pseudoSecond(i,i-2);
                    r2PSOB = pseudoSecond(i,i-1);
                end
            end
        end
    end
else
    error('Dimensions dont match!');
end
clear n m o p;

n = length(erLM);
m = length(erLPC);
o = length(erLF);
p = length(erLIF);

if n == m && m == o && o == p
    for i = 1:n
        if i == 1
            erLMA = erLM(i);
            erLPCA = erLPC(i);
            erLFA = erLF(i);
            erLIFA = erLIF(i);
        else
            if i == 2
                erLMB = erLM(i);
                erLPCB = erLPC(i);
                erLFB = erLF(i);
                erLIFB = erLIF(i);
            end
        end
    end
else
    error('Dimensions dont match!');
end

% Component A
if r2PFOA > r2PSOA 
    biggerKA = r2PFOA;
    adjustKA = {'Pseudo First Order'};
else
    if r2PSOA > r2PFOA
        biggerKA = r2PSOA;
        adjustKA = {'Pseudo Second Order'};
    end
end

% Component B
if r2PFOB > r2PSOB
    biggerKB = r2PFOB;
    adjustKB = {'Pseudo First Order'};
else
    if r2PSOB > r2PFOB
        biggerKB = r2PSOB;
        adjustKB = {'Pseudo Second Order'};
    end
end

if erLMA < erLPCA && erLMA < erLFA && erLMA < erLIFA
    biggerIA = erLMA;
    adjustIA = {'Langmuir Modified Model'};
else
    if erLPCA < erLMA && erLPCA < erLFA && erLPCA < erLIFA
        biggerIA = erLPCA;
        adjustIA = {'Langmuir Partially Competitive Multicomponent Model'};
    else
        if erLFA < erLMA && erLFA < erLPCA && erLFA < erLIFA
            biggerIA = erLFA;
            adjustIA = {'Freundlich Multicomponent Model'};
        else
            if erLIFA < erLMA && erLIFA < erLPCA && erLIFA < erLFA
                biggerIA = erLIFA;
                adjustIA = {'IAST-Freundlich Model'};
            end
        end
    end
end

if erLMB < erLPCB && erLMB < erLFB && erLMB < erLIFB
    biggerIB = erLMB;
    adjustIB = {'Langmuir Modified Model'};
else
    if erLPCB < erLMB && erLPCB < erLFB && erLPCB < erLIFB
        biggerIB = erLPCB;
        adjustIB = {'Langmuir Partially Competitive Multicomponent Model'};
    else
        if erLFB < erLMB && erLFB < erLPCB && erLFB < erLIFB
            biggerIB = erLFB;
            adjustIB = {'Freundlich Multicomponent Model'};
        else
            if erLIFB < erLMB && erLIFB < erLPCB && erLIFB < erLFB
                biggerIB = erLIFB;
                adjustIB = {'IAST-Freundlich Model'};
            end
        end
    end
end

biggerK = [biggerKA; biggerKB];
biggerI = [biggerIA; biggerIB];
adjustK = [adjustKA; adjustKB];
adjustI = [adjustIA; adjustIB];

n = length(biggerK);
for i = 1:n
    if i == 1
        if biggerK(i) == r2PFOA
            disp('Pseudo First Order Model');
            disp('fits best to the studied data');
            disp(['r^2 = ',num2str(r2PFOA)]);
            disp(['qe = ',num2str(qePFOA)]);
            disp(['k1 = ',num2str(k1A)]);
        else
            if biggerK(i) == r2PSOA
            disp('Pseudo Second Order Model');
            disp('fits best to the studied data');
            disp(['r^2 = ',num2str(r2PSOA)]);
            disp(['k = ',num2str(k2A)]);
            disp(['qe = ',num2str(qePSOA)]);
            end
        end
    else
        if i == 2
            if biggerK(i) == r2PFOB
                disp('Pseudo First Order Model');
                disp('fits best to the studied data');
                disp(['r^2 = ',num2str(r2PFOB)]);
                disp(['qe = ',num2str(qePFOB)]);
                disp(['k1 = ',num2str(k1B)]);
        else
            if biggerK(i) == r2PSOB
                disp('Pseudo Second Order Model');
                disp('fits best to the studied data');
                disp(['r^2 = ',num2str(r2PSOB)]);
                disp(['k = ',num2str(k2B)]);
                disp(['qe = ',num2str(qePSOB)]);
            end
            end
        end
    end
end

clear n;
n = length(biggerI);
for i = 1:n
    if i == 1
        if biggerI(i) == erLMA
            disp('Langmuir Modified Model');
            disp('fits better to the studied data');
            disp(['Average error (%) = ',num2str(erLMA)]);
        else
            if biggerI(i) == erLPCA
                disp('Langmuir Partially Competitive Multicomponent Model');
                disp('fits better to the studied data');
                disp(['Average error (%) = ',num2str(erLPCA)]);
            else
                if biggerI(i) == erLFA
                    disp('Freundlich Multicomponent Model');
                    disp('fits better to the studied data');
                    disp(['Average error (%) = ',num2str(erLFA)]);
                else
                    if biggerI(i) == erLIFA
                        disp('IAST-Freundlich Model');
                        disp('fits better to the studied data');
                        disp(['Average error (%) = ',num2str(erLIFA)]);
                    end
                end
            end
        end
    else
        if i == 2
            if biggerI(i) == erLMB
            disp('Langmuir Modified Model');
            disp('fits better to the studied data');
            disp(['Average error (%) = ',num2str(erLMB)]);
        else
            if biggerI(i) == erLPCB
                disp('Langmuir Partially Competitive Multicomponent Model');
                disp('fits better to the studied data');
                disp(['Average error (%) = ',num2str(erLPCB)]);
            else
                if biggerI(i) == erLFB
                    disp('Freundlich Multicomponent Model');
                    disp('fits better to the studied data');
                    disp(['Average error (%) = ',num2str(erLFB)]);
                else
                    if biggerI(i) == erLIFB
                        disp('IAST-Freundlich Model');
                        disp('fits better to the studied data');
                        disp(['Average error (%) = ',num2str(erLIFB)]);
                    end
                end
            end
            end
        end
    end
end

end

%% Choose the best adjust (single)

function [biggerK,biggerI,adjustK,adjustI]=bestAdjustS(pseudoFirst,pseudoSecond,langmuir,freundlich)

% pseudoFirst = [qePFO; k1; r2PFO];
% pseudoSecond = [qePSO; k2; r2PSO];
% langmuir = [qmax; kL; r2L];
% freundlich = [kF; n; r2F];

n = length(pseudoFirst);
m = length(pseudoSecond);
o = length(langmuir);
p = length(freundlich);
if n == m && m == o && o == p
    for i = 1:n
        if i == 1
            qePFO = pseudoFirst(i);
            qePSO = pseudoSecond(i);
            qmax = langmuir(i);
            kF = freundlich(i);
        else
            if i == 2
                k1 = pseudoFirst(i);
                k2 = pseudoSecond(i);
                kL = langmuir(i);
                n = langmuir(i);
            else
                if i == 3
                   r2PFO = pseudoFirst(i);
                   r2PSO = pseudoSecond(i);
                   r2L = langmuir(i);
                    r2F = freundlich(i);
                end
            end
        end
    end
else
    error('Dimensions dont match!');
end

biggerK = 0;
biggerI = 0;
adjustK = 0;
adjustI = 0;

if r2PFO > r2PSO 
   biggerK = r2PFO;
   adjustK = {'Pseudo First Order'};
else
   if r2PSO > r2PFO
       biggerK = r2PSO;
       adjustK = {'Pseudo Second Order'};
   end
end

if r2L > r2F
   biggerI = r2L;
   adjustI = {'Langmuir Model'};
else
   if r2F > r2L
       biggerI = r2F;
        adjustI = {'Freundlich Model'};
    end
end

if biggerK == r2PFO
    disp('Pseudo First Order Model');
    disp('fits best to the studied data');
    disp(['r^2 = ',num2str(r2PFO)]);
    disp(['qe = ',num2str(qePFO)]);
    disp(['k1 = ',num2str(k1)]);
else
    if biggerK == r2PSO
        disp('Pseudo Second Order Model');
        disp('fits best to the studied data');
        disp(['r^2 = ',num2str(r2PSO)]);
        disp(['k = ',num2str(k2)]);
        disp(['qe = ',num2str(qePSO)]);
    end
end

if biggerI == r2L
    disp('Langmuir Isotherm Model');
    disp('fits best to the studied data');
    disp(['r^2 = ',num2str(r2L)]);
    disp(['qmax = ',num2str(qmax)]);
    disp(['kL = ',num2str(kL)]);
else
    if biggerI == r2F
        disp('PFreundlich Isotherm Model');
        disp('fits best to the studied data');
        disp(['r^2 = ',num2str(r2F)]);
        disp(['kF= ',num2str(kF)]);
        disp(['n = ',num2str(n)]);
    end
end
end

%% Pseudo Second Order Model

function [qePSO,k2,r2PSO,PSO] = pseudoSecondOrder(txqt,x0)

t = txqt(:,1);
qt = txqt(:,2);

PSO = @(x,t) ((x(1)^2)*x(2).*t)./((x(1)*x(2).*t)+1);
% x0(1) = input('Initial Values: ');
% x0(2) = input('Initial Values: ');
% % x0 = [1,1];
par = lsqcurvefit(PSO,x0,t,qt);

qePSO = par(1);
k2 = par(2);

r2PSO = determinationCoefficient(t,qt,PSO,par);

end

%% Pseudo First Order Model

function [qePFO,k1,r2PFO,PFO] = pseudoFirstOrder(txqt,x0)

t = txqt(:,1);
qt = txqt(:,2);

PFO = @(x,t) x(1).*(1-exp(-x(2).*t));
% x0(1) = input('Initial Values: ');
% x0(2) = input('Initial Values: ');
% % x0 = [1,1];
par = lsqcurvefit(PFO,x0,t,qt);

qePFO = par(1);
k1 = par(2);

r2PFO = determinationCoefficient(t,qt,PFO,par);

end

%% IAST-Freundlch Model

function [cexqeAIF,cexqeBIF,erLIF] = iastFreunclich(freundlichSingle,cexqeA,cexqeB)

kFA = freundlichSingle(1,1);
kFB = freundlichSingle(1,2);
nA = freundlichSingle(2,1);
nB = freundlichSingle(2,2);

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

n = length(qeA);
m = length(qeB);

ceAIF = zeros(n,1);
ceBIF = zeros(m,1);

if n == m
    for i = 1:n
        ceAIF(i) = (qeA(i)/(qeA(i)+qeB(i)))*((nA*qeA(i)+nB*qeB(i))/(nA*kFA))^(nA);
        ceBIF(i) = (qeB(i)/(qeA(i)+qeB(i)))*((nA*qeA(i)+nB*qeB(i))/(nB*kFB))^(nB);
    end
else
    error('Dimension dont match!');
end

erA = avgerror(ceAIF,ceA);
erB = avgerror(ceBIF,ceB);
% [~,~,r2LA] = linRegression(ceAIF,qeA);
% [~,~,r2LB] = linRegression(ceBIF,qeB);
% [r2LA] = determinationCoefficient(ceAIF,qeA);
% [r2LB] = determinationCoefficient(ceBIF,qeB);

cexqeAIF = [ceAIF qeA];
cexqeBIF = [ceBIF qeB];
erLIF = [erA; erB];

end

%% Freundlich Multicomponent Model

function [cexqeAF,cexqeBF,erLF,a] = freundlichMulticomponent(freundlichSingle,cexqeA,cexqeB)

kFA = freundlichSingle(1,1);
kFB = freundlichSingle(1,2);
nA = freundlichSingle(2,1);
nB = freundlichSingle(2,2);

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

% [a,~,~] = linRegression(ceA,ceB);
% 
% a12 = a;
% a21 = (a12)^(-1);
% 
% clear a;

FA = @(a12,ceA) kFA.*ceA.*(ceA+a12.*ceB).^((1/nA)-1);
FB = @(a21,ceB) kFB.*ceB.*(a21.*ceA+ceB).^((1/nB)-1);

a12 = lsqcurvefit(FA,1,ceA,qeA);
a21 = lsqcurvefit(FB,1,ceB,qeB);

a = [a12 a21];

mA = (1/nA)-1;
mB = (1/nB)-1;

n = length(ceA);
m = length(ceB);

qeAF = zeros(n,1);
qeBF = zeros(m,1);

if n == m
    for i = 1:n
        qeAF(i) = kFA*ceA(i)*(ceA(i)+a12*ceB(i))^mA;
        qeBF(i) = kFB*ceB(i)*(a21*ceA(i)+ceB(i))^mB;
    end
else
    error('Dimensions dont match!');
end

erA = avgerror(qeAF,qeA);
erB = avgerror(qeBF,qeB);
% [~,~,r2LA] = linRegression(ceA,qeAF);
% [~,~,r2LB] = linRegression(ceB,qeBF);
% [r2LA] = determinationCoefficient(ceA,qeAF);
% [r2LB] = determinationCoefficient(ceB,qeBF);

cexqeAF = [ceA qeAF];
cexqeBF = [ceB qeBF];
erLF = [erA; erB];

end

%% Langmuir Partially Competitive Multicomponent Model

function [cexqeAPC,cexqeBPC,erLPC] = langmuirPartiallyCompetitive(langmuirSingle,cexqeA,cexqeB)

qmaxA = langmuirSingle(1,1);
qmaxB = langmuirSingle(1,2);
kLA = langmuirSingle(2,1);
kLB = langmuirSingle(2,2);

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

n = length(ceA);
m = length(ceB);

qeAPC = zeros(n,1);
qeBPC = zeros(m,1);

if m == n
    for i = 1:n
        qeAPC(i) = (((qmaxA-qmaxB)*kLA*ceA(i))/(1+kLA*ceA(i)))+((qmaxA*kLA*ceA(i))/(1+kLA*ceA(i)+kLB*ceB(i)));
        qeBPC(i) = (qmaxB*kLB*ceB(i))/(1+kLA*ceA(i)+kLB*ceB(i));
    end
else
    error('Dimensions dont match!');
end

erA = avgerror(qeAPC,qeA);
erB = avgerror(qeBPC,qeB);
% [~,~,r2LA] = linRegression(ceA,qeAPC);
% [~,~,r2LB] = linRegression(ceB,qeBPC);
% [r2LA] = determinationCoefficient(ceA,qeAPC);
% [r2LB] = determinationCoefficient(ceB,qeBPC);

cexqeAPC = [ceA qeAPC];
cexqeBPC = [ceB qeBPC];
erLPC = [erA; erB];

end

%% Modified Langmuir Model

function [cexqeAM,cexqeBM,erLM] = langmuirModified(langmuirSingle,cexqeA,cexqeB)

qmaxA = langmuirSingle(1,1);
qmaxB = langmuirSingle(1,2);
kLA = langmuirSingle(2,1);
kLB = langmuirSingle(2,2);

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

n = length(ceA);
m = length(ceB);

qeAM = zeros(n,1);
qeBM = zeros(m,1);

if n == m
    for i = 1:n
        qeAM(i) = (qmaxA*kLA*ceA(i))/(1+kLA*ceA(i)+kLB*ceB(i));
        qeBM(i) = (qmaxB*kLB*ceB(i))/(1+kLA*ceA(i)+kLB*ceB(i));
    end
else
    error('Dimensions dont match!');
end

erA = avgerror(qeAM,qeA);
erB = avgerror(qeBM,qeB);

% [~,~,r2LA] = linRegression(ceA,qeAM);
% [~,~,r2LB] = linRegression(ceB,qeBM);
% [r2LA] = determinationCoefficient(ceA,qeAM);
% [r2LB] = determinationCoefficient(ceB,qeBM);

cexqeAM = [ceA qeAM];
cexqeBM = [ceB qeBM];
erLM = [erA; erB];

end


%% Single Freundlich Model

function [kF,n,r2F]=freundlichModel(cexqe,x0)

[~,j]=size(cexqe);
if j~=2
    error('Just two columns are supported!');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

F = @(x,ce) (x(1).*ce.^(1/x(2)));
% x0(1) = input('Initial Values: ');
% x0(2) = input('Initial Values: ');
% x0 = [1 1];
par = lsqcurvefit(F,x0,ce,qe);

kF = par(1);
n = par(2);

r2F = determinationCoefficient(ce,qe,F,par);

end

%% Single Langmuir Model

function [qmax,kL,r2L]=langmuirModel(cexqe,x0)

[~,j]=size(cexqe);
if j~=2
    error('Just two columns are supported');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

L = @(x,ce) (x(1)*x(2).*ce)./(1+x(2).*ce);
% x0(1) = input('Initial Values: ');
% x0(2) = input('Initial Values: ');
% x0 = [1 1];
par = lsqcurvefit(L,x0,ce,qe);

qmax = par(1);
kL = par(2);

r2L = determinationCoefficient(ce,qe,L,par);

end

%% Average error

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

%% Linear Regression

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

%% Determination Coefficient

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

%% Preliminary Calculation - Multi v2.0

function [txqtA,txqtB,cexqeA,cexqeB,cexqeAS,cexqeBS,langmuirSingle,freundlichSingle] = preCalculationMulti(parK,parIS,parIM,V,concA,concB)

% to 'parK':
% j = 1 and 4 - time (min)
% j = 2 - Ce to component A (mg/L)
% j = 3 - mass to component A (mg)
% j = 5 - Ce to component B (mg/L)
% j = 6 - mass to component B (mg)

% to 'parIM':
% n = 1 and 3 - C0 (mg/L)
% n = 2 - Ce to component A (mg/L)
% n = 4 - Ce to component B (mg/L)
% n = 5 - mass (mg)

% to 'parISA':
% l = 1 - C0 (mg/L)
% l = 2 - Ce to component A (mg/L)
% l = 3 - mass to component A (mg)

% to 'parISB':
% l = 1 - C0 (mg/L)
% l = 2 - Ce to component B (mg/L)
% l = 3 - mass to component B (mg)

[i,~] = size(parK);
[m,~] = size(parIM);
[k,~] = size(parIS);
% [l,~] = size(parISB);

t = zeros(i,1);
cA = zeros(i,1);
cB = zeros(i,1);
massA = zeros(i,1);
massB = zeros(i,1);

for i = 1:i
    t(i) = parK(i,1);
    cA(i) = parK(i,2);
    cB(i) = parK(i,5);
    massA(i) = parK(i,3);
    massB(i) = parK(i,6);
end

ciA = zeros(k,1);
ceA = zeros(k,1);
mA = zeros(k,1);
ciB = zeros(k,1);
ceB = zeros(k,1);
mB = zeros(k,1);

for k = 1:k
    ciA(k) = parIS(k,1);
    ceA(k) = parIS(k,2);
    mA(k) = parIS(k,3);
    ciB(k) = parIS(k,4);
    ceB(k) = parIS(k,5);
    mB(k) = parIS(k,6);
end

ciAM = zeros(m,1);
ceAM = zeros(m,1);
ciBM = zeros(m,1);
ceBM = zeros(m,1);
mass = zeros(m,1);

for m = 1:m
    ciAM(m) = parIM(m,1);
    ceAM(m) = parIM(m,2);
    ciBM(m) = parIM(m,3);
    ceBM(m) = parIM(m,4);
    mass(m) = parIM(m,5);
end

% Kinetics

n = length(t);
qtA = zeros(n,1);
qtB = zeros(n,1);
for n = 1:n
    qtA(n) = ((concA-cA(n))*V)/(massA(n));
    qtB(n) = ((concB-cB(n))*V)/(massB(n));
end
txqtA = [t qtA];
txqtB = [t qtB];
% disp(txqtA);
% disp(txqtB);

% Isotherms

n = length(ceAM);
m = length(ceBM);
o = length(ceA);
u = length(ceB);
if n == m
    if o == u
        qeA = zeros(n,1);
        qeB = zeros(n,1);
        qeAS = zeros(o,1);
        qeBS = zeros(o,1);
        for n = 1:n
            % Multi-component
            qeA(n) = ((ciAM(n)-ceAM(n))*V)/(mass(n));
            qeB(n) = ((ciBM(n)-ceBM(n))*V)/(mass(n));
        end
        for n = 1:o
            % Single component
            qeAS(n) = ((ciA(n)-ceA(n))*V)/(mA(n));
            qeBS(n) = ((ciB(n)-ceB(n))*V)/(mB(n));
        end
        cexqeA = [ceAM qeA];
        cexqeB = [ceBM qeB];
        cexqeAS = [ceA qeAS];
        cexqeBS = [ceB qeBS];
%         disp(cexqeAS);
%         disp(cexqeBS);
%         disp(cexqeA);
%         disp(cexqeB);
    end
else
    error('Dimensions not equal');
end

disp('Langmuir Model - Component A');
x0(1) = xlsread('adsorptionData.xlsx','resultadoIsotermasMulti','C17');
x0(2) = xlsread('adsorptionData.xlsx','resultadoIsotermasMulti','D17');
[qmaxA,kLA,r2LA] = langmuirModel(cexqeAS,x0);
clear x0;
disp('Langmuir Model - Component B');
x0(1) = xlsread('adsorptionData.xlsx','resultadoIsotermasMulti','C18');
x0(2) = xlsread('adsorptionData.xlsx','resultadoIsotermasMulti','D18');
[qmaxB,kLB,r2LB] = langmuirModel(cexqeBS,x0);
clear x0;
disp('Freundlich Model - Component A');
x0(1) = xlsread('adsorptionData.xlsx','resultadoIsotermasMulti','C19');
x0(2) = xlsread('adsorptionData.xlsx','resultadoIsotermasMulti','D19');
[kFA,nA,r2FA] = freundlichModel(cexqeAS,x0);
clear x0;
disp('Freundlich Model - Component B');
x0(1) = xlsread('adsorptionData.xlsx','resultadoIsotermasMulti','C20');
x0(2) = xlsread('adsorptionData.xlsx','resultadoIsotermasMulti','D20');
[kFB,nB,r2FB] = freundlichModel(cexqeBS,x0);
clear x0;
langmuirSingle = [qmaxA qmaxB; kLA kLB; r2LA r2LB];
freundlichSingle = [kFA kFB; nA nB; r2FA r2FB];

end

%% Preliminary Calculation - Single v2.0

function [txqt,cexqe] = preCalculationSingle(parK,parIS,V,concA)

% to 'parK':
% j = 1 and 4 - time (min)
% j = 2 - Ce to compenent A (mg/L)
% j = 3 - mass to component A (mg)
% j = 5 - Ce to component B (mg/L)
% j = 6 - mass to component B (mg)

% to 'parIS':
% l = 1 and 4 - C0 (mg/L)
% l = 2 - Ce to component A (mg/L)
% l = 3 - mass to component A (mg)
% l = 5 - Ce to component B (mg/L)
% l = 6 - mass to component B (mg)

[i,~] = size(parK);
[k,~] = size(parIS);
    
for i = 1:i
    t(i) = parK(i,1);
    cA(i) = parK(i,2);
    massA(i) = parK(i,3);
end

for k = 1:k
    ci(k) = parIS(k,1);
    ceA(k) = parIS(k,2);
    mA(k) = parIS(k,3);
end

% Kinetics

n = length(t);
for n = 1:n
    qt(n) = ((concA-cA(n))*V)/(massA(n));
end
txqt = [t' qt'];

% Isotherms

n = length(mA);
m = length(ceA);
if n == m 
    for n = 1:n
        qe(n) = ((ci(n)-ceA(n))*V)/mA(n);
    end
    cexqe = [ceA' qe'];
else
    error('Dimensions not equal');
end
end


%% Preliminary Calculation - v1.0

function [txqt,txqtA,txqtB,cexqe,cexqeA,cexqeB,langmuirSingle,freundlichSingle] = preCalculation(parK,parIS,parIM,V,p,concA,concB)

% to 'parK':
% j = 1 and 4 - time (min)
% j = 2 - Ce to compenent A (mg/L)
% j = 3 - mass to component A (mg)
% j = 5 - Ce to component B (mg/L)
% j = 6 - mass to component B (mg)

% to 'parIM':
% n = 1 and 3 - C0 (mg/L)
% n = 2 - Ce to component A (mg/L)
% n = 4 - Ce to component B (mg/L)
% n = 5 - mass (mg)

% to 'parIS':
% l = 1 and 4 - C0 (mg/L)
% l = 2 - Ce to component A (mg/L)
% l = 3 - mass to component A (mg)
% l = 5 - Ce to component B (mg/L)
% l = 6 - mass to component B (mg)

disp('parK = ');
disp(parK);
disp('parIM = ');
disp(parIM);
disp('parIS = ');
disp(parIS);
disp(['concA = ',num2str(concA)]);
disp(['concB = ',num2str(concB)]);
disp(['V = ',num2str(V)]);

[i,~] = size(parK);
[m,~] = size(parIM);
[k,~] = size(parIS);
    
% disp('Random Debug Message');
% disp(['p = ',num2str(p)]);

if p == 1
    clear parIM;
    clear concB;
end

for i = 1:i
    t(i) = parK(i,1);
    cA(i) = parK(i,2);
    cB(i) = parK(i,5);
    massA(i) = parK(i,3);
    massB(i) = parK(i,6);
%     if p == 1
%         disp('Random Debug Message');
%         clear cB massB
%     end
end

if p == 1
    clear cB;
    clear massB;
end

% disp(t);
% disp(cA);
% disp(cB);
% disp(massA);
% disp(massB);


for k = 1:k
    ci(k) = parIS(k,1);
    ceA(k) = parIS(k,2);
    mA(k) = parIS(k,3);
    ceB(k) = parIS(k,5);
    mB(k) = parIS(k,6);
%     if p == 'S'
%         clear ceB mB
%     end
end

if p == 2
    for m = 1:m
        ceAM(m) = parIM(m,2);
        ceBM(m) = parIM(m,4);
        mass(m) = parIM(m,5);
    end
end

% V = V*(10^-3);

% Kinetics

if p == 1
    n = length(t);
%     disp(n);
    for n = 1:n
        qt(n) = ((concA-cA(n))*V)/(massA(n));
    end
%     disp(t);
%     disp(qt);
    txqt = [t' qt'];
%     disp(txqt);
else
    if p == 2
        n = length(t);
        for n = 1:n
            qtA(n) = ((concA-cA(n))*V)/(massA(n));
            qtB(n) = ((concB-cB(n))*V)/(massB(n));
        end
        txqtA = [t' qtA'];
        txqtB = [t' qtB'];
    end
end

% Isotherms

if p == 1
    clear ceB mB
    n = length(mA);
    m = length(ceA);
    if n == m 
        for n = 1:n
            qe(n) = ((ci(n)-ceA(n))*V)/mA(n);
        end
        cexqe = [ceA' qe'];
    else
        error('Dimensions not equal');
    end
else
    if p == 2
        n = length(ceAM);
        m = length(ceBM);
        o = length(ceA);
        u = length(ceB);
        if n == m && m == o && o == u
            for n = 1:n
                qeA(n) = ((ci(n)-ceAM(n))*V)/(mass(n));
                qeB(n) = ((ci(n)-ceBM(n))*V)/(mass(n));
                qeAS(n) = ((ci(n)-ceA(n))*V)/(mA(n));
                qeBS(n) = ((ci(n)-ceB(n))*V)/(mB(n));
            end
            cexqeA = [ceAM' qeA'];
            cexqeB = [ceBM' qeB'];
            cexqeAS = [ceA' qeAS'];
            cexqeBS = [ceB' qeBS'];
        else
            error('Dimensions not equal');
        end
        disp('Langmuir Model - Component A');
        [qmaxA,kLA,r2LA] = langmuirModel(cexqeAS);
        disp('Langmuir Model - Component B');
        [qmaxB,kLB,r2LB] = langmuirModel(cexqeBS);
        disp('Freundlich Model - Component A');
        [kFA,nA,r2FA] = freundlichModel(cexqeAS);
        disp('Freundlich Model - Component B');
        [kFB,nB,r2FB] = freundlichModel(cexqeBS);
        langmuirSingle = [qmaxA qmaxB; kLA kLB; r2LA r2LB];
        freundlichSingle = [kFA kFB; nA nB; r2FA r2FB];
    end
end
    
end

%% Loading Parameters

function [parK,parIM,parIS,V,concA,concB] = loadParameters()

parK = xlsread('adsorptionData.xlsx','dadosCinetica','R6:W10');
parIM = xlsread('adsorptionData.xlsx','dadosIsoterma','N33:R38');
parIS = xlsread('adsorptionData.xlsx','dadosIsoterma','O6:T11');
V = xlsread('adsorptionData.xlsx','dadosIsoterma','D3');
concA = xlsread('adsorptionData.xlsx','dadosCinetica','S2');
concB = xlsread('adsorptionData.xlsx','dadosCinetica','S3');
% disp(parIS);
% disp(parIM);
% disp(parK);

end