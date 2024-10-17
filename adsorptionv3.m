% Kinetics and Isotherms
% Raphael Gilioli Heineck

%% DEVELOPMENT %%
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
%            - Algorithm has to plot the adjusted data as well.
% 08/03/2021 - Graph plot was changed. best adjust no longer available
%            - Best adjust available, new equation
% 30/03/2021 - Average error function implemented, to use in multicomponent
%              models.

% Version 3.0
% 14/03/2021 - New main function implemented.
% 19/03/2021 - Linear Models available to.

%% MAIN FUNCTION %%

function adsorptionv3()

diary adsorption.txt;

close;
close all;
clear;
clearvars -global;
clc;
tic;

global options PFO PSO L F

disp('Kinetics and Isotherms - v3.0');
disp('Author: Raphael Gilioli Heineck, Chemical Engineer');
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','MaxIterations',1000,'FunctionTolerance',1e-99);

disp('Step 01 - Loading Parameters...');
[parK,parIM,parIS,V,concA,concB] = loadParameters();
disp('Loading Complete!');

disp('Options: ');
disp(' 1 - Nonlinear Single Component Kinetics');
disp(' 2 - Nonlinear Multi Component Kinetics');
disp(' 3 - Linear Single Component Kinetics ');
disp(' 4 - Linear Multi Component Kinetics ');
disp(' 5 - Nonlinear Single Component Isotherms');
disp(' 6 - Multi Component Isotherms');
disp(' 7 - Linear Single Component Isotherms');
disp(' 8 - Exit');
p = input('What do you wanna do? ');

while p == 1
    disp('Nonlinear Kinetics Single Component choosen');
    
    disp('Step 01 - Preliminary Calculations...');
    [txqt,~,~] = preCalculationKinetics(parK,V,concA,concB,p);
    disp('Step Complete');
    
    disp('Step 02 - Performing Adjusts...');
    [qePFO,k1,r2PFO,PFO] = pseudoFirstOrder(txqt);
    [qePSO,k2,r2PSO,PSO] = pseudoSecondOrder(txqt);
    pseudoFirst = [qePFO; k1; r2PFO];
    pseudoSecond = [qePSO; k2; r2PSO];
    disp('Step Complete');
    p = input('New choice: ');
end

while p == 2
    disp('Nonlinear Multi Component Kinetics');
    
    disp('Step 01 - Preliminary Calculations...');
    [~,txqtA,txqtB] = preCalculationKinetics(parK,V,concA,concB,p);
    disp('Step Complete');
    
    disp('Step 02 - Performing Adjusts...');
    [qePFOA,k1A,r2PFOA,~] = pseudoFirstOrder(txqtA);
    [qePFOB,k1B,r2PFOB,PFO] = pseudoFirstOrder(txqtB);
    [qePSOA, k2A,r2PSOA,~] = pseudoSecondOrder(txqtA);
    [qePSOB,k2B,r2PSOB,PSO] = pseudoSecondOrder(txqtB);
    pseudoFirst = [qePFOA qePFOB; k1A k1B; r2PFOA r2PFOB];
    pseudoSecond = [qePSOA qePSOB; k2A k2B; r2PSOA r2PSOB];
    disp('Step Complete');
    p = input('New choice: ');
end

while p == 3
    disp('Linear Single Component Kinetics');
    
    disp('Step 01 - Preliminary Calculations...');
    [txqt,~,~] = preCalculationKinetics(parK,V,concA,concB,p);
    disp('Step Complete');
    
    disp('Step 02 - Perfoming Adjusts...');
    [qePFOL,k1,r2PFOL] = pseudoFirstOrderLinear(txqt);
    [qePSOL,k2,r2PSOL] = pseudoSecondOrderLinear(txqt);
    pseudoFirstLinear = [qePFOL; k1; r2PFOL];
    pseudoSecondLinear = [qePSOL; k2; r2PSOL];
    disp('Step Complete');
    p = input('New Choice: ');
end

while p == 4
    disp('Linear Multi Component Kinetics');
    
    disp('Step 01 - Preliminary Calculations...');
    [~,txqtA,txqtB] = preCalculationKinetics(parK,V,concA,concB,p);
    disp('Step Complete');
    
    disp('Step 02 - Performing Adjusts...');
    [qePFOAL,k1A,r2PFOAL] = pseudoFirstOrderLinear(txqtA);
    [qePFOBL,k1B,r2PFOBL] = pseudoFirstOrderLinear(txqtB);
    [qePSOAL,k2A,r2PSOAL] = pseudoSecondOrderLinear(txqtA);
    [qePSOBL,k2B,r2PSOBL] = pseudoSecondOrderLinear(txqtB);
    pseudoFirstLinear = [qePFOAL qePFOBL; k1A k1B; r2PFOAL r2PFOBL];
    pseudoSecondLinear = [qePSOAL qePSOBL; k2A k2B; r2PSOAL r2PSOBL];
    disp('Step Complete');
    p = input('New Choice: ');
end

while p == 5
    disp('Nonlinear Single Component Isotherm');
    
    disp('Step 01 - Preliminary Calculations...');
    [cexqe,~,~,~,~,~,~] = preCalculationIsotherms(parIS,parIM,V,p);
    disp('Step Complete');
    
    disp('Step 02 - Performing Adjusts...');
    [qmax,kL,r2L,L] = lagmuirModel(cexqe);
    [kF,n,r2F,F] = freundlichModel(cexqe);
    langmuirSingle = [qmax; kL; r2L];
    freundlichSingle = [kF; n; r2F];
    disp('Step Complete');
    p = input('New choice: ');
end

while p == 6
    disp('Multi Component Isotherms');
    
    disp('Step 01 - Preliminary Calculations');
    [~,cexqeA,cexqeB,cexqeAS,cexqeBS,langmuirSingle,freundlichSingle] = preCalculationIsotherms(parIS,parIM,V,p);
    disp('Step Complete');
    
    disp('Step 02 - Performing Adjusts...');
    [cexqeAM,cexqeBM,erLM] = langmuirModified(langmuirSingle,cexqeA,cexqeB);
    [cexqeAPC,cexqeBPC,erLPC] = langmuirPartiallyCompetitive(langmuirSingle,cexqeA,cexqeB);
    [cexqeAF,cexqeBF,erLF] = freundlichMulticomponent(freundlichSingle,cexqeA,cexqeB);
    [cexqeAIF,cexqeBIF,erLIF] = iastFreunclich(freundlichSingle,cexqeA,cexqeB);
    langmuirM = [cexqeAM; cexqeBM];
    langmuirPC = [cexqeAPC; cexqeBPC];
    freundlichM = [cexqeAF; cexqeBF];
    IAST = [cexqeAIF; cexqeBIF];
    er = [erLM; erLPC; erLF; erLIF];
    disp('Step Complete');
    p = input('New choice: ');
end

while p == 7
    disp('Linear Single Component Isotherms');
    
    disp('Step 01 - Preliminary Calculation...');
    [cexqe,~,~,~,~,~,~] = preCalculationIsotherms(parIS,parIM,V,p);
    disp('Step Complete');
    
    disp('Step 02 - Performing Adjusts...');
    [qmax,kL,r2L] = lagmuirModelLinear(cexqe);
    [kF,n,r2F] = freundlichModelLinear(cexqe);
    langmuirSingleLinear = [qmax; kL; r2L];
    freundlichSingleLinear = [kF; n; r2F];
    disp('Step Complete');
    p = input('New choice: ');
end

if p == 8
    time = toc;
    disp(['Execution time = ',num2str(time),' seconds.']);
    diary off;
end
end

%% ISOTHERM MODELS %%

% IAST-Freundlch Model

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

cexqeAIF = [ceAIF' qeA];
cexqeBIF = [ceBIF' qeB];
erLIF = [erA; erB];

end

% Freundlich Multicomponent Model

function [cexqeAF,cexqeBF,erLF] = freundlichMulticomponent(freundlichSingle,cexqeA,cexqeB)

kFA = freundlichSingle(1,1);
kFB = freundlichSingle(1,2);
nA = freundlichSingle(2,1);
nB = freundlichSingle(2,2);

ceA = cexqeA(:,1);
qeA = cexqeA(:,2);
ceB = cexqeB(:,1);
qeB = cexqeB(:,2);

[a,~,~] = linRegression(ceA,ceB);

a12 = -(a);
a21 = 1/a12;

n = length(ceA);
m = length(ceB);
if n == m
    for i = 1:n
        qeAF(i) = kFA*ceA(i)*(ceA(i)+a12*ceB(i))^((1/nA)-1);
        qeBF(i) = kFB*ceB(i)*(ceA(i)*a21+ceB(i))^((1/nB)-1);
    end
else
    error('Dimensions dont match!');
end

erA = avgerror(qeAF,qeA);
erB = avgerror(qeBF,qeB);

cexqeAF = [ceA qeAF'];
cexqeBF = [ceB qeBF'];
erLF = [erA; erB];

end

% Langmuir Partially Competitive Multicomponent Model

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
if m == n
    for i = 1:n
        qeAPC(i) = (((qmaxA-qmaxB)*kLA*ceA(i))/1+kLA*ceA(i))+((qmaxA*kLA*ceA(i))/(1+kLA*ceA(i)+kLB*ceB(i)));
        qeBPC(i) = (qmaxB*kLB*ceB(i))/(1+kLA*ceA(i)+kLB*ceB(i));
    end
else
    error('Dimensions dont match!');
end

erA = avgerror(qeAPC,qeA);
erB = avgerror(qeBPC,qeB);

cexqeAPC = [ceA qeAPC'];
cexqeBPC = [ceB qeBPC'];
erLPC = [erA; erB];

end

% Modified Langmuir Model

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

cexqeAM = [ceA qeAM'];
cexqeBM = [ceB qeBM'];
erLM = [erA; erB];

end

% Single Freundlich Model Linear

function [kF,n,r2F] = freundlichModelLinear(cexqe)

[~,j] = size(cexqe);
if j ~= 2
    error('Just two columns are supported!');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

% lnqe = log (qe) -> y
% lnce = log (ce) -> x
% ln (kF) -> a
% 1/n -> b

n = length(ce);
m = length(qe);
for n = 1:n
    lnqe(n) = log(qe(n));
    lnce(n) = log(ce(n));
end

if m == n
    [a,b,r2F] = linRegression(lnqe,lnce);
end

kF = exp(b);
n = 1/a;

end

% Single Freundlich Model Nonlinear

function [kF,n,r2F,F]=freundlichModel(cexqe)

[~,j]=size(cexqe);
if j~=2
    error('Just two columns are supported!');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

F = @(x,ce) (x(1).*ce.^(1/x(2)));
x0 = [1 1];
par = lsqcurvefit(F,x0,ce,qe);

kF = par(1);
n = par(2);

r2F = determinationCoefficient(ce,qe,F,par);

end

% Single Langmuir Model Linear

function [qmax,kL,r2L] = lagmuirModelLinear(cexqe)

[~,j] = size(cexqe);
if j ~= 2
    error('Number of columns must be less than 2!');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

% ceqe = ce/qe -> y
% ce -> x
% 1/qmax -> a
% 1/kL qmax -> b

n = length(ce);
m = length(qe);
for n = 1:n
   ceqe(n) = ce(n)/qe(n);
end

if m == n
    [a,b,r2L] = linRegression(ce,ceqe);
end

qmax = 1/a;
kL = 1/(b*qmax);

end

% Single Langmuir Model Nonlinear

function [qmax,kL,r2L,L]=langmuirModel(cexqe)

[~,j]=size(cexqe);
if j~=2
    error('Just two columns are supported');
end

ce = cexqe(:,1);
qe = cexqe(:,2);

L = @(x,ce) (x(1)*x(2).*ce)./(1+x(2).*ce);
x0 = [1 1];
par = lsqcurvefit(L,x0,ce,qe);

qmax = par(1);
kL = par(2);

r2L = determinationCoefficient(ce,qe,L,par);

end

%% KINETIC MODELS %%

% Pseudo Second Order Model Nonlinear

function [qePSO,k2,r2PSO,PSO] = pseudoSecondOrder(txqt)

[~,j] = size(txqt);
if j ~= 2
    error('Number of columns must be less than 2!');
end

t = txqt(:,1);
qt = txqt(:,2);

PSO = @(x,t) ((x(1)^2)*x(2).*t)./((x(1)*x(2).*t)+1);
x0 = [1,1];
par = lsqcurvefit(PSO,x0,t,qt);

qePSO = par(1);
k2 = par(2);

r2PSO = determinationCoefficient(t,qt,PSO,par);

end

% Pseudo Second Order Model Linear

function [qePSOL,k2,r2PSOL] = pseudoSecondOrderLinear(txqt)

[~,j] = size(txqt);
if j ~= 2
    error('Number of columns must be less than 2!');
end

t = txqt(:,1);
qt = txqt(:,2);

% tqt = t/qt -> y
% t -> x
% 1/k2*qe^2 -> a
% 1/qe -> b

n = length(qt);
m = length(t);
for i = 1:n
    tqt = t(n)/qt(n);
end

if m == n
    [a,b,r2PSOL] = linRegression(t,tqt);
end

qePSOL = 1/a;
k2 = 1/b*(qePSOL^2);

end

% Pseudo First Order Model Nonlinear

function [qePFO,k1,r2PFO,PFO] = pseudoFirstOrder(txqt)

[~,j] = size(txqt);
if j ~= 2
    error('Number of columns must be less than 2!');
end

t = txqt(:,1);
qt = txqt(:,2);

PFO = @(x,t) x(1).*(1-exp(-x(2).*t));
x0 = [1,1];
par = lsqcurvefit(PFO,x0,t,qt);

qePFO = par(1);
k1 = par(2);

r2PFO = determinationCoefficient(t,qt,PFO,par);

end

% Pseudo First Order Model Linear

function [qePFOL,k1,r2PFOL] = pseudoFirstOrderLinear(txqt)

[~,j] = size(txqt);
if j ~= 2
    error('Number of columns must be less than 2!');
end

t = txqt(:,1);
qt = txqt(:,2);

n = length(qt);
qe = qt(n);
Qtn = extract(qt);
T = extract(t);
clear qt t
qt = Qtn;
t = T;

% lnqt = log(qe - qt) -> y
% lnqe = log(qe) -> b
% k1 -> a
% t -> x

n = length(qt);
m = length(t);
for n = 1:n
    lnqt(n) = log(qe - qt(n));
end

if m == n
    [a,b,r2PFOL] = linRegression(t,lnqt);
else
    error('Dimensions dont match');
end

k1 = -a;
qePFOL = exp(b);

end

%% AUXILIARY FUNCTIONS %%

% Extraction of vectors

function X = extract(x)

n = length(x);
for i = 1:n-1
    X(i) = x(i);
end

end

% Average error

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

% Linear Regression

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

% Determination of determination coefficient

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

%% PRELIMINARY CALCULATIONS %%

% Preliminary Calculation Isotherms

function [cexqe,cexqeA,cexqeB,cexqeAS,cexqeBS,langmuirSingle,freundlichSingle] = preCalculationIsotherms(parIS,parIM,V,p)

% to 'parIS':
% l = 1 and 4 - C0 (mg/L)
% l = 2 - Ce to component A (mg/L)
% l = 3 - mass to component A (mg)
% l = 5 - Ce to component B (mg/L)
% l = 6 - mass to component B (mg)

% to 'parIM':
% n = 1 and 3 - C0 (mg/L)
% n = 2 - Ce to component A (mg/L)
% n = 4 - Ce to component B (mg/L)
% n = 5 - mass (mg)

if p == 5 || 7
    ci = parIS(:,1);
    ceA = parIS(:,2);
    massA = parIS(:,3);
    
    n = length(massA);
    m = length(ceA);
    if n == m
        for n = 1:n
            qe(n) = ((ci(n)-ceA(n))*V)/massA(n);
        end
        cexqe = [ceA' qe'];
    else
        error('Dimension dont match');
    end
else
    if p == 6
        ciA = parIS(:,1);
        ceA = parIS(:,2);
        massA = parIS(:,3);
        ciB = parIS(:,4);
        ceB = parIS(:,5);
        massB = parIS(:,6);
        
        ciAM = parIM(:,1);
        ceAM = parIM(:,2);
        ciBM = parIM(:,3);
        ceBM = parIM(:,4);
        mass = parIM(:,5);
        
        n = length(ceAM);
        m = length(ceBM);
        o = length(ceA);
        u = length(ceB);
        
        if n == m && m == o && o == u
            for n = 1:n
                % Multi-component
                qeA(n) = ((ciAM(n)-ceAM(n))*V)/(mass(n));
                qeB(n) = ((ciBM(n)-ceBM(n))*V)/(mass(n));
                % Single component
                qeAS(n) = ((ciA(n)-ceA(n))*V)/(mA(n));
                qeBS(n) = ((ciB(n)-ceB(n))*V)/(mB(n));
            end
            cexqeA = [ceAM' qeA'];
            cexqeB = [ceBM' qeB'];
            cexqeAS = [ceA' qeAS'];
            cexqeBS = [ceB' qeBS'];
        else
            error('Dimension dont match');
        end
        
        [qmaxA,kLA,r2LA] = langmuirModel(cexqeAS);
        [qmaxB,kLB,r2LB] = langmuirModel(cexqeBS);
        [kFA,nA,r2FA] = freundlichModel(cexqeAS);
        [kFB,nB,r2FB] = freundlichModel(cexqeBS);
        langmuirSingle = [qmaxA qmaxB; kLA kLB; r2LA r2LB];
        freundlichSingle = [kFA kFB; nA nB; r2FA r2FB];
    end
end
end

% Preliminary Calculation Kinetics

function [txqt,txqtA,txqtB] = preCalculationKinetics(parK,V,concA,concB,p)

% to 'parK':
% j = 1 and 4 - time (min)
% j = 2 - Ce to compenent A (mg/L)
% j = 3 - mass to component A (mg)
% j = 5 - Ce to component B (mg/L)
% j = 6 - mass to component B (mg)


if p == 1 || 3
    t = parK(:,1);
    cA = parK(:,2);
    massA = parK(:,3);
    n = length(t);
    for n = 1:n
        qt(n) = ((concA-cA(n))*V)/massA(n);
    end
    txqt = [t' qt'];           

else
    if p == 2 || 4
        t = parK(:,1);
        cA = parK(:,2);
        cB = parK(:,5);
        massA = parK(:,3);
        massB = parK(:,6);
        n = length(t);
        for n =1:n
            qtA(n) = ((concA-cA(n))*V)/(massA(n));
            qtB(n) = ((concB-cB(n))*V)/(massB(n));
        end
        txqtA = [t' qtA'];
        txqtB = [t' qtB'];
    end
end
end

%% LOADING PARAMETERS %%

function [parK,parIM,parIS,V,concA,concB] = loadParameters()

parK = xlsread('adsorptionData.xlsx','dadosCinetica','R6:W11');
parIM = xlsread('adsorptionData.xlsx','dadosIsoterma','N20:R23');
parIS = xlsread('adsorptionData.xlsx','dadosIsoterma','O6:T9');
V = xlsread('adsorptionData.xlsx','dadosIsoterma','D3');
concA = xlsread('adsorptionData.xlsx','dadosCinetica','S2');
concB = xlsread('adsorptionData.xlsx','dadosCinetica','S3');
% disp(parIS);

end