% in this file we load simulation results and compute the stability
% condition

% IN THIS FILE LOOK ONLY AT HOUSEHOLDS LEARNING STABILITY

numSimulations = 100;

load rPi_Int;
load rR_Int;

% solve New Keynesian model steady state

ParameterValues = defineParameters(); 

[ SteadyStateValuesNK ] = solveNK_SteadyState( ParameterValues );

% calculate parameter values denoted in the learning tex file

[ ParameterValuesLearning ] = calculateLearningParams( SteadyStateValuesNK, ParameterValues );

numOfVars = 7;

syms eta Lst wst cst Rst alpha Xst rK theta beta Cw Hc Cx Hr Rd rPi rR kst Gk bkk akk akA rho bwk awk awA bmk amk amA bik aik aiA bink aink ainA 

% FIRMS
syms bkk_F akk_F akA_F bink_F aink_F ainA_F bmk_F amk_F amA_F

B  = sym('B',[3 3]);
B(:,:) = 0;
B(1,1) = 1;
B(2,1) = bkk;
B(2,2) = akk;
B(2,3) = akA;
B(3,3) = rho;

Ib = sym('B',[3 3]);
Ib(:,:) = 0;
Ib(1,1) = 1;
Ib(2,2) = 1;
Ib(3,3) = 1;

A = Gk * beta * (B')*((Ib-beta*B')^(-1));

% compute V
V      = sym('V',[numOfVars numOfVars]);
V(:,:) = 0;
% compute e
e      = sym('e',[numOfVars 1]);
e(:,1) = 0;
e(6,1) = 1;
% compute xi matrix
xiMatrix = sym('xiMatrix',[1 numOfVars]);
xiMatrix(1,:) = 0;
xiMatrix(1,1) = (1-alpha)*(kst^alpha)*(Lst^(-alpha));
xiMatrix(1,3) = -1;

V(1,1) = (eta-1)/Lst;
V(1,2) = -1/wst;
V(1,3) = 1/cst;
V(2,5) = -1;
V(2,6) = 1;
V(2,7) = -Rst;
V(3,1) = (1-alpha)/Lst;
V(3,4) = -1/Xst;
V(3,5) = -1/rK;
V(4,1) = -alpha/Lst;
V(4,2) = -1/wst;
V(4,4) = -1/Xst;
V(5,4) = (1-theta*beta)/Xst;
V(5,7) = theta/(1-theta);
V(6,2) = -Cw;
V(6,3) = Hc;
V(6,4) = -Cx;
V(6,6) = -Hr;
V(6,7) = Rd;
V(7,6) = 1;
V(7,7) = -Rst*(1-rR)*rPi;

% HOUSHOLDS
stabilityCondition = struct('capital',zeros(1,numSimulations),'wage',zeros(1,numSimulations),'interest',zeros(1,numSimulations),'inflation',zeros(1,numSimulations),'markup',zeros(1,numSimulations));

% FIRMS
stabilityCondition_F = struct('capital',zeros(1,numSimulations),'inflation',zeros(1,numSimulations),'markup',zeros(1,numSimulations));

for j = 1:numSimulations
    
    name      = sprintf('results_NK%d.mat',j); % set the output name
	inputFile = load(name);                    % load the output
    
    clear inputFile.M_
    
    [ REmatrix_A,~,~ ] = transition_matrix(inputFile.oo_);
    
    clear inputFile.oo_
    
    % households
    B_Capital     = define_B_matrix(bkk,akk,akA,rho);
    B_Wage        = define_B_matrix(bwk,awk,awA,rho);
    B_Markup      = define_B_matrix(bmk,amk,amA,rho);
    B_Interest    = define_B_matrix(bik,aik,aiA,rho);
    B_Inflation   = define_B_matrix(bink,aink,ainA,rho);
    
    % firms
    B_Capital_F   = define_B_matrix(bkk_F,akk_F,akA_F,rho);
    B_Inflation_F = define_B_matrix(bink_F,aink_F,ainA_F,rho);
    B_Markup_F    = define_B_matrix(bmk_F,amk_F,amA_F,rho);
    
    
    % define V matrix
    
    V      = sym('V',[numOfVars numOfVars]);
    V(:,:) = 0;
    % compute e
    e      = sym('e',[numOfVars 1]);
    e(:,1) = 0;
    e(6,1) = 1;
    % compute xi matrix
    xiMatrix = sym('xiMatrix',[1 numOfVars]);
    xiMatrix(1,:) = 0;
    xiMatrix(1,1) = (1-alpha)*(kst^alpha)*(Lst^(-alpha));
    xiMatrix(1,3) = -1;

    V(1,1) = (eta-1)/Lst;
    V(1,2) = -1/wst;
    V(1,3) = 1/cst;
    V(2,5) = -1;
    V(2,6) = 1;
    V(2,7) = -Rst;
    V(3,1) = (1-alpha)/Lst;
    V(3,4) = -1/Xst;
    V(3,5) = -1/rK;
    V(4,1) = -alpha/Lst;
    V(4,2) = -1/wst;
    V(4,4) = -1/Xst;
    V(5,4) = (1-theta*beta)/Xst;
    V(5,7) = theta/(1-theta);
    V(6,2) = -Cw;
    V(6,3) = Hc;
    V(6,4) = -Cx;
    V(6,6) = -Hr;
    V(6,7) = Rd;
    V(7,6) = 1;
    V(7,7) = -Rst*(1-rR)*rPi;
    
    % assign numerical values to symbols
    eta   = ParameterValues.eta;
    Lst   = SteadyStateValuesNK.L;
    wst   = SteadyStateValuesNK.w;
    cst   = SteadyStateValuesNK.c;
    Rst   = SteadyStateValuesNK.R;
    alpha = ParameterValues.alpha;
    Xst   = SteadyStateValuesNK.X;
    rK    = SteadyStateValuesNK.rk;
    theta = ParameterValues.theta;
    beta  = ParameterValues.beta;
    Cw    = ParameterValuesLearning.Cw;
    Hc    = (-ParameterValuesLearning.Cc)*(1+beta/((1-beta)^2));
    Cx    = ParameterValuesLearning.Cx;
    Hr    = (1/((1-beta)^2)*cst*(beta^2)*ParameterValuesLearning.Cc) - (beta/(1-beta))*(wst*Lst+(1-1/Xst)*(cst+1/beta)-cst)/SteadyStateValuesNK.R;
    Rd    = SteadyStateValuesNK.R*SteadyStateValuesNK.k;
    rPi   = rPi_Int(j); % we can choose this vale 
    rR    = rR_Int(j); % we can choose this vale
    kst   = SteadyStateValuesNK.k;
    rho   = ParameterValues.rho;
    Gk    = ParameterValuesLearning.Ck/(1-beta);
    
    % HOUSEHOLDS
    % capital estimates
    bkk   = SteadyStateValuesNK.k*(1-REmatrix_A(1,1));
    akk   = REmatrix_A(1,1);
    akA   = SteadyStateValuesNK.k*REmatrix_A(2,1);
    
    % wage estimates
    bwk = SteadyStateValuesNK.w*(1-REmatrix_A(1,2));
    awk = REmatrix_A(1,2)*SteadyStateValuesNK.w/SteadyStateValuesNK.k;
    awA = SteadyStateValuesNK.w*REmatrix_A(2,2);
    
    % interest rate estimates
    bik = SteadyStateValuesNK.R*(1-REmatrix_A(1,3));
    aik = REmatrix_A(1,3)*SteadyStateValuesNK.R/SteadyStateValuesNK.k;
    aiA = SteadyStateValuesNK.R*REmatrix_A(2,3);
    
    % inflation estimates
    bink = 1-REmatrix_A(1,4);
    aink = REmatrix_A(1,4)/SteadyStateValuesNK.k;
    ainA = REmatrix_A(2,4);
    
    % markup estimates
    bmk = SteadyStateValuesNK.X*(1-REmatrix_A(1,5));
    amk = REmatrix_A(1,5)*SteadyStateValuesNK.X/SteadyStateValuesNK.k;
    amA = SteadyStateValuesNK.X*REmatrix_A(2,5);
    
    % FIRMS
    % capital estimates
    bkk_F = SteadyStateValuesNK.k*(1-REmatrix_A(1,1));
    akk_F = REmatrix_A(1,1);
    akA_F = SteadyStateValuesNK.k*REmatrix_A(2,1);     
    
    % markup estimates
    bmk_F = SteadyStateValuesNK.X*(1-REmatrix_A(1,5));
    amk_F = REmatrix_A(1,5)*SteadyStateValuesNK.X/SteadyStateValuesNK.k;
    amA_F = SteadyStateValuesNK.X*REmatrix_A(2,5);
        
    % inflation estimates
    bink_F = 1-REmatrix_A(1,4);
    aink_F = REmatrix_A(1,4)/SteadyStateValuesNK.k;
    ainA_F = REmatrix_A(2,4);
    
    % define identity matrix
    Ib = sym('B',[3 3]);
    Ib(:,:) = 0;
    Ib(1,1) = 1;
    Ib(2,2) = 1;
    Ib(3,3) = 1;
    
    % HOUSEHOLDS
    A_Capital     = Gk * beta * (B_Capital')*((Ib-beta*B_Capital')^(-1));
    A_Inflation   = Gk * beta * (B_Inflation')*((Ib-beta*B_Inflation')^(-1));
    A_Interest    = Gk * beta * (B_Interest')*((Ib-beta*B_Interest')^(-1));
    A_Markup      = Gk * beta * (B_Markup')*((Ib-beta*B_Markup')^(-1));
    A_Wage        = Gk * beta * (B_Wage')*((Ib-beta*B_Wage')^(-1));
    
    % FIRMS
    A_Capital_F   = Gk * beta * (B_Capital_F')*((Ib-beta*B_Capital_F')^(-1));
    A_Markup_F    = Gk * beta * (B_Markup_F')*((Ib-beta*B_Markup_F')^(-1));
    A_Inflation_F = Gk * beta * (B_Inflation_F')*((Ib-beta*B_Inflation_F')^(-1));
    
    % evaluate matrices
    Ib          = eval(Ib);
    B_Capital   = eval(B_Capital);
    B_Inflation = eval(B_Inflation); 
    B_Wage      = eval(B_Wage);
    B_Interest  = eval(B_Interest);
    B_Markup    = eval(B_Markup);
    
    V           = eval(V);
    xiMatrix    = eval(xiMatrix);
    
    % HOUSEHOLDS
    A_Capital   = eval(A_Capital);
    A_Inflation = eval(A_Inflation);
    A_Interest  = eval(A_Interest);
    A_Markup    = eval(A_Markup);
    A_Wage      = eval(A_Wage);
    
    % FIRMS
    A_Capital_F   = eval(A_Capital_F);
    A_Markup_F    = eval(A_Markup_F);
    A_Inflation_F = eval(A_Inflation_F);
    
    % finally, compute H matrix
    
    % HOUSEHOLDS
    H_Capital   = (e'*(inv(V'))*(xiMatrix')*A_Capital - Ib);
    H_Inflation = (e'*(inv(V'))*(xiMatrix')*A_Inflation - Ib);
    H_Interest  = (e'*(inv(V'))*(xiMatrix')*A_Interest - Ib);
    H_Markup    = (e'*(inv(V'))*(xiMatrix')*A_Markup - Ib);
    H_Wage      = (e'*(inv(V'))*(xiMatrix')*A_Wage - Ib);
    
    % FIRMS
    F_Capital   = (e'*(inv(V'))*(xiMatrix')*A_Capital_F - Ib);
    F_Markup    = (e'*(inv(V'))*(xiMatrix')*A_Markup_F - Ib);
    F_Inflation = (e'*(inv(V'))*(xiMatrix')*A_Inflation_F - Ib);
    
    % HOUSEHOLDS
    stabilityCondition.capital(1,j)   = checkStability( eig(H_Capital) );
    stabilityCondition.inflation(1,j) = checkStability( eig(H_Inflation) );
    stabilityCondition.interest(1,j)  = checkStability( eig(H_Interest) );
    stabilityCondition.markup(1,j)    = checkStability( eig(H_Markup) );
    stabilityCondition.wage(1,j)      = checkStability( eig(H_Wage) );
    
    % FIRMS
    stabilityCondition_F.capital(1,j)   = checkStability( eig(F_Capital) );
    stabilityCondition_F.markup(1,j)    = checkStability( eig(F_Markup) );
    stabilityCondition_F.inflation(1,j) = checkStability( eig(F_Inflation) );
    
end

