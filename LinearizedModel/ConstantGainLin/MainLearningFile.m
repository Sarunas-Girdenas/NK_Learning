% in this file we test CG learning
tic

times = 100;

% define parameter values

ParameterValues = defineParameters(); 

% solve New Keynesian model steady state

[ SteadyStateValuesNK ] = solveNK_SteadyState( ParameterValues );

% calculate parameter values denoted in the learning tex file

[ ParameterValuesLearning ] = calculateLearningParams( SteadyStateValuesNK, ParameterValues );

% create containers to store firms and households learning data

% actual law of motion, that is, how variables evolve according to the
% model. Include all the variables

ActualLawOfMotion = struct('capital',zeros(1,times),'wage',zeros(1,times),'inflation',zeros(1,times), 'interestRate',zeros(1,times),'markup',zeros(1,times),'consumption',zeros(1,times),'labour',zeros(1,times),'capitalReturn',zeros(1,times),'A',zeros(1,times));

% Household_PLM       - households perceived law of motion 
% HouseholdParameters - households regression (learning) parameters

Household_PLM       = struct('capital',zeros(1,times),'wage',zeros(1,times),'inflation',zeros(1,times), 'interestRate',zeros(1,times),'markup',zeros(1,times),'A',zeros(1,times));
HouseholdParameters = struct('capital_Param',zeros(3,times),'wage_Param',zeros(3,times),'inflation_Param',zeros(3,times), 'interestRate_Param',zeros(3,times),'markup_Param',zeros(3,times));

% Firms_PLM       - firms perceived law of motion 
% FirmsParameters - firms regression (learning) parameters

Firms_PLM       = struct('capital',zeros(1,times),'inflation',zeros(1,times),'markup',zeros(1,times),'A',zeros(1,times));
FirmsParameters = struct('capital_Param',zeros(3,times),'inflation_Param',zeros(3,times),'markup_Param',zeros(3,times));

% generate shock

ActualLawOfMotion.A(1,1) = 0; % initial value of shock
ActualLawOfMotion.A(1,2) = 0.01; % upper bound sqrt(3)*0.01

for i = 3:times
    
    ActualLawOfMotion.A(1,i) = ParameterValues.rho * ActualLawOfMotion.A(1,i-1);    
    
end

% Rational Expectations Part
% load transition matrix from dynare (RE solution)

load REmatrix_A

% Households initial parameters in RE: capital, wage, interest rate,
% inflation and markup

HouseholdParameters.capital_Param(:,1)      = [ SteadyStateValuesNK.k*(1-REmatrix_A(1,1)) REmatrix_A(1,1) SteadyStateValuesNK.k*REmatrix_A(2,1) ]';
HouseholdParameters.wage_Param(:,1)         = [ SteadyStateValuesNK.w*(1-REmatrix_A(1,2)) REmatrix_A(1,2)*SteadyStateValuesNK.w/SteadyStateValuesNK.k SteadyStateValuesNK.w*REmatrix_A(2,2) ]';
HouseholdParameters.interestRate_Param(:,1) = [ SteadyStateValuesNK.R*(1-REmatrix_A(1,3)) REmatrix_A(1,3)*SteadyStateValuesNK.R/SteadyStateValuesNK.k SteadyStateValuesNK.R*REmatrix_A(2,3) ]';
HouseholdParameters.inflation_Param(:,1)    = [ 1-REmatrix_A(1,4) REmatrix_A(1,4)/SteadyStateValuesNK.k REmatrix_A(2,4) ]';
HouseholdParameters.markup_Param(:,1)       = [ SteadyStateValuesNK.X*(1-REmatrix_A(1,5)) REmatrix_A(1,5)*SteadyStateValuesNK.X/SteadyStateValuesNK.k SteadyStateValuesNK.X*REmatrix_A(2,5) ]';

% Firms initial parameters in RE: capital, inflation and markup

FirmsParameters.capital_Param(:,1)   = [ SteadyStateValuesNK.k*(1-REmatrix_A(1,1)) REmatrix_A(1,1) SteadyStateValuesNK.k*REmatrix_A(2,1) ]';
FirmsParameters.inflation_Param(:,1) = [ 1-REmatrix_A(1,4) REmatrix_A(1,4)/SteadyStateValuesNK.k REmatrix_A(2,4) ]';
FirmsParameters.markup_Param(:,1)    = [ SteadyStateValuesNK.X*(1-REmatrix_A(1,5)) REmatrix_A(1,5)*SteadyStateValuesNK.X/SteadyStateValuesNK.k SteadyStateValuesNK.X*REmatrix_A(2,5) ]';

% initialize variables in steady state at time 1.

ActualLawOfMotion.capital(1,1)       = SteadyStateValuesNK.k;
ActualLawOfMotion.wage(1,1)          = SteadyStateValuesNK.w;
ActualLawOfMotion.inflation(1,1)     = 1;
ActualLawOfMotion.interestRate(1,1)  = SteadyStateValuesNK.R;
ActualLawOfMotion.markup(1,1)        = SteadyStateValuesNK.X;
ActualLawOfMotion.consumption(1,1)   = SteadyStateValuesNK.c;
ActualLawOfMotion.labour(1,1)        = SteadyStateValuesNK.L;
ActualLawOfMotion.capitalReturn(1,1) = SteadyStateValuesNK.rk;

% Initialize Perceived Law of Motion (PLM) at steady state values
% Households
Household_PLM.capital(1,1)      = SteadyStateValuesNK.k;
Household_PLM.wage(1,1)         = SteadyStateValuesNK.w;
Household_PLM.interestRate(1,1) = SteadyStateValuesNK.R;
Household_PLM.inflation(1,1)    = 1;
Household_PLM.markup(1,1)       = SteadyStateValuesNK.X;
% Firms
Firms_PLM.capital(1,1)   = SteadyStateValuesNK.k;
Firms_PLM.inflation(1,1) = 1;
Firms_PLM.markup(1,1)    = SteadyStateValuesNK.X;

forecastPeriod = 100; % # of periods to compute forecast for Households and Firms

% load initial second moment matrix from RE (obtained from dynare)

load REvariance

initial_D_Matrix = [ 1 SteadyStateValuesNK.k 0; SteadyStateValuesNK.k SteadyStateValuesNK.k^2+REvariance(1,1) REvariance(1,9);0 REvariance(9,1) REvariance(9,9) ];

% initialize learning as a class. At this point we can change learning
% algorithm
        
% households
        
HH_Capital_Learning   = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.capital_Param(:,1), initial_D_Matrix, [ 1 Household_PLM.capital(1,1) ActualLawOfMotion.A(1,1) ]',ActualLawOfMotion.capital(1,1) );
HH_Wage_Learning      = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.wage_Param(:,1), initial_D_Matrix, [ 1 Household_PLM.capital(1,1) ActualLawOfMotion.A(1,1) ]',ActualLawOfMotion.wage(1,1) );
HH_Inflation_Learning = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.inflation_Param(:,1), initial_D_Matrix, [ 1 Household_PLM.capital(1,1) ActualLawOfMotion.A(1,1) ]',ActualLawOfMotion.inflation(1,1) );
HH_Interest_Learning  = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.interestRate_Param(:,1), initial_D_Matrix, [ 1 Household_PLM.capital(1,1) ActualLawOfMotion.A(1,1) ]',ActualLawOfMotion.interestRate(1,1) );
HH_Markup_Learning    = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.markup_Param(:,1), initial_D_Matrix, [ 1 Household_PLM.capital(1,1) ActualLawOfMotion.A(1,1) ]',ActualLawOfMotion.markup(1,1) );

% firms
        
FF_Capital_Learning   = CG_Learning(ParameterValues.Gamma_FF,FirmsParameters.capital_Param(:,1), initial_D_Matrix, [ 1 Firms_PLM.capital(1,1) ActualLawOfMotion.A(1,1) ]',ActualLawOfMotion.capital(1,1) );
FF_Inflation_Learning = CG_Learning(ParameterValues.Gamma_FF,FirmsParameters.inflation_Param(:,1), initial_D_Matrix, [ 1 Firms_PLM.capital(1,1) ActualLawOfMotion.A(1,1) ]',ActualLawOfMotion.inflation(1,1) );
FF_Markup_Learning    = CG_Learning(ParameterValues.Gamma_FF,FirmsParameters.markup_Param(:,1), initial_D_Matrix, [ 1 Firms_PLM.capital(1,1) ActualLawOfMotion.A(1,1) ]',ActualLawOfMotion.markup(1,1) );

% main learning loop

for t = 2:times
    
    % update the state variable (capital)
    
    ActualLawOfMotion.capital(1,t) = exp(ActualLawOfMotion.A(1,t))*ActualLawOfMotion.capital(1,t-1)^ParameterValues.alpha*ActualLawOfMotion.labour(1,t-1)^(1-ParameterValues.alpha)-ActualLawOfMotion.consumption(1,t-1)+(1-ParameterValues.delta)*ActualLawOfMotion.capital(1,t-1);
    
    % update capital in learning algorithms
    
    HH_Capital_Learning.variable = ActualLawOfMotion.capital(1,t);
    FF_Capital_Learning.variable = ActualLawOfMotion.capital(1,t);    
    
    % define zMatrix for households and firms with updated capital
    
    zMat = [1 ActualLawOfMotion.capital(1,t) ActualLawOfMotion.A(1,t)]';
    
    % update D matrix and parameters for both firms and households
    
    % households
    
    [ HouseholdParameters.capital_Param(:,t), HH_D_Out_Capital ]        = HH_Capital_Learning.do_CG_Learning();
    [ HouseholdParameters.wage_Param(:,t), HH_D_Out_Wage ]              = HH_Wage_Learning.do_CG_Learning();
    [ HouseholdParameters.inflation_Param(:,t), HH_D_Out_Inflation ]    = HH_Inflation_Learning.do_CG_Learning();
    [ HouseholdParameters.interestRate_Param(:,t), HH_D_Out_Interest ]  = HH_Interest_Learning.do_CG_Learning();
    [ HouseholdParameters.markup_Param(:,t), HH_D_Out_Markup ]          = HH_Markup_Learning.do_CG_Learning();
    
    % firms
    
    [ FirmsParameters.capital_Param(:,t), FF_D_Out_Capital ]     = FF_Capital_Learning.do_CG_Learning();
    [ FirmsParameters.inflation_Param(:,t), FF_D_Out_Inflation ] = FF_Inflation_Learning.do_CG_Learning();
    [ FirmsParameters.markup_Param(:,t), FF_D_Out_Markup ]       = FF_Markup_Learning.do_CG_Learning();
        
    % compute one step ahead forecast / PLM using updated parameters
    
    % households
    
    B_Households = [ 1 0 0; HouseholdParameters.capital_Param(:,t)'; 0 0 ParameterValues.rho ];
    
    tmpZ_HH = B_Households*zMat;
    
    Household_PLM.capital(1,t)      = tmpZ_HH(2,1);
    Household_PLM.A(1,t)            = tmpZ_HH(3,1);
    Household_PLM.wage(1,t)         = HouseholdParameters.wage_Param(:,t)'*tmpZ_HH;
    Household_PLM.interestRate(1,t) = HouseholdParameters.interestRate_Param(:,t)'*tmpZ_HH;
    Household_PLM.inflation(1,t)    = HouseholdParameters.inflation_Param(:,t)'*tmpZ_HH;
    Household_PLM.markup(1,t)       = HouseholdParameters.markup_Param(:,t)'*tmpZ_HH;
    
    B_Firms = [1 0 0; FirmsParameters.capital_Param(:,t)';0 0 ParameterValues.rho]; 
    
    tmpZ_FF = B_Firms*zMat;
    
    Firms_PLM.capital(1,t)   = tmpZ_FF(2,1);
    Firms_PLM.A(1,t)         = tmpZ_FF(3,1);
    Firms_PLM.inflation(1,t) = FirmsParameters.inflation_Param(:,t)'*tmpZ_FF;
    Firms_PLM.markup(1,t)    = FirmsParameters.markup_Param(:,t)'*tmpZ_FF;
    
    % now solve for actual law of motion (RE)
    
    [ E_K, E_S ] = Solve_Expecations( t, SteadyStateValuesNK, ParameterValuesLearning, ParameterValues, HouseholdParameters, FirmsParameters, ActualLawOfMotion );
     
    
 %   solution vector 
    
    b_vec = zeros(7,1);
    
    b_vec(2,1) = ( -ActualLawOfMotion.A(1,t) + (1-ParameterValues.alpha)*(ActualLawOfMotion.capital(1,t)/SteadyStateValuesNK.k-1) );
    
    b_vec(3,1) =  -ActualLawOfMotion.A(1,t) - ParameterValues.alpha*(ActualLawOfMotion.capital(1,t)/SteadyStateValuesNK.k-1);
    
    b_vec(5,1) = E_K;
    
    b_vec(6,1) = ParameterValuesLearning.Ca*ActualLawOfMotion.A(1,t) + ParameterValuesLearning.Ck*(ActualLawOfMotion.capital(1,t)-SteadyStateValuesNK.k) + SteadyStateValuesNK.R*(ActualLawOfMotion.capital(1,t-1)-SteadyStateValuesNK.k) + SteadyStateValuesNK.k*(ActualLawOfMotion.interestRate(1,t-1)-SteadyStateValuesNK.R)+ E_S;

    ALM_vec = ParameterValuesLearning.V_matrix*b_vec;
    
    ActualLawOfMotion.inflation(1,t) = ALM_vec(7,1) + 1;
    
    ActualLawOfMotion.interestRate(1,t) = ALM_vec(6,1) + SteadyStateValuesNK.R;
    
    ActualLawOfMotion.markup(1,t) = ALM_vec(4,1) +SteadyStateValuesNK.X;
    
    ActualLawOfMotion.capitalReturn(1,t) = ALM_vec(5,1) + SteadyStateValuesNK.rk;
    
    ActualLawOfMotion.labour(1,t) = ALM_vec(1,1) + SteadyStateValuesNK.L;
    
    ActualLawOfMotion.wage(1,t) =  ALM_vec(2,1) + SteadyStateValuesNK.w;
    
    ActualLawOfMotion.consumption(1,t) = ALM_vec(3,1) + SteadyStateValuesNK.c;
    
%     non-linear update

%     ActualLawOfMotion.inflation(1,t)     = SolveInflation2( ParameterValues, SteadyStateValuesNK, ParameterValuesLearning, t, ActualLawOfMotion, E_K, E_S );
%     ActualLawOfMotion.interestRate(1,t)  = SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl*ActualLawOfMotion.inflation(1,t)-SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl+SteadyStateValuesNK.R;
%     ActualLawOfMotion.markup(1,t)        = -ParameterValues.theta*SteadyStateValuesNK.X/(1-ParameterValues.theta)/(1-ParameterValues.theta*ParameterValues.beta)*ActualLawOfMotion.inflation(1,t)+ParameterValues.theta*SteadyStateValuesNK.X...
%    /(1-ParameterValues.theta)/(1-ParameterValues.theta*ParameterValues.beta)+SteadyStateValuesNK.X/(1-ParameterValues.theta*ParameterValues.beta)*E_K+SteadyStateValuesNK.X;
%     ActualLawOfMotion.capitalReturn(1,t) = ActualLawOfMotion.interestRate(1,t)/ActualLawOfMotion.inflation(1,t)-1+ParameterValues.delta;
%     ActualLawOfMotion.labour(1,t)        = ActualLawOfMotion.capital(1,t)*(ActualLawOfMotion.markup(1,t)*ActualLawOfMotion.capitalReturn(1,t)/ParameterValues.alpha/exp(ActualLawOfMotion.A(1,t)))^(1/(1-ParameterValues.alpha));
%     ActualLawOfMotion.wage(1,t)          = (1-ParameterValues.alpha)*exp(ActualLawOfMotion.A(1,t))*(ActualLawOfMotion.capital(1,t)/ActualLawOfMotion.labour(1,t))^ParameterValues.alpha/ActualLawOfMotion.markup(1,t);
%     ActualLawOfMotion.consumption(1,t)   = ActualLawOfMotion.wage(1,t)*ActualLawOfMotion.labour(1,t)^(1-ParameterValues.eta);
%        
    % Update learning algorithms for next iteration
    
    % households
        
    HH_Capital_Learning   = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.capital_Param(:,t), HH_D_Out_Capital, zMat,ActualLawOfMotion.capital(1,t) );
    HH_Wage_Learning      = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.wage_Param(:,t), HH_D_Out_Wage, zMat,ActualLawOfMotion.wage(1,t) );
    HH_Inflation_Learning = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.inflation_Param(:,t), HH_D_Out_Inflation, zMat,ActualLawOfMotion.inflation(1,t) );
    HH_Interest_Learning  = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.interestRate_Param(:,t), HH_D_Out_Interest, zMat,ActualLawOfMotion.interestRate(1,t) );
    HH_Markup_Learning    = CG_Learning(ParameterValues.Gamma_HH,HouseholdParameters.markup_Param(:,t), HH_D_Out_Markup, zMat,ActualLawOfMotion.markup(1,t) );

    % firms
        
    FF_Capital_Learning   = CG_Learning(ParameterValues.Gamma_FF,FirmsParameters.capital_Param(:,t), FF_D_Out_Capital, zMat,ActualLawOfMotion.capital(1,t) );
    FF_Inflation_Learning = CG_Learning(ParameterValues.Gamma_FF,FirmsParameters.inflation_Param(:,t), FF_D_Out_Inflation, zMat,ActualLawOfMotion.inflation(1,t) );
    FF_Markup_Learning    = CG_Learning(ParameterValues.Gamma_FF,FirmsParameters.markup_Param(:,t), FF_D_Out_Markup, zMat,ActualLawOfMotion.markup(1,t) );
    
        
end

toc

% figure
plot(ActualLawOfMotion.capital)
% hold
% plot(Household_PLM.capital,'r')
% figure
% plot(HouseholdParameters.capital_Param')
% figure
% plot(Household_PLM.inflation)





