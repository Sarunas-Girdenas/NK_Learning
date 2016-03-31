% in this file we check the stability stability condition
% we want to compute H = e((V')^(-1))*A - I

% we have to import steady state variables and parameters


% define parameter values

ParameterValues = defineParameters(); 

% solve New Keynesian model steady state

[ SteadyStateValuesNK ] = solveNK_SteadyState( ParameterValues );

% calculate parameter values denoted in the learning tex file

[ ParameterValuesLearning ] = calculateLearningParams( SteadyStateValuesNK, ParameterValues );

numOfVars = 7;

syms eta Lst wst cst Rst alpha Xst rK theta beta Cw Hc Cx Hr Rd rPi rR kst Gk bkk akk akA rho

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

% finally, compute H matrix



H = (e'*(inv(V'))*(xiMatrix')*A - Ib);

% assign numerical values to symbols

load REmatrix_A

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
rPi   = 1.1; % we can choose this vale 
rR    = 1.1; % we can choose this vale
kst   = SteadyStateValuesNK.k;
Gk    = ParameterValuesLearning.Ck/(1-beta);
bkk   = SteadyStateValuesNK.k*(1-REmatrix_A(1,1));
akk   = REmatrix_A(1,1);
akA   = SteadyStateValuesNK.k*REmatrix_A(2,1);
rho   = ParameterValues.rho;


a = subs(H);

% now compute eigenvalues of a

eigValues = eig(eval(a));

% now write loop to evaluate the stability
% first, compute pairs of parameters

% compute pairs and save the results
rPi_Int = linspace(0.2,1.25,10); 
rR_Int  = linspace(0.2,1.25,10);

pairs = zeros([ length(rPi_Int)^2 2 ]);

i = 0;
    
for val1 = rPi_Int
        
    for val2 = rR_Int
        
            i = i + 1;
            
            pairs(i,1) = val1;
            pairs(i,2) = val2;
    
    end
        
end

rPi_Int = pairs(:,1);
rR_Int  = pairs(:,2);

save rPi_Int rPi_Int
save rR_Int rR_Int


% main loop to compute eigenvalues


% if stabilityCondition(1,g) = 0 then it is stable

stabilityCondition = zeros(1,length(pairs));

for g = 1:length(pairs)
    
    rPi   = pairs(g,1);  % we can choose this vale 
    rR    = pairs(g,2);  % we can choose this vale
    
    a         = subs(H); % this line substitutes new values of rPi and rR
    eigValues = eig(eval(a));
    
    eigValues
    
    if sum(real(eigValues) > 0) ~= 0
        
        stabilityCondition(1,g) = 1;
        
    end
    
end
















