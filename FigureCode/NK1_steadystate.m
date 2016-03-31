% The code to solve the steady state values of the simple New Keynesian Model

function [ys,check]=NK1_steadystate(ys,exe)
global M_ lgy_
global  k_st L_st c_st w_st R_st rk_st X_st 
global  beta alpha delta theta epsilon eta rho
if isfield(M_,'param_nbr') == 1
NumberOfParameters = M_.param_nbr;
for i = 1:NumberOfParameters
  paramname = deblank(M_.param_names(i,:));
  eval([ paramname ' = M_.params(' int2str(i) ');']);
end
check = 0;
end

% Parameters

epsilon = 11; % elasticity of substitution between differentiated goods

alpha = 0.33;  % share of capital 

theta = 0.75;  % probability that price is fixed for retailers

delta = 0.025; % capital depreciation rate

beta = 0.95; % time discount factor

eta =2; % labour elasticity

rho = 0.9; % decaying parameter of TFP shock

% Steady state solution

R_st=1/beta;

X_st=epsilon/(epsilon-1);

rk_st=R_st-1+delta;

w_st=(1-alpha)/X_st*(alpha/X_st/rk_st)^(alpha/(1-alpha));

L_st=(w_st/((alpha/X_st/rk_st)^(alpha/(1-alpha))-delta*(alpha/X_st/rk_st)^(1/(1-alpha))))^(1/eta);

c_st=w_st/L_st^(eta-1);

k_st=(alpha/X_st/rk_st)^(1/(1-alpha))*L_st;


% Set log-linearized deviation initially equals zero

k=0;
L=0;
c=0;
w=0;
R=0;
X=0;
Infl=0;
rk=0;
A=0;
Y=0;

% This steady state variables have been checked match the dynare file.



for iter = 1:length(M_.params)
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

if isfield(M_,'param_nbr') == 1

if isfield(M_,'orig_endo_nbr') == 1
NumberOfEndogenousVariables = M_.orig_endo_nbr;
else
NumberOfEndogenousVariables = M_.endo_nbr;
end
ys = zeros(NumberOfEndogenousVariables,1);
for i = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(i,:));
  eval(['ys(' int2str(i) ') = ' varname ';']);
end
else
ys=zeros(length(lgy_),1);
for i = 1:length(lgy_)
    ys(i) = eval(lgy_(i,:));
end
check = 0;
end
end
