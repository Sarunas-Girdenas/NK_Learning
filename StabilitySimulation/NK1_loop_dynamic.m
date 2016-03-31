function [residual, g1, g2, g3] = NK1_loop_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
T19 = params(2)^(1-params(9));
lhs =params(1)*y(3);
rhs =params(1)*(1-params(10)+params(9)*params(1)^(params(9)-1)*T19)*y(1)+T19*params(1)^params(9)*y(11)+T19*(1-params(9))*params(1)^params(9)*y(8)-params(3)*y(9);
residual(1)= lhs-rhs;
lhs =y(8)*(params(13)-1);
rhs =y(4)-y(9);
residual(2)= lhs-rhs;
lhs =params(5)*(y(5)-y(6));
rhs =params(6)*y(10);
residual(3)= lhs-rhs;
lhs =(-y(9));
rhs =y(5)-y(13)-y(14);
residual(4)= lhs-rhs;
lhs =y(11)-y(7)+(1-params(9))*(y(8)-y(3));
rhs =y(10);
residual(5)= lhs-rhs;
lhs =y(11);
rhs =params(14)*y(2)+x(it_, 1);
residual(6)= lhs-rhs;
lhs =y(11)-y(7)+params(9)*(y(3)-y(8));
rhs =y(4);
residual(7)= lhs-rhs;
lhs =y(13)*params(8);
rhs =y(6)+y(7)*(1-params(8)*params(11))*(1-params(11))/params(11);
residual(8)= lhs-rhs;
lhs =y(5);
rhs =y(6)*params(15)*(1-params(16));
residual(9)= lhs-rhs;
lhs =y(12);
rhs =y(11)+y(3)*params(9)+(1-params(9))*y(8);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 15);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-(params(1)*(1-params(10)+params(9)*params(1)^(params(9)-1)*T19)));
  g1(1,3)=params(1);
  g1(1,8)=(-(T19*(1-params(9))*params(1)^params(9)));
  g1(1,9)=params(3);
  g1(1,11)=(-(T19*params(1)^params(9)));
  g1(2,4)=(-1);
  g1(2,8)=params(13)-1;
  g1(2,9)=1;
  g1(3,5)=params(5);
  g1(3,6)=(-params(5));
  g1(3,10)=(-params(6));
  g1(4,5)=(-1);
  g1(4,13)=1;
  g1(4,9)=(-1);
  g1(4,14)=1;
  g1(5,3)=(-(1-params(9)));
  g1(5,7)=(-1);
  g1(5,8)=1-params(9);
  g1(5,10)=(-1);
  g1(5,11)=1;
  g1(6,2)=(-params(14));
  g1(6,11)=1;
  g1(6,15)=(-1);
  g1(7,3)=params(9);
  g1(7,4)=(-1);
  g1(7,7)=(-1);
  g1(7,8)=(-params(9));
  g1(7,11)=1;
  g1(8,6)=(-1);
  g1(8,13)=params(8);
  g1(8,7)=(-((1-params(8)*params(11))*(1-params(11))/params(11)));
  g1(9,5)=1;
  g1(9,6)=(-(params(15)*(1-params(16))));
  g1(10,3)=(-params(9));
  g1(10,8)=(-(1-params(9)));
  g1(10,11)=(-1);
  g1(10,12)=1;
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,225);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,3375);
end
end
