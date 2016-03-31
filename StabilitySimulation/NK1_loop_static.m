function [residual, g1, g2] = NK1_loop_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 10, 1);

%
% Model equations
%

T19 = params(2)^(1-params(9));
lhs =params(1)*y(1);
rhs =y(1)*params(1)*(1-params(10)+params(9)*params(1)^(params(9)-1)*T19)+T19*params(1)^params(9)*y(9)+T19*(1-params(9))*params(1)^params(9)*y(6)-params(3)*y(7);
residual(1)= lhs-rhs;
lhs =y(6)*(params(13)-1);
rhs =y(2)-y(7);
residual(2)= lhs-rhs;
lhs =params(5)*(y(3)-y(4));
rhs =params(6)*y(8);
residual(3)= lhs-rhs;
lhs =(-y(7));
rhs =y(3)-y(4)-y(7);
residual(4)= lhs-rhs;
lhs =y(9)-y(5)+(1-params(9))*(y(6)-y(1));
rhs =y(8);
residual(5)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(14)+x(1);
residual(6)= lhs-rhs;
lhs =y(9)-y(5)+params(9)*(y(1)-y(6));
rhs =y(2);
residual(7)= lhs-rhs;
lhs =y(4)*params(8);
rhs =y(4)+y(5)*(1-params(8)*params(11))*(1-params(11))/params(11);
residual(8)= lhs-rhs;
lhs =y(3);
rhs =y(4)*params(15)*(1-params(16));
residual(9)= lhs-rhs;
lhs =y(10);
rhs =y(9)+y(1)*params(9)+(1-params(9))*y(6);
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

  %
  % Jacobian matrix
  %

  g1(1,1)=params(1)-params(1)*(1-params(10)+params(9)*params(1)^(params(9)-1)*T19);
  g1(1,6)=(-(T19*(1-params(9))*params(1)^params(9)));
  g1(1,7)=params(3);
  g1(1,9)=(-(T19*params(1)^params(9)));
  g1(2,2)=(-1);
  g1(2,6)=params(13)-1;
  g1(2,7)=1;
  g1(3,3)=params(5);
  g1(3,4)=(-params(5));
  g1(3,8)=(-params(6));
  g1(4,3)=(-1);
  g1(4,4)=1;
  g1(5,1)=(-(1-params(9)));
  g1(5,5)=(-1);
  g1(5,6)=1-params(9);
  g1(5,8)=(-1);
  g1(5,9)=1;
  g1(6,9)=1-params(14);
  g1(7,1)=params(9);
  g1(7,2)=(-1);
  g1(7,5)=(-1);
  g1(7,6)=(-params(9));
  g1(7,9)=1;
  g1(8,4)=params(8)-1;
  g1(8,5)=(-((1-params(8)*params(11))*(1-params(11))/params(11)));
  g1(9,3)=1;
  g1(9,4)=(-(params(15)*(1-params(16))));
  g1(10,1)=(-params(9));
  g1(10,6)=(-(1-params(9)));
  g1(10,9)=(-1);
  g1(10,10)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,100);
end
end
