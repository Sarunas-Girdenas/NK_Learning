clear variables

N = 300; % Initial Forecasting Horizon (See Binder and Pesaran, 1995)


msize     = 7;                   %  size of state space
rsize     = 1;                    %  number of shocks


%%%%% Parameters for the BGM Calibration %%%%%

beta      = 0.99;                     %  discount rate
sigma     = 1/1.3571;                 %  markup/subs elasticity
delta     = 0.025;                    % 10% rate of destruction of firms a year
lambda    = 0 ;                       % as in BGM;

nu    = 0.58333;
L     = 0.3;
fe    = 1;
f     = 0;

rho       = 0.979;                     %  persistence of the shock

%%%%% Composite Parameters %%%%%

%alpha_1 = (n*l)/L;
%alpha_2 = (f*n)/L;
%alpha_3 = (fe*ne)/L;

%%%%% which were calculated with SS file as %%%%%

alpha_1  = 0.7972;
alpha_2  = 0;
alpha_3  = 0.2028;


Lambda = alpha_1/(alpha_1 - (alpha_2*(sigma/(1-sigma))));


% Set Up the System      

% Chat * Z(t) = Ahat * Z(t-1) + Bhat * E(Z(t+1)|I(t))   %
%               + D1 * W(t) + D2 * E(W(t+1)|I(t))       %
% W(t) = R * W(t-1) + v(t)                              %

% Z MATRIX CAN BE INTERPRETED AS FOLLOWS                %
%                                                       %
% Row 1:      Labor Supply                              %
% Row 2:      Consumption                               %
% Row 3:      Firm Profits                              %
% Row 4:      Wages                                     %
% Row 5:      Firms                                     %
% Row 6:      New Firms                                 %
% Row 7:      TFP Shock                                 %


% Coefficient Matrices %


Chat = zeros(msize,msize); % 
Ahat = zeros(msize,msize); %
Bhat = zeros(msize,msize); % 

D1 = zeros(msize,1); % 
D2 = zeros(msize,1); % 

R = zeros(1,1);


% Filling Elements of the Matrices %

% Filling Elements of the Matrices                       %
% PERIOD T is Chat                                       %
% PERIOD T-1 is Ahat                                     %
% PERIOD T+1 is Bhat                                     %
% SHOCK is D                                             %

% EQUATION 1: Labor Supply

Chat(1,1) = 1;
Chat(1,2) = -alpha_1;
Chat(1,5) = -(alpha_2 + (alpha_1*((sigma-1)/sigma)));
Chat(1,6) = -alpha_3;
Chat(1,7) = 1;


% EQUATION 2: Firm Net Worth                                  %


Chat(2,2) = -Lambda;
Chat(2,3) = 1;
Chat(2,4) = Lambda-1;
Chat(2,5) = Lambda;
Chat(2,7) = 1-Lambda;


% EQUATION 3: Shares Euler                                         %

Chat(3,2)  = -1;
Chat(3,5)  = (1-sigma)/sigma;


Bhat(3,2)  = -1;
Bhat(3,3)  = 1-(beta*(1-delta));
Bhat(3,5)  = (beta*(1-delta))*((1-sigma)/sigma);

% EQUATION 4: Labor-Leisure                                         %

Chat(4,1) = nu*(L/(1-L));
Chat(4,2) = 1;
Chat(4,4) = -1;


% EQUATION 5: Free Entry                                           %

Chat(5,4) = -1;
Chat(5,5) = (1-sigma)/sigma;
Chat(5,7) = 1;

% EQUATION 6: Firm Dynamics                                           %

Chat(6,5) = 1/delta;

Ahat(6,5) = (1-delta)/delta;
Ahat(6,6) = 1;


% EQUATION 7: Shock Process                                           

Chat(7,7) = 1;

Ahat(7,7) = rho;

% SHOCK

D1(7,1) = 1; % a 1 percent shock

% Transform System to Canonical Form:


% x(t) = A * x(t-1) + B * E(x(t+1)|I(t)) 
%        + inv(Chat) * D1 * w(t) + inv(Chat) * D2 * E(w(t+1)|I(t))

 G = inv(Chat);
 B = inv(Chat)*Bhat;
 A = inv(Chat)*Ahat;


 dim1 = size(A,1);
 dim3 = size(R,1);

 
% Compute Matrix C Using Brute-Force Iterative Procedure


  C = eye(dim1);          % Initial Condition
  Y = eye(dim1);          % Initial Condition
  eps1 = 10^(-6);         % Convergence Criterion for F
  eps2 = 10^(-6);         % Convergence Criterion for C
  crit1 = 1; crit2 = 1;   % Initial Conditions
  iter = 0;

 
 while crit1 >= eps1 | crit2 >= eps2

    Yi = (eye(dim1)-B*C)\B;
    Ci = (eye(dim1)-B*C)\A;

    crit1 = max(max(abs(Yi-Y))); crit2 = max(max(abs(Ci-C)));
    C = Ci; Y = Yi;

    iter = iter+1;
    if iter > 100, 
       disp(' The brute-force iterative procedure did not converge after ')
       disp(' 100 iterations. See Binder and Pesaran (1995, 1997) for alternative ')
       disp(' algorithms to compute the matrix C. '),
    end
 end
 
 
 
 
 % Use Recursive Method of Binder and Pesaran (1995) to Compute the 
% Forward Part of the Solution - Determine N
 eps3 = 10^(-6);        % Convergence Criterion
 i = 0;
 aux3a = zeros(dim1,dim3);
 aux3b = zeros(dim1,dim3);
 while i <= N
    fp1 = Y^i/(eye(dim1)-B*C)/Chat*D1*R^i;
    fp2 = Y^i/(eye(dim1)-B*C)/Chat*D2*R^(i+1);
    aux3a = fp1+aux3a;
    aux3b = fp2+aux3b;
    i = i+1;
 end


 crit3 = max(max(abs(fp1+fp2)));

 while crit3 > eps3
    N = N+1;
    fp1 = F^N/(eye(dim1)-B*C)/chat*D1*R^N;
    fp2 = F^N/(eye(dim1)-B*C)/chat*D2*R^(N+1);
    aux3a = fp1+aux3a;
    aux3b = fp2+aux3b;
    crit3 = max(max(abs(fp1+fp2)));
 end
  H = aux3a+aux3b;


  z = H;
  for k = 2:1000
      z(:,k) = C*z(:,k-1);
  end
  
  
   % FIGURE TITLE: Impulse Responses to a TFP Shock
   
   Z=z';
   x=80;
   t=1:1:x;
   
   figure(1)
   
  
   subplot(2,3,1)
   plot(t,Z(1:x,2),'--','color',[0 0 0])
   hold on
   plot(t,Z(1:x,7),'-','color',[0 0 1])
   hold off
   grid on
   title('Consumption (C)');   
  
   subplot(2,3,2)
   plot(t,Z(1:x,6),'--','color',[0 0 0])
   hold on
   plot(t,Z(1:x,7),'-','color',[0 0 1])
   hold off
   grid on
   title('Entry (ne)');
   
   subplot(2,3,3)
   plot(t,Z(1:x,5),'--','color',[0 0 0])
   hold on
   plot(t,Z(1:x,7),'-','color',[0 0 1])
   hold off
   grid on
   title('Number of Products (n)');   

   subplot(2,3,4)
   plot(t,Z(1:x,3),'--','color',[0 0 0])
   hold on
   plot(t,Z(1:x,7),'-','color',[0 0 1])
   hold off
   grid on
   title('Profits (phi)');   
  
   subplot(2,3,5)
   plot(t,Z(1:x,4),'--','color',[0 0 0])
   hold on
   plot(t,Z(1:x,7),'-','color',[0 0 1])
   hold off
   grid on
   title('Wages (w)');   
  
   subplot(2,3,6)
   plot(t,Z(1:x,1),'--','color',[0 0 0])
   hold on
   plot(t,Z(1:x,7),'-','color',[0 0 1])
   hold off
   grid on
   title('Labor Supply (L)');   
  
     
  
   