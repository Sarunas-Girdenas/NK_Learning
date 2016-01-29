var k w R Infl X L c rk A Y;

varexo e_A;

parameters k_st L_st c_st w_st R_st rk_st X_st beta alpha delta theta epsilon eta rho r_Infl r_R;

r_Infl=0.27;

r_R=0.73;

model(linear);

//equation 1

k_st*k=(1-delta+alpha*k_st^(alpha-1)*L_st^(1-alpha))*k_st*k(-1)+k_st^alpha*L_st^(1-alpha)*A+(1-alpha)*k_st^alpha*L_st^(1-alpha)*L-c_st*c;

//equation 2

(eta-1)*L=w-c;

//equation 3

R_st*(R-Infl)=rk_st*rk;

//equation 4

-c=R-Infl(+1)-c(+1);

//equation 5

A-X+(1-alpha)*(L-k)=rk;

A=rho*A(-1)+e_A;

//equation 6

A-X+alpha*(k-L)=w;

// equation 7

beta*Infl(+1)=Infl+(1-theta*beta)*(1-theta)/theta*X;

//equation 8

R=r_Infl*(1-r_R)*Infl;

Y=A+alpha*k+(1-alpha)*L;

end;

shocks;

var e_A; // productivity shock

stderr 0.1;

end;


resid(1);

steady;

check;

stoch_simul(irf=200,order = 1);

conditional_variance_decomposition=1;
