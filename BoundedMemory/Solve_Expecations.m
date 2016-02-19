
function [ E_K, E_S ] = Solve_Expecations( t, SteadyStateValuesNK, ParameterValuesLearning, ParameterValues, HouseholdParameters, FirmsParameters, ActualLawOfMotion ) 
    
    % PURPOSE: given the state variables and estimated parameters at time t, it calculates the infinite horizon forecast of variables.
    % INPUT:   t - time period in which households and firms do forecasts
    %          SteadyStateValuesNK - steady state values of the  model
    %          ParameterValuesLearning - struct of learning parameters
    %          forecastPeriod - number of forecast periods
    %          ParameterValues - values of model parameters
    %          HouseholdParameters - households regression parameters
    %          FirmsParameters - firms regression parameters
    % OUTPUT:  E_K, E_S - infinite horizon forecast values from households and retailers
    % see paper for definition of E_S and E_K
    
    B1 = [ FirmsParameters.capital_Param(2,t) FirmsParameters.capital_Param(3,t);0 ParameterValues.rho ];   
    B2 = [ HouseholdParameters.capital_Param(2,t) HouseholdParameters.capital_Param(3,t);0 ParameterValues.rho ];    
    
    E_K  = 0;
    tmp1 = 0;
    tmp2 = 0;
    
    % first sum
    
    for i = 1:ParameterValues.forecastPeriod
        
        % changed B1^i to exp(i*log(B1)) to speed up the code
        
        if i == 1
            
            B1_i        = B1;
            betaTheta_i = ParameterValues.beta*ParameterValues.theta;
            
        else
            
            B1_i        = B1_i * B1;
            betaTheta_i = betaTheta_i * ParameterValues.beta*ParameterValues.theta;
            
        end
        
        tmp1 = [ FirmsParameters.inflation_Param(2,t) FirmsParameters.inflation_Param(3,t) ]*(B1_i)*[ ActualLawOfMotion.capital(1,t)-SteadyStateValuesNK.k ActualLawOfMotion.A(1,t)-SteadyStateValuesNK.A ]';
        tmp2 = [ FirmsParameters.markup_Param(2,t) FirmsParameters.markup_Param(3,t) ]*(B1_i)*[ ActualLawOfMotion.capital(1,t)-SteadyStateValuesNK.k ActualLawOfMotion.A(1,t)-SteadyStateValuesNK.A ]';
        E_K  = (betaTheta_i)*tmp1-(1-ParameterValues.theta*ParameterValues.beta)/SteadyStateValuesNK.X*(betaTheta_i)*tmp2+E_K;
  
    end

    E_S  = 0;
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    tmp4 = 0;
    tmp5 = 0;
    
    % second sum    
    
    for i = 1:ParameterValues.forecastPeriod
        
        if i == 1
            
            B2_i   = B2;
            B1_i   = B1;
            beta_i = ParameterValues.beta;
            
        else
            
            B2_i   = B2_i * B2;
            B1_i   = B1_i * B1;
            beta_i = beta_i * ParameterValues.beta;
            
        end

        tmp1 = (B2_i)*[ActualLawOfMotion.capital(1,t)-SteadyStateValuesNK.k ActualLawOfMotion.A(1,t)-SteadyStateValuesNK.A]';
        tmp2 = [HouseholdParameters.interestRate_Param(2,t) HouseholdParameters.interestRate_Param(3,t)]*(B2_i)*[ActualLawOfMotion.capital(1,t)-SteadyStateValuesNK.k ActualLawOfMotion.A(1,t)-SteadyStateValuesNK.A]';
        tmp3 = [HouseholdParameters.inflation_Param(2,t) HouseholdParameters.inflation_Param(3,t)]*(B2_i)*[ActualLawOfMotion.capital(1,t)-SteadyStateValuesNK.k ActualLawOfMotion.A(1,t)-SteadyStateValuesNK.A]';
        tmp4 = [HouseholdParameters.wage_Param(2,t) HouseholdParameters.wage_Param(3,t)]*(B2_i)*[ActualLawOfMotion.capital(1,t)-SteadyStateValuesNK.k ActualLawOfMotion.A(1,t)-SteadyStateValuesNK.A]';
        tmp5 = [HouseholdParameters.markup_Param(2,t) HouseholdParameters.markup_Param(3,t)]*(B1_i)*[ActualLawOfMotion.capital(1,t)-SteadyStateValuesNK.k ActualLawOfMotion.A(1,t)-SteadyStateValuesNK.A]';
        E_S  = ParameterValuesLearning.Ck/(1-ParameterValues.beta)*(beta_i)*tmp1(1,1)+ParameterValuesLearning.Ca/(1-ParameterValues.beta)*(beta_i)*tmp1(2,1)+(SteadyStateValuesNK.c*ParameterValuesLearning.Cc*ParameterValues.beta/(1-ParameterValues.beta)^2-SteadyStateValuesNK.M*ParameterValues.beta/SteadyStateValuesNK.R/(1-ParameterValues.beta))*(beta_i)*tmp2+(SteadyStateValuesNK.M/(1-ParameterValues.beta)-SteadyStateValuesNK.c*ParameterValuesLearning.Cc/(1-ParameterValues.beta)^2)*(beta_i)*tmp3+ParameterValuesLearning.Cw/(1-ParameterValues.beta)*(beta_i)*tmp4+ParameterValuesLearning.Cx/(1-ParameterValues.beta)*(beta_i)*tmp5+E_S;
     
    end
    

end

    


