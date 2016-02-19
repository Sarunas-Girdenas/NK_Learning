function [ Final_Inflation ] = SolveInflation2( ParameterValues, SteadyStateValuesNK, ParameterValuesLearning, t, ActualLawOfMotion, E_K, E_S )
    
    % PURPOSE: solve the inflation given the state variable and infinite horizon expectations
    % INPUT: state variables
    % OUTPUT: inflation

    function [ InflationOut ] = calculateInflation ( inflation )
        
        if t == 1
        
            R_minus = SteadyStateValuesNK.R;
            k_minus = SteadyStateValuesNK.k;
        
            else
        
            R_minus = ActualLawOfMotion.interestRate(1,t-1);
            k_minus = ActualLawOfMotion.capital(1,t-1);
        
        end
        
        A = ActualLawOfMotion.A(1,t);

        k = ActualLawOfMotion.capital(1,t);

        R = SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl*inflation-SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl+SteadyStateValuesNK.R;

        X = -ParameterValues.theta*SteadyStateValuesNK.X/(1-ParameterValues.theta)/(1-ParameterValues.theta*ParameterValues.beta)*inflation+ParameterValues.theta*SteadyStateValuesNK.X...
            /(1-ParameterValues.theta)/(1-ParameterValues.theta*ParameterValues.beta)+SteadyStateValuesNK.X/(1-ParameterValues.theta*ParameterValues.beta)*E_K+SteadyStateValuesNK.X;

        denominator = (SteadyStateValuesNK.R-SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl+(SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl-1+ParameterValues.delta)*inflation);

        First_term = ParameterValuesLearning.Hc*(1-ParameterValues.alpha)*exp(A)*k^(1-ParameterValues.eta)*(ParameterValues.alpha*exp(A))^((ParameterValues.eta+ParameterValues.alpha-1)/(1-ParameterValues.alpha))*(inflation)^((ParameterValues.eta+ParameterValues.alpha-1)/(1-ParameterValues.alpha));

        Second_term = ParameterValuesLearning.Cw*(1-ParameterValues.alpha)*exp(A)*(ParameterValues.alpha*exp(A))^(ParameterValues.alpha/(1-ParameterValues.alpha))*(inflation)^(ParameterValues.alpha/(1-ParameterValues.alpha))*denominator^((ParameterValues.eta-1)/(1-ParameterValues.alpha))*X^((ParameterValues.eta-1)/(1-ParameterValues.alpha));

        Third_term = ParameterValuesLearning.Cx*X*X^(ParameterValues.eta/(1-ParameterValues.alpha))*denominator^((ParameterValues.eta+ParameterValues.alpha-1)/(1-ParameterValues.alpha));

        Fourth_term = (-ParameterValuesLearning.Ca*exp(A)-ParameterValuesLearning.Ck*k+SteadyStateValuesNK.R*SteadyStateValuesNK.k*inflation-ParameterValuesLearning.HR*R-SteadyStateValuesNK.k*R_minus-SteadyStateValuesNK.R*k_minus)*X^(ParameterValues.eta/(1-ParameterValues.alpha))*denominator^((ParameterValues.eta+ParameterValues.alpha-1)/(1-ParameterValues.alpha));

        Constant_term = (ParameterValuesLearning.Hc*SteadyStateValuesNK.c-ParameterValuesLearning.Cw*SteadyStateValuesNK.w-ParameterValuesLearning.Cx*SteadyStateValuesNK.X-ParameterValuesLearning.Ca-ParameterValuesLearning.Ck*SteadyStateValuesNK.k-ParameterValuesLearning.HR*SteadyStateValuesNK.R-SteadyStateValuesNK.R*SteadyStateValuesNK.k+E_S)*denominator^((ParameterValues.eta+ParameterValues.alpha-1)/(1-ParameterValues.alpha))*X^(ParameterValues.eta/(1-ParameterValues.alpha));

        InflationOut = First_term - Second_term - Third_term + Fourth_term - Constant_term;


    end


  Final_Inflation = fzero(@calculateInflation, 1);


end