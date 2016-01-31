function [ Final_Inflation ] = SolveInflation2_PAR( ParameterValues, SteadyStateValuesNK, ParameterValuesLearning, t, storeResults, E_K, E_S,j,h )
    
    % PURPOSE: solve the inflation given the state variable and infinite horizon expectations
    % INPUT: state variables
    % OUTPUT: inflation

    function [ InflationOut ] = calculateInflation ( inflation )
        
        if t == 1
        
            R_minus = SteadyStateValuesNK.R;
            k_minus = SteadyStateValuesNK.k;
        
            else
        
            R_minus = storeResults{j}.Actual{h}.interestRate(1,t-1);
            k_minus = storeResults{j}.Actual{h}.capital(1,t-1);
        
        end
        
        A = storeResults{j}.Actual{h}.A(1,t);
        k = storeResults{j}.Actual{h}.capital(1,t);
        
        R = SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl*inflation-SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl+SteadyStateValuesNK.R;
        
        X = -ParameterValues.theta*SteadyStateValuesNK.X/(1-ParameterValues.theta)/(1-ParameterValues.theta*ParameterValues.beta)*inflation+ParameterValues.theta*SteadyStateValuesNK.X...
            /(1-ParameterValues.theta)/(1-ParameterValues.theta*ParameterValues.beta)+SteadyStateValuesNK.X/(1-ParameterValues.theta*ParameterValues.beta)*E_K{j}{h}+SteadyStateValuesNK.X;
        
        rk = (SteadyStateValuesNK.R-SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl)/inflation+SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl-1+ParameterValues.delta;
    
        w = (1-ParameterValues.alpha)*exp(A)/(ParameterValues.alpha*exp(A))^(-ParameterValues.alpha/(1-ParameterValues.alpha))*X^(-1/(1-ParameterValues.alpha))*rk^(-ParameterValues.alpha/(1-ParameterValues.alpha));
    
        c = (1-ParameterValues.alpha)*exp(A)*k^(1-ParameterValues.eta)/(ParameterValues.alpha*exp(A))^((1-ParameterValues.eta-ParameterValues.alpha)/(1-ParameterValues.alpha))...
            *X^(-ParameterValues.eta/(1-ParameterValues.alpha))*rk^((1-ParameterValues.eta-ParameterValues.alpha)/(1-ParameterValues.alpha));

        InflationOut = ParameterValuesLearning.Cw*(w-SteadyStateValuesNK.w)+ParameterValuesLearning.Cx*(X-SteadyStateValuesNK.X)+ParameterValuesLearning.Ca*A+ParameterValuesLearning.Ck*(k-SteadyStateValuesNK.k)...
            -SteadyStateValuesNK.R*SteadyStateValuesNK.k*(inflation-1)+ParameterValuesLearning.HR*(R-SteadyStateValuesNK.R)+SteadyStateValuesNK.k*(R_minus-SteadyStateValuesNK.R)+SteadyStateValuesNK.R*(k_minus...
            -SteadyStateValuesNK.k)+E_S{j}{h}-ParameterValuesLearning.Hc*(c-SteadyStateValuesNK.c);

    end


    Final_Inflation = fzero(@calculateInflation, 1);


end