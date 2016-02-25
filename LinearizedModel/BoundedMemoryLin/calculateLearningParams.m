function [ ParameterValuesLearning ] = calculateLearningParams( SteadyStateValuesNK, ParameterValues )

    % PURPOSE: calculate parameter values for learning procedure
    % INPUT:   SteadyStateValuesNK - steady state values of NK model
    %          ParameterValues - ParameterValues used in NK model
    % OUTPUT:  ParameterValuesLearning - SteadyStateValuesNKValuesNK for learning procedure
    
    % for precise parameter definition see the paper
    
    ParameterValuesLearning = struct('Cc',0,'Cw',0,'Cx',0,'Ca',0,'Ck',0,'Hc',0,'HR',0, 'V_matrix',zeros(7,7));
    
    ParameterValuesLearning.Cc = -1/(ParameterValues.eta-1)*SteadyStateValuesNK.w^(ParameterValues.eta/(ParameterValues.eta-1))*SteadyStateValuesNK.c^(-ParameterValues.eta/(ParameterValues.eta-1))...
        -(1-ParameterValues.alpha)/(ParameterValues.eta-1)*(1-1/SteadyStateValuesNK.X)*SteadyStateValuesNK.k^ParameterValues.alpha*SteadyStateValuesNK.w^((1-ParameterValues.alpha)/(ParameterValues.eta-1))...
        *SteadyStateValuesNK.c^((ParameterValues.alpha-ParameterValues.eta)/(ParameterValues.eta-1))-1;

    ParameterValuesLearning.Cw = ParameterValues.eta/(ParameterValues.eta-1)*SteadyStateValuesNK.w^(1/(ParameterValues.eta-1))*SteadyStateValuesNK.c^(-1/(ParameterValues.eta-1))...
        +(1-ParameterValues.alpha)/(ParameterValues.eta-1)*(1-1/SteadyStateValuesNK.X)*SteadyStateValuesNK.k^ParameterValues.alpha*...
        SteadyStateValuesNK.w^((2-ParameterValues.alpha-ParameterValues.eta)/(ParameterValues.eta-1))*SteadyStateValuesNK.c^(-(1-ParameterValues.alpha)/(ParameterValues.eta-1));

    ParameterValuesLearning.Cx = 1/SteadyStateValuesNK.X^2*SteadyStateValuesNK.k^ParameterValues.alpha*SteadyStateValuesNK.w^((1-ParameterValues.alpha)/(ParameterValues.eta-1))*SteadyStateValuesNK.c^(-(1-ParameterValues.alpha)/(ParameterValues.eta-1));

    ParameterValuesLearning.Ca = (1-1/SteadyStateValuesNK.X)*SteadyStateValuesNK.k^ParameterValues.alpha*SteadyStateValuesNK.w^((1-ParameterValues.alpha)/(ParameterValues.eta-1))*SteadyStateValuesNK.c^(-(1-ParameterValues.alpha)/(ParameterValues.eta-1));

    ParameterValuesLearning.Ck = ParameterValues.alpha*(1-1/SteadyStateValuesNK.X)*SteadyStateValuesNK.k^(ParameterValues.alpha-1)*SteadyStateValuesNK.w^((1-ParameterValues.alpha)/(ParameterValues.eta-1))*SteadyStateValuesNK.c^(-(1-ParameterValues.alpha)/(ParameterValues.eta-1));

    ParameterValuesLearning.Hc = -ParameterValuesLearning.Cc*(1+ParameterValues.beta/(1-ParameterValues.beta)^2);

    ParameterValuesLearning.HR = ParameterValuesLearning.Cc*SteadyStateValuesNK.c*ParameterValues.beta^2/(1-ParameterValues.beta)^2-SteadyStateValuesNK.M/SteadyStateValuesNK.R*ParameterValues.beta/(1-ParameterValues.beta);

    V_matrix(1,1) =  (ParameterValues.eta-1)/SteadyStateValuesNK.L;
    
    V_matrix(1,2) =  -1/SteadyStateValuesNK.w;
    
    V_matrix(1,3) =  1/SteadyStateValuesNK.c;
    
    V_matrix(2,1) =  (1-ParameterValues.alpha)/SteadyStateValuesNK.L;
    
    V_matrix(2,4) =  -1/SteadyStateValuesNK.X;
    
    V_matrix(2,5) =  -1/SteadyStateValuesNK.rk;
    
    V_matrix(3,1) =  -ParameterValues.alpha/SteadyStateValuesNK.L;
    
    V_matrix(3,2) =  -1/SteadyStateValuesNK.w;
    
    V_matrix(3,4) =  -1/SteadyStateValuesNK.X;
    
    V_matrix(4,5) =  -1;
    
    V_matrix(4,6) =  1;
    
    V_matrix(4,7) =  -SteadyStateValuesNK.R;
    
    V_matrix(5,4) =  (1-ParameterValues.theta*ParameterValues.beta)/SteadyStateValuesNK.X;
    
    V_matrix(5,7) = (ParameterValues.theta)/(1-ParameterValues.theta);
    
    V_matrix(6,2) =  -ParameterValuesLearning.Cw;
    
    V_matrix(6,3) =  ParameterValuesLearning.Hc;
    
    V_matrix(6,4) =  -ParameterValuesLearning.Cx;
    
    V_matrix(6,6) = -ParameterValuesLearning.HR;
    
    V_matrix(6,7) = SteadyStateValuesNK.R*SteadyStateValuesNK.k;
    
    V_matrix(7,6) = 1;
    
    V_matrix(7,7) = -SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl; 
    
    ParameterValuesLearning.V_matrix = pinv(V_matrix);

end