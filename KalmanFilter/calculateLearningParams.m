function [ ParameterValuesLearning ] = calculateLearningParams( SteadyStateValuesNK, ParameterValues )

    % PURPOSE: calculate parameter values for learning procedure
    % INPUT:   SteadyStateValuesNK - steady state values of NK model
    %          ParameterValues - ParameterValues used in NK model
    % OUTPUT:  ParameterValuesLearning - SteadyStateValuesNKValuesNK for learning procedure
    
    % for precise parameter definition see the paper
    
    ParameterValuesLearning = struct('Cc',0,'Cw',0,'Cx',0,'Ca',0,'Ck',0,'Hc',0,'HR',0);
    
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



end