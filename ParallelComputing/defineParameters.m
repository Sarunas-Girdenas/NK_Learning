function [ ParametersOut ] = defineParameters()

    % PURPOSE: define struct of parameters
    
    % INPUT: None
    % OUTPUT: Struct of Parameters
    

    ParametersOut = struct('epsilon',11,'alpha',0.33,'theta',0.75,'delta',0.025,'beta',0.95,'eta',2,'r_R',0.73,'r_Infl',0.27,'gamma',0.1,'rho',0.9,'forecastPeriod',100);

end