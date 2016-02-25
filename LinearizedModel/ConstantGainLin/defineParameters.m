function [ ParametersOut ] = defineParameters()

    % PURPOSE: define struct of parameters
    
    % INPUT: None
    % OUTPUT: Struct of Parameters
    
    Gamma_Household = 0.8;
    
    Gamma_Firm = 0.8;

    ParametersOut = struct('epsilon',11,'alpha',0.33,'theta',0.75,'delta',0.025,'beta',0.95,'eta',2,'r_R',0.73,'r_Infl',0.27,'Gamma_HH',Gamma_Household,'Gamma_FF',Gamma_Firm,'rho',0.9,'forecastPeriod',100);

end