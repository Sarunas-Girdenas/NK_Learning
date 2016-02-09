function [ SteadyStateValues ] = solveNK_SteadyState( Parameters )
    
    % PURPOSE: solve for New Keynesian model steady state
    % INPUT:   Parameters - struct of parameters
    % OUTPUT:  SteadyStateValues - struct of steady state values
    
    SteadyStateValues = struct('R',0,'X',0,'rk',0,'w',0,'L',0,'c',0,'k',0,'M',0,'A',0);
    
    % steady state solution of the model
    
    SteadyStateValues.R  = 1/Parameters.beta;
    SteadyStateValues.X  = Parameters.epsilon/(Parameters.epsilon-1);
    SteadyStateValues.rk = SteadyStateValues.R-1+Parameters.delta;
    SteadyStateValues.w  = (1-Parameters.alpha)/SteadyStateValues.X*(Parameters.alpha/SteadyStateValues.X/SteadyStateValues.rk)^(Parameters.alpha/(1-Parameters.alpha));
    SteadyStateValues.L  = (SteadyStateValues.w/((Parameters.alpha/SteadyStateValues.X/SteadyStateValues.rk)^(Parameters.alpha/(1-Parameters.alpha))-Parameters.delta*(Parameters.alpha/SteadyStateValues.X/SteadyStateValues.rk)^(1/(1-Parameters.alpha))))^(1/Parameters.eta);
    SteadyStateValues.c  = SteadyStateValues.w/SteadyStateValues.L^(Parameters.eta-1);
    SteadyStateValues.k  = (Parameters.alpha/SteadyStateValues.X/SteadyStateValues.rk)^(1/(1-Parameters.alpha))*SteadyStateValues.L;
    SteadyStateValues.M  = SteadyStateValues.k-SteadyStateValues.R*SteadyStateValues.k;
    SteadyStateValues.A  = 0;
    
end