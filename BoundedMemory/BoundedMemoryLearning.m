classdef BoundedMemoryLearning < handle
    
    % Purpose: do Bounded Memory learning
    
    % Input: windowLength, current values of variables (two independent
    % variables and one dependent variable)
    
    % Output: parameter estimates
    
    
    properties
        
        % things we need to instatiate the class
        
        windowLength;      % constant gain parameter
        curVarOne;         % current value of the first variable
        curVarTwo;         % current value of the second variable
        firstInterval;     % interval of the first variable, x1
        secondInterval;    % interval of the second variable, x2
        thirdInterval;     % dependent variable, y
        curVarThree;       % new obs of dependent variable
        initialParameters; % initial parameters we need in case algorithm returns NaNs, zeros or Inf
        
    end
    
    methods
        
        function obj = BoundedMemoryLearning(windowLength,initialParameters,initialValue,initialValue2)
            
            % class constructor
            
            obj.windowLength      = windowLength;
            obj.initialParameters = initialParameters;
            obj.firstInterval     = initialValue2; % initial value of capital (regression)
            obj.secondInterval    = zeros(1,1);    % value of the shock (A)
            obj.thirdInterval     = initialValue;  % initial value of variable (y in regression)
            
        end
        
        function obj = UpdateIntervals(obj)
            
            % update intervals, the length of both intervals has to be the
            % same
            
            if length(obj.firstInterval) < obj.windowLength
                
                obj.firstInterval  = [ obj.firstInterval obj.curVarOne ];
                obj.secondInterval = [ obj.secondInterval obj.curVarTwo ];
                obj.thirdInterval  = [ obj.thirdInterval obj.curVarThree ];
                
            else
                
                % update intervals
                % x1
                
                obj.firstInterval = obj.firstInterval(1,2:end);
                obj.firstInterval = [ obj.firstInterval obj.curVarOne ];
                
                % x2
                
                obj.secondInterval = obj.secondInterval(1,2:end);
                obj.secondInterval = [ obj.secondInterval obj.curVarTwo ];              
                
                % y
                
                obj.thirdInterval = obj.thirdInterval(1,2:end);
                obj.thirdInterval = [ obj.thirdInterval obj.curVarThree ]; 
                
            end
            
        end
        
        
        function [ paramsOut ] = do_BM_Learning(obj)
            
            % this function does OLS on the given intervals
            % paramsOut consists of three parameters:
            % paramsOut = [intercept varOne varTwo]
            % those parameters for three variables (intercept and two
            % variables)
            
            X = [ ones(size(obj.firstInterval')) obj.firstInterval' obj.secondInterval' ];

            paramsOut = inv((X'*X))*X'*obj.thirdInterval';
            
            % add some exceptions, double check if those are consistent
            % with the model
            
            
            
            if ( isinf(paramsOut(1,1)) ) || ( isnan(paramsOut(1,1)) )
                
                paramsOut(1,1) = obj.initialParameters(1,1);
                
            end
                
            if ( isinf(paramsOut(2,1)) ) || ( isnan(paramsOut(2,1)) )
                
                paramsOut(2,1) = obj.initialParameters(2,1);
                
            end
                
            if ( isinf(paramsOut(3,1)) ) || ( isnan(paramsOut(3,1)) )
                
                paramsOut(3,1) = obj.initialParameters(3,1);
                
            end
            

        end
        
        
    end
    
    
end
























