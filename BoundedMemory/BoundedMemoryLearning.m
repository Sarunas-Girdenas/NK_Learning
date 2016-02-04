classdef BoundedMemoryLearning < handle
    
    % Purpose: do Bounded Memory learning
    
    % Input: windowLength, current values of variables (two independent
    % variables and one dependent variable)
    
    % Output: parameter estimates
    
    
    properties
        
        % things we need to instatiate the class
        
        windowLength;   % constant gain parameter
        curVarOne;      % current value of the first variable
        curVarTwo;      % current value of the second variable
        firstInterval;  % interval of the first variable, x1
        secondInterval; % interval of the second variable, x2
        thirdInterval;  % dependent variable, y
        curVarThree;    % new obs of dependent variable
        
    end
    
    methods
        
        function obj = BoundedMemoryLearning(windowLength)
            
            % class constructor
            
            obj.windowLength   = windowLength;
            obj.firstInterval  = zeros(1,1);
            obj.secondInterval = zeros(1,1);
            obj.thirdInterval  = zeros(1,1);
            
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
            
            paramsOut = zeros(1,3);
            
            % parameter for the first variable
            
            paramsOut(1,2) = sum((obj.firstInterval - mean(obj.firstInterval))*(obj.thirdInterval - mean(obj.thirdInterval))') / sum((obj.thirdInterval - mean(obj.thirdInterval)));
            
            % parameter for the second variable
            
            paramsOut(1,3) = sum((obj.secondInterval - mean(obj.secondInterval))*(obj.thirdInterval - mean(obj.thirdInterval))') / sum((obj.thirdInterval - mean(obj.thirdInterval)));
            
            % parameter for intercept
            
            paramsOut(1,1) = mean(obj.thirdInterval) - paramsOut(1,2)*mean(obj.firstInterval) - paramsOut(1,2)*mean(obj.secondInterval);
            
            % add some exceptions, double check if those are consistent
            % with the model !!!!!
            % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
            if ( isinf(paramsOut(1,1)) ) || ( paramsOut(1,1) == 0 ) || ( isnan(paramsOut(1,1)) )
                
                paramsOut(1,1) = 0.01*randn;
                
            end
                
            if ( isinf(paramsOut(1,2)) ) || ( paramsOut(1,2) == 0 ) || ( isnan(paramsOut(1,2)) )
                
                paramsOut(1,2) = 0.01*randn;
                
            end
                
            if ( isinf(paramsOut(1,3)) ) || ( paramsOut(1,3) == 0 ) || ( isnan(paramsOut(1,3)) )
                
                paramsOut(1,3) = 0.01*randn;
                
            end
            
            % make it consistent with the existing struct in the rest of
            % the code
            
            paramsOut = paramsOut';

        end
        
        
    end
    
    
end
























