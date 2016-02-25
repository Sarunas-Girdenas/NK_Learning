classdef BoundedMemoryNeuralNetwork < handle
    
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
        
    end
    
    methods
        
        function obj = BoundedMemoryNeuralNetwork(windowLength,initialValue,initialValue2)
            
            % class constructor
            
            obj.windowLength      = windowLength;
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
        
        
        function [ paramsOut ] = do_BM_Perceptron(obj)
            
            % this function uses Linear Perceptron (Neural Network) to
            % approximate OLS
            
            X = [ ones(size(obj.firstInterval')) obj.firstInterval' obj.secondInterval' ];
            
            % normalize the data
            
            Y = obj.thirdInterval';
            
            % do linear perceptron now
            
            w = zeros(1,3);
            
            threshold = 1e-7; % threshold of convergence
            
            hOld    = 1;
            H       = 0;
            counter = 0;
            
            eta = 0.005;
            
            while abs(H - hOld) > threshold
                
                % stochastic gradient descent
                
                counter = counter + 1;
                hOld    = H;
                grad_t  = zeros(1,3); % initial parameters
                
                for h = 1:length(Y)
                    
                    x_t = X(h,:);
                    y_t = Y(h);
                    
                    % compute hypothesis
                    
                    H = w*x_t' - y_t;
                    grad_t = grad_t + 2*H*x_t;%*exp(-w*x_t')/((1+exp(-w*x_t'))^(2));
                                     
                end
                
                % Bold Driver Algorithm
                
                if H - hOld > 10e-4
                    
                    eta = eta*0.4;
                    
                else
                    
                    eta = eta*1.04;
                    
                end
                
                % End of Bold Driver
                      
                w = w - eta*grad_t;
                
                if counter >= 10e7
                    
                    break
                    
                end
                
            end
            
            paramsOut = w';       

        end
        
        
    end
    
    
end
























