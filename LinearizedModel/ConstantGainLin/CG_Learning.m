classdef CG_Learning < handle
    
    % Purpose: do CG learning
    
    % Input: gainParameter, previous estimates of params,
    %        moments matrix, new data point
    
    % Output: parameter estimates, new value of variable
    
    
    properties
        
        % things we need to instatiate the class
        
        gainParameter;       % constant gain parameter
        previousParameters;  % estimate of previous parameters
        D_Matrix;            % moments matrix for estimation
        zMat;                % variables that are used for estimation
        variable;            % that is the previous value of variable we are estimating
        
    end
    
    methods
        
        function obj = CG_Learning(gainParam,previousParameters,D_Matrix,zMat,variable)
            
            % class constructor
            
            obj.gainParameter      = gainParam;
            obj.previousParameters = previousParameters;
            obj.D_Matrix           = D_Matrix;
            obj.zMat               = zMat;
            obj.variable           = variable;
            
        end
        
        function [ paramsOut, D_Out] = do_CG_Learning(obj)
           
            % estimate parameters using CG
            
            paramsOut = obj.previousParameters + obj.gainParameter * inv(obj.D_Matrix)*obj.zMat *( obj.variable - obj.previousParameters' * obj.zMat );
            
            % update moments matrix
            
            D_Out = obj.D_Matrix + obj.gainParameter * (obj.zMat * obj.zMat' - obj.D_Matrix);
            
        end
        
        
    end
    
    
end
























