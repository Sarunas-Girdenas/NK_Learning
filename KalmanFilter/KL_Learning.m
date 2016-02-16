classdef KL_Learning < handle
    
    % Purpose: do CG learning
    
    % Input: gainParameter, previous estimates of params,
    %        moments matrix, new data point
    
    % Output: parameter estimates, new value of variable
    
    
    properties
        
        % things we need to instatiate the class
        
        previousParameters;  % estimate of previous parameters
        P_MatrixOld;         % moments matrix for estimation
        zMat;                % variables that are used for estimation
        variable;            % that is the previous value of variable we are estimating
        H_matrix;
        
    end
    
    methods
        
        function obj = KL_Learning(previousParameters,P_MatrixOld,H_matrix,zMat,variable)
            
            % class constructor
            
            obj.previousParameters = previousParameters;
            obj.P_MatrixOld        = P_MatrixOld;
            obj.zMat               = zMat;
            obj.variable           = variable;
            obj.H_matrix           = H_matrix;
            
        end
        
        function [ paramsOut, D_Out] = do_KL_Learning(obj)
            
            % predict
            
            newParameters = obj.previousParameters; % here we assume that a=1
            P_matrixNew   = P_matrixOld; % again, assuming that a=1
            
            KalmanGain = P_matrixNew/()
             
           
            
        end
        
        
    end
    
    
end
























