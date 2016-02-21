classdef Kalman_Learning < handle
    
    % Purpose: do Kalman Filter learning
    
    % Input: 
    
    % Output: 
    
    
    properties
        
        % things we need to instatiate the class
        
        Previous_beta;        %  Estimated parameters, 3 by 1 vector.
        Prior_beta;
        Posterior_beta;
        
        Previous_P;           %  Estimated variance of the parameters beta
        Prior_P;
        Posterior_P;        
       
        Q_mat;                %  Variance ajustment matrix, need to initialized carefully
        
        H_vec;                %  Previous state variables, 1 by 3 vector
        
        K_vec;                %  Kalman gain, 3 by 1 vector
        S_mat;                %  Intermediate output of deriving Kalman gain
        r_parameter;          %  Adjustment parameter in Kalman gain, scalar, need to be initialized carefully
        
        variable;             % actual value of capital.
         
        
    end
    
    methods
        
        function obj = Kalman_Learning(Previous_beta,Previous_P,H_vec,Q_mat,r_parameter,variable)
            
            % class constructor
            
            obj.Previous_beta       = Previous_beta;
            obj.Previous_P          = Previous_P;
            obj.H_vec               = H_vec;
            obj.Q_mat               = Q_mat;
            obj.r_parameter         = r_parameter;
            obj.variable            = variable;
            
        end
        
        function Prior_beta_Output = Predict_Kalman_Learning(obj)
            
            % predict
            
            obj.Prior_beta = obj.Previous_beta; % here we assume that a=1
            
            Prior_beta_Output = obj.Prior_beta;
            
            obj.Prior_P    = obj.Previous_P + obj.Q_mat; % again, assuming that a=1
     
        end
        
        function Posterior_beta_Output = Update_Kalman_Learning(obj)
        
        % update
        
        % update the actual law of motion of capital before Updating Kalman
        % gain
        
        Prediction_error = obj.variable - obj.H_vec*obj.Prior_beta; % H vector should be updated in the main loop 
        
        obj.S_mat = obj.H_vec*obj.Prior_P*obj.H_vec' + obj.r_parameter;
        
        obj.K_vec = obj.Prior_P*obj.H_vec' / obj.S_mat;
        
        obj.Posterior_beta = obj.Prior_beta + obj.K_vec * Prediction_error;
        
        obj.Posterior_P = (eye(3)- obj.K_vec*obj.H_vec)*obj.Prior_P;
        
        Posterior_beta_Output = obj.Posterior_beta;
        
        obj.Previous_beta = obj.Posterior_beta;
        
        obj.Previous_P = obj.Posterior_P;
        
        
        end
        
        
    end
    
end
























