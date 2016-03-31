function [ stability ] = checkStability(eigValues)
    %Purpose: check stability of system by summing real parts of
    %eigenvalues
    % Input: eigValues: vector of eigenvalues
    % Output: 1 or 0 dependent on the sum of eigvalues real part
    
    if sum(real(eigValues) > 0) ~= 0
        
        stability = 1;
        
    else
        
        stability = 0;
        
    end

end