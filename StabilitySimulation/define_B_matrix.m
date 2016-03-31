function [B_Out] = define_B_matrix(bkk,akk,akA,rho)
    
    % Purpose: define symbolic matrix B for symbolic computations
    % Input: bkk, akk, akA - initial values of coefficients
    %        rho - autoregressive parameter from shock
    % Output: symbolic matrix B
    % NOTE: B matrix is assumed to be 3 by 3
    
    B_Out      = sym('B',[3 3]);
    B_Out(:,:) = 0;
    B_Out(1,1) = 1;
    B_Out(2,1) = bkk;
    B_Out(2,2) = akk;
    B_Out(2,3) = akA;
    B_Out(3,3) = rho;


end