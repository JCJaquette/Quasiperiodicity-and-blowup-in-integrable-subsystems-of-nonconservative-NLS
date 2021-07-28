function [Z1, Z2 ] = Z_bound(C_n_j,N)
%Z_bound -- Computes the Z_1 constant and Z_2(r) polynomial from Lemma 3.6
%           This file uses interval arithmetic. 
%  Input
%     C_n_j -   Cell array of (interval) coefficients C_{n,j}. 
%               See ComputeCoeff_intval for data structure. 
%     N     -   The Galkerin projection. 
%  Outputs (both single intervals)   
%     Z1 and Z2 (the linear coefficient of the Z_2 bound) 

% The Z2 bound 
    Z2 = intval(2)/(N+1)^2;

% The Z1 bound 

% We divide each c_{n,j} by ( n^2 +2n(N+1) - j )
    for n = 1:N
        j = (1:n^2);
        divider = intval(  n^2+2*n*(N+1)- j);
        C_n_j{n} = C_n_j{n}./divider;
    end
% We sum up everything and multiply by 4.
    Kell_1 = intval(0)*(1:N);
    for n=1:N
        Kell_1(n) = sum(abs( C_n_j{n}));
    end
    
    Z1 = 4*sum(Kell_1);

end

