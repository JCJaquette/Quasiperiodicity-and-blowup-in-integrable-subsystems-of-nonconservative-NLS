function [Y0] = Y_bound(C_n_j,N)
%Y_bound -- Computes the Y_0 constant from Lemma 3.6
%           This file uses interval arithmetic. 
%  Input
%     C_n_j -   Cell array of (interval) coefficients C_{n,j}. 
%               See ComputeCoeff_intval for data structure. 
%     N     -   The Galkerin projection. 
%  Outputs 
%     Y1    -   single interval

    
%   We define the sequence b from the paper, 
%     summing the norms of the time-coefficients 
    C_n_Norms = intval(0)*(1:N);
    for n = 1:N
        C_n_Norms(n) = sum(abs(C_n_j{n}));
    end
%     'Approx norm'
%     sum(C_n_Norms)
    
%     We compute the Cauchy product (b*b)_n,    for  1<=n<=2N
    out = [0 quadratic_cauchy_product( [C_n_Norms],[C_n_Norms])' ];
%     We compute (b*b)_n / (n-1)
    dividor = intval( (1:(2*N)) - 1 );
    out = out ./ dividor;
%     Sum from N+1 to 2N
    Y0 = sum( out(N+1:end));
    
end

