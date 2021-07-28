% This script computes (WITH interval arithemetic) the space-time
% Fourier coefficients of monochromatic initial data to the 1d quadratic
% NLS:   i u_t = d_xx u + u^2
% 
% Used to create the data file "Coeff_110_intval.mat" 

clear
close all

% Computational Parameters

N = 110;  % N=110 takes ~ 60 min.
A = intval(3);

% We create array C_n_j to store the Fourier coefficients
% In each of the slots 1 <= n <= N, it stores a vector of length n^2.  
C_n_j = cell(N,1);

% Note: This program runs rather slowly (especially compared to the
% non-interval version). The program is not terribly optimized/vectorized
% to run quickly in MATLAB. This could be improved in future versions. 

% Define the initial value for the first mode 
C_n_j{1} = [A];

tic
% We compute the sequence c_{n,j} according to Lemma 2.6, for 2<=n<=N
for n = 2:N
    n
% We compute the product (c^2)_{n,j} = sum_{k=1}^{n-1} c_{k} * c_{n-k}
%   The variable C_n(j) here corresponds with (c^2)_{n,j}
    C_n = intval(0)*(1:n^2);

%   Since the product is symmetric, we can just compute half the sum and
%   multiply by two, with the exception that if n is even we only 
%   add c_{n/2} * c_{n/2} once.  This reduces the computation by ~1/2. 
    
    ODD = mod(n,2);
    for k=1:  floor( (n-1)/2 )
        nKn = quadratic_cauchy_product([C_n_j{k}],[C_n_j{n-k}])';
        nKn_length = length(nKn);
        C_n(n:nKn_length+1) = C_n(n:nKn_length+1) + 2*nKn(n-1:nKn_length);
    end
    
    if ODD == 0
        k = n/2;
        nKn = quadratic_cauchy_product([C_n_j{k}],[C_n_j{n-k}])';            
        nKn_length = length(nKn);
        C_n(n:nKn_length+1) = C_n(n:nKn_length+1) + nKn(n-1:nKn_length);
    end
    
%   We apply the operator K to C_n = (c^2)_n
    dividor = intval(n*n- (1:(n*n)));
    dividor(n*n)=intval(1);             % This prevents dividing by zero; note though that C_n(n*n)=0;
    C_n = C_n./dividor;
    
% Apply the operator L    
    C_n(n*n) = -sum(fliplr(C_n));
% Store the Fourier coefficients
    C_n_j{n}=C_n;

end% for n = 2:N
toc
