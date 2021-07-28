% This script computes (WITHOUT interval arithemetic) the space-time
% Fourier coefficients of monochromatic initial data to the 1d quadratic
% NLS:   i u_t = d_xx u + u^2. 
% 
% Used to create the data file "Coeff_300.mat"

clear 
close all

% Computational Parameters
N = 300; % N=300 takes ~ 70 sec.
A = 1;

% We create array C_n_j to store the Fourier coefficients
% In each of the slots 1 <= n <= N, it stores a vector of length n^2. 
C_n_j = cell(N,1);

C_n_j{1} = [A];
tic
for n = 2:N
    n
% We compute the product (c^2)_{n,j} = sum_{k=1}^{n-1} c_{k} * c_{n-k}
%   The variable C_n(j) here corresponds with (c^2)_{n,j-1} 

% NOTE the indexing is different from the intval version; 
%       we include j=0, but then remove it when we store the data.
    C_n = 0*(0:n^2);

    for k=1:n-1         
        nKn = conv([0 C_n_j{k}],[0 C_n_j{n-k}]);
        nKn_length = length(nKn);
        C_n(1:nKn_length) = C_n(1:nKn_length) + nKn;
    end
    
% Remove the zero mode. 
    C_n = C_n(2:end);
% Apply the operator K to c^2.
    dividor = n*n- (1:(n*n));
    dividor(n*n)=1;             % This prevents dividing by zero; note though that C_n(n*n)=0;
    C_n = C_n./dividor;
    
% Apply the operator L
    C_n(n*n) = -sum(C_n);
% Store the Fourier coefficients    
    C_n_j{n}=C_n;
    
end% for n = 2:N 
toc

