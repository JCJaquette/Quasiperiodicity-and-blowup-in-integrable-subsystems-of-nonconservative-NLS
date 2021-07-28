% (Attempts to) produce a computer-assisted-proof of the existence of 
% a periodic orbit 

% We load the space-time Fourier coefficents. 
% This will also work with other outputs from "ComputeCoeff_intval".
load('Coeff_110_intval.mat')

% We compute the Y bounds and the Z bounds. 
Y0      = Y_bound(C_n_j,N);
[Z1,Z2] = Z_bound(C_n_j,N);

% The corresponding radii polynomial is a quadratic. We find its roots.
zeros_of_p = roots(mid([Z2,-1+Z1,Y0]));

% We define r_0 to be 1 pct larger than the minimal root, then rounded up
% to the nearest integer, assumed real. 
r_0 = intval( ceil(real( min(zeros_of_p))*1.01  ));

display(['r_0       = ',num2str(mid(r_0))])


% We compute the value of the radii polynomial
p_of_r0 = -r_0 + Y0 + Z1*r_0 + Z2*r_0^2 ;

display(['P(r_0)    = ',num2str(mid(p_of_r0 ))])

% If the radii polynomial is negative, then the computer assisted proof is
% successful!
if p_of_r0 < 0 
    disp('Success!')
else
    disp('Failure ')
end
