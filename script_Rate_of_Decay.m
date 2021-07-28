% % %  
% This script loads the non-rigorous computation of Fourier coefficents 
% from 'Coeff_300', where it is assumed that "A=1". 
% It then make a scatter plot, and calculate the rate of decay, from which
% we get the critical value $A*$.

% To see the residual error, on the figure which is produced 
% go to Tools->Basic Fitting, 
% then select 'Linear', and then under 'Error Estimation (Residuals)',
% select a Plot Style


clear
close all
load('Coeff_300.mat')
format long

%%% Computational Parameters

trunc_m = 100;  % lower truncation of 'n' for fitting data
trunc_M = 300;  % upper truncation of 'n' for fitting data
trunc_M = min([trunc_M,N]);


% = = = = = = = = =

% We calculate the ell^1, ell^2, and ell^infty norms, 
% by summing the c_n_j across j 
ell_1 = 0*(1:N);
ell_2 = 0*(1:N);
ell_infty = 0*(1:N);
for n=1:N
    ell_1(n) = sum(abs( C_n_j{n}));
    ell_2(n) = sqrt(sum(abs( C_n_j{n}.^2  )));    
    ell_infty(n) = max(abs( C_n_j{n}));
end

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

% For each n (between trunc_m & trunc_M) we plot the log_10 of these norms

figure
hold on 
scatter(trunc_m:trunc_M,log(ell_1(trunc_m:trunc_M))/log(10))
scatter(trunc_m:trunc_M,log(ell_2(trunc_m:trunc_M))/log(10))
scatter(trunc_m:trunc_M,log(ell_infty(trunc_m:trunc_M))/log(10))
hold off

leg = legend('$\ell_{1}$','$\ell_{2}$','$\ell_{\infty}$');
leg.Interpreter = 'latex';
xlabel('n');
ylabby=ylabel('$ \log_{10} \| c_{n,\cdot}\|$');
ylabby.Interpreter = 'latex';


% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

%  We now compute an estimate for A^*.

%  Recall that we've defined ell_1, ell_2, ell_\infty such that 
%    ell_1(n) = \sum_{n<=j<=n^2} | c_{n,j| 
%  	 ell_2(n) = ...

%  Using a linear regression, we find a line of best fit (y = mx +b) where:
%       x = n,    y(x) = log ( ell_1(n) )

mdl_1 = fitlm((trunc_m:trunc_M),log(ell_1(trunc_m:trunc_M)))
mdl_2 = fitlm((trunc_m:trunc_M),log(ell_2(trunc_m:trunc_M)))
mdl_infty = fitlm((trunc_m:trunc_M),log(ell_infty(trunc_m:trunc_M)))

% Below is the predicted slope coefficients
predicted_slope_1 = mdl_1.Coefficients{2,1};
predicted_slope_2 = mdl_2.Coefficients{2,1};
predicted_slope_infty = mdl_infty.Coefficients{2,1};

% We transform the equation so we get the predictor 
%       ell_1(n) ~  B e^{ m n} = B  ( exp(m) )^n
% Since A=1 in our numerical data, we get A* = 1/exp(m). 
disp([newline,'Estimates for A*'])
A_star_1 = exp(-predicted_slope_1)
A_star_2 = exp(-predicted_slope_2)
A_star_infty = exp(-predicted_slope_infty)

% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

disp([newline,'Goodness of fit'])
disp('  R^2 for ell_1')
disp(mdl_1.Rsquared)

disp('  R^2 for ell_2')
disp(mdl_2.Rsquared)

disp('  R^2 for ell_infty')
disp(mdl_infty.Rsquared)
