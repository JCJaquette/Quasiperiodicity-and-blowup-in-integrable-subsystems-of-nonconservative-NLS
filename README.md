# Quasiperiodicity-and-blowup-in-integrable-subsystems-of-nonconservative-NLS
This repository contains the MATLAB codes associated with the paper: 

"Quasiperiodicity and blowup in integrable subsystems 
of nonconservative nonlinear Schroedinger equations".

by Jonathan Jaquette
___________________________________________________________________________
SUMMARY

This code computes the space-time Fourier coefficients of monochromatic 
initial data to the 1d quadratic NLS:   
    i u_t = d_xx u + u^2
where x is taken in the 1-torus R / 2 \pi Z. 

In particular, this means that the variable omega in the paper is taken to 
equal 1, and as such does not appear as a variable in the code. 
The initial data considered in this software is 
    u_0(x) := A e^{i x}
for a complex scalar A. 

With the Fourier coefficients computed, this code is set up to 
(1) Produce a computer-assisted-proof of the existence of a periodic orbit
(2) Generate figure displaying the periodic orbits
(3) Estimate the critical value of A* between where the Fourier 
    coefficients grow/decay exponentially. 

These code for part (1) requires the INTLAB interval arithmetic library: 
http://www.ti3.tu-harburg.de/rump/intlab/

===========================================================================

_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
Data Files 

Coeff_110_intval.mat    -- The Fourier coefficients c_{n,j} computed  
                            USING interval arithmetic, 
                            for A=3 & 1<=n<=110

Coeff_300.mat           -- The Fourier coefficients c_{n,j} computed  
                            NOT USING interval arithmetic, 
                            for A=1 & 1<=n<=300

(Note: By the rescaling in Theorem 3.1, the coefficients for one value of A 
can be easily used to calculate the coefficients for another value of A. 
This is utilized in the Script_ProduceFigures. However the computation for 
the rescaling will introduce some rounding error.)

_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
Script files

ComputeCoeff        - Used to create the data file "Coeff_110_intval.mat"
ComputeCoeff_intval - Used to create the data file "Coeff_300.mat"


script_CAP            - Produce a computer-assisted-proof of the existence
                        of a periodic orbit, with the data Coeff_110_intval
script_ProduceFigures - Generates figures displaying the periodic orbits
script_Rate_of_Decay  - Estimate the critical value of A* between where the
                        Fourier coefficients grow/decay exponentially.

_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
Auxiliary Functions 

quadratic_cauchy_product ... 
        -- Computes the Cauchy product of two sequences (vectors) with 
            possibly different lengths
Y_bound -- Computes the Y_0 constant from Lemma 3.6
Z_bound -- Computes the Z_1 constant and Z_2(r) polynomial from Lemma 3.6


===========================================================================

Copyright (C) 2021  Jonathan Jaquette.

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
