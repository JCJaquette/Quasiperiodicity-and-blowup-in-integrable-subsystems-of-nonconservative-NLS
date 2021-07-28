% This script will generate figures displaying the periodic orbits having 
% monochromatic initial data, that is u_0(x) := A e^{i x} . 

% Both x and t are taken in the 1-torus R / 2 \pi Z. 

clear 
close all

%--------------------------%
% Computational Parameters %
%--------------------------%

% The different initial amplitudes
A_list = [ 1 2 3];

SAVE_FILE       = 1; % Option to save figures as a png
CLOSE_FIGURE    = 1; % Matlab figures can take up a lot of memory, so we can close the figures after we save them

DotsPerInch = 300; % Figure Resolution


% Set the spatial-temporal resolution to graph with:
    % Discretization of time ( from 1:2^m)
    m = 11; 
    % Discretization of space ( from 1:2^m_space)
    m_space =11;


%%%%%%%%%%%%%%%%%%%%%%%%

% We load the Fourier coefficients computed for A=1,
load('Coeff_300.mat')
% For every initial amplitude A, we rescale the coefficients to match. 
% We fix the coefficients before rescaling everything. This mitigates
% potential problems from repeated rescalings.
C_n_j_old = C_n_j;

Old_A = C_n_j{1}(1); % We expect this to be '1'. Other values are also valid


% We define the X(time) Y(space) coordinates for plotting later.
X = (1:(2^m))*( 2*pi/ (2^m));
Y = (1:(2^m_space))*( 2*pi/ (2^m_space));

% We iterate through all the initial amplitudes
for New_A = A_list
%     This rescales the Fourier Coefficients
    rescaling = New_A/Old_A;
    for n = 1:N
        C_n_j{n} = (rescaling^n)*C_n_j_old{n};            
    end 

%     This stores the value of the function to-be evaluated on a grid
    space_time = zeros(2^m_space,2^m);

%     We sum up the contribution from each Fourier mode with fft, padding 
%     to obtain a sufficient graphical resolution

    for n = 1:N % Spatial modes from c_{n,j}
        n
        padding = 0*(1:2^m);                % initialize the # of time modes
        local_space = 0*(1:2^m_space);      % initialize the # of space modes
        local_space(1+n)=2^m_space;         % add single mode corresponding to e^{inx}
        local_space = ifft(local_space);    % convert (w/ inverse fft) to grid space 

        % We create the function, a_n(t)
        % First we create a vector of length 2^m, with coefficients c_{n,j}         
        if length(C_n_j{n}) < 2^m           % Max time modes is less than 2^m, We need to pad
            padding(2:length(C_n_j{n})+1) = C_n_j{n};
        else                                % Max time modes is more than 2^m; We need to truncate time modes
            padding(2:end) = C_n_j{n}(1:2^m-1);
        end
        time = 2^m*ifft(padding);           % convert (w/ inverse fft) to grid space 
        % The contribution of a_n(t)*e^{inx} is given by the exterior
        % product:  local_space'*time 
        % We add to the total function u(x,t) evaluated on the grid. 
        space_time = space_time + local_space'*time;
    end% n = 1:N
    
%--------------------------%
%       Graphing           %
%--------------------------%


%     Now that we have stored the function value, we graph it
    figure
    set(gcf, 'Position',  [100, 100, 1200, 250])

    %     The real component!
    subplot(1,2,1)
    surf(X,Y,real(space_time),'EdgeColor','none')
%     title('Real')
    
    xlim([0,2*pi])
    X_Ticks_1st_fig = linspace(0,2*pi,9);
    xticks(X_Ticks_1st_fig )
    xticklabels({'0',' ','\pi/2',' ','\pi',' ','3\pi/2',' ','2\pi'})
    xlabel('t')
    ylim([0,2*pi])
    ylabel('x')        
    Y_Ticks_1st_fig = linspace(0,2*pi,9);
    yticks(Y_Ticks_1st_fig )
    yticklabels({'0',' ','\pi/2',' ','\pi',' ','3\pi/2',' ','2\pi'})
    
    colormap jet
    view(2) % Looking down on the surface plot
    colorbar
    
    %------The imaginary component!-------%
    subplot(1,2,2)
    s= surf(X,Y,imag(space_time),'EdgeColor','none')
%     title('Imaginary')
    
    xlim([0,2*pi])
    X_Ticks_1st_fig = linspace(0,2*pi,9);
    xticks(X_Ticks_1st_fig )
    xticklabels({'0',' ','\pi/2',' ','\pi',' ','3\pi/2',' ','2\pi'})
    xlabel('t')
    ylim([0,2*pi])
    ylabel('x')        
    Y_Ticks_1st_fig = linspace(0,2*pi,9);
    yticks(Y_Ticks_1st_fig )
    yticklabels({'0',' ','\pi/2',' ','\pi',' ','3\pi/2',' ','2\pi'})
    
    colormap jet
    colorbar
    view(2) % Looking down on the surface plot

    
%-------------------%
%      Saving       %
%-------------------%

%     We save the plot (if so desired)
    title_str = ['Plot_', num2str(New_A), '.png'];

    if (SAVE_FILE)
        exportgraphics(gcf,title_str,'Resolution',DotsPerInch);
    end
    if CLOSE_FIGURE
        close;
    end    
    
end