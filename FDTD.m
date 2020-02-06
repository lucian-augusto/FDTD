% =============== - Created by Lucian Augusto - ===============
% =============== - February, 2020 - ===============
clear all
close all
clc

%% Fundamental Constants
epsilon_0 = 8.85e-12; % Simple medium permittivity
mu_0 = 4*pi*1e-7; % Simple medium permeability
eta_0 = sqrt(mu_0/epsilon_0); % Simple medium impedance
c = 1/sqrt(epsilon_0*mu_0); % Speed of light in the free space

%% Grid Parameters
deltaZ = 1e-9; % Size of the step on the z-direction
deltaT = deltaZ/(2*c); % Size of time step (based on the Courant Condition)
Nz = 1000; % Number of steps in the z-direction
Nt = 3000; % Number of time steps
z = [1:Nz]*deltaZ; % The dimentions of z considered in the simulation
t = [1:Nt]*deltaT; % The time interval considered in the simulaiton

%% Parameters of the Electromagnetic Field Source
sourcePosition = 200; % Initial position of the EM Field source
wavelength = 500e-9; % Wavelength of of the field emitted by the source
omega = 2 * pi * c/wavelength; % Angular frequency of the wave

%% Variable Initialization
Ex(1:Nz) = 0; % Initializing the Electric Field in the space grid
Hy(1:Nz) = 0; % Initializing the Magnetic Field in the space grid 

% Creating the buffers that will receive the values at the limits of the 
% space gird to avoid the total reflection on the edges
ExBuffer1 = 0;
ExBuffer2 = 0;
HyBuffer1 = 0;
HyBuffer2 = 0;

%% Simulation Kernel (Application of Yee's Algorithim)

for n = 2:Nt
    % ------------- Electric Field Equations/Simulation -------------
    % Discrete form of Maxwell's Equations for calculating the
    % Electric Field
    Ex(2:Nz) = Ex(2:Nz) - deltaT/(epsilon_0*deltaZ) *...
        (Hy(2:Nz) - Hy(1:Nz-1));
    
    % Modeling the Electromagnetic Field Source (Sinusoidal)
    Ex(sourcePosition) = sin(omega*t(n))/2 + Ex(sourcePosition);
    
    % Fixing the total reflection on the left side of the simulation area
    % (beginning of the space grid) 
    Ex(1) = ExBuffer1;
    ExBuffer1 = ExBuffer2;
    ExBuffer2 = Ex(2);
    
    % ------------- Magnetic Field Equations/Simulation -------------
    % Discrete form of Maxwell's Equations for calculating the
    % Magnetic Field
    Hy(1:Nz-1) = Hy(1:Nz-1) - deltaT/(mu_0*deltaZ) *...
        (Ex(2:Nz) - Ex(1:Nz-1)); 
    
    %Sinusoidal Magnetic Field Source (forces unidirectional propagation)
    Hy(sourcePosition) = sin(omega*t(n))/eta_0/2 + Hy(sourcePosition);
    
    % Fixing the total reflection on the right side of the simulation area
    % (ending of the space grid) 
    Hy(Nz) = HyBuffer1;
    HyBuffer1 = HyBuffer2;
    HyBuffer2 = Hy(Nz-1);
    
    % Plotting the Fields
    subplot(2,1,1), plot(z,Ex), axis([z(1) z(Nz) -1.5 1.5]),...
        xlabel('z (nm)'), ylabel('Ex');
    subplot(2,1,2), plot(z,Hy), axis([z(1) z(Nz) -1.5/eta_0 1.5/eta_0]),...
        xlabel('z (nm)'), ylabel('Hy');
    drawnow
        
end
