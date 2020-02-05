% =============== - Created by Lucian Augusto - ===============
% =============== - December, 2019 - ===============

clear all
close all
clc

%% Fundamental Constants
epsilon_0 = 8.85e-12; % Simple medium permitivity
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
sourcePosition = 200;
wavelength = 500e-9;
omega = 2 * pi * c/wavelength; % Angular frequency of the wave
pulse_centre = 200*deltaT;
pulse_width = 50*deltaT;

%% Variable Initialisation
Ex(1:Nz) = 0;
Hy(1:Nz) = 0;

%% Simulation Kernel (Application of Yee's Algorithim)

for n = 2:Nt
    % Discrete form of Maxwell's Equations
    Ex(2:Nz) = Ex(2:Nz) - deltaT/(epsilon_0*deltaZ) *...
        (Hy(2:Nz) - Hy(1:Nz-1));
    
    % Modeling the Electromagnetic Field Source (Sinusoidal)
    Ex(sourcePosition) = sin(omega*t(n))/2 + Ex(sourcePosition);

    Hy(1:Nz-1) = Hy(1:Nz-1) - deltaT/(mu_0*deltaZ) *...
        (Ex(2:Nz) - Ex(1:Nz-1)); 
    
    Hy(sourcePosition) = sin(omega*t(n))/eta_0/2 + Hy(sourcePosition); %Sinusoidal Magnetic Field Source (forces unidirectional propagation)
    
    % Plotting the Fields
    subplot(2,1,1), plot(z,Ex), axis([z(1) z(Nz) -1.5 1.5]),...
        xlabel('z (nm)'), ylabel('Ex');
    subplot(2,1,2), plot(z,Hy), axis([z(1) z(Nz) -1.5/eta_0 1.5/eta_0]),...
        xlabel('z (nm)'), ylabel('Hy');
    drawnow
        
end
