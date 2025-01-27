%% Computes Rayleigh, Prandtl, Stefan numbers at 10C and 20C given a characteristic length scale.

clear all; close all;

%% General constants
g = 9.8; %m/s^2, gravity
T_0 = 0; %C, melting temperature
% My pinnacles
L = 0.08 / .8; %m, characteristic length scale. Initial height: 8cm.
% Rollover
% L = .05;

%% 10C
beta = 8.8 * 1e-5; %1/C, thermal expansion coefficient
T_inf = 10; %C, far-field temperature
nu = 1.3065*1e-6; %m^2/s, kinematic viscosity
kappa_T = 1.38*1e-7; %m^2/s, thermal diffusivity
c_p = 4.1955*1e3; %J / kg K, specific heat capacity
latent_heat = 3.35*1e5; %J / kg, latent heat of fusion

%% 20C
% beta = 2.07 * 1e-4;
% T_inf = 20;
% nu = 1.0035*1e-6;
% kappa_T = 1.43*1e-7;
% c_p = 4.1844*1e3;
% latent_heat = 2.4535*1e6;

%% 0C
% beta = -.000050; %1/C at 0C
% nu = 1.789*1e-6; %m^2/s at 0C
% kappa_T = .132*1e-6; %m^2/s at 0C



Ra = g * beta * (T_inf - T_0) * L^3 / (nu * kappa_T)
Pr = nu / kappa_T
St = c_p * (T_inf - T_0) / latent_heat
