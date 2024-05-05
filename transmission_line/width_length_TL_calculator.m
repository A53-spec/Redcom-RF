close all
clear
clc

%% INPUT PARAMETERS
e_r = 4.2; %relative permittivity of the dielectric
Z0 = 50; % characteristic impedance
d = 1.58; % heigth of the dielectric in mm
phase_shift = 270; % phase shift required in degree
frequency = 10e9; %frequency of interest for the phase shift


%% CALCULATIONS
A = Z0/60 *sqrt((e_r+1)/2) + (e_r-1)/(e_r+1)*(0.23 + 0.11/e_r);
B = 377*pi/(2*Z0*sqrt(e_r));
W_over_d = 8*exp(A)/(exp(2*A)-2);

disp("W and l are in mm : ");
if W_over_d < 2
  W = W_over_d * d
elseif W_over_d > 2
  W =  d * (2/pi * (B - 1 -log(2*B-1) + (e_r-1)/(2*e_r)*(log(B-1) + 0.39 - 0.61/e_r)))
end

k_0 = 2*pi.*frequency ./ (3*10^8);
e_e = (e_r+1)/2 + (e_r-1)/2 * 1/(sqrt(1+12*d/W));
l = (phase_shift*pi/180)/(sqrt(e_e)*k_0) *10^3