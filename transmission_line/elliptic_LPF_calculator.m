clear all
clc

%% INPUT PARAMETERS
fc = 3e9; %cutoff frequency of filter
Z0 = 50; % Input impedance of filter
Z_0C = 14; %Characteristic impedance of capacitors
Z_0L = 93; %Characteristic impedance of inductance

e_r = 10.8; %relative permittivity of the dielectric
d = 1.27; % heigth of the dielectric in mm


%g_L = [0.8214 0.3892 1.1880 0.7413 1.117];
%        g1      g'2     g3    g'4   g5
g_L = [0.7422 0.3313 1.2276 0.6260 1.1413];
%g_C = [1.0840 0.9077 1.1360];
%        g2      g4      g6
g_C = [1.1189 0.9746 1.0273];

L_i = 1/(2*pi*fc)*Z0*g_L;
C_i = 1/(2*pi*fc)*1/Z0 *g_C;

W_50 = width_Z_calculator (Z0, e_r, d);
W_C = width_Z_calculator (Z_0C, e_r, d);
W_L = width_Z_calculator (Z_0L, e_r, d);

e_e_C = (e_r+1)/2 + (e_r-1)/2 * 1/(sqrt(1+12*d/W_C));
e_e_L = (e_r+1)/2 + (e_r-1)/2 * 1/(sqrt(1+12*d/W_L));

guided_wave_C = 3*10^8/(fc*sqrt(e_e_C));
guided_wave_L = 3*10^8/(fc*sqrt(e_e_L));

l_Li = guided_wave_L/(2*pi) * asin(2*pi*fc*L_i/Z_0L)*10^3
l_Ci = guided_wave_C/(2*pi) * asin(2*pi*fc*C_i*Z_0C)*10^3




