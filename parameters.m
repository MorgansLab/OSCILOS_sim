%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TARS Simulink library for the simulation of thermoacoustic instabilities %
%in gas turbin combustors. README.pdf contains a short introduction to the%
%library.                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters for Simulink simulation
s = tf('s');

%% General Rijke tube constants
%Gas constants
gamma = 1.4;
R = 287.058;
cp = R / (1 - 1/gamma);

%geometric constants
xu = 0.8;
xd = 1.0;
xr = 0.5;
    
%Flame parameters
nu = 0.8; %combustion efficiency
phi = 0.7; %mean equivalence ratio
H_st = 47.161 * 10^6; %calorific value of ethene
FAR = 0.0678; %stoicheiometric fuel-to-air ratio
f = phi * FAR; %mean fuel-to-air ratio

%Mean upstream conditions
M1 = 0.11;
p1 = 101325;
T1 = 293.15;
rho1 = p1 / (R * T1);
c1 = sqrt(gamma * R * T1);
u1 = M1 * c1;

%Mean heat release per unit area
Q = nu * phi * H_st * FAR / (1 + phi * FAR) * rho1 * u1;

%Coefficients of quadratic solved for u2
a = - gamma / (gamma - 1) * rho1 * u1 + 0.5 * rho1 * u1;
b = gamma / (gamma - 1) * rho1 * u1^2 + gamma / (gamma - 1) * p1;
c = - gamma / (gamma - 1) * p1 * u1 - 0.5 * rho1 * u1^3 - Q;

u2 = (- b + sqrt(b^2 - 4 * a * c)) / (2 * a);

%Use u2 to calculate mean downstream conditions
p2 = p1 - rho1 * u1 * (u2 - u1);
rho2 = rho1 * u1 / u2;
T2 = p2 / (rho2 * R);
c2 = sqrt(gamma * R * T2);
M2 = u2 / c2;

%Reflection coefficients
upstream = menu('Choose upstream boundary condition','Open','Closed', ...
    'Choked (Marble&Candel)', 'Choked (Stow)');
if upstream == 1
    Ru = -1;
elseif upstream == 2
    Ru = 1;
elseif upstream == 3
    Ru = (1-M1)/(1+M1);
elseif upstream == 4
    Ru = (1-gamma*M1+(gamma-1)*M1^2)/(1+gamma*M1+(gamma-1)*M1^2);
end

downstream = menu('Choose downstream boundary condition','Open', ...
    'Closed','Choked (Marble&Candel)');
if downstream == 1
    Rd = -1;
    Rs = 0;
elseif downstream == 2
    Rd = 1;
    Rs = 0;
elseif downstream == 3
    Rd = (1 - 0.5 * (gamma - 1) * M2) / (1 + 0.5 * (gamma - 1) * M2);
    Rs = -((gamma * p2 / cp) * 0.5 * M2) / (1 + 0.5 * (gamma - 1) * M2);
end

%time delays
tau_u = 2 * xu / (c1 * (1 - M1^2));
tau_d = 2 * xd / (c2 * (1 - M2^2));
tau_s = xd / (u2 * (1 - M2));
tau_r1 = xr / (c2 - u2);
tau_r2 = xr / (c2 + u2);

%matrix coefficients
X11 = 1 + M1 * (rho1 * c1) / (rho2 * c2);
X12 = -1 + M1 * (2 - u2 / u1) - M1^2 * (1 - u2 / u1);
X13 = 0;
X21 = (c2 / c1) * (1 + gamma * M2) / (gamma - 1) + M1 * M2 * rho1 / rho2;
X22 = (1 - gamma * M1) / (gamma - 1) + M1^2 - M1^2 * (1 - M1) * 0.5 * ((u2 / u1)^2 - 1);
X23 = 0;
X31 = (1 + M2) / c2;
X32 = (1 - M1) / c1;
X33 = -u2 * rho2 / cp;
X = [X11, X12, X13; X21, X22, X23; X31, X32, X33];

Y11 = (M1 * (rho1 * c1) / (rho2 * c2) - 1) * Rd;
Y12 = (1 + M1 * (2 - u2 / u1) + M1^2 * (1 - u2 / u1)) * Ru;
Y13 = (M2 - 1) * Rs;
Y21 = ((c2 / c1) * (1 - gamma * M2) / (gamma - 1) + M1 * M2 * rho1 / rho2) * Rd;
Y22 = ((1 + gamma * M1) / (gamma - 1) + M1^2 - M1^2 * (1 + M1) * 0.5 * ((u2 / u1)^2 - 1)) * Ru;
Y23 = ((c2 / c1) * (1 - gamma * M2) / (gamma - 1) + M1 * M2 * rho1 / rho2) * Rs;
Y31 = ((1 - M2) / c2) * Rd;
Y32 = ((1 + M1) / c1) * Ru;
Y33 = ((1 - M2) / c2) * Rs;
Y = [Y11, Y12, Y13; Y21, Y22, Y23; Y31, Y32, Y33];

%% Constants for Simulink models
%time step
dt = 10^-5;

%Process Matrix for calculation of states B and C
M = inv(X) * Y;
M11 = M(1,1);
M12 = M(1,2);
M13 = M(1,3);
M21 = M(2,1);
M22 = M(2,2);
M23 = M(2,3);
M31 = M(3,1);
M32 = M(3,2);
M33 = M(3,3);

%Matrix characterising effect of flame on B and C
N = inv(X) * [0; Q/c1; 0];
N1 = N(1);
N2 = N(2);
N3 = N(3);

%% Flame models
%Parameters for flame model related to flow velocity fluctuations
nf = 1; %flame interaction index (strength of interaction between Q and u')
tau_f = 0.001; %time delay
tau_1 = 0.01 / (2*pi); %reciprocal of corner frequency (related to size of flame holder)

%Parameters for flame model related to equivalence ratio fluctuations
nphi = 1; %flame interaction index (strength of interaction between Q and phi')
tau_phi = 0.001; %time delay
tau_1 = 0.01 / (2*pi); %reciprocal of corner frequency (related to size of flame holder)

%% Valve models
%first order model
tau_v = 1 / (2*pi * 1000);
%second order model
tau_v1 = 1 / (2*pi * 1000);
tau_v2 = 1 / (2*pi * 800);

W_ac = 1 / (tau_v1 * s + 1) * 1 / (tau_v2 * s + 1);

%% Controllers
%Lag compensator
lag_beta = 3;
lag_omega = 134.83 * 2 * pi;
lag_a = lag_omega / sqrt(lag_beta);
lag_k = -0.2;

K_lag = lag_k * ((s + lag_beta * lag_a) / (s + lag_a))^2;

%Lead compensator
lead_beta = 6.2;
lead_omega = 134.83 * 2 * pi;
lead_a = lead_omega / sqrt(lead_beta);
lead_k = 3.0;

%% H infinity loop shaping

%Model nonlinear transfer function
H_phi = nphi * exp(-tau_phi * s) / (tau_1 * s + 1);
Z = inv(X - Y * [exp(-tau_d * s), 0, 0; 0, exp(-tau_u * s), 0; 0, 0, exp(-tau_s * s)]);
G = 1/(rho1*u1*c1^2) * (Ru*exp(-tau_u*s)-1) * [0, 1, 0] * Z * [0; Q; 0];
F = 1/(p1*c1) * (exp(-xr/(c2+u2)*s) + Rd*exp((xr/(c2-u2)-tau_d)*s)) ...
    * [1, 0, 0] * Z * [0; Q; 0];

TF = W_ac * feedback(H_phi, G) * F;

%Obtain reduced order TF for controller synthesis
w_est = logspace(log10(300), log10(2700), 30000);

resp = freqresp(TF, w_est);
sys = frd(resp, w_est);
reduced_TF = fitfrd(sys, 4);

[Kinfinity,CL,GAM,INFO]=ncfsyn(reduced_TF);

%% Adaptive controller
zc = 170;