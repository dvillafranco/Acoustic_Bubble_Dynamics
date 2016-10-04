function Rdot2 = Rayleigh(t,y)
rho = 998;        % kg/m^3
P_a = 0  ;        % Acoustic Pressure amplitude (Pascals)
w = 2*pi*20*10^3; % Acoustic Driving Frequency (2*pi*20 kHz)
Pv= 2.33*10^3;    % Vapor Pressure (Water/Air) Boundary               
P_inf_e=1*10^5;   % Pressure far away from the bubble at equil(Pa)
k = 1.4;          % Polytropic Exponent of gas (air adiabitc)
sigma = 0.072;    % Surface Tension for Air/Water
nu  = 0.001;      % Shear(dynamic) liquid viscosity for water(Pas)
Ro  = 1*10^-5;    % Equlibrium Bubble Radius in meters
Rdot = y(2);
R = y(1);
% Pressure in the gas @ equlibrium
P_ge = P_inf_e - Pv + (2*sigma/Ro);

% Pressure far away from the bubble
P_inf = P_inf_e + P_a.*sin(w*t);
P = P_ge*(Ro/R)^(3*k) + Pv - P_inf - (2*sigma/R) - (4*nu*Rdot/R);
Rdot2 = [y(2); (P/rho - (3/2)*y(2).^2)*(1./y(1))];