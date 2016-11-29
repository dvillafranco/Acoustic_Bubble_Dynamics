function Rdot2 = Rayleigh(t,y)
global P sigma k rho P_a nu Ro p_stat Pv f_use w P_inf_e
%f = 100*10^3:10:100*10^6;
rho = 998;            % kg/m^3
P_a = 7*10^4  ;            % Acoustic Pressure amplitude (Pascals)
w   = 2*pi*f_use;   % Acoustic Driving Frequency (2*pi*20 kHz)
Pv  = 2.33*10^3;      % Vapor Pressure (Water/Air) Boundary               
P_inf_e = 1*10^5;     % Pressure far away from the bubble at equil ~ Pa
k   = 1.33;            % Polytropic Exponent of gas (air adiabitc approx.)
sigma = 0.0725;        % Surface Tension for Air/Water
nu  = 0.001;          % Shear (dynamic) liquid viscosity for water in Pa*s
Ro  = 10*10^-6;        % Equlibrium Bubble Radius in meters
P_c = 1*10^5;
Rdot = y(2);
R = y(1);
p_stat = P_inf_e;
% Pressure in the gas @ equlibrium
P_ge = P_inf_e - Pv + (2*sigma/Ro);

% Pressure far away from the bubble
 P_inf =  P_inf_e - P_a.*sin(w*t);
% P_inf = P_inf_e + P_c*heaviside(t-0.0002);
P = P_ge*(Ro/R)^(3*k) + Pv - P_inf - (2*sigma/R) - (4*nu*Rdot/R);

Rdot2 = [y(2); (P/rho - (3/2)*y(2).^2)*(1./y(1))];