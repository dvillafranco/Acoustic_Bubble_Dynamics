function dy = Runge_Kutta(t,y)

rho = 998;            % kg/m^3
P_a = 0  ;       % Acoustic Pressure amplitude (Pascals)
w   = 2*pi*20*10^3;   % Acoustic Driving Frequency (2*pi*20 kHz)
Pv  = 2.33*10^3;      % Vapor Pressure (Water/Air) Boundary               
P_inf_e = 1*10^5;     % Pressure far away from the bubble at equil ~ Pa
k   = 1.4;            % Polytropic Exponent of gas (air adiabitc approx.)
sigma = 0.072;        % Surface Tension for Air/Water
nu  = 0.001;          % Shear (dynamic) liquid viscosity for water in Pa*s
Ro  = 1*10^-5;          % Equlibrium Bubble Radius in meters
dy = zeros(2,1);
dy(1) = y(2);
dy(2) = -P_a/(rho*Ro).*sin(w*t)-((4*nu)/(rho*Ro^2))*y(2)-(1/(rho*Ro^2))*...
    (3*k*(P_inf_e - Pv + (2*sigma/Ro)) - (2*sigma/Ro))*y(1);