% This script will model the Rayleigh Plesset Equations
% Author: Dorien Villafranco
% Department of Mechanical Engineering, Boston University
% Requirements for ME721: Acoustic Bubble Dynamics ~ HW1
clear
clc
close all
%% Define Variables
%  rho = 998;            % kg/m^3
% P_a = 1*10^4  ;       % Acoustic Pressure amplitude (Pascals)
%  w   = 2*pi*20*10^3;   % Acoustic Driving Frequency (2*pi*20 kHz)
% Pv  = 2.33*10^3;      % Vapor Pressure (Water/Air) Boundary               
%  P_inf_e = 1*10^5;     % Pressure far away from the bubble at equil ~ Pa
%  k   = 1.4;            % Polytropic Exponent of gas (air adiabitc approx.)
% sigma = 0.072;        % Surface Tension for Air/Water
% nu  = 0.001;          % Shear (dynamic) liquid viscosity for water in Pa*s
Ro  = 1*10^-5;          % Equlibrium Bubble Radius in meters
% %tspan= [0 0.0001];   % Time vector in seconds
% 
t0 = 0;
tf = 0.001;
y0 = [1.5*Ro, 0]';
[t,y] = ode23(@Rayleigh2,[t0,tf],y0);
t1 = t;
% tplot = t-t0;

figure(1)
hp = plot(t,y(:,1));
grid on
set(hp,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('Time (s)')
ylabel('Bubble Radius (m)')

% Question 3
% Lauterborn's Figure 2 was digitized and points importeded into data set
% format long
% LAT = dlmread('Lauterborn_Figure2.txt');
% t_lat = LAT(:,1);
% R_Ro_lat = LAT(:,2);
% figure(2)
% hp2 = plot(tplot,y(:,1)/Ro,'x');
% hold on;
% hp3 = plot(t_lat,R_Ro_lat);
% leg1 = legend('Computed Results','Lauterborn Results');
% set(leg1,'FontSize',13)
% xlim([0 17*10^-6])
% ylim([0 2])
% grid on
% set(hp2,'LineWidth',2)
% set(hp3,'LineWidth',2)
% set(gca,'FontSize',13)
% xlabel('Time (s)')
% ylabel('R/Ro')

% Question 4
clear t
tf = 0.001;
[T,Y] = ode45(@Runge_Kutta, [t0,tf],y0);
figure(3)
hp4 = plot(T,Y(:,1)+Ro);
set(gca,'FontSize',12)
set(hp4,'LineWidth',2)
xlabel('Time (s)')
ylabel('Bubble Radius (m)')
grid on


% Quesiton 5 
 rho = 998;            % kg/m^3
 P_a = 1*10^2  ;       % Acoustic Pressure amplitude (Pascals)
 w   = 2*pi*20*10^3;   % Acoustic Driving Frequency (2*pi*20 kHz)
 Pv  = 2.33*10^3;      % Vapor Pressure (Water/Air) Boundary               
 P_inf_e = 1*10^5;     % Pressure far away from the bubble at equil ~ Pa
 k   = 1.4;            % Polytropic Exponent of gas (air adiabitc approx.)
 sigma = 0.072;        % Surface Tension for Air/Water
 nu  = 0.001;
 tspan= [0:0.000001:0.001];   % Time vector in seconds

w_sq = (3*k*P_inf_e)/(rho*Ro^2);
damp = (2*nu)/(sqrt(rho*Ro^2*3*k*P_inf_e));
M = ((1-w.^2/w_sq)).^2 + (2*damp*w.^2/w_sq)^(-0.5);
phi = -atan((2*damp*w/w_sq.^0.5)/(1-(w.^2/w_sq)));
y_approx = -P_a/(rho*Ro*w_sq)*M*sin(w*tspan+phi);

figure(4)
plot(tspan,y_approx)

figure(5)
hp10 = plot(t1,y(:,1),'x');
hold on;
hp20 = plot(T,Y(:,1)+Ro);

hp30 = plot(tspan,y_approx+Ro);
xlim([0 0.001])
set(gca,'FontSize',14)
grid on
leg2 = legend('Equation 1','Equation 3','Equation 5');
xlabel('Time (s)')
ylabel('Bubble Radius (m)')
set(leg2,'FontSize',14)
set(hp10,'LineWidth',1)
set(hp20,'LineWidth',1)
set(hp30,'LineWidth',2)
