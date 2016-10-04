% This script will model the Rayleigh Plesset Equations
% Author: Dorien Villafranco
% Department of Mechanical Engineering, Boston University
% Requirements for ME721: Acoustic Bubble Dynamics ~ HW1
clear
clc

t0 = 0;
tf = 0.001;
y0 = [Ro, 0]';
[t,y] = ode23(@Rayleigh,[t0,tf],y0);

figure(1)
hp = plot(t,y(:,1));
grid on
set(hp,'LineWidth',2)
set(gca,'FontSize',12)
xlabel('Time (s)')
ylabel('Bubble Radius (m)')
