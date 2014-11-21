%% Euler's bench - a test of first-order integration schemes
%
% This script integrates a mass-on-a-spring system using different integration
% schemes (all first order) and compares accuracy and stability.

% run controls
clear
close all
clc
physunits on
si = setUnits;

% Uncomment next 3 lines to run in single precision (I think)
% clear all
% physunits off
% si = structfun(@single, setUnits, 'uniformoutput', false);

%% Setup
% Made-up spring system
m = 1*si.kg;
k = 1*si.newton/si.m;
A = 10*si.m;
w = sqrt(k/m);
P = 2*pi/w;

% Common integration parameters
nb_steps = 100;
t_target = 2*P;
t = linspace(0*si.s,t_target,nb_steps)';
dt = t(2)-t(1);

%% Forward Euler
x_fe(1) = A;
v_fe(1) = 0*si.m/si.s;
for j=1:nb_steps-1
    x_fe(j+1) = x_fe(j) + v_fe(j)*dt; %#ok<*SAGROW>
    v_fe(j+1) = v_fe(j) - k/m*x_fe(j)*dt;
end
T_fe = 0.5*m*v_fe.^2;
V_fe = 0.5*k*x_fe.^2;
E_fe = T_fe + V_fe;

%% Backward Euler
x_be(1) = A;
v_be(1) = 0*si.m/si.s;
fac = (1 + k/m*dt^2)^-1;
for j=1:nb_steps-1
    x_be(j+1) = fac*(x_be(j) + v_be(j)*dt);
    v_be(j+1) = fac*(v_be(j) - k/m*x_be(j)*dt);
end
T_be = 0.5*m*v_be.^2;
V_be = 0.5*k*x_be.^2;
E_be = T_be + V_be;

%% Symplectic Euler
x_se(1) = A;
v_se(1) = 0*si.m/si.s;
for j=1:nb_steps-1
    v_se(j+1) = v_se(j) - k/m*x_se(j)*dt;
    x_se(j+1) = x_se(j) + v_se(j+1)*dt;
end
T_se = 0.5*m*v_se.^2;
V_se = 0.5*k*x_se.^2;
E_se = T_se + V_se;

%% Analytic
x_an = A*cos(w*t);
v_an = -w*A*sin(w*t);
T_an = 0.5*m*v_an.^2;
V_an = 0.5*k*x_an.^2;
E_an = T_an + V_an;

%% Normalize
t_norm = t/P;
x_fe_norm = x_fe/A;
x_be_norm = x_be/A;
x_se_norm = x_se/A;
x_an_norm = x_an/A;
v_fe_norm = v_fe/w/A;
v_be_norm = v_be/w/A;
v_se_norm = v_se/w/A;
v_an_norm = v_an/w/A;
E_fe_norm = E_fe/(0.5*k*A^2);
E_be_norm = E_be/(0.5*k*A^2);
E_se_norm = E_se/(0.5*k*A^2);
E_an_norm = E_an/(0.5*k*A^2);

%% Position plot
figure;
lh = [];
ah = axes;
hold(ah,'all')
lh(end+1) = plot(t_norm,x_fe_norm,'linewidth',2,'displayname','forward euler');
lh(end+1) = plot(t_norm,x_be_norm,'linewidth',2,'displayname','backward euler');
lh(end+1) = plot(t_norm,x_se_norm,'linewidth',2,'displayname','symplectic euler');
lh(end+1) = plot(t_norm,x_an_norm,'k--','linewidth',1,'displayname','reality'); %#ok<NASGU>
legend(ah,'location','nw')
xlabel('Normalized time')
ylabel('Normalized position')

%% Energy plot
figure;
lh = [];
ah = axes;
hold(ah,'all')
lh(end+1) = plot(t_norm,E_fe_norm,'linewidth',2,'displayname','forward euler');
lh(end+1) = plot(t_norm,E_be_norm,'linewidth',2,'displayname','backward euler');
lh(end+1) = plot(t_norm,E_se_norm,'linewidth',2,'displayname','symplectic euler');
lh(end+1) = plot(t_norm,E_an_norm,'k--','linewidth',1,'displayname','reality');
legend(ah,'location','nw')
xlabel('Normalized time')
ylabel('Normalized energy')