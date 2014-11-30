%% Euler's bench - a test of first-order integration schemes
%
% This script integrates a mass-on-a-spring system using different integration
% schemes (all first order) and compares accuracy and stability.

% run controls
clear
close all
clc
physunits off
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
nb_steps = 8999;
t_target = 35.9960*si.s;
% nb_steps = 100;
% t_target = P;
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

%% Half-step
x_hs(1) = A;
v_hs(1) = 0*si.m/si.s - k/m*x_hs(1)*dt/2; % v_1/2
for j=1:nb_steps-1
    v_hs(j+1) = v_hs(j) - k/m*x_hs(j)*dt; % v_j+1/2
    x_hs(j+1) = x_hs(j) + v_hs(j+1)*dt;
end
T_hs = 0.5*m*v_hs.^2;
V_hs = 0.5*k*x_hs.^2;
E_hs = T_hs + V_hs;

%% Mid-point
x_mp(1) = A;
v_mp(1) = 0*si.m/si.s;
for j=1:nb_steps-1
    v_mp(j+1) = v_mp(j) - k/m*x_mp(j)*dt;
    x_mp(j+1) = x_mp(j) + 0.5*(v_mp(j) + v_mp(j+1))*dt;
end
T_mp = 0.5*m*v_mp.^2;
V_mp = 0.5*k*x_mp.^2;
E_mp = T_mp + V_mp;

%% Leapfrog
x_lf(1) = A;
v_lf(1) = 0*si.m/si.s;
for j=1:nb_steps-1
    x_jph = x_lf(j) + 0.5*v_lf(j)*dt;
    v_lf(j+1) = v_lf(j) - k/m*x_jph*dt;
    x_lf(j+1) = x_jph + 0.5*v_lf(j+1)*dt;
end
T_lf = 0.5*m*v_lf.^2;
V_lf = 0.5*k*x_lf.^2;
E_lf = T_lf + V_lf;

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

%% PhysX
out_file = '/Verify.out';
raw = importdata(out_file,' ',headcount(out_file));
t_px = raw.data(:,1);
x_px = raw.data(:,2);
v_px = raw.data(:,3);
T_px = 0.5*v_px.^2;
V_px = 0.5*x_px.^2;
E_px = T_px + V_px;

%% Analytic
x_an = A*cos(w*t);
v_an = -w*A*sin(w*t);
T_an = 0.5*m*v_an.^2;
V_an = 0.5*k*x_an.^2;
E_an = T_an + V_an;

%% Normalize
t_norm = t/P;
t_px_norm = t_px/(2*pi);
x_fe_norm = x_fe/A;
x_be_norm = x_be/A;
x_hs_norm = x_hs/A;
x_mp_norm = x_mp/A;
x_lf_norm = x_lf/A;
x_se_norm = x_se/A;
x_an_norm = x_an/A;
x_px_norm = x_px/10;
v_fe_norm = v_fe/w/A;
v_be_norm = v_be/w/A;
v_hs_norm = v_hs/w/A;
v_mp_norm = v_mp/w/A;
v_lf_norm = v_lf/w/A;
v_se_norm = v_se/w/A;
v_an_norm = v_an/w/A;
v_px_norm = v_px/10;
E_fe_norm = E_fe/(0.5*k*A^2);
E_be_norm = E_be/(0.5*k*A^2);
E_hs_norm = E_hs/(0.5*k*A^2);
E_mp_norm = E_mp/(0.5*k*A^2);
E_lf_norm = E_lf/(0.5*k*A^2);
E_se_norm = E_se/(0.5*k*A^2);
E_an_norm = E_an/(0.5*k*A^2);
E_px_norm = E_px/50;


%% Position plot
figure;
lh = [];
ah = axes;
hold(ah,'all')
lh(end+1) = plot(t_norm,x_fe_norm,'linewidth',2,'displayname','forward euler');
lh(end+1) = plot(t_norm,x_be_norm,'linewidth',2,'displayname','backward euler');
lh(end+1) = plot(t_norm,x_se_norm,'linewidth',2,'displayname','symplectic euler');
lh(end+1) = plot(t_norm,x_hs_norm,'linewidth',2,'displayname','half-step');
lh(end+1) = plot(t_norm,x_mp_norm,'linewidth',2,'displayname','mid-point');
lh(end+1) = plot(t_norm,x_lf_norm,'linewidth',2,'displayname','leapfrog');
lh(end+1) = plot(t_px_norm,x_px_norm,'linewidth',2,'displayname','physx');
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
lh(end+1) = plot(t_norm,E_hs_norm,'linewidth',2,'displayname','half-step');
lh(end+1) = plot(t_norm,E_mp_norm,'linewidth',2,'displayname','mid-point');
lh(end+1) = plot(t_norm,E_lf_norm,'linewidth',2,'displayname','leapfrog');
lh(end+1) = plot(t_px_norm,E_px_norm,'linewidth',2,'displayname','physx');
lh(end+1) = plot(t_norm,E_an_norm,'k--','linewidth',1,'displayname','reality');
legend(ah,'location','nw')
xlabel('Normalized time')
ylabel('Normalized energy')