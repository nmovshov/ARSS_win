%% SPRING_TEST - a test of the PhysX engine forward integration
%
% A mass of 1 cu is placed as x=(10,0,0) cu and is subject to a force F=-kx. With
% different choices of dt, we want to look at the time evolution of x and at the
% conservation of 1/2*m*v^2 + 1/2*k*x^2

clear
close all
clc
fh = [];
ah = [];
lh = [];

%% Load data from file
out_file = '../Verify.out';
raw = importdata(out_file,' ',headcount(out_file));
t = raw.data(1:end,1);
x = raw.data(1:end,2);
v = raw.data(1:end,3);
k = 1;
dt = t(2) - t(1);

%% Plot position vs. time
fh(end+1) = figure;
ah(end+1) = axes;
hold(ah(end),'all')
lh(end+1) = plot(t,x,'linewidth',2);
set(lh(end),'displayname','PhysX 3.2')
lh(end+1) = plot(t,10*cos(1/k*t),'k--');
set(lh(end),'displayname','reality','linewidth',2)
title(['Position vs. time (dt = ',num2str(dt),')'])
xlabel('t [CU]')
ylabel('x [CU]')
legend(ah(end),'show','location','se')
grid(ah(end),'on')

%% Plot energy vs. time
fh(end+1) = figure;
ah(end+1) = axes;
hold(ah(end),'all')
lh(end+1) = plot(t,0.5*v.^2 + 0.5*k*x.^2,'linewidth',2);
set(lh(end),'displayname','PhysX 3.2')
lh(end+1) = plot([t(1) t(end)],[50*k 50*k],'k--');
set(lh(end),'displayname','reality','linewidth',2)
title(['Energy vs. time (dt = ',num2str(dt),')'])
xlabel('t [CU]')
ylabel('E [CU]')
legend(ah(end),'show','location','se')
grid(ah(end),'on')