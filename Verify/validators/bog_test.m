%% BOG_TEST - a test of the PhysX engine forward integration
%
% A ball of radius 1 cu was dropped from a height of 6 cu. Ground plane and
% constant gravity (g=10 cu) implemented by PhysX promitives.

clear
close all
clc

%% Load data from file
out_file = '../Verify.out';
raw = importdata(out_file,' ',4);
t = raw.data(1:end,1);
y = raw.data(1:end,2);
v = raw.data(1:end,3);
g = 10;
dt = t(2) - t(1);

%% Plot of fractional delta energy vs. time
E = 0.5*v.^2 + g*y;
dE = (E - E(1))/E(1);
fh1 = figure;
ah1 = axes;
hold(ah1,'all')
lh1 = plot(t,dE);
title('Fractional delta energy vs time')
xlabel('t [CU]')
ylabel('E [CU]')

%% Plot position vs time
fh2 = figure;
ah2 = axes;
hold(ah2,'all')
lh2 = plot(t,y,'linewidth',2);
lh3 = plot(t,6 - 0.5*g*t.^2,'k--');
set(lh2,'displayname','PhysX 3.2')
set(lh3,'displayname','reality')
title('Position vs time')
xlabel('t [CU]')
ylabel('y [CU]')
legend(ah2,'show')

%% Plot v vs t
fh3 = figure;
ah3 = axes;
hold(ah3,'all')
lh4 = plot(t,v,'linewidth',2);
lh5 = plot(t,-g*t,'k--');
set(lh4,'displayname','PhysX 3.2')
set(lh5,'displayname','reality')
title('Velocity vs time')
xlabel('t [CU]')
ylabel('v [CU]')
legend(ah3,'show')
