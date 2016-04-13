%% Check SGP orbit accuracy
clear
close all
clc
filename = '../../SGP.out';

%% Load orbit log
try
    raw = importdata(filename,' ',headcount(filename));
catch
    [filename, pathname] = uigetfile('*.out');
    filename = fullfile(pathname, filename);
    raw = importdata(filename,' ',headcount(filename));
end
disp(raw.textdata);
t = raw.data(:,1);
x = raw.data(:,2);
y = raw.data(:,3);
z = raw.data(:,4);

%% Calculate some measured diagnostics
r = hypot(x,y);
zMaxRelative = max(abs(z))/min(r);
[q, ind] = min(r);
q_p = [x(ind), y(ind)];
[Q, ind] = max(r);
Q_p = [x(ind), y(ind)];
ecc = (Q - q)/(Q + q);

%% Plot orbit in lab frame
plot(x,y)
hold
plot(q_p(1), q_p(2), 'r+')
plot(Q_p(1), Q_p(2), 'rx')
axis equal
figure
plot(t,r)
