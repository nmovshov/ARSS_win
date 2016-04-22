%% Check SGP orbit accuracy
clear
close all
clc
filename = 'test.orb';

%% Load orbit log
try
    raw = importdata(filename,' ',headcount(filename));
catch
    [filename, pathname] = uigetfile('*.orb');
    filename = fullfile(pathname, filename);
    raw = importdata(filename,' ',headcount(filename));
end
disp(raw.textdata);
t = raw.data(:,1);
x = raw.data(:,2);
y = raw.data(:,3);

%% Calculate some measured diagnostics
r = hypot(x,y);
[q, ind] = min(r);
q_p = [x(ind), y(ind)];
[Q, ind] = max(r);
Q_p = [x(ind), y(ind)];
ecc = (Q - q)/(Q + q);

%% Plot orbit in lab frame
plot(x,y)
hold
plot(q_p(1), q_p(2), 'r+')
plot(0,0,'ko');
axis equal
%figure
%plot(t,r)
