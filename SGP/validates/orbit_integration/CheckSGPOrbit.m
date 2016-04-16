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
orbit_type = cell2mat(textscan(raw.textdata{3},'%*[^:]%*c%d'));
t = raw.data(:,1);
x = raw.data(:,2);
y = raw.data(:,3);
z = raw.data(:,4);

%% Calculate some measured diagnostics
r = hypot(x,y);
zMaxRelative = max(abs(z))/min(r);
[q, ind] = min(r);
q_p = [x(ind), y(ind)];
if orbit_type == 1
    [Q, ind] = max(r);
    Q_p = [x(ind), y(ind)];
end
if orbit_type == 1
    ecc = (Q - q)/(Q + q);
else
    vinf = cell2mat(textscan(raw.textdata{5},'%*[^V]%*[^=]%*c%f'));
    bigG = cell2mat(textscan(raw.textdata{9},'%*[^=]%*c%f'));
    bigM = cell2mat(textscan(raw.textdata{4},'%*[^=]%*c%f'));
    ecc = 1 + q*vinf^2/(bigG*bigM);
end

%% Plot orbit in lab frame
plot(x,y)
hold
plot(q_p(1), q_p(2), 'r+')
plot(0,0,'ko');
if orbit_type == 1, plot(Q_p(1), Q_p(2), 'rx'); end
axis equal
%figure
%plot(t,r)
