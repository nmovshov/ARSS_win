%% SI1D_TEST - a test of the PhysX engine forward integration
%
% Two solid spheres were set at rest with their centers 8 cu (code units) apart on
% the x axis. Their separation and relative velocity was saved at each time step.
% The timescale of the dynamics in cu is ~30. Runs were made with succesively
% smaller time steps, starting at 1 cu. The time step is included in the output
% file's name as well as in the header.

clear
close all
clc

%% Load data from files and do initial conversions
out_files = dir('*.out');
G = 6.674e-2;
M = 12.5664;
for k = 1:numel(out_files)
    raw(k) = importdata(out_files(k).name,' ',4); %#ok<*SAGROW>
    s = textscan(raw(k).textdata{2},'%[^=]%c%f');
    data(k).dt = s{end};
    data(k).t = raw(k).data(:,1);
    data(k).R = raw(k).data(:,2);
    data(k).V = raw(k).data(:,3);
    data(k).K = 0.5*M/2*data(k).V.^2;
    data(k).U = -G*M^2./data(k).R;
    data(k).E = data(k).K + data(k).U;
end

%% Plot of total energy vs. time
fh = figure;
ah = axes;
hold(ah,'all');
for k = 1:numel(data)
    lh = plot(data(k).t,data(k).E);
    set(lh,'DisplayName',['dt = ',num2str(data(k).dt), ' cu']);
end
% Add kinetic and potential envelope and a line of fixed energy
[~,k] = min([data.dt]);
plot(data(k).t,data(k).K,'k--')
plot(data(k).t,data(k).U,'k--')
plot([0,data(k).t(end)],[-G*M^2/8,-G*M^2/8],'k--')
legend(ah,'hide')