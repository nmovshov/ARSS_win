%% SI1D_TEST - a test of the PhysX engine forward integration
%
% Two unit spheres were set at rest with their centers 8 cu (code units) apart
% on the x axis. Their separation and relative velocity was saved at each time
% step. Runs were made with succesively smaller time steps, starting at 1 cu.
% The time step is included in the output file's name as well as in the header,
% which also contains the code units and G and rho in cu.

clear
close all
clc

%% Load data from files and do initial conversions
out_files = dir('*.out');
for k = 1:numel(out_files)
    raw = importdata(out_files(k).name,' ',headcount(out_files(k).name));
    rho = cell2mat(textscan(raw.textdata{4},'%*[^=]%*c%f'));
    bigG = cell2mat(textscan(raw.textdata{8},'%*[^=]%*c%f'));
    dt = cell2mat(textscan(raw.textdata{6},'%*[^=]%*c%f'));
    M = 4*pi/3*rho;
    data(k).dt = dt; %#ok<*SAGROW>
    data(k).t = raw.data(:,1);
    data(k).R = raw.data(:,2);
    data(k).V = raw.data(:,3);
    data(k).K = 0.5*M/2*data(k).V.^2; % system kinetic energy
    data(k).U = -bigG*M^2./data(k).R; % system potential energy
    data(k).E = data(k).K + data(k).U; % system total energy
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
plot([0,data(k).t(end)],[-bigG*M^2/8,-bigG*M^2/8],'k--')

%% Touch up and anotate
xlabel('time [cu]')
ylabel('E [cu]')
set(ah, 'box', 'on')
% Show only items with a display name
mask = cellfun(@isempty, get(ah.Children,'DisplayName'));
display_list = (ah.Children(~mask));
lh = legend(ah, display_list, 'location', 'nw');
