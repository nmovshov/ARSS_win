%% SGP_ENERGY_TEST - a test of the PhysX engine forward integration
%
% A self-gravitating pile is set up in loose configuration and allowed to
% collapse. Measured quantities are system kinetic and potential energy. For
% this test the PhysX params were manually modified to approximate perfect
% bouncing and zero friction.
clear
clc

%% Load data from files and do initial conversions
out_files = dir('*.out');
for k = 1:numel(out_files)
    raw = importdata(out_files(k).name,' ',headcount(out_files(k).name));
    %bigG = cell2mat(textscan(raw.textdata{5},'%*[^=]%*c%f'));
    dt = cell2mat(textscan(raw.textdata{5},'%*[^=]%*c%f'));
    data(k).dt = dt; %#ok<*SAGROW>
    data(k).t = raw.data(:,1);
    data(k).K = raw.data(:,6); % system kinetic energy
    data(k).U = raw.data(:,5); % system potential energy
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
plot(data(k).t,data(k).K,'r--','linewidth',2)
plot(data(k).t,data(k).U,'g--','linewidth',2)
plot([0,data(k).t(end)],[data(k).E(1), data(k).E(1)],'k--')

%% Touch up and anotate
xlabel('time [cu]')
ylabel('E [cu]')
set(ah, 'box', 'on')
ah.XLim = [0, min(arrayfun(@(x)max(x.t),data))];
% Show only items with a display name
mask = cellfun(@isempty, get(ah.Children,'DisplayName'));
display_list = (ah.Children(~mask));
lh = legend(ah, display_list, 'location', 'nw');
