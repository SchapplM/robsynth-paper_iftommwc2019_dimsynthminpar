% Erstelle Pseudo-Datenreihen für Veranschaulichungsbild der Optimierung
% (Vergleich mehrerer Roboterstrukturen)

clc
clear

% Zufallsgenerator zurücksetzen
rng default

% Exportiere Bilder in den Ordner dieses Skriptes
this_path = fileparts( mfilename('fullpath') );
if isempty(this_path)
  exportpath = pwd;
else
  exportpath = this_path;
end
%% Init
n = 20;
gen = (1:n)';
data = NaN(n,4);

data(:,1) = 5+ 5*exp(-gen / 10);
data(:,1) = data(:,1) + rand(n,1);

data(:,2) = 3+ 18*exp(-gen / 6);
data(:,2) = data(:,2) + 2*rand(n,1).*exp(-gen / 5);

data(:,3) = 8+ 5*exp(-gen / 15);
data(:,3) = data(:,3) + 10*rand(n,1).*exp(-gen / 10);

data(:,4) = 10+ 3*exp(-gen / 7);
data(:,4) = data(:,4) + 6*rand(n,1).*exp(-gen / 3);
data(end,:) = NaN;
%% Bild
figHandle = figure(1);clf;
hold on;
plot(gen, data(:,1), 'bo');
plot(gen, data(:,2), 'gx');
plot(gen, data(:,3), 'rs');
plot(gen, data(:,4), 'kv');

axhdl = gca;
for i = 1:length(axhdl(:))
  % Rand vollständig zeichnen
  set(axhdl(i), 'Box', 'on')
end

% Schriftgröße und -art setzen
set_font_fontsize(gcf,'Times',12)

% Weißer Hintergrund
set(gcf,'color','w');

xlabel('Iterations');
ylabel('Fitness Fcn');
set(gca, 'xticklabel', {});
set(gca, 'yticklabel', {});
ylim([3,20])

set_size_plot_subplot(figHandle, ...
  5, 3.5, ...
  axhdl, ...
  0.14, 0.01, 0.03, 0.15, ... % l r u d
  0.14, 0.1) % x y

Filebasename_res = 'Vgl_Optimierung';
export_fig(fullfile(exportpath, [Filebasename_res, '.pdf']));