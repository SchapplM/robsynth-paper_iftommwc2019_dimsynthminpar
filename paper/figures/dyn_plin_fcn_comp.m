% Bild: Vergleich der Rechenzeit für InvDyn-Direktaufruf mit
% Regressor-Multiplikation
% 
% Vorher ausführen: dimsynth_timing_example.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-03
% (C) Institut für Mechatronische Systeme, Universität Hannover

clc
clear
this_path = fileparts(which('dyn_plin_fcn_comp.m'));
Rob_Names = {'S4RRPR1', 'S6RRRRRR10', 'S7RRRRRRR1'};

roblegnames = {'SCARA', 'Ind.Rob.', 'LWR'};
methodnames = {'InvDyn', 'RegMult'};

% Schalter für zwei Einzelne Bilder oder ein kombiniertes Bild
combinedfig = true;
%% Bild: Vergleich der Rechenzeit über die Anzahl der Dyn.-Eval.
axhdl = NaN(1,2);

if combinedfig
  figure(3);clf;
  axhdl(1)=subplot(1,2,1);hold on;
else
  figure(1);clf;hold on;
end
nt = 1000;
linhdl = NaN(6,1);
legnames = {};
markers = {'o', 's', 'v'};
for rr = 1:length(Rob_Names)
  Name = Rob_Names{rr};
  resfile = fullfile(this_path, '..', '..', 'calc', sprintf('%s_timing_example_s0.mat',Name));
  res = load(resfile);
  i = find(res.nt_val == nt);
  set(gca, 'ColorOrderIndex', rr);
  linhdl(2*rr-1)=plot(res.nid_val, 1e3*mean(res.T2m_ges(:,:, i),2)./res.nid_val', sprintf('%s--',markers{rr}));
  set(gca, 'ColorOrderIndex', rr);
  linhdl(2*rr)=plot(res.nid_val, 1e3*mean(res.T4m_gesn(:,:, i),2)./res.nid_val', sprintf('%s-',markers{rr}));
  for kk = 1:2
    legnames = {legnames{:}, sprintf('%s: %s', roblegnames{rr}, methodnames{kk})}; %#ok<CCAT>
  end
end
xlabel('# of dynamics iterations', 'interpreter', 'latex');
ylabel(sprintf('Time for $N_\\mathrm{T}{=}%d$ dyn. eval. in ms',nt), 'interpreter', 'latex');
xlim([0, 30])
ylim([0,2])
grid on;

set(1, 'WindowStyle', 'normal');
format = {'r', 'o', '--', 7; ...
          'r', 'o', '-', 9; ...
          'k', 's', '--', 5; ...
          'k', 's', '-', 12; ...
          'b', 'v',  '--', 11; ...
          'b', 'v', '-', 13};
leghdl = line_format_publication(linhdl, format, legnames);
legend(leghdl, legnames);

if ~combinedfig
  figure_format_publication(gca);
  set_size_plot_subplot(1, ...
    12.2, 6, ... 
    gca, ...
    0.08, 0.01, 0.02, 0.18, ... % l, r, u, d
    0.1, 0.08); % dx, dy
  export_fig(1, fullfile(this_path, 'time_cmp_per_traj.fig'));
  export_fig(1, fullfile(this_path, 'time_cmp_per_traj.pdf'));
  fprintf('Bild "time_cmp_per_traj" generiert\n');
end

%% Bild: Vergleich des Zeitverhältnisses bei verschiedener Traj.-Länge
Name = Rob_Names{3};
resfile = fullfile(this_path, '..', '..', 'calc', sprintf('%s_timing_example_s0.mat',Name));
res = load(resfile);
legnames = {};
if combinedfig
  axhdl(2) = subplot(1,2,2);hold on;
else
  figure(2);clf;hold on;
end
nt_val_sel = [10, 50, 500, 1000];
linhdl = NaN(3,1);
k=0;
for i = 1:length(res.nt_val)
  if ~any(res.nt_val(i) == nt_val_sel)
    continue
  end
  k=k+1;
  linhdl(k)=plot(res.nid_val, 100*mean(res.T4m_gesn(:,:, i),2)./mean(res.T2m_ges(:,:, i),2));
  legnames = {legnames{:}, sprintf('$N_\\mathrm{T}$=%d', res.nt_val(i))}; %#ok<CCAT>
end
xlabel('# of dynamics iterations', 'interpreter', 'latex');

ylabel('$T_\mathrm{RegMult}$ / $T_\mathrm{InvDyn}$ in \% ', 'interpreter', 'latex');
legend(legnames, 'interpreter', 'LaTex');
xlim([1, 48]);
ylim([0, 150]);
yticks(20:20:140)
grid on
set(2, 'WindowStyle', 'normal');
format = {[0 80 155 ]/255, 'd', '-', 12; ...
          [231 123 41 ]/255, 'x', '-', 7; ...
          [200 211 23 ]/255, '+', '-', 5; ...
          'k', '', '--', 0};

leghdl = line_format_publication(linhdl, format, legnames);
legend(leghdl, legnames);
if ~combinedfig
  figure_format_publication(gca);
  set_size_plot_subplot(2, ...
    12.2, 6, ... 
    gca, ...
    0.09, 0.01, 0.02, 0.14, ... % l, r, u, d
    0.1, 0.08); % dx, dy
  export_fig(2, fullfile(this_path, sprintf('%s_time_ratio.fig', Name)));
  export_fig(2, fullfile(this_path, sprintf('%s_time_ratio.pdf', Name)));
  fprintf('Bild "time_ratio" generiert\n');
end

%% Gemeinsames Bild formatieren und speichern
if combinedfig
  figure_format_publication(axhdl);
  set_size_plot_subplot(3, ...
    12.2, 5.5, ...  % IFToMM-Springer: 12.2cm Breite, damit keine overfull hbox
    axhdl, ...
    0.08, 0.01, 0.02, 0.14, ... % l, r, u, d
    0.1, 0.08); % dx, dy
  export_fig(3, fullfile(this_path, 'dyn_timing_eval.fig'));
  export_fig(3, fullfile(this_path, 'dyn_timing_eval.pdf'));
  fprintf('Bild mit Zeitvergleich und Verhältnisvergleich generiert\n');
end