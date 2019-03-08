% Beispiel für die Struktursynthese mit Entwurfsoptimierung der Antriebe
% Zeitmessung für die Ausnutzung der parameterlinearen Form
% 

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-03
% (C) Institut für Mechatronische Systeme, Universität Hannover

clc
clear
warning off % für IK-Warnungen

%% Benutzereingaben
short = false; % true: Nur wenige Wiederholungen, damit das Bild schnell gezeichnet wird
K_ra = 5; % Faktor für Änderung der Dynamikparameter mit Zufallszahlen. Auf Null setzen für Test der Dynamik-Berechnung

Rob_Names = {'S7RRRRRRR1', 'S6RRRRRR10', 'S4RRPR1'};

%% Schleife über alle Roboter
for rr = 1:length(Rob_Names)
  %% Initialisierung

  % Typ des seriellen Roboters auswählen
  SName=Rob_Names{rr};

  % Instanz der Roboterklasse erstellen
  RS = serroblib_create_robot_class(SName);
  RS.DynPar.mode = 4; % Nehme Minimalparameterform für Dynamik
  RS.fill_fcn_handles(true, true); % nehme mex-Funktionen für die IK

  RS.gen_testsettings(true, true);

  %% Beispiel-Trajektorie
  % Keine Trajektorie berechnen. Einfach Zufallswerte nehmen
  Q = rand(1e5, RS.NQJ);
  QD = rand(1e5, RS.NQJ);
  QDD = rand(1e5, RS.NQJ);
  %% Anfangs-Parameter generieren
  rSges = rand(RS.NQJ+1,3);
  mges = rand(RS.NQJ+1,1);
  Icges =  rand(RS.NQJ+1,6);
  RS.update_dynpar1(mges, rSges, Icges);
  MPV_rob = RS.DynPar.mpv;
  RS.gravity = rand(3,1);

  %% Zeit für Dynamikberechnung messen

  if short
    n_statwdh = 2; % Wiederholungen für Statistik
    nid_val = [10, 25, 50];
    nt_val = [100, 200];
  else
    n_statwdh = 25; % Wiederholungen für Statistik
    nid_val = [0:50, 75, 100, 150, 200]; % , 500, 1000
    nt_val = [5, 10, 50, 100, 200, 500, 750, 1000]; % , 5000, 10000
  end
  Ts3m_ges = NaN(length(nt_val), n_statwdh);
  T2m_ges = NaN(length(nid_val), n_statwdh, length(nt_val));
  T4m_ges = NaN(length(nid_val), n_statwdh, length(nt_val));

  for iis = 0:n_statwdh % Schleife über Wiederholungen für Statistik.
    % Starte mit nicht-gezählter Iteration Null, damit
    % Initialisierungsprozesse nicht die Zeitmessung verfälschen.
    fprintf('%s: Starte Iteration %d/%d: %d Traj.-längen, %d ID-Iterationen.\n', ...
      SName, iis, n_statwdh, length(nt_val), length(nid_val));
    for ii_nt = 1:length(nt_val) % Schleife für unterschiedlich lange Traj.
      t0 = tic();
      % Anzahl-Trajpunkte und Begrenzung der Trajektorie
      nt = nt_val(ii_nt);
      Q_ii = Q(1:nt,:);
      QD_ii = QD(1:nt,:);
      QDD_ii = QDD(1:nt,:);

      %% Berechne Regressor-Matrix für ganze Trajektorie
      % Regressormatrix-Vektoren mit Funktion erzeugen
      % ... mit kompilierter Funktion
      tic();
      RMV_ges = RS.invdynregmat_traj(Q_ii, QD_ii, QDD_ii);
      if iis > 0
        Ts3m_ges(ii_nt,iis) = toc();   
      end
      % Iterationen der inversen Dynamik simulieren
      for ii_id = 1:length(nid_val) % Schleife über verschiedene Anzahl Dyn.-Iterationen
        nd = nid_val(ii_id); % Anzahl Dynamik-Iterationen
        %% Schleife für Maßsynthese
        % Simuliere die Optimierung mit Beiden ID-Verfahren
        TAU2 = NaN(nt, RS.NQJ);
        TAU3 = NaN(nt, RS.NQJ);
        TAU4 = NaN(nt, RS.NQJ);

        %% Berechnung mit Inverse-Dynamik Funktionsaufruf
        t2m = tic();
        for di = 1:nd
          MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
          RS.update_dynpar_mpv(MPV_new);
          TAU2 = RS.invdyn2_traj(Q_ii, QD_ii, QDD_ii);
        end
        if iis > 0
          T2m_ges(ii_id,iis, ii_nt) = toc(t2m);
        end
        %% Berechnung mit Regressormatrix-Vektor als Funktion
        % Berechnung mit Funktionsaufruf (mex)
        t4m = tic();
        for di = 1:nd
          MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
          RS.update_dynpar_mpv(MPV_new);
          TAU4 = RS.invdyn3_traj(RMV_ges);
        end
        if iis > 0
          T4m_ges(ii_id,iis, ii_nt) = toc(t4m);
        end
      end
      fprintf('\tTraj.-längen Iteration %d/%d beendet (%d Bahnpunkte). Dauer: %1.1fs\n', ...
        ii_nt, length(nt_val), nt_val(ii_nt), toc(t0));
    end
    %% Inhaltliche Prüfung der Ergebnisse
    if iis == 1 && K_ra == 0
      % Prüfe, ob Ergebnisse übereinstimmen
      test = TAU2 - TAU4;
      if max(abs(test(:))) > 1e-10
        error('Inverse Dynamik stimmt nicht zwischen Berechnungsmethoden');
      end
    end
  end
  %% Netto-Zeiten berechnen (Dauer für Regressor-Berechnung mit einbeziehen)
  % Nehme die beste Regressor-Implementierung Ts3m_ges
  T4m_gesn = NaN(size(T4m_ges));
  for i = 1:length(nt_val)
    % Egal wie viele Dynamik-Iterationen gemacht werden (Dimension 2 von
    % T3m_ges), die Zeit für die Berechnung des Regressors ist immer gleich
    T4m_gesn(:,:,i) = T4m_ges(:,:,i) + repmat(Ts3m_ges(i,:), length(nid_val), 1);
  end

  %% Vergleich der Zeitauswertung
  figure(rr);clf;
  for i = 1:length(nt_val)
    n_op = nt_val(i) * nid_val';
    % Gesamtdauer der Trajektorie
    subplot(length(nt_val), 2, sprc2no(length(nt_val), 2, i,1)); hold on; grid on;
    % Mittelwert kennzeichnen
    plot(nid_val, 1e3*mean(T2m_ges(:,:, i),2), 'r-o');
    plot(nid_val, 1e3*mean(T4m_gesn(:,:, i),2), 'b-v');
    % Spannbreite der Daten kennzeichnen
    for k = 1:length(nid_val)
      plot([1;1]*nid_val(k), 1e3*minmax(T2m_ges(k, :, i))', 'r-');
      plot([1;1]*nid_val(k), 1e3*minmax(T4m_gesn(k, :, i))', 'b--');
    end

    title(sprintf('nt val = %d', nt_val(i)));
    ylabel('Zeit (Traj.) in ms');
    legend({'InvDyn-Funktion', 'Regressor-Vektor'});
    if i == length(nt_val)
      xlabel('Anzahl Dynamik-Iterationen');
    end

    % Dauer einer einzelnen Dynamik-Auswertung
    subplot(length(nt_val), 2, sprc2no(length(nt_val), 2, i,2)); hold on; grid on;
    % Mittelwert
    plot(nid_val, 1e6*mean(T2m_ges(:,:, i),2)./n_op, 'r-o');
    plot(nid_val, 1e6*mean(T4m_gesn(:,:, i),2)./n_op, 'b-x');

    % Spannbreite der Daten kennzeichnen
    for k = 1:length(nid_val)
      plot([1;1]*nid_val(k), 1e6*minmax(T2m_ges(k, :, i)/n_op(k))', 'r-');
      plot([1;1]*nid_val(k), 1e6*minmax(T4m_gesn(k, :, i)/n_op(k))', 'b-');
    end
    title(sprintf('nt val = %d', nt_val(i)));
    ylabel('Zeit (pro Zeitschritt) in us');
    legend({'InvDyn-Funktion', 'Regressor-Vektor'});
  end
  linkxaxes
  set(rr, 'Name', sprintf('%s_Zeiten_Vgl', SName), 'NumberTitle', 'off');
  
  % Prozentuale Verringerung
  figure(10+rr);clf;hold on;
  legtxts = cell(length(nt_val),1);
  for i = 1:length(nt_val)
    % Gesamtdauer der Trajektorie
    % Mittelwert kennzeichnen
    plot(nid_val, 100*mean(T4m_gesn(:,:, i),2)./mean(T2m_ges(:,:, i),2));
    legtxts{i} = sprintf('nt=%d', nt_val(i));
  end
  xlabel('Anzahl Dynamik-Iterationen');
  ylabel('T_{Reg} / T_{InvDyn} in %');
  legend(legtxts);
  grid on;
  title(sprintf('%s: Zeitverhältnis der Berechnungsarten', SName));
  set(10+rr, 'Name', sprintf('%s_Zeitverhältnis', SName), 'NumberTitle', 'off');
  %% Ergebnisse speichern
  resdir = fileparts(which('dimsynth_timing_example.m'));
  resfile = fullfile(resdir, sprintf('%s_timing_example_s%d.mat', SName, short));
  save(resfile, 'Ts3m_ges', 'T2m_ges', 'T4m_ges', 'T4m_gesn', 'nt_val', 'nid_val');
  saveas(rr,    fullfile(resdir, sprintf('%s_time_comparison_s%d.fig', SName, short)));
  saveas(10+rr, fullfile(resdir, sprintf('%s_time_ratio_s%d.fig', SName, short)));
end