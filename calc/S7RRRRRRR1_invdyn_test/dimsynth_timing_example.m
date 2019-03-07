% Beispiel für die Struktursynthese mit Entwurfsoptimierung der Antriebe
% Zeitmessung für die Ausnutzung der parameterlinearen Form
% 

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-03
% (C) Institut für Mechatronische Systeme, Universität Hannover

clc
clear
warning off % für IK-Warnungen

%% Benutzereingaben
debug = false; % Alternative Varianten durchrechnen (zeitaufwändiger)
traj = false; % Trajektorie oder Zufallswerte für Gelenkkonfigurationen
short = false; % true: Nur wenige Wiederholungen, damit das Bild schnell gezeichnet wird
K_ra = 5; % Faktor für Änderung der Dynamikparameter mit Zufallszahlen. Auf Null setzen für Test der Dynamik-Berechnung
%% Initialisierung

% Typ des seriellen Roboters auswählen
SName='S7RRRRRRR1';

% Instanz der Roboterklasse erstellen
RS = serroblib_create_robot_class(SName);
RS.I_EE = logical([1 1 1 0 0 1]); % nur 4FG in InvKin betrachten
RS.DynPar.mode = 4; % Nehme Minimalparameterform für Dynamik
RS.fill_fcn_handles(true, true); % nehme mex-Funktionen für die IK
%% Anfangsannahme für Kinematik-Parameter
kin = false;
if kin
  RName='S7RRRRRRR1_KUKA1';
  RS_tmp = serroblib_create_robot_class(SName, RName);
  pkin_start = RS_tmp.pkin;
  RS.update_mdh(pkin_start);
  RS.qref = RS_tmp.qref;
  RS.qlim = RS_tmp.qlim;
  RS.qDlim = RS_tmp.qDlim;
else
  RS.gen_testsettings(true, true);
end
%% Beispiel-Trajektorie

if traj
  % Grundstellung
  q0 = RS.qref+rand(RS.NQJ,1)*1e-2;
  T_E = RS.fkineEE(q0);
  x0 = RS.t2x(T_E);

  % Start in Grundstellung
  k=1; XE = x0';
  % Fahrt in Startstellung
  k=k+1; XE(k,:) = XE(k-1,:) + [ -0.1,0.0,0, 0,0,0];
  % Beginn Quadrat
  d1=0.1;
  % k=k+1; XE(k,:) = XE(k-1,:) + [0,-d1/2,0, 0,0,0];
  % k=k+1; XE(k,:) = XE(k-1,:) + [-d1,0,0  0,0,0];
  % k=k+1; XE(k,:) = XE(k-1,:) + [0,d1,0, 0,0,0];
  % k=k+1; XE(k,:) = XE(k-1,:) + [d1,0,0  0,0,0];
  % % TODO: Sinnvoll weitermachen

  [X,XD,XDD,T] = traj_trapez2_multipoint(XE, 1.0, 100e-3, 10e-3, 1e-3, 0.100);

  figure(2);clf;
  s_plot = struct( 'ks', [1:RS.NJ, RS.NJ+2], 'straight', 0);
  hold on;
  grid on;
  xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');view(3);
  RS.plot( q0, s_plot );
  title(sprintf('Grundstellung: %s', SName));
  plot3(X(:,1), X(:,2), X(:,3));

  % Inverse Kinematik
  [Q,QD,QDD] = RS.invkin2_traj(X,XD,XDD,T,q0);


  % IK-Fehler
  Phi_ges = NaN(length(T), 6);
  for i = 1:length(T)
    Phi_i = RS.constr2(Q(i,:)', X(i,:)');
    Phi_ges(i,:) = Phi_i;
    if max(abs(Phi_i)) > 1e-3
      % Eigentlich ist für Dynamik auch vollkommen egal, ob die IK
      % konvergiert ist ...
  %     error('IK stimmt nicht');
    end
  end

  figure(8);clf;
  for k = 1:6
      subplot(4,6,sprc2no(4,6,1,k));hold on;
      plot(T, X(:,k));
      grid on;
      if k < 4, ylabel(sprintf('%s_E in m', char(119+k)));
      else,     ylabel(sprintf('\\phi_{E,1} in rad')); end
      xlabel('t / s');
      subplot(4,6,sprc2no(4,6,2,k));hold on;
      plot(T, XD(:,k));
      xlabel('t / s');
      if k < 4, ylabel(sprintf('D%s_E in m/s', char(119+k)));
      else,     ylabel(sprintf('D\\phi_{E,1} in rad/s')); end
      grid on;
      subplot(4,6,sprc2no(4,6,3,k));hold on;
      plot(T, XDD(:,k));
      xlabel('t / s');
      if k < 4, ylabel(sprintf('DD%s_E in m/s²', char(119+k)));
      else,     ylabel(sprintf('DD\\phi_{E,1} in rad/s²')); end
      grid on;
      subplot(4,6,sprc2no(4,6,4,k));hold on;
      plot(T, Phi_ges(:,k));
      if k < 4, ylabel(sprintf('\\Delta %s in m', char(119+k)));
      else,     ylabel(sprintf('\\Delta \\phi_{1} in rad')); end
      grid on;
  end
  figure(9);clf;
  for k = 1:RS.NQJ
    subplot(RS.NQJ,3,sprc2no(RS.NQJ,3,k,1));hold on;
    plot(T, Q(:,k));
    ylabel(sprintf('q_%d',k));
    subplot(RS.NQJ,3,sprc2no(RS.NQJ,3,k,2));hold on;
    plot(T, QD(:,k));
    ylabel(sprintf('qD_%d',k));
    subplot(RS.NQJ,3,sprc2no(RS.NQJ,3,k,3));hold on;
    plot(T, QDD(:,k));
    ylabel(sprintf('qDD_%d',k));
  end


  % Trajektorie vervielfachen
  for k = 1:6
    Q = [Q; Q];
    QD = [QD; QD];
    QDD = [QDD; QDD];
  end
else
  % Keine Trajektorie berechnen. Einfach Zufallswerte nehmen
  Q = rand(1e5, RS.NQJ);
  QD = rand(1e5, RS.NQJ);
  QDD = rand(1e5, RS.NQJ);
end
%% Anfangs-Parameter generieren
rSges = rand(RS.NQJ+1,3);
mges = rand(RS.NQJ+1,1);
Icges =  rand(RS.NQJ+1,6);
RS.update_dynpar1(mges, rSges, Icges);
MPV_rob = RS.DynPar.mpv;
RS.gravity = rand(3,1);

%% Init

matlabfcn2mex({'S7RRRRRRR1_invdyn_mdp_traj'});
matlabfcn2mex({'S7RRRRRRR1_invdyn_mdp_traj2'});
matlabfcn2mex({'S7RRRRRRR1_invdyn_mdp_traj3'});
matlabfcn2mex({'S7RRRRRRR1_invdynJ_fixb_regmin_traj'});

pkin = RS.pkin;
g = RS.gravity;

%% Zeit für Dynamikberechnung messen

if short
  n_statwdh = 2; % Wiederholungen für Statistik
  nid_val = [10, 25, 50];
  nt_val = [100, 200];
else
  n_statwdh = 5; % Wiederholungen für Statistik
  nid_val = [0:50, 75, 100, 150, 200, 500, 1000];
  nt_val = [200, 500, 750, 1000, 5000, 10000];
end
Ts1_ges = NaN(length(nt_val), n_statwdh); % Regressor in Matlab-Skript
Ts2_ges = Ts1_ges; % Regressormatrix als Vektor in Matlab
Ts3_ges = Ts1_ges; Ts3m_ges = Ts3_ges;% Matrix-Vektor als Funktion
T1_ges = NaN(length(nid_val), n_statwdh, length(nt_val)); 
T2_ges = Ts1_ges; T2m_ges = Ts1_ges;
T3_ges = Ts1_ges; T3m_ges = Ts1_ges;
T4_ges = Ts1_ges; T4m_ges = Ts1_ges;


for iis = 1:n_statwdh % Schleife über Wiederholungen für Statistik
  for ii_nt = 1:length(nt_val) % Schleife für unterschiedlich lange Traj.
    % Anzahl-Trajpunkte und Begrenzung der Trajektorie
    nt = nt_val(ii_nt);
    Q_ii = Q(1:nt,:);
    QD_ii = QD(1:nt,:);
    QDD_ii = QDD(1:nt,:);

    %% Berechne Regressor-Matrix für ganze Trajektorie
    % Regressor-Matrizen stapeln
    t0 = tic;
    RM_ges = NaN(nt*RS.NQJ, length(MPV_rob) );
    for ii = 1:nt
      q = Q_ii(ii,:)';
      qD = QD_ii(ii,:)';
      qDD = QDD_ii(ii,:)';
      [t,tr] = RS.invdyn(q, qD, qDD);
      RM_ges((ii-1)*RS.NQJ+1:ii*RS.NQJ,:) = tr;
    end
    Ts1_ges(ii_nt,iis) = toc(t0);
  
    % Regressor-Matrix als Vektor aufbauen und diese Vektoren stapeln
    t0 = tic();
    rmv_dummy = S7RRRRRRR1_regmat2regmatvector( rand(RS.NQJ, length(RS.DynPar.mpv)) );
    RMV_ges = NaN(nt, length(rmv_dummy));
    for ii = 1:nt
      q = Q_ii(ii,:)';
      qD = QD_ii(ii,:)';
      qDD = QDD_ii(ii,:)';
      [t,tr] = RS.invdyn(q, qD, qDD);
      RMV_ges(ii,:) = S7RRRRRRR1_regmat2regmatvector(tr);
    end
    Ts2_ges(ii_nt,iis) = toc(t0);
    
    % Regressormatrix-Vektoren mit Funktion erzeugen
    tic();
    RMV_ges2 = S7RRRRRRR1_invdynJ_fixb_regmin_traj(Q_ii, QD_ii, QDD_ii, g, pkin);
    Ts3_ges(ii_nt,iis) = toc();
    
    % ... mit kompilierter Funktion
    tic();
    RMV_ges3 = S7RRRRRRR1_invdynJ_fixb_regmin_traj_mex(Q_ii, QD_ii, QDD_ii, g, pkin);
    Ts3m_ges(ii_nt,iis) = toc();   

    % Iterationen der inversen Dynamik simulieren
    for ii_id = 1:length(nid_val) % Schleife über verschiedene Anzahl Dyn.-Iterationen
      nd = nid_val(ii_id); % Anzahl Dynamik-Iterationen


      %% Schleife für Maßsynthese
      % Simuliere die Optimierung mit Beiden ID-Verfahren
      TAU1 = NaN(nt, RS.NQJ);
      TAU2 = NaN(nt, RS.NQJ);
      TAU3 = NaN(nt, RS.NQJ);
      TAU4 = NaN(nt, RS.NQJ);

      %% Berechnung mit Regressor-Matrix (Matlab Skript)
      if debug
        t1 = tic;
        for di = 1:nd
          MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
          for ii = 1:nt
            tr = RM_ges((ii-1)*RS.NQJ+1:ii*RS.NQJ,:);
            TAU1(ii,:) = tr * MPV_new;
          end
        end
        T1_ges(ii_id,iis, ii_nt) = toc(t1);
      end
      %% Berechnung mit Inverse-Dynamik Funktionsaufruf
      if debug
        t2 = tic;
        for di = 1:nd
          MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
          TAU2 = S7RRRRRRR1_invdyn_mdp_traj(Q_ii, QD_ii, QDD_ii, g, pkin, MPV_new);
        end
        T2_ges(ii_id,iis, ii_nt) = toc(t2);
      end

      t2m = tic;
      for di = 1:nd
        MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
        TAU2 = S7RRRRRRR1_invdyn_mdp_traj_mex(Q_ii, QD_ii, QDD_ii, g, pkin, MPV_new);
      end
      T2m_ges(ii_id,iis, ii_nt) = toc(t2m);

      %% Berechnung mit Regressor-Matrix, Funktion
      if debug
        t0 = tic;
        for di = 1:nd
          MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
          TAU3 = S7RRRRRRR1_invdyn_mdp_traj2(RM_ges, RS.NQJ, MPV_new);
        end
        T3_ges(ii_id,iis, ii_nt) = toc(t0);
      end

      t0 = tic;
      for di = 1:nd
        MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
        TAU3 = S7RRRRRRR1_invdyn_mdp_traj2_mex(RM_ges, RS.NQJ, MPV_new);
      end
      T3m_ges(ii_id,iis, ii_nt) = toc(t0);
      
      %% Berechnung mit Regressormatrix-Vektor als Funktion
      if debug
        t0 = tic;
        for di = 1:nd
          MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
          TAU4 = S7RRRRRRR1_invdyn_mdp_traj3(RMV_ges, RS.NQJ, MPV_new);
        end
        T4_ges(ii_id,iis, ii_nt) = toc(t0);
      end

      % Berechnung mit Funktionsaufruf (mex)
      t4m = tic;
      for di = 1:nd
        MPV_new = MPV_rob + K_ra*rand(length(MPV_rob), 1);
        TAU4 = S7RRRRRRR1_invdyn_mdp_traj3_mex(RMV_ges, RS.NQJ, MPV_new);
      end
      T4m_ges(ii_id,iis, ii_nt) = toc(t4m);
    end
  end
  %% Inhaltliche Prüfung der Ergebnisse
  if iis == 1 && K_ra == 0
    % Prüfe, ob Ergebnisse übereinstimmen
    test2 = TAU1 - TAU2;
    test3 = TAU1 - TAU3;
    test4 = TAU1 - TAU4;
    if max(abs([test2(:); test3(:); test4(:)])) > 1e-10
      error('Inverse Dynamik stimmt nicht zwischen Berechnungsmethoden');
    end
  end
end
%% Netto-Zeiten berechnen (Dauer für Regressor-Berechnung mit einbeziehen)
% Nehme die beste Regressor-Implementierung Ts3m_ges
T3m_gesn = NaN(size(T3m_ges));
T4m_gesn = NaN(size(T4m_ges));
for i = 1:length(nt_val)
  % Egal wie viele Dynamik-Iterationen gemacht werden (Dimension 2 von
  % T3m_ges), die Zeit für die Berechnung des Regressors ist immer gleich
  T3m_gesn(:,:,i) = T3m_ges(:,:,i) + repmat(Ts3m_ges(i,:), length(nid_val), 1);
  T4m_gesn(:,:,i) = T4m_ges(:,:,i) + repmat(Ts3m_ges(i,:), length(nid_val), 1);
end

%% Vergleich der Zeitauswertung
figure(1);clf;
for i = 1:length(nt_val)
  n_op = nt_val(i) * nid_val';
  % Gesamtdauer der Trajektorie
  subplot(length(nt_val), 2, sprc2no(length(nt_val), 2, i,1)); hold on; grid on;
  % Mittelwert kennzeichnen
  plot(nid_val, 1e3*mean(T2m_ges(:,:, i),2), 'r-o');
  plot(nid_val, 1e3*mean(T3m_gesn(:,:, i),2), 'g-x');
  plot(nid_val, 1e3*mean(T4m_gesn(:,:, i),2), 'b-v');
  % Spannbreite der Daten kennzeichnen
  for k = 1:length(nid_val)
    plot([1;1]*nid_val(k), 1e3*minmax(T2m_ges(k, :, i))', 'r-');
    plot([1;1]*nid_val(k), 1e3*minmax(T3m_gesn(k, :, i))', 'g--');
    plot([1;1]*nid_val(k), 1e3*minmax(T4m_gesn(k, :, i))', 'b--');
  end

  title(sprintf('nt val = %d', nt_val(i)));
  ylabel('Zeit (Traj.) in ms');
  legend({'InvDyn-Funktion', 'Regressor-Matrix', 'Regressor-Vektor'});
  if i == length(nt_val)
    xlabel('Anzahl Dynamik-Iterationen');
  end
  
  % Dauer einer einzelnen Dynamik-Auswertung
  subplot(length(nt_val), 2, sprc2no(length(nt_val), 2, i,2)); hold on; grid on;
  % Mittelwert
  plot(nid_val, 1e6*mean(T2m_ges(:,:, i),2)./n_op, 'r-o');
  plot(nid_val, 1e6*mean(T3m_gesn(:,:, i),2)./n_op, 'g-x');
  plot(nid_val, 1e6*mean(T4m_gesn(:,:, i),2)./n_op, 'b-x');
  
  % Spannbreite der Daten kennzeichnen
  for k = 1:length(nid_val)
    plot([1;1]*nid_val(k), 1e6*minmax(T2m_ges(k, :, i)/n_op(k))', 'r-');
    plot([1;1]*nid_val(k), 1e6*minmax(T3m_gesn(k, :, i)/n_op(k))', 'g-');
    plot([1;1]*nid_val(k), 1e6*minmax(T4m_gesn(k, :, i)/n_op(k))', 'b-');
  end
  title(sprintf('nt val = %d', nt_val(i)));
  ylabel('Zeit (pro Zeitschritt) in us');
  legend({'InvDyn-Funktion', 'Regressor-Matrix', 'Regressor-Vektor'});
end
linkxaxes

% Vergleche Zeit für Regressor-Berechnung
if debug
  figure(2);clf;
  subplot(2,1,1);hold on;
  plot(nt_val,mean(Ts1_ges, 2), 's-')
  plot(nt_val,mean(Ts2_ges, 2), 'v-')
  plot(nt_val,mean(Ts3_ges, 2), 'o-')
  plot(nt_val,mean(Ts3m_ges,2), 'x-')
  legend({'Klasse Matrix', 'Klasse Vektor', 'Funktion', 'Funktion mex'});
  xlabel('Anzahl Traj.-Punkte');
  ylabel('Zeit für Traj. [s]');
  grid on;
  subplot(2,1,2);hold on;
  % Mittelwert
  plot(nt_val,1e3*mean(Ts1_ges, 2)./nt_val', 's-r')
  plot(nt_val,1e3*mean(Ts2_ges, 2)./nt_val', 'v-g')
  plot(nt_val,1e3*mean(Ts3_ges, 2)./nt_val', 'o-b')
  plot(nt_val,1e3*mean(Ts3m_ges,2)./nt_val', 'x-m')
  % Spannbreite
  for k = 1:length(nt_val)
    plot(nt_val(k)*[1;1],1e3*minmax(Ts1_ges (k,:))/nt_val(k), '-r')
    plot(nt_val(k)*[1;1],1e3*minmax(Ts2_ges (k,:))/nt_val(k), '-g')
    plot(nt_val(k)*[1;1],1e3*minmax(Ts3_ges (k,:))/nt_val(k), '-b')
    plot(nt_val(k)*[1;1],1e3*minmax(Ts3m_ges(k,:))/nt_val(k), '-m')
  end
  legend({'Klasse Matrix', 'Klasse Vektor', 'Funktion', 'Funktion mex'});
  xlabel('Anzahl Traj.-Punkte');
  ylabel('Zeit für Traj.-Pkt. [ms]');
  grid on;
end

% Prozentuale Verringerung
figure(3);clf;hold on;
legtxts = cell(length(nt_val),1);
for i = 1:length(nt_val)
  n_op = nt_val(i) * nid_val';
  % Gesamtdauer der Trajektorie
  % Mittelwert kennzeichnen
  plot(nid_val, 100*mean(T4m_gesn(:,:, i),2)./mean(T2m_ges(:,:, i),2));
  legtxts{i} = sprintf('nt=%d', nt_val(i));
end
xlabel('Anzahl Dynamik-Iterationen');
ylabel('T_{Reg} / T_{InvDyn} in %');
legend(legtxts);
grid on;

%% Ergebnisse speichern
resdir = fileparts(which('dimsynth_timing_example.m'));
resfile = fullfile(resdir, sprintf('timing_example_d%d_s%d.mat', debug, short));
Ts1_ges = NaN(length(nt_val), n_statwdh); % Regressor in Matlab-Skript
Ts2_ges = Ts1_ges; % Regressormatrix als Vektor in Matlab
Ts3_ges = Ts1_ges; Ts3m_ges = Ts3_ges;% Matrix-Vektor als Funktion
T1_ges = NaN(length(nid_val), n_statwdh, length(nt_val)); 
T2_ges = Ts1_ges; T2m_ges = Ts1_ges;
T3_ges = Ts1_ges; T3m_ges = Ts1_ges;
T4_ges = Ts1_ges; T4m_ges = Ts1_ges;
save(resfile, 'Ts1_ges', 'Ts2_ges', 'Ts3_ges', 'Ts3m_ges', ...
  'T1_ges', 'T2_ges', 'T2m_ges', 'T3_ges', 'T3m_ges', 'T4_ges', 'T4m_ges');