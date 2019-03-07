% Inverse Dynamik für komplette Trajektorie für
% S4RRPR1
%
% Eingabe:
% Q
%   Trajektorie von Gelenkpositionen (Lösung der IK)
% QD
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD
%   Trajektorie von Gelenkbeschleunigungen
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR1_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU

function TAU = S7RRRRRRR1_invdyn_mdp_traj(Q, QD, QDD, g, pkin, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,7]),
%$cgargs  coder.newtype('double',[inf,7]),
%$cgargs  coder.newtype('double',[inf,7]),
%$cgargs  zeros(3,1), zeros(4,1), zeros(45,1)}

%% Initialisierung
% Vorbelegung der Ausgabe
TAU = NaN(size(Q));

%% Iterative Berechnung der gesamten Trajektorie
for k = 1:size(Q,1)
  TAU(k,:) = S7RRRRRRR1_invdynJ_fixb_mdp_slag_vp(Q(k,:)', QD(k,:)', QDD(k,:)', g, pkin, MDP);
end