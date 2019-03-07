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
%   see S7RRRRRRR1_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU

function TAU = S7RRRRRRR1_invdyn_mdp_traj3(RMV_ges, NQJ, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,191]), 7, zeros(45,1)}

%% Initialisierung
% Vorbelegung der Ausgabe
TAU = NaN(size(RMV_ges,1), NQJ);

for ii = 1:size(RMV_ges,1)
  TAU(ii,:) = S7RRRRRRR1_invdynJ_fixb_mdp2_slag_vp(RMV_ges(ii,:), MDP);
end