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

function TAU = S7RRRRRRR1_invdyn_mdp_traj2(RM_ges, NQJ, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,45]), 6, zeros(45,1)}

%% Initialisierung
% Vorbelegung der Ausgabe
TAU = NaN(size(RM_ges,1)/NQJ, NQJ);

for ii = 1:size(RM_ges,1)/NQJ
  tr = RM_ges((ii-1)*NQJ+1:ii*NQJ,:);
  TAU(ii,:) = tr * MDP;
end