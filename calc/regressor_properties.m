% Teste Eigenschaften der Regressormatrix für ausgewählte Roboter

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-11
% (C) Institut für mechatronische Systeme, Universität Hannover

clear
clc

%% Init 1
addpath('/home/schappler/IMES/REPO/maple/codeexport/S6RRRRRR10/testfcn');
addpath('/home/schappler/IMES/REPO/maple/codeexport/S7RRRRRRR1/testfcn');
addpath('/home/schappler/IMES/REPO/maple/codeexport/kuka6dof/matlabfcn');
addpath('/home/schappler/IMES/REPO/maple/codeexport/kuka6dof/testfcn');
%% Benutzereingabe

% robot_name = 'S4RRPR1';
% robot_name = 'S6RRRRRR10';
% robot_name = 'S7RRRRRRR1';
robot_name = 'kuka6dof';

%% Init 2
serroblib_addtopath({robot_name});
eval(sprintf('TSS = %s_varpar_testfunctions_parameter()', robot_name));
pkin = TSS.pkin;

%% Teste Strukturelle Eigenschaften
  
eval(sprintf('func_regmin = @%s_invdynJ_fixb_regmin_slag_vp;', robot_name));
eval(sprintf('func_reg2 = @%s_invdynJ_fixb_reg2_slag_vp;', robot_name));
NQJ = TSS.NQJ;
NPV2 = TSS.NJ*10;
tmp=func_regmin(zeros(NQJ,1),zeros(NQJ,1),zeros(NQJ,1),zeros(3,1),pkin);
NMPV = size(tmp,2);

tau_regmin_mask = zeros(NQJ,NMPV);
tau_reg2_mask = zeros(NQJ,NPV2);

for i = 1:size(TSS.Q,1)
  q = TSS.Q(i,:)';
  qD = TSS.QD(i,:)';
  qDD = TSS.QDD(i,:)';
  g_base = TSS.G(i,:)';

  
  tau_regmin = func_regmin(q, qD, qDD, g_base, pkin);
  tau_reg2 = func_reg2(q, qD, qDD, g_base, pkin);
  
  tau_regmin_mask = tau_regmin_mask + abs(tau_regmin);
  tau_reg2_mask = tau_reg2_mask + abs(tau_reg2);
end
I_reg2 = tau_reg2_mask ~= 0;
I_regmin = tau_regmin_mask ~= 0;
% Alles links von einer "eins" in dem Index-Vektor muss auch eine Eins
% sein (strukturell). Wenn dort eine "Null" steht, liegt es an der
% Besonderheit des Systems
for i = 1:size(I_reg2,1)
  II_set_i = find(I_reg2(i,:));
  first1 = II_set_i(1);
  I_reg2(i, first1:end) = true;
  % das gleiche für Minimalparameter-Regressor
  IImin_set_i = find(I_regmin(i,:));
  first1min = IImin_set_i(1);
  I_regmin(i, first1min:end) = true;
end
figure(1);clf;hold on;grid on;
title(sprintf('%s: Belegung der MinPar-Regressormatrix', robot_name));
II = 1:size(I_regmin,2);
for i = 1:size(I_regmin,1)
  plot(II(I_regmin(i,:)==1), i, 'kx');
  if any(I_regmin(i,:)==0)
    plot(II(I_regmin(i,:)==0), i, 'ro');
  end
end
xlabel('Eintrag #'); ylabel('Achse #');
num_reg2_zero = sum( tau_reg2_mask(:) == 0 );

num_regmin_zero = sum( tau_regmin_mask(:) == 0 )