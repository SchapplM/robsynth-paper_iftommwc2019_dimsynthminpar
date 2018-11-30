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

num_reg2_zero = sum( tau_reg2_mask(:) == 0 );

num_regmin_zero = sum( tau_regmin_mask(:) == 0 )