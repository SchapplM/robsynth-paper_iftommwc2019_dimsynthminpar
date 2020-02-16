% Teste Eigenschaften der Regressormatrix für ausgewählte Roboter
% Daten für Tabelle 1 im Paper (Zeile E)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-11
% (C) Institut für mechatronische Systeme, Universität Hannover

clear
clc

%% Benutzereingabe

Rob_Names = {'S7RRRRRRR1', 'S6RRRRRR10V2', 'S4RRPR1'};
%% Schleife über alle Roboter
for rr = 1:length(Rob_Names)
  robot_name = Rob_Names{rr};
  %% Init 2
  RS = serroblib_create_robot_class(robot_name);
  TSS = RS.gen_testsettings(true, true);
  % serroblib_addtopath({robot_name});
  % eval(sprintf('TSS = %s_varpar_testfunctions_parameter();', robot_name));
  pkin = RS.pkin;

  %% Teste Strukturelle Eigenschaften

  eval(sprintf('func_regmin = @%s_invdynJ_fixb_regmin_slag_vp;', robot_name));
  eval(sprintf('func_reg2 = @%s_invdynJ_fixb_reg2_slag_vp;', robot_name));
  NQJ = RS.NQJ;
  NPV2 = RS.NJ*10;
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
  figure(rr);clf;hold on;grid on;
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

  num_regmin_zero = sum( tau_regmin_mask(:) == 0 );
  fprintf('Robot %s: %d/%d zeros in regressor matrix\n', robot_name, num_regmin_zero, length(tau_regmin_mask(:)));
end