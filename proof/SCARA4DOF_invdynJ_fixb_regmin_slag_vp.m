% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% SCARA4DOF
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-04-13 17:19
% Revision: 28dffc800b4c7a157726d217f92010bbd0fa4111
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = SCARA4DOF_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
assert(isa(qJ,'double') && isreal(qJ) && all(size(qJ) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] double');
assert(isa(qJD,'double') && isreal(qJD) && all(size(qJD) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] double');
assert(isa(qJDD,'double') && isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] double');
assert(isa(g,'double') && isreal(g) && all(size(g) == [3 1]), ...
  'SCARA4DOF_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] double');
assert(isa(pkin,'double') && isreal(pkin) && all(size(pkin) == [2 1]), ...
  'SCARA4DOF_invdynJ_fixb_regmin_slag_vp: pkin has to be [2x1] double');

%% Variable Initialization
qJ1s = qJ(1);
qJ2s = qJ(2);
qJ3s = qJ(3);
qJ4s = qJ(4);

qJD1s = qJD(1);
qJD2s = qJD(2);
qJD3s = qJD(3);
qJD4s = qJD(4);

qJDD1s = qJDD(1);
qJDD2s = qJDD(2);
qJDD3s = qJDD(3);
qJDD4s = qJDD(4);

g1 = g(1);
g2 = g(2);
g3 = g(3);

L1 = pkin(1);
L2 = pkin(2);

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-04-13 17:18:31
% EndTime: 2018-04-13 17:18:32
% DurationCPUTime: 0.26s
% Computational Cost: add. (227->78), mult. (350->107), div. (0->0), fcn. (204->8), ass. (0->51)
t21 = qJ1s + qJ2s;
t17 = sin(t21);
t18 = cos(t21);
t59 = g1 * t17 - g2 * t18;
t20 = qJD1s + qJD2s;
t26 = cos(qJ2s);
t50 = qJD1s * t26;
t7 = -L1 * t50 - t20 * L2;
t58 = L1 * t7;
t44 = -qJD2s - qJD4s;
t16 = qJD1s - t44;
t56 = t16 * t7;
t55 = t26 * L1;
t22 = sin(qJ4s);
t23 = sin(qJ2s);
t54 = t22 * t23;
t25 = cos(qJ4s);
t53 = t23 * t25;
t52 = g1 * t18 + g2 * t17;
t51 = qJD1s * t23;
t49 = qJD2s * t26;
t48 = qJD4s * t22;
t47 = qJD4s * t25;
t46 = qJDD1s * t23;
t19 = qJDD1s + qJDD2s;
t15 = qJDD4s + t19;
t45 = -t15 - qJDD1s;
t43 = L1 * t51;
t28 = L1 ^ 2;
t42 = t28 * t50;
t41 = t22 * t49;
t40 = t25 * t49;
t39 = t22 * t46;
t38 = t25 * t46;
t37 = qJD1s * (-qJD2s + t20);
t36 = qJD2s * (-qJD1s - t20);
t35 = t44 * t16;
t14 = qJDD1s * t55;
t34 = t14 + t59;
t33 = (-qJD1s - t16) * t49;
t32 = t22 * t15 + t16 * t47;
t31 = -t25 * t15 + t16 * t48;
t2 = -t19 * L2 + qJD2s * t43 - t14;
t3 = t17 * t22 - t18 * t25;
t4 = -t17 * t25 - t18 * t22;
t30 = -g1 * t3 - g2 * t4 + t22 * t2 + t43 * t48 + t7 * t47;
t29 = -g1 * t4 + g2 * t3 - t25 * t2 + t7 * t48;
t27 = cos(qJ1s);
t24 = sin(qJ1s);
t13 = -L2 - t55;
t1 = [qJDD1s g1 * t24 - g2 * t27 g1 * t27 + g2 * t24 t19 (t19 * t26 + t23 * t36) * L1 + t34 ((-qJDD1s - t19) * t23 + t26 * t36) * L1 + t52 t2 * t13 - g1 * (-t24 * L1 - t17 * L2) - g2 * (t27 * L1 + t18 * L2) + (t28 * t46 + (0.2e1 * t42 + t58) * qJD2s) * t23 t15 t31 * t13 + (t22 * t33 + (t45 * t22 + (-qJD1s * qJD4s + t35) * t25) * t23) * L1 + t29 t32 * t13 + (t25 * t33 + (-t22 * t35 + t45 * t25) * t23) * L1 + t30; 0 0 0 t19 t23 * L1 * t37 + t34 (t26 * t37 - t46) * L1 + t52 (-t42 - t58) * t51 + (-t2 + t59) * L2 t15 -t31 * L2 + (-t39 + (-t23 * t47 - t41 - (-t22 * t26 - t53) * t16) * qJD1s) * L1 + t29 -t32 * L2 + (-t38 + (-t40 + (t25 * t26 - t54) * t16) * qJD1s) * L1 + t30; 0 0 0 0 0 0 qJDD3s + g3 0 0 0; 0 0 0 0 0 0 0 t15 -t22 * t56 + (-t39 + (-t41 + (-qJD4s + t16) * t53) * qJD1s) * L1 + t29 -t25 * t56 + (-t38 + (-t16 * t54 - t40) * qJD1s) * L1 + t30;];
tau_reg  = t1 ;
