% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:06
% EndTime: 2019-03-08 18:35:06
% DurationCPUTime: 0.22s
% Computational Cost: add. (385->92), mult. (662->125), div. (0->0), fcn. (386->12), ass. (0->59)
t42 = sin(qJ(2));
t60 = qJDD(1) * t42;
t45 = cos(qJ(2));
t61 = qJD(1) * t45;
t71 = pkin(1) * (qJD(2) * t61 + t60);
t38 = qJ(1) + qJ(2);
t34 = sin(t38);
t35 = cos(t38);
t70 = g(1) * t34 - g(2) * t35;
t37 = qJD(1) + qJD(2);
t33 = qJD(4) + t37;
t69 = qJD(4) - t33;
t68 = pkin(1) * t42;
t39 = sin(pkin(7));
t67 = pkin(2) * t39;
t65 = t45 * pkin(1);
t64 = t39 * t42;
t40 = cos(pkin(7));
t63 = t40 * t42;
t62 = g(1) * t35 + g(2) * t34;
t29 = pkin(7) + qJ(4) + t38;
t23 = sin(t29);
t24 = cos(t29);
t41 = sin(qJ(4));
t18 = pkin(1) * t61 + t37 * pkin(2);
t58 = qJD(1) * t68;
t8 = t39 * t18 + t40 * t58;
t59 = qJD(4) * t41 * t8 + g(1) * t24 + g(2) * t23;
t36 = qJDD(1) + qJDD(2);
t57 = qJD(1) * (-qJD(2) + t37);
t56 = qJD(2) * (-qJD(1) - t37);
t30 = pkin(2) + t65;
t55 = -pkin(1) * t64 + t40 * t30;
t31 = qJDD(1) * t65;
t54 = t31 + t70;
t44 = cos(qJ(4));
t7 = t40 * t18 - t39 * t58;
t5 = t37 * pkin(3) + t7;
t53 = -t41 * t5 - t44 * t8;
t25 = t40 * pkin(2) + pkin(3);
t52 = t41 * t25 + t44 * t67;
t51 = t44 * t25 - t41 * t67;
t50 = pkin(1) * (-t39 * t45 - t63);
t49 = pkin(1) * (t40 * t45 - t64);
t10 = t36 * pkin(2) - qJD(2) * t58 + t31;
t3 = t40 * t10 - t39 * t71;
t2 = t36 * pkin(3) + t3;
t4 = t39 * t10 + t40 * t71;
t48 = g(1) * t23 - g(2) * t24 + t44 * t2 - t41 * t4;
t46 = cos(qJ(1));
t43 = sin(qJ(1));
t32 = qJDD(4) + t36;
t16 = pkin(1) * t63 + t39 * t30;
t15 = qJD(2) * t49;
t14 = qJD(1) * t49;
t13 = qJD(2) * t50;
t12 = qJD(1) * t50;
t11 = pkin(3) + t55;
t1 = [qJDD(1), g(1) * t43 - g(2) * t46, g(1) * t46 + g(2) * t43, t36 (t36 * t45 + t42 * t56) * pkin(1) + t54 ((-qJDD(1) - t36) * t42 + t45 * t56) * pkin(1) + t62, t4 * t16 + t8 * t15 + t3 * t55 + t7 * t13 - g(1) * (-t43 * pkin(1) - pkin(2) * t34) - g(2) * (t46 * pkin(1) + pkin(2) * t35) t32 (t44 * t13 - t41 * t15) * t33 + (t44 * t11 - t41 * t16) * t32 + ((-t41 * t11 - t44 * t16) * t33 + t53) * qJD(4) + t48 (-(-qJD(4) * t16 + t13) * t33 - t11 * t32 - t2) * t41 + (-(qJD(4) * t11 + t15) * t33 - t16 * t32 - t4 - qJD(4) * t5) * t44 + t59; 0, 0, 0, t36, t57 * t68 + t54 (t45 * t57 - t60) * pkin(1) + t62, -t7 * t12 - t8 * t14 + (t3 * t40 + t39 * t4 + t70) * pkin(2), t32, t51 * t32 - (t44 * t12 - t41 * t14) * t33 + (-t52 * t33 + t53) * qJD(4) + t48, -t52 * t32 - t44 * t4 - t41 * t2 + (t41 * t12 + t44 * t14) * t33 + (-t51 * t33 - t44 * t5) * qJD(4) + t59; 0, 0, 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t32, t69 * t53 + t48 (-t8 * t33 - t2) * t41 + (-t69 * t5 - t4) * t44 + t59;];
tau_reg  = t1;
