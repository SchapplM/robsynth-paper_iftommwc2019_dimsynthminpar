% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:02
% EndTime: 2019-03-08 18:35:03
% DurationCPUTime: 0.39s
% Computational Cost: add. (718->114), mult. (1296->148), div. (0->0), fcn. (762->14), ass. (0->85)
t65 = cos(pkin(7));
t67 = sin(qJ(2));
t89 = qJDD(1) * t67;
t70 = cos(qJ(2));
t90 = qJD(1) * t70;
t77 = pkin(1) * (qJD(2) * t90 + t89);
t108 = t65 * t77;
t63 = qJ(1) + qJ(2);
t57 = sin(t63);
t58 = cos(t63);
t107 = g(1) * t57 - g(2) * t58;
t60 = qJDD(1) + qJDD(2);
t102 = t60 * pkin(2);
t99 = t70 * pkin(1);
t53 = qJDD(1) * t99;
t106 = pkin(1) * t67;
t88 = qJD(1) * t106;
t22 = -qJD(2) * t88 + t102 + t53;
t64 = sin(pkin(7));
t97 = t64 * t22;
t14 = t97 + t108;
t61 = qJD(1) + qJD(2);
t33 = pkin(1) * t90 + t61 * pkin(2);
t17 = t65 * t33 - t64 * t88;
t15 = t61 * pkin(3) + t17;
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t13 = t65 * t22 - t64 * t77;
t8 = t60 * pkin(3) + t13;
t18 = t64 * t33 + t65 * t88;
t94 = t66 * t18;
t1 = (qJD(4) * t15 + t14) * t69 - qJD(4) * t94 + t66 * t8;
t105 = pkin(2) * t57;
t104 = pkin(2) * t64;
t101 = t65 * pkin(2);
t68 = sin(qJ(1));
t100 = t68 * pkin(1);
t95 = t65 * t67;
t80 = pkin(1) * (-t64 * t70 - t95);
t24 = qJD(1) * t80;
t96 = t64 * t67;
t79 = pkin(1) * (t65 * t70 - t96);
t26 = qJD(1) * t79;
t46 = pkin(3) + t101;
t31 = t69 * t104 + t66 * t46;
t98 = -t31 * qJD(4) - t69 * t24 + t66 * t26;
t29 = -t66 * t104 + t69 * t46;
t93 = t29 * qJD(4) - t66 * t24 - t69 * t26;
t92 = g(1) * t58 + g(2) * t57;
t50 = pkin(2) * t58;
t71 = cos(qJ(1));
t91 = t71 * pkin(1) + t50;
t56 = pkin(7) + t63;
t85 = qJD(1) * (-qJD(2) + t61);
t84 = qJD(2) * (-qJD(1) - t61);
t52 = pkin(2) + t99;
t28 = -pkin(1) * t96 + t65 * t52;
t83 = t53 + t107;
t44 = sin(t56);
t82 = -pkin(3) * t44 - t105;
t81 = g(1) * t68 - g(2) * t71;
t6 = t66 * t15 + t69 * t18;
t23 = pkin(3) + t28;
t30 = pkin(1) * t95 + t64 * t52;
t11 = t69 * t23 - t66 * t30;
t12 = t66 * t23 + t69 * t30;
t2 = -t6 * qJD(4) - t66 * t14 + t69 * t8;
t51 = qJ(4) + t56;
t42 = sin(t51);
t43 = cos(t51);
t76 = g(1) * t43 + g(2) * t42 - t1;
t45 = cos(t56);
t75 = g(1) * t45 + g(2) * t44 - t108;
t74 = g(1) * t42 - g(2) * t43 + t2;
t73 = g(1) * t44 - g(2) * t45 + t13;
t62 = qJDD(3) - g(3);
t55 = qJD(4) + t61;
t54 = qJDD(4) + t60;
t41 = pkin(3) * t45;
t27 = qJD(2) * t79;
t25 = qJD(2) * t80;
t5 = t69 * t15 - t94;
t4 = -t12 * qJD(4) + t69 * t25 - t66 * t27;
t3 = t11 * qJD(4) + t66 * t25 + t69 * t27;
t7 = [0, 0, 0, 0, 0, qJDD(1), t81, g(1) * t71 + g(2) * t68, 0, 0, 0, 0, 0, 0, 0, t60 (t60 * t70 + t67 * t84) * pkin(1) + t83 ((-qJDD(1) - t60) * t67 + t70 * t84) * pkin(1) + t92, 0 (t81 + (t67 ^ 2 + t70 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t60, t25 * t61 + t28 * t60 + t73, -t27 * t61 - t30 * t60 + t75 - t97, 0, t14 * t30 + t18 * t27 + t13 * t28 + t17 * t25 - g(1) * (-t100 - t105) - g(2) * t91, 0, 0, 0, 0, 0, t54, t11 * t54 + t4 * t55 + t74, -t12 * t54 - t3 * t55 + t76, 0, t1 * t12 + t6 * t3 + t2 * t11 + t5 * t4 - g(1) * (t82 - t100) - g(2) * (t41 + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t85 * t106 + t83 (t70 * t85 - t89) * pkin(1) + t92, 0, 0, 0, 0, 0, 0, 0, t60, t60 * t101 - t24 * t61 + t73, t26 * t61 + (-t22 - t102) * t64 + t75, 0, -t17 * t24 - t18 * t26 + (t13 * t65 + t14 * t64 + t107) * pkin(2), 0, 0, 0, 0, 0, t54, t29 * t54 + t98 * t55 + t74, -t31 * t54 - t93 * t55 + t76, 0, t1 * t31 + t2 * t29 - g(1) * t82 - g(2) * (t41 + t50) + t93 * t6 + t98 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t6 * t55 + t74, t5 * t55 + t76, 0, 0;];
tau_reg  = t7;
