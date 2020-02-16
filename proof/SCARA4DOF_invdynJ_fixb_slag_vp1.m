% Calculate vector of inverse dynamics joint torques for
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
% m_num_mdh [5x1]
%   mass of all robot links (including the base)
% rSges_num_mdh [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges_num_mdh [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-04-13 17:19
% Revision: 28dffc800b4c7a157726d217f92010bbd0fa4111
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = SCARA4DOF_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m_num, rSges_num_mdh, Icges_num_mdh)
%% Coder Information
%#codegen
assert(isa(qJ,'double') && isreal(qJ) && all(size(qJ) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp1: qJ has to be [4x1] double');
assert(isa(qJD,'double') && isreal(qJD) && all(size(qJD) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp1: qJD has to be [4x1] double');
assert(isa(qJDD,'double') && isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] double');
assert(isa(g,'double') && isreal(g) && all(size(g) == [3 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp1: g has to be [3x1] double');
assert(isa(pkin,'double') && isreal(pkin) && all(size(pkin) == [2 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp1: pkin has to be [2x1] double');
assert(isa(m_num,'double') && isreal(m_num) && all(size(m_num) == [5 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp1: m_num has to be [5x1] double'); 
assert(isa(rSges_num_mdh,'double') && isreal(rSges_num_mdh) && all(size(rSges_num_mdh) == [5,3]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp1: rSges_num_mdh has to be [5x3] double');
assert(isa(Icges_num_mdh,'double') && isreal(Icges_num_mdh) && all(size(Icges_num_mdh) == [5 6]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp1: Icges_num_mdh has to be [5x6] double'); 

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

M0 = m_num(1);
M1 = m_num(2);
M2 = m_num(3);
M3 = m_num(4);
M4 = m_num(5);

SX0 = rSges_num_mdh(1,1);
SY0 = rSges_num_mdh(1,2);
SZ0 = rSges_num_mdh(1,3);
SX1 = rSges_num_mdh(2,1);
SY1 = rSges_num_mdh(2,2);
SZ1 = rSges_num_mdh(2,3);
SX2 = rSges_num_mdh(3,1);
SY2 = rSges_num_mdh(3,2);
SZ2 = rSges_num_mdh(3,3);
SX3 = rSges_num_mdh(4,1);
SY3 = rSges_num_mdh(4,2);
SZ3 = rSges_num_mdh(4,3);
SX4 = rSges_num_mdh(5,1);
SY4 = rSges_num_mdh(5,2);
SZ4 = rSges_num_mdh(5,3);

XXC0 = Icges_num_mdh(1,1);
XYC0 = Icges_num_mdh(1,4);
XZC0 = Icges_num_mdh(1,5);
YYC0 = Icges_num_mdh(1,2);
YZC0 = Icges_num_mdh(1,6);
ZZC0 = Icges_num_mdh(1,3);
XXC1 = Icges_num_mdh(2,1);
XYC1 = Icges_num_mdh(2,4);
XZC1 = Icges_num_mdh(2,5);
YYC1 = Icges_num_mdh(2,2);
YZC1 = Icges_num_mdh(2,6);
ZZC1 = Icges_num_mdh(2,3);
XXC2 = Icges_num_mdh(3,1);
XYC2 = Icges_num_mdh(3,4);
XZC2 = Icges_num_mdh(3,5);
YYC2 = Icges_num_mdh(3,2);
YZC2 = Icges_num_mdh(3,6);
ZZC2 = Icges_num_mdh(3,3);
XXC3 = Icges_num_mdh(4,1);
XYC3 = Icges_num_mdh(4,4);
XZC3 = Icges_num_mdh(4,5);
YYC3 = Icges_num_mdh(4,2);
YZC3 = Icges_num_mdh(4,6);
ZZC3 = Icges_num_mdh(4,3);
XXC4 = Icges_num_mdh(5,1);
XYC4 = Icges_num_mdh(5,4);
XZC4 = Icges_num_mdh(5,5);
YYC4 = Icges_num_mdh(5,2);
YZC4 = Icges_num_mdh(5,6);
ZZC4 = Icges_num_mdh(5,3);

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-04-13 17:18:18
% EndTime: 2018-04-13 17:18:20
% DurationCPUTime: 0.80s
% Computational Cost: add. (1198->115), mult. (1176->128), div. (0->0), fcn. (792->6), ass. (0->73)
t73 = sin(qJ1s);
t74 = cos(qJ1s);
t75 = qJD1s ^ 2;
t119 = (-qJDD1s * t73 - t74 * t75) * L1 - g1;
t71 = qJ1s + qJ2s;
t65 = sin(t71);
t66 = cos(t71);
t70 = qJD1s + qJD2s;
t68 = t70 ^ 2;
t69 = qJDD1s + qJDD2s;
t118 = (-t65 * t69 - t66 * t68) * L2 + t119;
t72 = sin(qJ4s);
t99 = cos(qJ4s);
t38 = t65 * t72 - t66 * t99;
t39 = t65 * t99 + t66 * t72;
t18 = t39 * SX4 - t38 * SY4;
t64 = qJD4s + t70;
t16 = t38 * t64;
t17 = t64 * t39;
t4 = -SX4 * t16 - SY4 * t17;
t63 = qJDD4s + t69;
t117 = -t18 * t63 - t4 * t64 + t118;
t102 = L1 * t73;
t67 = t74 * L1;
t115 = qJDD1s * t67 - t75 * t102 - g2;
t114 = t18 * t64;
t101 = L2 * t65;
t113 = -t18 - t101;
t59 = t66 * L2;
t112 = -t68 * t101 + t69 * t59 + t115;
t94 = L1 * qJD1s;
t89 = t73 * t94;
t41 = SX2 * t65 + SY2 * t66;
t98 = t41 * t70;
t23 = -t89 - t98;
t56 = t66 * SY3;
t40 = SX3 * t65 - t56;
t42 = t66 * SX3 + t65 * SY3;
t111 = -t69 * t40 - t42 * t68 + t118;
t95 = t66 * t70;
t97 = t65 * t70;
t28 = SX2 * t95 - SY2 * t97;
t110 = -t28 * t70 - t41 * t69 + t119;
t20 = -t38 * SX4 - t39 * SY4;
t3 = -t17 * SX4 + t16 * SY4;
t109 = t20 * t63 + t3 * t64 + t112;
t49 = SY3 * t95;
t108 = t69 * t42 + t70 * (-SX3 * t97 + t49) + t112;
t43 = t66 * SX2 - SY2 * t65;
t107 = t43 * t69 - t70 * t98 + t115;
t106 = t3 + t114;
t103 = t59 + t42;
t105 = t103 * t70;
t104 = t70 * t40 + t49;
t100 = -L2 - SX3;
t91 = L2 * t95;
t54 = ZZC4 * t63;
t90 = t54 + (ZZC2 + ZZC3) * t69;
t14 = t59 + t20;
t88 = t74 * t94;
t29 = t100 * t65 + t56;
t51 = SX1 * t74 - SY1 * t73;
t50 = SX1 * t73 + SY1 * t74;
t83 = -t20 * t64 - t91;
t81 = -L2 * t97 - t89;
t79 = -t4 - t91;
t21 = -t89 + (-t40 - t101) * t70;
t22 = t88 + t105;
t76 = (t21 * t100 * t66 + (-t21 * SY3 + t22 * t100) * t65) * t70;
t24 = t43 * t70 + t88;
t8 = -t83 + t88;
t7 = t81 - t114;
t1 = [ZZC1 * qJDD1s + t90 + (t107 * (t43 + t67) + t110 * (-t41 - t102) + (-t28 - t88 + t24) * t23) * M2 + (g1 * t50 - g2 * t51 + (t50 ^ 2 + t51 ^ 2) * qJDD1s) * M1 + (t7 * (t79 - t88) + t109 * (t67 + t14) + (t7 + t106) * t8 + t117 * (t113 - t102)) * M4 + (-t21 * t88 + t76 + t108 * (t67 + t103) + t111 * (t29 - t102) + (t21 - t81 - t89 + t104) * t22) * M3; t90 + ((-t83 + t79) * t7 + t109 * t14 + t106 * t8 + t117 * t113) * M4 + (t21 * t105 + t76 + t108 * t103 + t111 * t29 + (t101 * t70 + t104) * t22) * M3 + (-t23 * t28 - t24 * t98 + (t23 * t70 + t107) * t43 + (t24 * t70 - t110) * t41) * M2; (M3 + M4) * (g3 + qJDD3s); t54 + (t8 * t3 - t7 * t4 + (t7 * t64 + t109) * t20 + (t8 * t64 - t117) * t18) * M4;];
tau  = t1 ;
