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
% mrSges_num_mdh [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges_num_mdh [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-04-13 17:19
% Revision: 28dffc800b4c7a157726d217f92010bbd0fa4111
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = SCARA4DOF_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m_num, mrSges_num_mdh, Ifges_num_mdh)
%% Coder Information
%#codegen
assert(isa(qJ,'double') && isreal(qJ) && all(size(qJ) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp2: qJ has to be [4x1] double');
assert(isa(qJD,'double') && isreal(qJD) && all(size(qJD) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp2: qJD has to be [4x1] double');
assert(isa(qJDD,'double') && isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] double');
assert(isa(g,'double') && isreal(g) && all(size(g) == [3 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp2: g has to be [3x1] double');
assert(isa(pkin,'double') && isreal(pkin) && all(size(pkin) == [2 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp2: pkin has to be [2x1] double');
assert(isa(m_num,'double') && isreal(m_num) && all(size(m_num) == [5 1]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp2: m_num has to be [5x1] double'); 
assert(isa(mrSges_num_mdh,'double') && isreal(mrSges_num_mdh) && all(size(mrSges_num_mdh) == [5,3]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp2: mrSges_num_mdh has to be [5x3] double');
assert(isa(Ifges_num_mdh,'double') && isreal(Ifges_num_mdh) && all(size(Ifges_num_mdh) == [5 6]), ...
  'SCARA4DOF_invdynJ_fixb_slag_vp2: Ifges_num_mdh has to be [5x6] double'); 

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

MX0 = mrSges_num_mdh(1,1);
MY0 = mrSges_num_mdh(1,2);
MZ0 = mrSges_num_mdh(1,3);
MX1 = mrSges_num_mdh(2,1);
MY1 = mrSges_num_mdh(2,2);
MZ1 = mrSges_num_mdh(2,3);
MX2 = mrSges_num_mdh(3,1);
MY2 = mrSges_num_mdh(3,2);
MZ2 = mrSges_num_mdh(3,3);
MX3 = mrSges_num_mdh(4,1);
MY3 = mrSges_num_mdh(4,2);
MZ3 = mrSges_num_mdh(4,3);
MX4 = mrSges_num_mdh(5,1);
MY4 = mrSges_num_mdh(5,2);
MZ4 = mrSges_num_mdh(5,3);

XX0 = Ifges_num_mdh(1,1);
XY0 = Ifges_num_mdh(1,4);
XZ0 = Ifges_num_mdh(1,5);
YY0 = Ifges_num_mdh(1,2);
YZ0 = Ifges_num_mdh(1,6);
ZZ0 = Ifges_num_mdh(1,3);
XX1 = Ifges_num_mdh(2,1);
XY1 = Ifges_num_mdh(2,4);
XZ1 = Ifges_num_mdh(2,5);
YY1 = Ifges_num_mdh(2,2);
YZ1 = Ifges_num_mdh(2,6);
ZZ1 = Ifges_num_mdh(2,3);
XX2 = Ifges_num_mdh(3,1);
XY2 = Ifges_num_mdh(3,4);
XZ2 = Ifges_num_mdh(3,5);
YY2 = Ifges_num_mdh(3,2);
YZ2 = Ifges_num_mdh(3,6);
ZZ2 = Ifges_num_mdh(3,3);
XX3 = Ifges_num_mdh(4,1);
XY3 = Ifges_num_mdh(4,4);
XZ3 = Ifges_num_mdh(4,5);
YY3 = Ifges_num_mdh(4,2);
YZ3 = Ifges_num_mdh(4,6);
ZZ3 = Ifges_num_mdh(4,3);
XX4 = Ifges_num_mdh(5,1);
XY4 = Ifges_num_mdh(5,4);
XZ4 = Ifges_num_mdh(5,5);
YY4 = Ifges_num_mdh(5,2);
YZ4 = Ifges_num_mdh(5,6);
ZZ4 = Ifges_num_mdh(5,3);

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-04-13 17:18:18
% EndTime: 2018-04-13 17:18:20
% DurationCPUTime: 0.54s
% Computational Cost: add. (402->102), mult. (743->139), div. (0->0), fcn. (352->8), ass. (0->58)
t64 = L1 * qJD1s;
t43 = sin(qJ2s);
t79 = L1 * t43;
t46 = cos(qJ2s);
t78 = t46 * t64;
t40 = qJD1s + qJD2s;
t66 = MY2 - MY3;
t67 = MX2 + MX3;
t77 = t40 * (t67 * t43 + t66 * t46);
t76 = M3 + M4;
t59 = qJD1s * qJD2s;
t23 = (-qJDD1s * t46 + t43 * t59) * L1;
t73 = t46 * L1;
t24 = (qJDD1s * t43 + t46 * t59) * L1;
t42 = sin(qJ4s);
t45 = cos(qJ4s);
t25 = -t40 * L2 - t78;
t58 = t43 * t64;
t7 = -t42 * t25 + t45 * t58;
t39 = qJDD1s + qJDD2s;
t8 = -t39 * L2 + t23;
t3 = -t7 * qJD4s - t42 * t24 - t45 * t8;
t35 = qJDD4s + t39;
t72 = t3 * MX4 + ZZ4 * t35;
t41 = qJ1s + qJ2s;
t37 = sin(t41);
t38 = cos(t41);
t18 = t37 * t42 - t38 * t45;
t19 = -t37 * t45 - t38 * t42;
t71 = t19 * MX4 + t18 * MY4;
t36 = qJD4s + t40;
t70 = t36 * MX4;
t69 = t42 * t43;
t68 = t43 * t45;
t65 = -t18 * MX4 + t19 * MY4;
t62 = qJD4s * t42;
t61 = qJD4s * t45;
t57 = L1 * M2 + MX1;
t56 = g1 * t37 - g2 * t38;
t55 = -t42 * t46 - t68;
t54 = t45 * t46 - t69;
t6 = -t45 * t25 - t42 * t58;
t52 = t66 * t37 - t67 * t38 - t65;
t51 = (t25 + t78) * t79;
t50 = t67 * t37 + t66 * t38 - t71;
t49 = -t23 * MX2 - t8 * MX3 + ZZ2 * t39 - t66 * t24 + t72;
t47 = cos(qJ1s);
t44 = sin(qJ1s);
t33 = -L2 - t73;
t17 = L1 * t68 - t42 * t33;
t16 = -L1 * t69 - t45 * t33;
t15 = t54 * t64;
t14 = t55 * t64;
t13 = t24 * t79;
t5 = t33 * t62 + (t55 * qJD2s - t43 * t61) * L1;
t4 = -t33 * t61 + (t54 * qJD2s - t43 * t62) * L1;
t2 = t6 * qJD4s + t45 * t24 - t42 * t8;
t1 = [t49 + M3 * (qJD2s * t51 + t8 * t33 + t13) + (t47 * MY1 + t57 * t44 - t76 * (-t44 * L1 - t37 * L2) + t50) * g1 + (t44 * MY1 - t57 * t47 - t76 * (t47 * L1 + t38 * L2) + t52) * g2 + (-t33 * MX3 + ZZ3 + (MX2 * t46 - t66 * t43) * L1) * t39 - qJD2s * L1 * t77 + M4 * (t3 * t16 + t2 * t17 + t7 * t4 + t6 * t5) + M2 * (-t23 * t73 + t13) + (-t17 * t35 - t4 * t36 - t2) * MY4 + (t16 * t35 + t5 * t36) * MX4 + ZZ1 * qJDD1s; -t14 * t70 + ZZ3 * t39 + (t15 * t36 - t2) * MY4 - M4 * (t6 * t14 + t7 * t15) - M3 * qJD1s * t51 + t52 * g2 + t50 * g1 + t64 * t77 + (t39 * MX3 + (-t42 * t35 - t36 * t61) * MY4 + (t45 * t35 - t36 * t62) * MX4 + (t2 * t42 + t3 * t45 - t6 * t62 + t7 * t61 + t56) * M4 + (t56 - t8) * M3) * L2 + t49; t76 * (g3 + qJDD3s); t7 * t70 - g1 * t71 - g2 * t65 + (t6 * t36 - t2) * MY4 + t72;];
tau  = t1 ;
