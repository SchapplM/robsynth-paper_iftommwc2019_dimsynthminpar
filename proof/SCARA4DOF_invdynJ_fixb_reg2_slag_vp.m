% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-04-13 17:19
% Revision: 28dffc800b4c7a157726d217f92010bbd0fa4111
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = SCARA4DOF_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
assert(isa(qJ,'double') && isreal(qJ) && all(size(qJ) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] double');
assert(isa(qJD,'double') && isreal(qJD) && all(size(qJD) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] double');
assert(isa(qJDD,'double') && isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'SCARA4DOF_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] double');
assert(isa(g,'double') && isreal(g) && all(size(g) == [3 1]), ...
  'SCARA4DOF_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] double');
assert(isa(pkin,'double') && isreal(pkin) && all(size(pkin) == [2 1]), ...
  'SCARA4DOF_invdynJ_fixb_reg2_slag_vp: pkin has to be [2x1] double');

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
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-04-13 17:18:40
% EndTime: 2018-04-13 17:18:41
% DurationCPUTime: 0.30s
% Computational Cost: add. (371->92), mult. (643->124), div. (0->0), fcn. (356->8), ass. (0->71)
t44 = cos(qJ4s);
t42 = sin(qJ2s);
t71 = qJDD1s * t42;
t45 = cos(qJ2s);
t74 = qJD2s * t45;
t51 = (qJD1s * t74 + t71) * L1;
t86 = t44 * t51;
t37 = qJDD1s + qJDD2s;
t38 = qJD1s + qJD2s;
t40 = qJ1s + qJ2s;
t35 = sin(t40);
t36 = cos(t40);
t63 = -g1 * t36 - g2 * t35 + t51;
t85 = (t37 * t42 + t38 * t74) * L1 + t63;
t76 = qJD1s * t45;
t67 = L1 * t76;
t18 = -t38 * L2 - t67;
t84 = L1 * t18;
t83 = g2 * t36;
t82 = t35 * L2;
t33 = t37 * L2;
t81 = t45 * L1;
t41 = sin(qJ4s);
t80 = t41 * t42;
t79 = t42 * t44;
t78 = L1 * qJD1s;
t77 = qJD1s * t42;
t75 = qJD2s * t42;
t73 = qJD4s * t41;
t72 = qJD4s * t44;
t47 = L1 ^ 2;
t70 = qJDD1s * t47;
t68 = L1 * t77;
t62 = t41 * t68;
t31 = qJDD1s * t81;
t8 = qJD2s * t68 - t31 - t33;
t69 = -qJD4s * t62 - t18 * t72 - t41 * t8;
t66 = t47 * t76;
t65 = t42 * t72;
t61 = g1 * t35 + t31 - t83;
t43 = sin(qJ1s);
t46 = cos(qJ1s);
t60 = -g1 * (-t43 * L1 - t82) - g2 * (t46 * L1 + t36 * L2);
t59 = g1 * t43 - g2 * t46;
t58 = (-qJD1s - t38) * t75;
t57 = -t41 * t45 - t79;
t56 = t44 * t45 - t80;
t13 = t35 * t41 - t36 * t44;
t14 = -t35 * t44 - t36 * t41;
t53 = -g1 * t13 - g2 * t14 - t69;
t52 = t38 * t67 - t63;
t50 = (-qJD2s + t38) * t68 + t61;
t49 = t53 - t86;
t2 = -t44 * t8 + t18 * t73 + (-t41 * t71 + (-t41 * t74 - t65) * qJD1s) * L1;
t48 = -g1 * t14 + g2 * t13 + t2;
t39 = qJDD3s + g3;
t34 = qJD4s + t38;
t32 = qJDD4s + t37;
t29 = t42 ^ 2 * t70;
t28 = -L2 - t81;
t24 = g1 * t82;
t12 = L1 * t79 - t41 * t28;
t11 = -L1 * t80 - t44 * t28;
t10 = t56 * t78;
t9 = t57 * t78;
t6 = -t41 * t18 + t44 * t68;
t5 = -t44 * t18 - t62;
t4 = t28 * t73 + (t57 * qJD2s - t65) * L1;
t3 = -t28 * t72 + (t56 * qJD2s - t42 * t73) * L1;
t1 = t69 + t86;
t7 = [0 0 0 0 0 qJDD1s t59 g1 * t46 + g2 * t43 0 0 0 0 0 0 0 t37 (t37 * t45 + t58) * L1 + t61 -t85 0 t45 ^ 2 * t70 + t59 * L1 + t29 0 0 0 0 0 t37 L1 * t58 - t28 * t37 + t33 + t61 t85 0 t8 * t28 + t29 + (0.2e1 * t66 + t84) * t75 + t60 0 0 0 0 0 t32 t11 * t32 + t4 * t34 + t48 -t12 * t32 - t3 * t34 + t49 0 t1 * t12 + t2 * t11 + t6 * t3 + t5 * t4 + t60; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t37 t50 t52 0 0 0 0 0 0 0 t37 0.2e1 * t33 + t50 -t52 0 t24 + (-t66 - t84) * t77 + (-t8 - t83) * L2 0 0 0 0 0 t32 -t9 * t34 + (t32 * t44 - t34 * t73) * L2 + t48 -t41 * L2 * t32 + t10 * t34 + (-L2 * qJD4s * t34 - t51) * t44 + t53 0 -t6 * t10 - t5 * t9 + t24 + (-t83 + t1 * t41 + t2 * t44 + (-t41 * t5 + t44 * t6) * qJD4s) * L2; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t39 0 0 0 0 0 0 0 0 0 t39; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t32 t6 * t34 + t48 t5 * t34 + t49 0 0;];
tau_reg  = t7 ;
