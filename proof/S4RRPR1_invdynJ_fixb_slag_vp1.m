% Calculate vector of inverse dynamics joint torques for
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
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:59
% EndTime: 2019-03-08 18:35:00
% DurationCPUTime: 0.64s
% Computational Cost: add. (1318->117), mult. (942->131), div. (0->0), fcn. (448->8), ass. (0->73)
t67 = sin(qJ(1));
t68 = cos(qJ(1));
t69 = qJD(1) ^ 2;
t106 = (-qJDD(1) * t67 - t68 * t69) * pkin(1) - g(1);
t62 = t68 * pkin(1);
t93 = pkin(1) * t67;
t105 = qJDD(1) * t62 - t69 * t93 - g(2);
t66 = qJ(1) + qJ(2);
t60 = sin(t66);
t61 = cos(t66);
t65 = qJD(1) + qJD(2);
t63 = t65 ^ 2;
t64 = qJDD(1) + qJDD(2);
t104 = (-t60 * t64 - t61 * t63) * pkin(2) + t106;
t52 = pkin(2) * t61;
t92 = pkin(2) * t60;
t103 = t64 * t52 - t63 * t92 + t105;
t59 = pkin(7) + t66;
t53 = qJ(4) + t59;
t46 = sin(t53);
t47 = cos(t53);
t27 = rSges(5,1) * t46 + rSges(5,2) * t47;
t58 = qJD(4) + t65;
t19 = t58 * t27;
t84 = pkin(1) * qJD(1);
t80 = t67 * t84;
t34 = rSges(3,1) * t60 + rSges(3,2) * t61;
t88 = t34 * t65;
t20 = -t80 - t88;
t50 = sin(t59);
t51 = cos(t59);
t31 = rSges(4,1) * t50 + rSges(4,2) * t51;
t90 = rSges(4,2) * t50;
t38 = t65 * t90;
t44 = t51 * rSges(4,1);
t101 = -t64 * t31 - t65 * (t44 * t65 - t38) + t104;
t49 = t61 * rSges(3,1);
t85 = t60 * t65;
t26 = -rSges(3,2) * t85 + t49 * t65;
t100 = -t26 * t65 - t34 * t64 + t106;
t43 = t47 * rSges(5,1);
t28 = -rSges(5,2) * t46 + t43;
t57 = qJDD(4) + t64;
t99 = -t19 * t58 + t28 * t57 + (-t50 * t63 + t51 * t64) * pkin(3) + t103;
t32 = t44 - t90;
t98 = -t31 * t63 + t64 * t32 + t103;
t35 = -rSges(3,2) * t60 + t49;
t97 = t35 * t64 - t65 * t88 + t105;
t87 = t46 * t58;
t33 = rSges(5,2) * t87;
t86 = t47 * t58;
t16 = rSges(5,1) * t86 - t33;
t96 = t16 * t58 + t27 * t57 - (-t50 * t64 - t51 * t63) * pkin(3) - t104;
t94 = t32 + t52;
t76 = pkin(3) * t51 + t52;
t77 = -pkin(3) * t50 - t92;
t7 = t65 * t77 - t19 - t80;
t89 = t28 * t58;
t48 = Icges(5,3) * t57;
t82 = t48 + (Icges(3,3) + Icges(4,3)) * t64;
t81 = t68 * t84;
t22 = -t31 - t92;
t42 = rSges(2,1) * t68 - rSges(2,2) * t67;
t41 = rSges(2,1) * t67 + rSges(2,2) * t68;
t14 = t28 + t76;
t13 = -t27 + t77;
t11 = t22 * t65 - t80;
t12 = t65 * t94 + t81;
t70 = (t11 * (-t44 - t52) + t12 * t22) * t65;
t24 = t65 * t31;
t21 = t35 * t65 + t81;
t8 = t65 * t76 + t81 + t89;
t1 = [Icges(2,3) * qJDD(1) + t82 + (t97 * (t35 + t62) + t100 * (-t34 - t93) + (-t26 - t81 + t21) * t20) * m(3) + (g(1) * t41 - g(2) * t42 + (t41 ^ 2 + t42 ^ 2) * qJDD(1)) * m(2) + (t7 * (-t16 - t81) + (-t7 * t76 + t77 * t8) * t65 - t96 * (t13 - t93) + t8 * (-rSges(5,1) * t87 - rSges(5,2) * t86 - t80) + t99 * (t14 + t62)) * m(5) + (t11 * (t38 - t81) + t70 + t98 * (t94 + t62) + t101 * (t22 - t93) + (pkin(2) * t85 + t11 + t24) * t12) * m(4); t82 + (t12 * t24 - (-t11 * t94 - t12 * t92) * t65 + t11 * t38 + t70 + t98 * t94 + t101 * t22) * m(4) + (-t20 * t26 - t21 * t88 + (t20 * t65 + t97) * t35 + (t21 * t65 - t100) * t34) * m(3) + ((-t43 * t58 + t33 + t89) * t7 + t99 * t14 - t96 * t13) * m(5); (m(4) + m(5)) * (-g(3) + qJDD(3)); t48 + (-t7 * t16 + (t58 * t7 + t99) * t28 + t96 * t27) * m(5);];
tau  = t1;
