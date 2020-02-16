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
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:54
% EndTime: 2019-03-08 18:34:55
% DurationCPUTime: 0.50s
% Computational Cost: add. (748->133), mult. (1380->176), div. (0->0), fcn. (738->14), ass. (0->75)
t90 = m(4) + m(5);
t63 = sin(pkin(7));
t69 = cos(qJ(2));
t64 = cos(pkin(7));
t66 = sin(qJ(2));
t82 = t64 * t66;
t75 = pkin(1) * (-t63 * t69 - t82);
t24 = qJD(1) * t75;
t83 = t63 * t66;
t74 = pkin(1) * (t64 * t69 - t83);
t26 = qJD(1) * t74;
t87 = pkin(2) * t64;
t47 = pkin(3) + t87;
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t88 = pkin(2) * t63;
t31 = t47 * t65 + t68 * t88;
t92 = -t31 * qJD(4) - t24 * t68 + t26 * t65;
t29 = t47 * t68 - t65 * t88;
t91 = -t29 * qJD(4) + t24 * t65 + t26 * t68;
t89 = pkin(1) * t69;
t62 = qJ(1) + qJ(2);
t58 = cos(t62);
t48 = pkin(2) * t58;
t56 = pkin(7) + t62;
t49 = qJ(4) + t56;
t41 = sin(t49);
t86 = g(1) * t41;
t67 = sin(qJ(1));
t85 = g(1) * t67;
t80 = pkin(1) * qJD(1);
t77 = t66 * t80;
t34 = -qJD(2) * t77 + qJDD(1) * t89;
t60 = qJDD(1) + qJDD(2);
t22 = pkin(2) * t60 + t34;
t79 = qJD(2) * t69;
t35 = (qJD(1) * t79 + qJDD(1) * t66) * pkin(1);
t15 = t22 * t63 + t35 * t64;
t61 = qJD(1) + qJD(2);
t36 = pkin(2) * t61 + t69 * t80;
t17 = t64 * t36 - t63 * t77;
t16 = pkin(3) * t61 + t17;
t18 = t36 * t63 + t64 * t77;
t7 = t16 * t65 + t18 * t68;
t14 = t64 * t22 - t35 * t63;
t8 = pkin(3) * t60 + t14;
t3 = -t7 * qJD(4) - t15 * t65 + t68 * t8;
t54 = qJDD(4) + t60;
t84 = t3 * mrSges(5,1) + Ifges(5,3) * t54;
t70 = cos(qJ(1));
t81 = t70 * pkin(1) + t48;
t50 = pkin(2) + t89;
t28 = -pkin(1) * t83 + t64 * t50;
t6 = t16 * t68 - t18 * t65;
t23 = pkin(3) + t28;
t30 = pkin(1) * t82 + t50 * t63;
t11 = t23 * t68 - t30 * t65;
t12 = t23 * t65 + t30 * t68;
t2 = t6 * qJD(4) + t15 * t68 + t65 * t8;
t42 = cos(t49);
t76 = g(1) * t42 + g(2) * t41 - t2;
t38 = t42 * mrSges(5,1);
t45 = sin(t56);
t46 = cos(t56);
t57 = sin(t62);
t73 = -t58 * mrSges(3,1) - t46 * mrSges(4,1) + t57 * mrSges(3,2) + t45 * mrSges(4,2) - t38;
t72 = t34 * mrSges(3,1) + t14 * mrSges(4,1) - t35 * mrSges(3,2) + t84 + (Ifges(3,3) + Ifges(4,3)) * t60;
t71 = t58 * mrSges(3,2) + t46 * mrSges(4,2) + (m(5) * pkin(3) + mrSges(4,1)) * t45 + (t90 * pkin(2) + mrSges(3,1)) * t57;
t55 = qJD(4) + t61;
t40 = pkin(3) * t46;
t27 = qJD(2) * t74;
t25 = qJD(2) * t75;
t5 = -t12 * qJD(4) + t25 * t68 - t27 * t65;
t4 = t11 * qJD(4) + t25 * t65 + t27 * t68;
t1 = [(-t12 * t54 - t4 * t55 - t2) * mrSges(5,2) + (-t27 * t61 - t30 * t60 - t15) * mrSges(4,2) + (t11 * t54 + t5 * t55) * mrSges(5,1) + (t25 * t61 + t28 * t60) * mrSges(4,1) + t72 + (-t70 * mrSges(2,1) + t67 * mrSges(2,2) - m(4) * t81 - m(5) * (t40 + t81) + t41 * mrSges(5,2) + t73) * g(2) + m(5) * (t11 * t3 + t12 * t2 + t4 * t7 + t5 * t6) + m(4) * (t14 * t28 + t15 * t30 + t17 * t25 + t18 * t27) + (t67 * mrSges(2,1) + t41 * mrSges(5,1) + t70 * mrSges(2,2) + t42 * mrSges(5,2) + t71) * g(1) + Ifges(2,3) * qJDD(1) + (t90 * t85 + (-t60 * t66 - t61 * t79) * mrSges(3,2) + (-qJD(2) * t61 * t66 + t60 * t69) * mrSges(3,1) + (-g(2) * t70 + t34 * t69 + t35 * t66 + t85) * m(3)) * pkin(1); t72 + (-t31 * t54 + t91 * t55 + t76) * mrSges(5,2) + (-t24 * mrSges(4,1) + t26 * mrSges(4,2) + (mrSges(3,1) * t66 + mrSges(3,2) * t69) * t80) * t61 + (t29 * t54 + t92 * t55 + t86) * mrSges(5,1) + (-t60 * t88 - t15) * mrSges(4,2) + t73 * g(2) + t71 * g(1) + t60 * mrSges(4,1) * t87 + ((-t40 - t48) * g(2) + t2 * t31 + t3 * t29 - t91 * t7 + t92 * t6) * m(5) + ((t14 * t64 + t15 * t63) * pkin(2) - t48 * g(2) - t17 * t24 - t18 * t26) * m(4); t90 * (-g(3) + qJDD(3)); -g(2) * t38 + (t55 * t7 + t86) * mrSges(5,1) + (t55 * t6 + t76) * mrSges(5,2) + t84;];
tau  = t1;
