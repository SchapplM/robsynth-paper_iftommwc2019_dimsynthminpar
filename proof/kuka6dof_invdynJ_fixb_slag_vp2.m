% Calculate vector of inverse dynamics joint torques for
% kuka6dof
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-08-07 17:57
% Revision: 4fc915f170f947983f72d2b29299dd86a0cf3511
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = kuka6dof_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'kuka6dof_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'kuka6dof_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'kuka6dof_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'kuka6dof_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'kuka6dof_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'kuka6dof_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'kuka6dof_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'kuka6dof_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-08-07 17:44:34
% EndTime: 2018-08-07 17:45:39
% DurationCPUTime: 34.45s
% Computational Cost: add. (20850->999), mult. (44080->1373), div. (0->0), fcn. (34308->14), ass. (0->456)
t380 = sin(qJ(3));
t386 = cos(qJ(3));
t387 = cos(qJ(2));
t520 = qJD(1) * t387;
t381 = sin(qJ(2));
t521 = qJD(1) * t381;
t310 = t380 * t521 - t386 * t520;
t321 = t380 * t387 + t381 * t386;
t311 = t321 * qJD(1);
t243 = -pkin(4) * t310 + pkin(5) * t311;
t379 = sin(qJ(4));
t385 = cos(qJ(4));
t616 = pkin(3) * qJD(2);
t504 = t380 * t616;
t213 = t243 * t379 - t385 * t504;
t378 = sin(qJ(5));
t384 = cos(qJ(5));
t519 = qJD(2) * t386;
t503 = pkin(3) * t519;
t510 = qJD(5) * t385;
t511 = qJD(5) * t384;
t515 = qJD(4) * t384;
t701 = -t384 * t213 + t378 * t503 + pkin(4) * t511 + (-t378 * t510 - t379 * t515) * pkin(5);
t528 = t384 * t385;
t223 = t310 * t378 - t311 * t528;
t512 = qJD(5) * t379;
t477 = t378 * t512;
t759 = (t223 + t477) * pkin(6);
t505 = pkin(3) * t521;
t227 = t243 + t505;
t615 = pkin(3) * qJD(3);
t502 = t380 * t615;
t746 = -t227 * t385 + t379 * t502;
t629 = pkin(3) * t387;
t364 = pkin(2) + t629;
t340 = t364 * qJD(1);
t217 = pkin(4) * t311 + pkin(5) * t310 - t340;
t375 = qJD(2) + qJD(3);
t333 = pkin(5) * t375 + t503;
t181 = t217 * t379 + t333 * t385;
t332 = -pkin(4) * t375 + t504;
t161 = t181 * t378 + t332 * t384;
t518 = qJD(3) * t386;
t317 = (qJD(2) * t518 + qJDD(2) * t380) * pkin(3);
t374 = qJDD(2) + qJDD(3);
t293 = -pkin(4) * t374 + t317;
t508 = qJD(1) * qJD(2);
t324 = qJDD(1) * t387 - t381 * t508;
t325 = qJDD(1) * t381 + t387 * t508;
t533 = t380 * t381;
t418 = -t386 * t387 + t533;
t208 = qJD(1) * qJD(3) * t418 - t324 * t380 - t325 * t386;
t403 = t321 * qJD(3);
t209 = -qJD(1) * t403 + t324 * t386 - t325 * t380;
t302 = -qJDD(1) * pkin(2) - pkin(3) * t324;
t120 = -pkin(4) * t208 - pkin(5) * t209 + t302;
t630 = pkin(3) * t386;
t316 = -qJD(2) * t502 + qJDD(2) * t630;
t292 = pkin(5) * t374 + t316;
t180 = t385 * t217 - t333 * t379;
t517 = qJD(4) * t180;
t70 = t120 * t379 + t292 * t385 + t517;
t34 = -qJD(5) * t161 - t293 * t378 + t384 * t70;
t420 = t310 * t385 - t375 * t379;
t149 = qJD(4) * t420 - t209 * t379 + t385 * t374;
t146 = qJDD(5) + t149;
t260 = t310 * t379 + t385 * t375;
t148 = qJD(4) * t260 + t209 * t385 + t374 * t379;
t304 = qJD(4) + t311;
t201 = -t378 * t304 - t384 * t420;
t205 = qJDD(4) - t208;
t79 = -qJD(5) * t201 - t148 * t378 - t205 * t384;
t75 = qJDD(6) - t79;
t662 = t75 / 0.2e1;
t257 = qJD(5) + t260;
t377 = sin(qJ(6));
t383 = cos(qJ(6));
t156 = t201 * t383 + t257 * t377;
t200 = -t304 * t384 + t378 * t420;
t78 = qJD(5) * t200 + t148 * t384 - t205 * t378;
t30 = -qJD(6) * t156 + t146 * t383 - t377 * t78;
t669 = t30 / 0.2e1;
t155 = -t201 * t377 + t257 * t383;
t29 = qJD(6) * t155 + t146 * t377 + t383 * t78;
t670 = t29 / 0.2e1;
t24 = pkin(6) * t146 + t34;
t71 = -qJD(4) * t181 + t385 * t120 - t292 * t379;
t25 = -pkin(6) * t78 + t71;
t121 = -pkin(6) * t201 + t180;
t162 = t181 * t384 - t378 * t332;
t125 = pkin(6) * t257 + t162;
t59 = t121 * t383 - t125 * t377;
t3 = qJD(6) * t59 + t24 * t383 + t25 * t377;
t60 = t121 * t377 + t125 * t383;
t4 = -qJD(6) * t60 - t24 * t377 + t25 * t383;
t684 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t8 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t75;
t679 = t684 + Ifges(7,5) * t670 + Ifges(7,6) * t669 + Ifges(7,3) * t662 - t78 * Ifges(6,4) / 0.2e1 - t79 * Ifges(6,2) / 0.2e1 - t146 * Ifges(6,6) / 0.2e1 + t8 / 0.2e1;
t758 = -t34 * mrSges(6,3) + t679;
t376 = qJ(2) + qJ(3);
t369 = sin(t376);
t370 = cos(t376);
t291 = -t369 * t378 + t370 * t528;
t538 = t378 * t385;
t407 = t369 * t384 + t370 * t538;
t443 = mrSges(7,1) * t383 - mrSges(7,2) * t377;
t415 = mrSges(6,1) + t443;
t674 = pkin(6) * m(7);
t496 = mrSges(7,3) + t674;
t459 = mrSges(6,2) - t496;
t757 = t291 * t415 - t407 * t459;
t756 = -mrSges(6,3) - mrSges(5,2);
t755 = mrSges(7,1) * t377;
t442 = mrSges(7,2) * t383 + t755;
t414 = -mrSges(6,3) - t442;
t754 = mrSges(5,2) - t414;
t212 = t385 * t243 + t379 * t504;
t623 = pkin(6) * t384;
t487 = -pkin(5) - t623;
t514 = qJD(4) * t385;
t753 = t487 * t514 - t212 + t759;
t363 = pkin(5) + t630;
t462 = -t363 - t623;
t752 = t462 * t514 + t746 + t759;
t516 = qJD(4) * t379;
t500 = pkin(6) * t516;
t552 = t311 * t379;
t751 = pkin(6) * t552 + t500 - t701;
t581 = t257 * Ifges(6,3);
t584 = t201 * Ifges(6,5);
t585 = t200 * Ifges(6,6);
t113 = t581 + t584 + t585;
t577 = t304 * Ifges(5,6);
t580 = t260 * Ifges(5,2);
t607 = Ifges(5,4) * t420;
t166 = t577 + t580 - t607;
t707 = t166 + t113;
t632 = pkin(3) * t380;
t362 = -pkin(4) + t632;
t501 = pkin(3) * t518;
t400 = t363 * t510 + t501;
t457 = t385 * t502;
t698 = -t363 * t516 - t457;
t191 = -t400 * t378 + (-qJD(5) * t362 + t698) * t384;
t750 = (pkin(6) * t311 + t227 * t384) * t379 - t191 + t500;
t536 = t379 * t383;
t173 = -t223 * t377 + t311 * t536;
t455 = qJD(6) + t515;
t456 = -qJD(6) * t384 - qJD(4);
t216 = t456 * t536 + (-t385 * t455 + t477) * t377;
t749 = t173 - t216;
t541 = t377 * t379;
t174 = t223 * t383 + t311 * t541;
t513 = qJD(5) * t378;
t529 = t383 * t385;
t215 = t455 * t529 + (t377 * t456 - t383 * t513) * t379;
t748 = t174 - t215;
t323 = t378 * pkin(4) + pkin(5) * t528;
t482 = t378 * t516;
t747 = -pkin(5) * t482 + qJD(5) * t323 - t213 * t378 - t384 * t503;
t745 = t516 + t552;
t198 = qJD(6) - t200;
t648 = t198 / 0.2e1;
t651 = t156 / 0.2e1;
t653 = t155 / 0.2e1;
t683 = t59 * mrSges(7,1) - t60 * mrSges(7,2);
t681 = -t180 * mrSges(6,1) - t683;
t744 = Ifges(7,5) * t651 + Ifges(7,6) * t653 + Ifges(7,3) * t648 - t681;
t649 = -t198 / 0.2e1;
t652 = -t156 / 0.2e1;
t654 = -t155 / 0.2e1;
t592 = Ifges(7,3) * t198;
t593 = Ifges(7,6) * t155;
t596 = Ifges(7,5) * t156;
t61 = t592 + t593 + t596;
t667 = -t61 / 0.2e1;
t743 = Ifges(7,5) * t652 + Ifges(7,6) * t654 + Ifges(7,3) * t649 + t667 - t683;
t742 = t332 * mrSges(5,1) + t162 * mrSges(6,2) - t181 * mrSges(5,3);
t741 = (mrSges(4,2) - mrSges(5,3)) * t370;
t738 = t148 / 0.2e1;
t737 = t149 / 0.2e1;
t736 = t205 / 0.2e1;
t735 = -t381 / 0.2e1;
t733 = t181 * mrSges(5,2);
t718 = mrSges(6,3) * t201;
t710 = mrSges(6,1) * t257 + mrSges(7,1) * t155 - mrSges(7,2) * t156 - t718;
t525 = mrSges(5,1) * t304 - mrSges(6,1) * t200 + mrSges(6,2) * t201 + mrSges(5,3) * t420;
t35 = t378 * t70 - t332 * t513 + (qJD(5) * t181 + t293) * t384;
t575 = t35 * t378;
t408 = t161 * t511 + t575;
t732 = -t162 * t513 + t34 * t384 + t408;
t571 = t385 * t70;
t573 = t379 * t71;
t730 = t571 - t573;
t447 = mrSges(3,1) * t381 + mrSges(3,2) * t387;
t453 = -pkin(4) * t370 - pkin(5) * t369;
t613 = mrSges(5,1) * t385;
t631 = pkin(3) * t381;
t659 = -m(6) - m(5);
t694 = -m(7) + t659;
t728 = (t453 - t631) * t694 + m(4) * t631 + mrSges(5,3) * t369 + t447 + (mrSges(4,1) + t613) * t370;
t727 = t161 * mrSges(6,1) - t580 / 0.2e1 - t577 / 0.2e1 + t742;
t102 = -mrSges(5,2) * t205 + mrSges(5,3) * t149;
t486 = t180 * t514;
t409 = -t486 - t573;
t570 = -mrSges(5,1) * t205 + mrSges(6,1) * t79 - mrSges(6,2) * t78 + mrSges(5,3) * t148;
t726 = m(6) * t409 + t385 * t102 + t570 * t379;
t617 = -mrSges(6,1) * t146 - mrSges(7,1) * t30 + mrSges(7,2) * t29 + mrSges(6,3) * t78;
t282 = -t362 * t378 + t363 * t528;
t373 = t385 * pkin(6);
t268 = t282 + t373;
t305 = t462 * t379;
t206 = -t268 * t377 + t305 * t383;
t720 = qJD(6) * t206 + t377 * t752 - t383 * t750;
t207 = t268 * t383 + t305 * t377;
t719 = -qJD(6) * t207 + t377 * t750 + t383 * t752;
t717 = t180 * mrSges(6,2);
t715 = t332 * mrSges(5,2);
t714 = -t378 * t496 - t384 * t443 - mrSges(5,1);
t713 = -mrSges(5,2) - t442;
t307 = t373 + t323;
t326 = t487 * t379;
t239 = -t307 * t377 + t326 * t383;
t712 = qJD(6) * t239 + t377 * t753 - t383 * t751;
t240 = t307 * t383 + t326 * t377;
t711 = -qJD(6) * t240 + t377 * t751 + t383 * t753;
t708 = t385 * t525;
t706 = t746 * t180;
t554 = t260 * t384;
t177 = -t377 * t554 + t383 * t420;
t509 = qJD(6) * t378;
t396 = t377 * t511 + t383 * t509;
t705 = t177 - t396;
t530 = t383 * t384;
t178 = t260 * t530 + t377 * t420;
t397 = -t377 * t509 + t383 * t511;
t704 = t178 + t397;
t703 = -m(4) * t340 + mrSges(4,1) * t311 - mrSges(4,2) * t310;
t611 = mrSges(4,3) * t310;
t702 = -mrSges(4,1) * t375 - mrSges(5,1) * t260 - mrSges(5,2) * t420 - t611;
t498 = pkin(4) * m(7) + mrSges(4,1);
t417 = t498 + t613;
t497 = pkin(5) * m(7) + mrSges(5,3);
t700 = t369 * t497 + t370 * t417 + t659 * t453;
t322 = -t384 * pkin(4) + pkin(5) * t538;
t699 = t161 * t747 + t322 * t35;
t561 = t181 * t379;
t562 = t180 * t385;
t695 = qJD(4) * (-t561 - t562);
t609 = Ifges(3,4) * t381;
t693 = pkin(2) * t447 + (Ifges(3,1) * t387 - t609) * t735;
t426 = t377 * t60 + t383 * t59;
t428 = Ifges(7,5) * t383 - Ifges(7,6) * t377;
t599 = Ifges(7,4) * t383;
t433 = -Ifges(7,2) * t377 + t599;
t600 = Ifges(7,4) * t377;
t438 = Ifges(7,1) * t383 - t600;
t150 = Ifges(7,4) * t155;
t63 = Ifges(7,1) * t156 + Ifges(7,5) * t198 + t150;
t572 = t383 * t63;
t601 = Ifges(7,4) * t156;
t62 = Ifges(7,2) * t155 + Ifges(7,6) * t198 + t601;
t574 = t377 * t62;
t690 = -t426 * mrSges(7,3) + t433 * t653 + t438 * t651 + t428 * t648 - t574 / 0.2e1 + t572 / 0.2e1;
t688 = t71 * mrSges(5,1) - t70 * mrSges(5,2) + Ifges(5,5) * t148 + Ifges(5,6) * t149 + Ifges(5,3) * t205;
t544 = t370 * t384;
t288 = t369 * t538 - t544;
t545 = t370 * t378;
t289 = -t369 * t528 - t545;
t687 = -t289 * mrSges(6,1) - (-t289 * t377 + t369 * t536) * mrSges(7,2) - (t289 * t383 + t369 * t541) * mrSges(7,1) + t756 * t369 * t379 + t741 + (-mrSges(6,2) + mrSges(7,3)) * t288;
t388 = cos(qJ(1));
t534 = t379 * t388;
t546 = t369 * t388;
t686 = -t370 * t534 * t754 - mrSges(4,2) * t546 + t388 * t757;
t382 = sin(qJ(1));
t493 = t370 * t536;
t532 = t382 * t379;
t547 = t369 * t382;
t685 = -mrSges(4,2) * t547 + (-t755 + t756) * t370 * t532 + (-mrSges(7,2) * t493 + t757) * t382;
t668 = Ifges(5,1) * t738 + Ifges(5,4) * t737 + Ifges(5,5) * t736;
t594 = Ifges(6,6) * t257;
t604 = Ifges(6,4) * t201;
t114 = Ifges(6,2) * t200 + t594 + t604;
t197 = Ifges(6,4) * t200;
t597 = Ifges(6,5) * t257;
t115 = Ifges(6,1) * t201 + t197 + t597;
t165 = -Ifges(5,5) * t420 + Ifges(5,6) * t260 + Ifges(5,3) * t304;
t256 = Ifges(5,4) * t260;
t578 = t304 * Ifges(5,5);
t579 = t420 * Ifges(5,1);
t167 = t256 + t578 - t579;
t143 = Ifges(6,3) * t146;
t73 = Ifges(6,6) * t79;
t74 = Ifges(6,5) * t78;
t19 = t74 + t73 + t143;
t222 = -t310 * t384 - t311 * t538;
t608 = Ifges(4,4) * t310;
t228 = -Ifges(4,2) * t311 + Ifges(4,6) * t375 - t608;
t303 = Ifges(4,4) * t311;
t229 = -Ifges(4,1) * t310 + Ifges(4,5) * t375 - t303;
t313 = t377 * t385 + t379 * t530;
t481 = t378 * t514;
t398 = t379 * t511 + t481;
t399 = t384 * t514 - t477;
t535 = t379 * t384;
t406 = t377 * t535 - t529;
t429 = -Ifges(6,5) * t378 - Ifges(6,6) * t384;
t430 = Ifges(6,5) * t384 - Ifges(6,6) * t378;
t431 = Ifges(5,5) * t385 - Ifges(5,6) * t379;
t603 = Ifges(6,4) * t378;
t434 = -Ifges(6,2) * t384 - t603;
t602 = Ifges(6,4) * t384;
t435 = -Ifges(6,2) * t378 + t602;
t605 = Ifges(5,4) * t385;
t436 = -Ifges(5,2) * t379 + t605;
t439 = -Ifges(6,1) * t378 - t602;
t440 = Ifges(6,1) * t384 - t603;
t606 = Ifges(5,4) * t379;
t441 = Ifges(5,1) * t385 - t606;
t444 = mrSges(6,1) * t378 + mrSges(6,2) * t384;
t446 = mrSges(5,1) * t379 + mrSges(5,2) * t385;
t465 = t511 / 0.2e1;
t466 = -t511 / 0.2e1;
t467 = -t513 / 0.2e1;
t540 = t378 * t379;
t551 = t311 * t385;
t582 = t205 * Ifges(5,6);
t588 = t149 * Ifges(5,2);
t590 = t148 * Ifges(5,4);
t57 = t582 + t588 + t590;
t637 = -t310 / 0.2e1;
t638 = -t304 / 0.2e1;
t640 = t420 / 0.2e1;
t641 = -t260 / 0.2e1;
t642 = t257 / 0.2e1;
t643 = -t257 / 0.2e1;
t644 = t201 / 0.2e1;
t645 = -t201 / 0.2e1;
t646 = t200 / 0.2e1;
t647 = -t200 / 0.2e1;
t650 = -t167 / 0.2e1;
t655 = t146 / 0.2e1;
t656 = -t115 / 0.2e1;
t657 = t114 / 0.2e1;
t658 = -t114 / 0.2e1;
t660 = t79 / 0.2e1;
t661 = t78 / 0.2e1;
t663 = t63 / 0.2e1;
t664 = -t63 / 0.2e1;
t665 = t62 / 0.2e1;
t666 = -t62 / 0.2e1;
t21 = t78 * Ifges(6,1) + t79 * Ifges(6,4) + t146 * Ifges(6,5);
t671 = t21 / 0.2e1;
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t30 + Ifges(7,5) * t75;
t673 = t10 / 0.2e1;
t676 = Ifges(7,4) * t670 + Ifges(7,2) * t669 + Ifges(7,6) * t662;
t680 = (t115 * t384 + t378 * t61 + t167) * t514 / 0.2e1 + (t57 + t19) * t385 / 0.2e1 + (-Ifges(4,1) * t311 + t165 + t608) * t310 / 0.2e1 + (t260 * t436 + t304 * t431 - t420 * t441) * qJD(4) / 0.2e1 + (-mrSges(6,1) * t385 + mrSges(7,1) * t406 + mrSges(7,2) * t313 + mrSges(6,3) * t535) * t35 + (Ifges(7,5) * t313 - Ifges(7,6) * t406) * t662 + (Ifges(7,4) * t313 - Ifges(7,2) * t406) * t669 - t406 * t676 + (Ifges(7,1) * t313 - Ifges(7,4) * t406) * t670 + (Ifges(4,2) * t310 + t229 - t303) * t311 / 0.2e1 + (mrSges(5,1) * t310 - mrSges(6,1) * t222 - mrSges(6,2) * t223 - mrSges(5,3) * t551) * t180 + t340 * (-mrSges(4,1) * t310 - mrSges(4,2) * t311) + (-Ifges(5,3) * t310 - t311 * t431) * t638 + (-Ifges(5,5) * t310 - t311 * t441) * t640 + (-Ifges(5,6) * t310 - t311 * t436) * t641 - t375 * (-Ifges(4,5) * t311 + Ifges(4,6) * t310) / 0.2e1 + t758 * t540 + (-Ifges(6,4) * t645 - Ifges(6,2) * t647 - Ifges(6,6) * t643 + t657 + t743) * t222 + t744 * t398 + ((t222 - t398) * mrSges(6,3) + t745 * mrSges(6,2)) * t162 + ((-t223 + t399) * mrSges(6,3) - t748 * mrSges(7,2) + t749 * mrSges(7,1) + t745 * mrSges(6,1)) * t161 + (-t3 * t406 - t313 * t4 + t59 * t748 - t60 * t749) * mrSges(7,3) + t707 * (-t552 / 0.2e1 - t516 / 0.2e1) + (t180 * t399 - t34 * t385) * mrSges(6,2) + (Ifges(6,5) * t223 + Ifges(6,3) * t552) * t643 + (Ifges(7,5) * t215 + Ifges(7,6) * t216) * t648 + (Ifges(7,5) * t174 + Ifges(7,6) * t173) * t649 + (Ifges(7,4) * t215 + Ifges(7,2) * t216) * t653 + (Ifges(7,4) * t174 + Ifges(7,2) * t173) * t654 + (Ifges(6,4) * t223 + Ifges(6,6) * t552) * t647 - t310 * t733 + t293 * (mrSges(5,2) * t379 - t613) + t223 * t656 + t481 * t658 + (Ifges(6,6) * t385 + t379 * t435) * t660 + (Ifges(6,5) * t385 + t379 * t440) * t661 + t215 * t663 + t174 * t664 + t216 * t665 + t173 * t666 + t379 * t668 + t535 * t671 + t304 * t332 * t446 + t228 * t637 + (t429 * t512 + (-Ifges(6,3) * t379 + t385 * t430) * qJD(4)) * t642 + (t439 * t512 + (-Ifges(6,5) * t379 + t385 * t440) * qJD(4)) * t644 + (t434 * t512 + (-Ifges(6,6) * t379 + t385 * t435) * qJD(4)) * t646 - t551 * t650 + (Ifges(6,3) * t385 + t379 * t430) * t655 + t444 * t573 + t313 * t673 + Ifges(4,3) * t374 + (Ifges(7,1) * t215 + Ifges(7,4) * t216) * t651 + (Ifges(7,1) * t174 + Ifges(7,4) * t173) * t652 - t316 * mrSges(4,2) - t317 * mrSges(4,1) + Ifges(4,6) * t208 + Ifges(4,5) * t209 + t379 * t61 * t465 + t379 * t114 * t466 + t379 * t115 * t467 + (Ifges(5,5) * t379 + Ifges(5,6) * t385) * t736 + (Ifges(5,2) * t385 + t606) * t737 + (Ifges(5,1) * t379 + t605) * t738 + (-t181 * t552 + t730) * mrSges(5,3) + (Ifges(6,1) * t223 + Ifges(6,5) * t552) * t645;
t675 = pkin(2) * m(3);
t639 = -t420 / 0.2e1;
t633 = t387 / 0.2e1;
t626 = pkin(6) * t200;
t624 = pkin(6) * t288;
t621 = g(3) * t370;
t620 = t34 * mrSges(6,2);
t359 = t370 * pkin(5);
t614 = mrSges(4,1) * t369;
t610 = mrSges(4,3) * t311;
t598 = Ifges(3,5) * t387;
t595 = Ifges(3,6) * t381;
t252 = t375 * t533 + (-t518 - t519) * t387;
t253 = -qJD(2) * t321 - t403;
t171 = -pkin(4) * t252 - pkin(5) * t253 + t381 * t616;
t244 = pkin(4) * t321 + pkin(5) * t418 - t364;
t558 = t244 * t385;
t569 = t171 * t562 + t71 * t558;
t566 = t161 * t378;
t564 = t180 * t378;
t563 = t180 * t384;
t560 = t227 * t379;
t557 = t253 * t379;
t556 = t253 * t385;
t555 = t260 * t378;
t543 = t370 * t388;
t542 = t377 * t378;
t539 = t378 * t383;
t210 = -mrSges(5,2) * t304 + mrSges(5,3) * t260;
t537 = t379 * t210;
t531 = t382 * t385;
t526 = t385 * t388;
t523 = pkin(5) * t543 + t388 * t364;
t506 = -mrSges(3,3) - mrSges(4,3) + mrSges(2,2);
t499 = -mrSges(2,1) - t675;
t495 = t227 * t540;
t485 = t181 * t516;
t473 = t378 * t657;
t463 = -pkin(4) * t369 + t359;
t299 = t369 * t531 + t534;
t461 = t299 * t384 + t382 * t545;
t460 = t299 * t378 - t382 * t544;
t454 = t161 * t171 * t540;
t451 = g(1) * t382 - g(2) * t388;
t450 = -t418 * t510 - t252;
t116 = -t450 * t378 + (-qJD(5) * t321 + t418 * t516 + t556) * t384;
t449 = -qJD(6) * t379 * t418 - t116;
t448 = pkin(6) * t418 + t244 * t384;
t437 = Ifges(7,1) * t377 + t599;
t432 = Ifges(7,2) * t383 + t600;
t427 = Ifges(7,5) * t377 + Ifges(7,6) * t383;
t192 = -t362 * t513 - t363 * t482 - t378 * t457 + t384 * t400;
t281 = t362 * t384 + t363 * t538;
t425 = t161 * t192 + t281 * t35;
t422 = t162 * t384 + t566;
t238 = -t321 * t378 - t418 * t528;
t168 = -pkin(6) * t238 + t558;
t189 = t448 * t379;
t106 = t168 * t383 - t189 * t377;
t107 = t168 * t377 + t189 * t383;
t157 = -mrSges(6,2) * t257 + mrSges(6,3) * t200;
t395 = t157 * t384 - t378 * t710 + t210;
t394 = -qJD(6) * t238 + t418 * t514 - t557;
t301 = t369 * t526 - t532;
t250 = -t301 * t384 - t378 * t543;
t393 = -g(1) * t250 + g(2) * t461 - g(3) * t291 + t3 * t383 - t377 * t4;
t366 = Ifges(3,4) * t520;
t346 = pkin(4) * t547;
t339 = -mrSges(3,1) * t387 + mrSges(3,2) * t381;
t309 = Ifges(3,1) * t521 + Ifges(3,5) * qJD(2) + t366;
t308 = Ifges(3,6) * qJD(2) + (t387 * Ifges(3,2) + t609) * qJD(1);
t300 = t369 * t534 + t531;
t298 = t369 * t532 - t526;
t277 = -mrSges(4,2) * t375 - t610;
t249 = -t301 * t378 + t384 * t543;
t194 = t250 * t383 + t300 * t377;
t193 = -t250 * t377 + t300 * t383;
t186 = t238 * t383 + t418 * t541;
t185 = -t238 * t377 + t418 * t536;
t142 = -pkin(6) * t554 - t181;
t138 = pkin(6) * t420 + t563;
t126 = t180 * t212;
t109 = -t161 * t383 - t377 * t626;
t108 = t161 * t377 - t383 * t626;
t104 = mrSges(7,1) * t198 - mrSges(7,3) * t156;
t103 = -mrSges(7,2) * t198 + mrSges(7,3) * t155;
t82 = -mrSges(5,1) * t149 + mrSges(5,2) * t148;
t81 = t138 * t383 + t142 * t377;
t80 = -t138 * t377 + t142 * t383;
t67 = t448 * t514 + (-pkin(6) * t253 + t171 * t384 - t244 * t513) * t379;
t51 = -pkin(6) * t116 + t171 * t385 - t244 * t516;
t48 = t377 * t449 + t383 * t394;
t47 = t377 * t394 - t383 * t449;
t41 = -mrSges(6,2) * t146 + mrSges(6,3) * t79;
t15 = -mrSges(7,2) * t75 + mrSges(7,3) * t30;
t14 = mrSges(7,1) * t75 - mrSges(7,3) * t29;
t13 = -qJD(6) * t107 - t377 * t67 + t383 * t51;
t12 = qJD(6) * t106 + t377 * t51 + t383 * t67;
t1 = [t556 * t715 + (-mrSges(6,1) * t250 + t301 * mrSges(5,1) - m(7) * t523 - t194 * mrSges(7,1) - t193 * mrSges(7,2) + t756 * t300 + t659 * (-pkin(4) * t546 + t523) + t459 * t249 + t506 * t382 + (-m(4) * t364 + t498 * t369 + t499 + t741) * t388) * g(2) + (Ifges(7,5) * t186 + Ifges(7,6) * t185) * t662 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t648 + (Ifges(7,4) * t186 + Ifges(7,2) * t185) * t669 + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t653 + (t71 * mrSges(6,1) - Ifges(6,4) * t661 - Ifges(6,2) * t660 - Ifges(6,6) * t655 + t758) * (t321 * t384 - t418 * t538) + (-t162 * mrSges(6,3) + t658 - Ifges(6,6) * t642 - Ifges(6,4) * t644 - Ifges(6,2) * t646 + t61 / 0.2e1 + t744) * (t253 * t538 - t321 * t513 + t384 * t450 + t418 * t482) + t304 * (Ifges(5,5) * t556 - Ifges(5,3) * t252) / 0.2e1 + (Ifges(7,1) * t186 + Ifges(7,4) * t185) * t670 + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t651 + (Ifges(5,1) * t556 - Ifges(5,5) * t252) * t639 + (t185 * t3 - t186 * t4 - t47 * t59 + t48 * t60) * mrSges(7,3) + (t302 * mrSges(4,1) - t316 * mrSges(4,3) - Ifges(4,4) * t209 - Ifges(4,2) * t208 - Ifges(4,6) * t374 + t688) * t321 + (t161 * mrSges(6,3) + t717 + Ifges(6,5) * t642 + Ifges(6,1) * t644 + Ifges(6,4) * t646 + t115 / 0.2e1) * t116 + (mrSges(3,1) * t324 - mrSges(3,2) * t325 + (-t339 + t675) * qJDD(1)) * pkin(2) + (-m(4) * t302 + mrSges(4,1) * t208 - mrSges(4,2) * t209) * t364 + (mrSges(6,2) * t71 + mrSges(6,3) * t35 + Ifges(6,1) * t661 + Ifges(6,4) * t660 + Ifges(6,5) * t655 + t671) * t238 + m(6) * (t162 * t171 * t535 + t454 + t569) + m(5) * (t171 * t561 + t569) + t260 * (Ifges(5,4) * t556 - Ifges(5,6) * t252) / 0.2e1 + t167 * t556 / 0.2e1 + (t309 * t633 + t308 * t735 + (t598 / 0.2e1 - t595 / 0.2e1) * qJD(2) + ((Ifges(3,4) * t387 - Ifges(3,2) * t381) * t633 - t693) * qJD(1) + (t703 * t381 + (t252 * t386 + t253 * t380) * mrSges(4,3)) * pkin(3)) * qJD(2) + (-mrSges(5,1) * t252 - mrSges(5,3) * t556) * t180 + t47 * t663 + t48 * t665 + (Ifges(4,1) * t253 + Ifges(4,4) * t252) * t637 + (mrSges(3,1) * t451 + Ifges(3,4) * t325 + Ifges(3,2) * t324 + Ifges(3,6) * qJDD(2)) * t387 + (-mrSges(3,2) * t451 + Ifges(3,1) * t325 + Ifges(3,4) * t324 + Ifges(3,5) * qJDD(2)) * t381 + (t379 * t395 + t708) * t171 + m(7) * (t106 * t4 + t107 * t3 + t12 * t60 + t13 * t59 + t454) + (-t707 / 0.2e1 - Ifges(5,4) * t639 - Ifges(6,5) * t644 - Ifges(6,6) * t646 - Ifges(6,3) * t642 + t727) * t557 - (t302 * mrSges(4,2) + t317 * mrSges(4,3) + Ifges(4,1) * t209 + Ifges(4,4) * t208 + Ifges(4,5) * t374 + (t293 * mrSges(5,2) - t71 * mrSges(5,3) + 0.2e1 * t668) * t385 + (-t19 / 0.2e1 - t57 / 0.2e1 + t293 * mrSges(5,1) - t582 / 0.2e1 - t590 / 0.2e1 - t588 / 0.2e1 - t73 / 0.2e1 + t35 * mrSges(6,1) + t620 - t143 / 0.2e1 - t74 / 0.2e1 - t70 * mrSges(5,3)) * t379 + ((-t715 - t578 / 0.2e1 - t256 / 0.2e1 + t579 / 0.2e1 + t650 + t180 * mrSges(5,3)) * t379 + (t607 / 0.2e1 - t166 / 0.2e1 - t113 / 0.2e1 - t581 / 0.2e1 - t584 / 0.2e1 - t585 / 0.2e1 + t727) * t385) * qJD(4)) * t418 + t186 * t673 + t185 * t676 + t375 * (Ifges(4,5) * t253 + Ifges(4,6) * t252) / 0.2e1 + Ifges(2,3) * qJDD(1) - t340 * (-mrSges(4,1) * t252 + mrSges(4,2) * t253) - t311 * (Ifges(4,4) * t253 + Ifges(4,2) * t252) / 0.2e1 + t252 * t228 / 0.2e1 - t252 * t165 / 0.2e1 + t253 * t229 / 0.2e1 + t35 * (-mrSges(7,1) * t185 + mrSges(7,2) * t186) + t161 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + t107 * t15 + t12 * t103 + t13 * t104 + t106 * t14 + (-m(7) * t346 - mrSges(5,1) * t299 - t415 * t461 + t754 * t298 + t659 * (t346 + (-t364 - t359) * t382) + t459 * t460 + t506 * t388 + (-t614 + (m(7) + m(4)) * t364 + (-mrSges(4,2) + t497) * t370 - t499) * t382) * g(1) + t252 * t733 + (((m(5) * t181 + m(6) * t422 + m(7) * t566 + t395) * qJD(4) - t570) * t385 + (t384 * t41 + t102 + t617 * t378 - t525 * qJD(4) + (-t378 * t157 - t384 * t710) * qJD(5) + m(7) * t408 + m(6) * (-t517 + t732) + m(5) * (t70 - t517)) * t379) * t244; (-t485 - t486) * mrSges(5,3) - (-Ifges(3,2) * t521 + t309 + t366) * t520 / 0.2e1 + t525 * (-t363 * t514 + t746) + (-t181 * t560 + t293 * t362 + (-t181 * t380 * t385 + t332 * t386) * t615 + t706) * m(5) + (t162 * t191 + t282 * t34 - t422 * t560 + t425 + t706) * m(6) + (-t227 * t535 + t191) * t157 - t227 * t537 - t503 * t611 + t698 * t210 - t710 * (-t495 + t192) - (-t595 + t598) * t508 / 0.2e1 - (mrSges(4,1) * t374 - mrSges(4,3) * t209) * t632 + t680 + t719 * t104 + t720 * t103 + (-t161 * t495 + t206 * t4 + t207 * t3 + t59 * t719 + t60 * t720 + t425) * m(7) + t617 * t281 + t702 * t501 - t703 * t505 - t277 * t502 + m(4) * (t316 * t386 + t317 * t380) * pkin(3) + (-mrSges(4,2) * t374 + mrSges(4,3) * t208) * t630 + t308 * t521 / 0.2e1 + t362 * t82 + Ifges(3,6) * t324 + Ifges(3,5) * t325 + t282 * t41 + t504 * t610 + t206 * t14 + t207 * t15 + Ifges(3,3) * qJDD(2) + t693 * qJD(1) ^ 2 + (-m(4) * t629 + m(7) * t624 + t369 * t613 + t339 + t614 + t694 * (t463 + t629) + t687) * g(3) + (t382 * t728 + t685) * g(2) + (t388 * t728 + t686) * g(1) + ((t695 + t730) * m(5) + t726) * t363; ((t277 + t610) * t380 + (-t611 - t702) * t386) * t616 - t710 * t747 + (t382 * t700 + t685) * g(2) + (t388 * t700 + t686) * g(1) + (t162 * t701 + t323 * t34 - t126 + t699) * m(6) + ((-t359 + t624) * g(3) + t239 * t4 + t240 * t3 + t712 * t60 + t711 * t59 + t699) * m(7) + t701 * t157 + (m(5) * (t409 - t485 + t571) + (-t537 - t708) * qJD(4) + t726) * pkin(5) + t680 - m(5) * (t181 * t213 + t332 * t503 + t126) + (t417 * t369 + t659 * t463 + t687) * g(3) + t712 * t103 + (-m(5) * t293 - t82) * pkin(4) + t711 * t104 - t525 * t212 + t617 * t322 + t323 * t41 + t239 * t14 + t240 * t15 - t213 * t210 + mrSges(5,3) * t695; (Ifges(7,5) * t178 + Ifges(7,6) * t177) * t649 + (-m(7) * t564 + mrSges(6,1) * t420 + mrSges(7,1) * t705 - mrSges(7,2) * t704) * t161 + (t427 * t648 + t432 * t653 + t437 * t651) * t509 - m(7) * (t59 * t80 + t60 * t81) + (-g(1) * t300 - g(2) * t298 + t379 * t621 + t71) * (mrSges(6,1) * t384 - mrSges(6,2) * t378) + (t377 * t63 + t383 * t62) * t509 / 0.2e1 - (t200 * t435 + t201 * t440 + t257 * t430) * qJD(5) / 0.2e1 - (-Ifges(6,5) * t645 - Ifges(5,2) * t641 - Ifges(5,6) * t638 - Ifges(6,6) * t647 - Ifges(6,3) * t643 - t742) * t420 + t743 * t555 + (Ifges(7,4) * t178 + Ifges(7,2) * t177) * t654 + (t256 + t167) * t641 + (-t21 / 0.2e1 - t428 * t662 - t433 * t669 - t438 * t670 + (t3 * t377 + t383 * t4 + (-t377 * t59 + t383 * t60) * qJD(6)) * t674) * t378 + (t103 * t396 + t104 * t397 + t14 * t539 + t15 * t542) * pkin(6) + t607 * t640 + (mrSges(7,1) * t313 - mrSges(7,2) * t406 + t496 * t540 + t446) * t621 + t710 * t564 - t683 * t513 + t426 * t511 * t674 + (t714 * t300 + t713 * t301) * g(1) + (t298 * t714 + t299 * t713) * g(2) + (Ifges(5,1) * t640 + Ifges(5,5) * t638 + t430 * t643 + t435 * t647 + t440 * t645 + t473 - t715) * t260 - t442 * t575 + (t3 * t542 + t4 * t539 + t59 * t704 - t60 * t705) * mrSges(7,3) + t707 * t639 + t688 - t10 * t539 / 0.2e1 + t554 * t656 + t434 * t660 + t439 * t661 + t178 * t664 + t177 * t666 + (t473 + (-Ifges(7,3) * t378 - t384 * t428) * t648 + (-Ifges(7,5) * t378 - t384 * t438) * t651 + (-Ifges(7,6) * t378 - t384 * t433) * t653) * qJD(5) + ((-t444 + mrSges(5,3)) * t260 - t444 * qJD(5) - m(6) * (-t181 + t422) - t210) * t180 + t525 * t181 - t157 * t563 + t429 * t655 + t679 * t384 + (t572 + t115) * t466 + t61 * t467 + (Ifges(7,1) * t178 + Ifges(7,4) * t177) * t652 + t542 * t676 - t81 * t103 - t80 * t104 + (-g(1) * t301 - g(2) * t299 - t161 * t554 + t162 * t555 + t385 * t621 - t732) * mrSges(6,3) + t465 * t574; t377 * t673 + (t594 / 0.2e1 - t592 / 0.2e1 - t593 / 0.2e1 - t596 / 0.2e1 + t667 + t657 + t604 / 0.2e1 + t681) * t201 + t393 * mrSges(7,3) + (m(7) * t393 - t377 * t14 + t383 * t15) * pkin(6) + t19 + (t161 * t442 + (-m(7) * t426 - t377 * t103 - t383 * t104) * pkin(6) + t690) * qJD(6) + (-t597 / 0.2e1 - t197 / 0.2e1 - t717 + t656 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t201 + t414 * t161 - t690) * t200 + t427 * t662 + t432 * t669 + t437 * t670 + (-mrSges(6,2) * t461 - t415 * t460) * g(2) + t161 * t157 - t108 * t104 - t109 * t103 - m(7) * (t108 * t59 + t109 * t60) + (mrSges(6,2) * t250 + t249 * t415) * g(1) - t415 * t35 + (mrSges(6,2) * t291 + t407 * t415) * g(3) + t383 * t676 - t620 + (-m(7) * t161 + t710 + t718) * t162; -t161 * (mrSges(7,1) * t156 + mrSges(7,2) * t155) + (Ifges(7,1) * t155 - t601) * t652 + t62 * t651 + (Ifges(7,5) * t155 - Ifges(7,6) * t156) * t649 - t59 * t103 + t60 * t104 - g(1) * (mrSges(7,1) * t193 - mrSges(7,2) * t194) - g(2) * ((t298 * t383 + t377 * t461) * mrSges(7,1) + (-t298 * t377 + t383 * t461) * mrSges(7,2)) - g(3) * ((-t291 * t377 - t493) * mrSges(7,1) + (-t291 * t383 + t370 * t541) * mrSges(7,2)) + (t155 * t59 + t156 * t60) * mrSges(7,3) + t8 + (-Ifges(7,2) * t156 + t150 + t63) * t654 + t684;];
tau  = t1;
