% Calculate vector of inverse dynamics joint torques for
% S3PPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% MDP [5x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S3PPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-28 11:56
% Revision: 2bf3b907e1213de0593c9d1d0a7eb98ef6ddbfca (2019-02-28)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S7RRRRRRR1_invdynJ_fixb_mdp2_slag_vp(RMV, MDP)
%% Coder Information
%#codegen
% %$cgargs {zeros(191,1), zeros(5,1)}
% assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
%   'S3PPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [5x6] (double)'); 
MDP1 = MDP(1);
MDP2 = MDP(2);
MDP3 = MDP(3);
MDP4 = MDP(4);
MDP5 = MDP(5);
MDP6 = MDP(6);
MDP7 = MDP(7);
MDP8 = MDP(8);
MDP9 = MDP(9);
MDP10 = MDP(10);
MDP11 = MDP(11);
MDP12 = MDP(12);
MDP13 = MDP(13);
MDP14 = MDP(14);
MDP15 = MDP(15);
MDP16 = MDP(16);
MDP17 = MDP(17);
MDP18 = MDP(18);
MDP19 = MDP(19);
MDP20 = MDP(20);
MDP21 = MDP(21);
MDP22 = MDP(22);
MDP23 = MDP(23);
MDP24 = MDP(24);
MDP25 = MDP(25);
MDP26 = MDP(26);
MDP27 = MDP(27);
MDP28 = MDP(28);
MDP29 = MDP(29);
MDP30 = MDP(30);
MDP31 = MDP(31);
MDP32 = MDP(32);
MDP33 = MDP(33);
MDP34 = MDP(34);
MDP35 = MDP(35);
MDP36 = MDP(36);
MDP37 = MDP(37);
MDP38 = MDP(38);
MDP39 = MDP(39);
MDP40 = MDP(40);
MDP41 = MDP(41);
MDP42 = MDP(42);
MDP43 = MDP(43);
MDP44 = MDP(44);
MDP45 = MDP(45);
tau = [RMV(185) * MDP45 + RMV(178) * MDP44 + RMV(164) * MDP42 + RMV(171) * MDP43 + RMV(150) * MDP40 + RMV(157) * MDP41 + RMV(137) * MDP38 + RMV(143) * MDP39 + RMV(125) * MDP36 + RMV(131) * MDP37 + RMV(119) * MDP35 + RMV(107) * MDP33 + RMV(113) * MDP34 + RMV(96) * MDP31 + RMV(101) * MDP32 + RMV(81) * MDP28 + RMV(86) * MDP29 + RMV(91) * MDP30 + RMV(71) * MDP26 + RMV(76) * MDP27 + RMV(58) * MDP23 + RMV(62) * MDP24 + RMV(66) * MDP25 + RMV(42) * MDP19 + RMV(46) * MDP20 + RMV(50) * MDP21 + RMV(54) * MDP22 + RMV(26) * MDP14 + RMV(29) * MDP15 + RMV(32) * MDP16 + RMV(35) * MDP17 + RMV(38) * MDP18 + RMV(13) * MDP9 + RMV(15) * MDP10 + RMV(17) * MDP11 + RMV(20) * MDP12 + RMV(23) * MDP13 + RMV(1) * MDP1 + RMV(2) * MDP2 + RMV(3) * MDP3 + RMV(4) * MDP4 + RMV(6) * MDP5 + RMV(8) * MDP6 + RMV(10) * MDP7; RMV(186) * MDP45 + RMV(179) * MDP44 + RMV(165) * MDP42 + RMV(172) * MDP43 + RMV(151) * MDP40 + RMV(158) * MDP41 + RMV(138) * MDP38 + RMV(144) * MDP39 + RMV(126) * MDP36 + RMV(132) * MDP37 + RMV(120) * MDP35 + RMV(108) * MDP33 + RMV(114) * MDP34 + RMV(97) * MDP31 + RMV(102) * MDP32 + RMV(82) * MDP28 + RMV(87) * MDP29 + RMV(92) * MDP30 + RMV(72) * MDP26 + RMV(77) * MDP27 + RMV(59) * MDP23 + RMV(63) * MDP24 + RMV(67) * MDP25 + RMV(43) * MDP19 + RMV(47) * MDP20 + RMV(51) * MDP21 + RMV(55) * MDP22 + RMV(27) * MDP14 + RMV(30) * MDP15 + RMV(33) * MDP16 + RMV(36) * MDP17 + RMV(39) * MDP18 + RMV(11) * MDP7 + RMV(12) * MDP8 + RMV(14) * MDP9 + RMV(16) * MDP10 + RMV(18) * MDP11 + RMV(21) * MDP12 + RMV(24) * MDP13 + RMV(5) * MDP4 + RMV(7) * MDP5 + RMV(9) * MDP6; RMV(180) * MDP44 + RMV(187) * MDP45 + RMV(166) * MDP42 + RMV(173) * MDP43 + RMV(152) * MDP40 + RMV(159) * MDP41 + RMV(139) * MDP38 + RMV(145) * MDP39 + RMV(127) * MDP36 + RMV(133) * MDP37 + RMV(121) * MDP35 + RMV(109) * MDP33 + RMV(115) * MDP34 + RMV(98) * MDP31 + RMV(103) * MDP32 + RMV(83) * MDP28 + RMV(88) * MDP29 + RMV(93) * MDP30 + RMV(68) * MDP25 + RMV(73) * MDP26 + RMV(78) * MDP27 + RMV(56) * MDP22 + RMV(60) * MDP23 + RMV(64) * MDP24 + RMV(44) * MDP19 + RMV(48) * MDP20 + RMV(52) * MDP21 + RMV(28) * MDP14 + RMV(31) * MDP15 + RMV(34) * MDP16 + RMV(37) * MDP17 + RMV(40) * MDP18 + RMV(19) * MDP11 + RMV(22) * MDP12 + RMV(25) * MDP13; RMV(181) * MDP44 + RMV(188) * MDP45 + RMV(167) * MDP42 + RMV(174) * MDP43 + RMV(153) * MDP40 + RMV(160) * MDP41 + RMV(140) * MDP38 + RMV(146) * MDP39 + RMV(128) * MDP36 + RMV(134) * MDP37 + RMV(116) * MDP34 + RMV(122) * MDP35 + RMV(110) * MDP33 + RMV(94) * MDP30 + RMV(99) * MDP31 + RMV(104) * MDP32 + RMV(79) * MDP27 + RMV(84) * MDP28 + RMV(89) * MDP29 + RMV(69) * MDP25 + RMV(74) * MDP26 + RMV(57) * MDP22 + RMV(61) * MDP23 + RMV(65) * MDP24 + RMV(41) * MDP18 + RMV(45) * MDP19 + RMV(49) * MDP20 + RMV(53) * MDP21; RMV(70) * MDP25 + RMV(75) * MDP26 + RMV(80) * MDP27 + RMV(85) * MDP28 + RMV(90) * MDP29 + RMV(95) * MDP30 + RMV(100) * MDP31 + RMV(105) * MDP32 + RMV(111) * MDP33 + RMV(117) * MDP34 + RMV(123) * MDP35 + RMV(129) * MDP36 + RMV(135) * MDP37 + RMV(141) * MDP38 + RMV(147) * MDP39 + RMV(154) * MDP40 + RMV(161) * MDP41 + RMV(168) * MDP42 + RMV(175) * MDP43 + RMV(182) * MDP44 + RMV(189) * MDP45; RMV(106) * MDP32 + RMV(112) * MDP33 + RMV(118) * MDP34 + RMV(124) * MDP35 + RMV(130) * MDP36 + RMV(136) * MDP37 + RMV(142) * MDP38 + RMV(148) * MDP39 + RMV(155) * MDP40 + RMV(162) * MDP41 + RMV(169) * MDP42 + RMV(176) * MDP43 + RMV(183) * MDP44 + RMV(190) * MDP45; RMV(149) * MDP39 + RMV(156) * MDP40 + RMV(163) * MDP41 + RMV(170) * MDP42 + RMV(177) * MDP43 + RMV(184) * MDP44 + RMV(191) * MDP45;];