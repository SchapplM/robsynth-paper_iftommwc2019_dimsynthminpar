# Review 1 (25.02.)

This paper presents a method to reduce the cost of computations repeated in the optimal design including structural and dimensional analysis of serial robots. Focusing on the linearity of the inertial parameters in the equation of motion, cost reduction of computations for changing the values of inertial parameters was achieved. Demonstrative examples for serial robots show the effectiveness of the proposed method.

The reviewer suggests the authors to revise the paper according to the following points:
* [_] If the proposed method is applicable only to serial robots, this fact and the reason of this limitation should be discussed or mentioned in the earlier place than the end of chapter 3.
  * Es geht prinzipiell auch für PKM, aber der volle Vereinfachungsgrad (über die Dreieckform) geht dort nicht. Weiter ausführen.

* [_] It is hard to understand the correspondence of Eq. (4) to \\phi_{tau} and \\phi_{w}.
  * Beispiel nennen: Industrieroboter?
  * Dimensionen hinschreiben?

# Review 2 (26.02.)

This paper proposes an improved method for the calculation of robots inverse dynamic models. The approach is based on the linearization of dynamics equations using regressor matrix. The paper proposes the exploitation of the linear formulation of dynamics problem to a combined procedure of structural-dimensional synthesis. The efficiency of the linear dynamics equations was demonstrated through 3 examples of serial robots

* [x] Page 4, Line 2: which are listed in [15] ?
  * Bezieht sich sprachlich auf "are". Verweis auf [15] sprachlich vorziehen. Sieht auch blöd aus, soll ja nur ein kurzer Verweis werden. Lasse das "are".
* [N] Page 4, Line 9: incomplete sentence
  * Dritter Punkt der Aufzählung verbindet zwei Sätze. Das wurde nicht richtig erkannt. Von Kollegen prüfen lassen.
* [x] Page 7, Section 5: There is no any details about robots’ dynamic equations, this part needs more explanation. 
  * Hinweis: Matlab-Implementierung, einzelne Funktionsdatei, die mex-kompilierbar ist und für Echtzeit-Implementierung (SPS, Simulink-basierter Prüfstand) geeignet ist. Automatisch generierter Code ist sehr umfangreich, daher keine Darstellung (und weil außerhalb des Fokus des Beitrags)
  * Mehr Details zum Rechenweg

# Review 3 (28.02.)

The article presents a method to reduce computational effort in dynamic evaluations by exploding the properties of the dynamics regressor form in design optimization. A method not fully exploited by the community.
The topic is very interesting and application to some robots are described in results.

* [_] Any way some important comparative results are missing. Some simulation times and comparisons should be included showing the real advantage of the proposed method.
  * Simulative Beispielrechnung geht auch ohne implementierte Maßsynthese. Beispiel zeigen.
