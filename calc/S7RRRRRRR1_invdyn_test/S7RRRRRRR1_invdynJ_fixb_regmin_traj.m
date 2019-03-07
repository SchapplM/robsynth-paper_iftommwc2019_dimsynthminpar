function RMV_Traj = S7RRRRRRR1_invdynJ_fixb_regmin_traj(Q, QD, QDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,7]),
%$cgargs  coder.newtype('double',[inf,7]),
%$cgargs  coder.newtype('double',[inf,7]),
%$cgargs  zeros(3,1), zeros(4,1)}

%% Trajektorie der Regressor-Vektoren aufbauen
RMV_Traj = NaN(size(Q,1), 191);

for ii = 1:size(Q,1)
  RMV_Traj(ii,:) = S7RRRRRRR1_regmat2regmatvector( ...
    S7RRRRRRR1_invdynJ_fixb_regmin_slag_vp(Q(ii,:)', QD(ii,:)', QDD(ii,:)', g, pkin) );
end