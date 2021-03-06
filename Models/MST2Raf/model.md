% latex table generated in R 3.5.1 by xtable 1.8-2 package
% Wed Aug 29 10:42:34 2018
\begin{table}[ht]
\centering
\begin{tabular}{ll}
  \hline
a & b \\ 
  \hline
\$d/dt[ pRaf1]\$ & \$ -1*((k1*RasGTP*pRaf1)/(Km1 + pRaf1)) +1*(  (k2a*LATS1a*Raf1)/(Km2a + Raf1) +   (k2b*Kin*Raf1)/(Km2b + Raf1)) -1*(  -(k5r*pMpR) + k5f*pMST2*pRaf1) -1*(  -(k6r*MpR) + k6f*MST2*pRaf1) \$ \\ 
  \$d/dt[ Raf1a]\$ & \$ 1*( (V3*(1 + (Fa*ppERK)/Ka)*Raf1)/  ((1 + ppERK/Ka)*(Km3 + Raf1))) -1*( (V4*Raf1a)/(Km4 + Raf1a)) -1*(  k16af*MEK*Raf1a - k16ar*RaMk) +1*( k16b*RaMk) -1*(  k18af*pMEK*Raf1a - k18ar*RapMk) +1*( k18b*RapMk) \$ \\ 
  \$d/dt[ Raf1]\$ & \$ 1*((k1*RasGTP*pRaf1)/(Km1 + pRaf1)) -1*(  (k2a*LATS1a*Raf1)/(Km2a + Raf1) +   (k2b*Kin*Raf1)/(Km2b + Raf1)) -1*( (V3*(1 + (Fa*ppERK)/Ka)*Raf1)/  ((1 + ppERK/Ka)*(Km3 + Raf1))) +1*( (V4*Raf1a)/(Km4 + Raf1a)) \$ \\ 
  \$d/dt[ pMST2]\$ & \$ -1*(  -(k5r*pMpR) + k5f*pMST2*pRaf1) +1*( (AKTa*k7*(1 + Kact*RasGTP)*MST2)/  (Km7 + MST2)) -1*( (k8*PP2A*pMST2)/(Km8 + pMST2)) \$ \\ 
  \$d/dt[ MST2]\$ & \$ -1*(  -(k6r*MpR) + k6f*MST2*pRaf1) -1*( (AKTa*k7*(1 + Kact*RasGTP)*MST2)/  (Km7 + MST2)) +1*( (k8*PP2A*pMST2)/(Km8 + pMST2)) -2*( k9*MST2*MST2) +1*(  (k10*PP2A*MST2a)/(Km10 + MST2a)) -1*( -(k12r*MF1A) + k12f*RASSF1A*MST2) \$ \\ 
  \$d/dt[ MST2a]\$ & \$ 2*( k9*MST2*MST2) -1*(  (k10*PP2A*MST2a)/(Km10 + MST2a)) -1*( -(k11r*MaF1A) +   k11f*RASSF1A*MST2a) \$ \\ 
  \$d/dt[ MaF1A]\$ & \$ 1*( -(k11r*MaF1A) +   k11f*RASSF1A*MST2a) +1*(  (V13*MF1A)/(Km13 + MF1A)) \$ \\ 
  \$d/dt[ MF1A]\$ & \$ 1*( -(k12r*MF1A) + k12f*RASSF1A*MST2) -1*(  (V13*MF1A)/(Km13 + MF1A)) \$ \\ 
  \$d/dt[ RASSF1A]\$ & \$ -1*( -(k11r*MaF1A) +   k11f*RASSF1A*MST2a) -1*( -(k12r*MF1A) + k12f*RASSF1A*MST2) \$ \\ 
  \$d/dt[ LATS1]\$ & \$ -1*(  (k14b*LATS1*MaF1A)/(Km14b + LATS1) + (k14a*LATS1*MST2a)/   (Km14a + LATS1)) +1*( (V15*LATS1a)/(Km15 + LATS1a)) \$ \\ 
  \$d/dt[ LATS1a]\$ & \$ 1*(  (k14b*LATS1*MaF1A)/(Km14b + LATS1) + (k14a*LATS1*MST2a)/   (Km14a + LATS1)) -1*( (V15*LATS1a)/(Km15 + LATS1a)) \$ \\ 
  \$d/dt[ MEK]\$ & \$ -1*(  k16af*MEK*Raf1a - k16ar*RaMk) +1*(  (V17*pMEK)/(Km17 + pMEK + (Km17*ppMEK)/Km19)) \$ \\ 
  \$d/dt[ pMEK]\$ & \$ 1*( k16b*RaMk) -1*(  (V17*pMEK)/(Km17 + pMEK + (Km17*ppMEK)/Km19)) -1*(  k18af*pMEK*Raf1a - k18ar*RapMk) +1*(  (V19*ppMEK)/(Km19 + (Km19*pMEK)/Km17 + ppMEK)) \$ \\ 
  \$d/dt[ ppMEK]\$ & \$ 1*( k18b*RapMk) -1*(  (V19*ppMEK)/(Km19 + (Km19*pMEK)/Km17 + ppMEK)) \$ \\ 
  \$d/dt[ ERK]\$ & \$ -1*(  (k20*ERK*ppMEK)/(Km20 + ERK + (Km20*pERK)/Km22)) +1*(  (V21*pERK)/(Km21 + (Km21*ERK)/Ki + pERK + (Km21*ppERK)/Km23)) \$ \\ 
  \$d/dt[ pERK]\$ & \$ 1*(  (k20*ERK*ppMEK)/(Km20 + ERK + (Km20*pERK)/Km22)) -1*(  (V21*pERK)/(Km21 + (Km21*ERK)/Ki + pERK + (Km21*ppERK)/Km23)) -1*(  (k22*pERK*ppMEK)/(Km22 + (Km22*ERK)/Km20 + pERK)) +1*(  (V23*ppERK)/(Km23 + (Km23*ERK)/Ki + (Km23*pERK)/Km21 + ppERK)) \$ \\ 
  \$d/dt[ ppERK]\$ & \$ 1*(  (k22*pERK*ppMEK)/(Km22 + (Km22*ERK)/Km20 + pERK)) -1*(  (V23*ppERK)/(Km23 + (Km23*ERK)/Ki + (Km23*pERK)/Km21 + ppERK)) \$ \\ 
  \$d/dt[ pMpR]\$ & \$ 1*(  -(k5r*pMpR) + k5f*pMST2*pRaf1) \$ \\ 
  \$d/dt[ MpR]\$ & \$ 1*(  -(k6r*MpR) + k6f*MST2*pRaf1) \$ \\ 
  \$d/dt[ RaMk]\$ & \$ 1*(  k16af*MEK*Raf1a - k16ar*RaMk) -1*( k16b*RaMk) \$ \\ 
  \$d/dt[ RapMk]\$ & \$ 1*(  k18af*pMEK*Raf1a - k18ar*RapMk) -1*( k18b*RapMk) \$ \\ 
   \hline
\end{tabular}
\end{table}

