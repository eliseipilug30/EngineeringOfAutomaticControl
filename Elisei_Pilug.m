n = 5.22;

%% IDENTIFICAREA BUCLEI DE REGLARE A DEBITULUI - Partea 1

% a) referitor la transportul melcat

%amplificatorul de putere
Kap = 20 + 0.25*n
tauap = 10^(-2)

%motorul de antrenare
K1 = 0.25 + 0.01*n
K2 = 5 + 0.1*n
Tm1 = 0.05 + 2*0.001*n
Tm2 = 0.5 + 0.01*n

%tahogeneratorul de masurare a turatiei
KTomega = 0.1
TTomega = 0.01

%transportul melcat TM
TQ1 = 5
TTM = 5
KTM = 0.12 + 0.01 * n
TB = 60

% b) referitor la transportul cu cupe TC 
K = 0.9
T = 10
omegasteaTC = 1
taum = 60/omegasteaTC

% c) doza gavimetrica cu adaptor
TG = 2
KG = 0.16

%% IDENTIFICAREA BUCLEI DE REGLARE A TEMPERATURII 


% a) venitul pneumatic UP si convertorul electropneumatic
Kpg = 0.01
KceKv= 0.025
Tv = 4

% b) cuptorul C
Kc = 200 + 2*n
KtetaT = 0.8
KtetaZ = 0.3
TtetaZ = 120
Tc = 600 + 5*n
tauc = 0.15 * Tc
TtetaT = 100 + 2*n
taut = 0.05*TtetaT

% d) traductoarele de temperatura
 Ktetam = 0.16
 Ttetam = 4

 Ktetac = 0.1
 Ttetac = 16


 %% Partea 2 - Calculul regulatoarelor prin metoda repartitiei poli-zerouri

 % a) se neglijeaza perturbatiile
 Mr = 0

 % b) se aproximeaza timpul mort sub forma
 timp_mort = 1/(1+tauap)

 % c) se transfigurează  sistemul la forma reacţiei negative unitare
 TMstelat = Tm2 + TTomega
 KMstelat = KTomega * K2
 
 H = tf(KMstelat , [TMstelat , 1]) 

 %% 1.2 Calculul regulatorului HR1(s) pt cazul sistemului echivalent de ordinul doi necorectat

 % a) se impune urmatorul set de performante

  Estpimpus = 0;
  suprareglajimpus = 0.15;
  trimpus = 1.2
  Wb <= 12


% b) Parametrii sistemului echivalent de ordinul doi
tzita = abs(log(suprareglajimpus))/sqrt((log(suprareglajimpus))^2 + pi^2)
Wn = 4/(trimpus*tzita)
 
Ho2 = tf(Wn^2 , [1 2*tzita*Wn Wn^2])

% c) Verificarile impuse in aceasta situatie sunt:

tr = 4/ (tzita*Wn)
Wb = Wn*(1-2*tzita^2+(2-4*tzita^2+4*tzita^4)^(1/2))^(1/2)
Estp = 1-evalfr(Ho2,0)
Estvinit = 2*tzita/Wn
step(Ho2); title('Raspuns la treapta sistem inchis in cazul cu HR1') % => suprareglaj de 15
t = 0:0.1:10;
figure
lsim(Ho2,t,t);  title('Raspuns la rampa sistem inchis in cazul HR1')
% d) determinarea analitica a regulatorului 

numf = [Kap*K1*KMstelat];
denf = [tauap*Tm1*TMstelat tauap*TMstelat+tauap*Tm1+Tm1*TMstelat tauap+Tm1+TMstelat 1];
Hf = tf(numf, denf)

% bucla directa este:
Hd = tf(Wn^2 , [1 2*tzita*Wn 0 ])

HR1 = Hd*(1/Hf)

%d2)
TMdublustelat = tauap + TMstelat
num1 = [Tm1*TMdublustelat*Wn/(2*tzita) (Tm1+TMdublustelat)*Wn/(2*tzita) Wn/(2*tzita) ];
den1 = [Kap*K1*KMstelat/(2*tzita*Wn) Kap*K1*KMstelat 0];
HR1modif = tf(num1, den1)
Hdmodif = HR1modif*Hf
H02modif = feedback(Hdmodif , 1)

step(H02modif,Ho2); title('Raspuns la treapta sistem inchis in cazul cu HR1secund')
t = 0:0.1:10;
figure
lsim(H02modif,t,t); title('Raspuns la rampa sistem inchis in cazul HR1secund')
hold on;
lsim(Ho2,t,t);


%d1) alegem:
tzitaaprox = 0.278;
Wnaprox = 3.2;
% se simplifica cu polul din TMstelat*s+1;
numaprox = [tauap*Tm1 tauap+Tm1 1];
denaprox = [1 0];
HR1aprox = tf(Wnaprox/(tzitaaprox*2)*numaprox , Kap*KMstelat*K1*denaprox)
Hdaprox = HR1aprox*Hf
H02aprox = feedback(Hdaprox , 1)
figure
step(H02aprox)
hold on
step(Ho2); title('Raspuns la treapta sistem inchis in cazul cu HR1prim')
t = 0:0.1:10;
figure
lsim(H02aprox,t,t); title('Raspuns la rampa sistem inchis in cazul HR1prim')
hold on;
lsim(Ho2,t,t);

% e) 

zpk(HR1modif) 
% T1 = 0.06 , T2 = 0.57 , TN = 0.14, UT = 1.009
 

%% 1.3 Calcul HR2 pentru cazul sistemului de ordinul doi corectat (Corectia cu Dipol)

% a) setul de performante:
%Estp = 0;
suprareglajstelat = 0.1
trstelat = 1;
%Wb <= 12;
 Estv = 0.05
deltasuprareglajc = 0.05

suprareglaj2 = suprareglajstelat-deltasuprareglajc
tzita2 = abs(log(suprareglaj2))/sqrt((log(suprareglaj2))^2 + pi^2)
Wn2 = 4/trstelat*tzita2



Ho3 = tf(Wn2^2 , [1 2*tzita2*Wn2 Wn2^2])
Hd3 = Ho3 / (1-Ho3)

Estp3 = 1-evalfr(Ho3,0)
step(Ho3) % => suprareglaj de 5%
tr3 = 4/ (tzita2*Wn2) % nu se indeplineste aceasta conditie
Wb3 = Wn2*(1-2*tzita2^2+(2-4*tzita2^2+4*tzita2^4)^(1/2))^(1/2)  %mai mica decat Wbimpus
Estv1 = 2*tzita2/Wn2 %nu se respecta, = 0.5

tzita3 = tzita2;
Wn3 = Wb3 / sqrt((1-2*tzita2^2)*sqrt(2-4*tzita2^2+4*tzita2^4))

%b)
pc = deltasuprareglajc / ((2*tzita3/Wn3) - Estv)
zc = pc / (1 + deltasuprareglajc)

Hcorectie = tf([pc zc*pc],[zc pc*zc])
Hoc = series(Ho3 , Hcorectie)

figure
step(Hoc) % => tr=2.37
t = 0:0.1:10;
figure
lsim(Hoc,t,t); 

%c) 
tv = 4/(tzita3*Wn3) %se respecta


%d)
a = 2*tzita3*Wn3 - pc
b = 2*tzita3*Wn3 - Wn3^2*deltasuprareglajc

num3 = [1 a b 0];
den3 = 1;
L = tf(num3,den3)

beta1 = -0.5298;
beta2 = 17.5349;


numr2 = conv([Wn3^2+deltasuprareglajc*Wn3^2 zc*Wn2^2+zc*Wn3^2*deltasuprareglajc] , [tauap*Tm1*TMstelat tauap*TMstelat+tauap*Tm1+TMstelat*Tm1 tauap+Tm1+TMstelat 1])
denr2 = conv([1 beta1 0],[1 beta2]);
HR2 = tf(numr2 , Kap*KMstelat*K1*denr2)

% d1) simplificari
num_d1 = [conv([0.01 1.008431 0.8431] , [0.56 1])];
den_d1 = [1 -0.5294 0];
HR2_d1 = tf(1680.0133*1.05*num_d1 , Kap*K1*KMstelat*den_d1)

% d2) simplificari
num_d2 = [conv([1/0.8431 1] , [0.57 1])];
den_d2 = [-1/0.5294 1];
HR2_d2 = tf(-280.94*num_d2 , Kap*K1*KMstelat*den_d2)


% d3) simplificari
num_d3 = conv([1.186 1],[-1.31 1]);
HR2_d3 = tf(-79.34*num_d3,1)


% e) identificarea parametrilor regulatorului




%% 1.4 CALCULUL REGULATORULUI HR3 - CORECTIA CU ZERO SI POL DE BALAST

% a) setul de performante:
%Estpimpus = 0;
suprareglaj_init = 0.1;
tr_init = 1;
%Wbimpus <= 12;
 Estv_stelat = 0.05;


 % b) parametrii sistemului inchis
 tzita_stelat =  abs(log(suprareglaj_init))/sqrt((log(suprareglaj_init))^2 + pi^2)
 tzita_prim = 1.1*tzita_stelat
 Wn_prim = 4/(tr_init*tzita_prim)

 cv_stelat = 1/Estv_stelat

 zc1 = Wn_prim/(2*tzita_prim)
 pc1 = 3.25*Wn_prim

 % c) verificarea performantelor impuse


Ho4 = tf(Wn_prim^2 , [1 2*tzita_prim*Wn_prim Wn_prim^2])
Hd4 = Ho4 / (1-Ho4)

 Hcorectie_1 = tf([pc1 zc1*pc1],[zc1 pc1*zc1])
Hoc_1 = series(Ho4 , Hcorectie_1)

figure
step(Hoc_1) % => settling time = 0.838 ,  overhoot = 21.5%

Estp4 = 1-evalfr(Ho4,0)
Wb4 = Wn_prim*(1-2*tzita_prim^2+(2-4*tzita_prim^2+4*tzita_prim^4)^(1/2))^(1/2)
Estv4 = 2*tzita_prim/Wn_prim


%% CALCULUL REGULATOARELOR PRIN METODE FRECVENTIALE 

Tf = tauap + Tm1
Kf = Kap * K1

Tstelat = TMstelat/KMstelat 

% 2.2 Determinarea lui VR a unui regulator P

% a) setul de performante

sigma_stea = 0.15;
tr_stea = 1;
cv_stea = 5;
Wb_stea = 10;

% b) calculul factorului de amplificare VR

%partea fixata:
den_VR = [Tf*Tstelat Tstelat 0]
Hf_VR = tf(Kf , den_VR)

Wf_det = 14.28;
tzita_det = abs(log(sigma_stea))/sqrt((log(sigma_stea))^2 + pi^2)
A = 1/(4*(tzita_det^2))
VR = 9.7724

bode(Hf_VR) 
Wt_det = 5.66 %din bode

% c) verificari ale urmatoarelor performante
Wn_det = 2*tzita_det*Wt_det
tr_det = 4/(tzita_det*Wn_det)
Wb_det = Wt_det
sigma_det = 0.15
cv_det = Wn_det/(2*tzita_det)

HR_det = VR ; 
Hd_det = HR_det*Hf_VR
figure
Ho_det = Hd_det/(1+Hd_det)
% Ho_det = tf(Wn_det^2 , [1 2*tzita_det*Wn_det Wn_det^2])
step(Ho_det); title('Regulator P')
t = 0:0.1:10;
figure
lsim(Ho_det,t,t);  title('Regulator P')




% 2.3 Determinarea parametrilor unui regulator PI

% a) set de performante
sigma_stea1 = 0.075;
tr_stea1 = 1;
cv_stea1 = 10;
Wb_stea1 = 10;

% b) calculul parametrilor VR , Tz , Tp

%partea fixata:
den_VR1 = [Tf*Tstelat Tstelat 0]
Hf_VR1 = tf(Kf , den_VR1)

Wf_det1 = 14.28;
bode(Hf_VR1)
Wt_det1 = 5.66;
tzita_det1 = abs(log(sigma_stea1))/sqrt((log(sigma_stea1))^2 + pi^2)
A1  = 1/(4*(tzita_det1^2))
VR1 = 14.7126
Wn_det1 = 2*tzita_det1*Wt_det1
cv_det1 = Wn_det1/(2*tzita_det1)


Wz = 0.1*Wt_det1
Wp = (cv_det1/cv_stea1)*Wz

Tz = 1/Wz
Tp = 1/Wp

HPI = tf(VR1*[Tz 1],[Tp 1])
HdPI = HPI*Hf_VR1
Ho_det1 = HdPI/(1+HdPI)
bode(Ho_det1) % => Wb = 138 la -3dB
step(Ho_det1); title('Regulator PI')
figure
t = 0:0.016:10;
lsim(Ho_det1,t,t); title('Regulator PI')





% 2.3 Determinarea parametrilor unui regulator PD

% a) set de performante
sigma_stea2 = 0.1;
tr_stea2 = 0.5;
cv_stea2 = 5;
Wb_stea2 = 10;

% b) calculul parametrilor VRt , taud , taun

%partea fixata:
den_VR2 = [Tf*Tstelat Tstelat 0]
Hf_VR2 = tf(Kf , den_VR2)

Wf_det2 = 14.28;
bode(Hf_VR2)
Wt1_det2 = 5.66;
tzita_det2 = abs(log(sigma_stea2))/sqrt((log(sigma_stea2))^2 + pi^2)
A2 = 1/(4*(tzita_det2^2))
VRt = 12.6911
tr_det2 = 2/(tzita_det2^2*Wt1_det2)

Wt2_det2 = Wt1_det2 * (tr_det2/tr_stea2)
taud = Tf
taun = taud * (tr_stea2/tr_det2)

HPD = tf(VRt*[taud 1],[taun 1])

HdPD = HPD*Hf_VR2
Ho_det2 = HdPD/(1+HdPD)
bode(Ho_det2) % => Wb = 69.7 la -3dB
figure
step(Ho_det2); title('Regulator PD')
figure
t = 0:0.016:10;
lsim(Ho_det2,t,t); title('Regulator PD')



% 2.5 Determinarea unui regulator PID 

% a) set de performante
sigma_stea3 = 0.075
tr_stea3 = 0.5
Wb_stea3 = 10
cv_stea3 = 12


Wt1_det3 = 5.66;
tzita_det3 = abs(log(sigma_stea3))/sqrt((log(sigma_stea3))^2 + pi^2)
A3 = 1/(4*(tzita_det3^2))
tr_det3 = 2/(tzita_det3^2*Wt1_det3)
Wt2_det3 = Wt1_det3 * (tr_det3/tr_stea3)
VR_det3 = 14.71
Wn_det3 = 2*tzita_det3*Wt1_det3
cv_det3 = Wn_det3/(2*tzita_det3)


 Wz_det3 = 0.1*Wt2_det3
 Wp_det3 = Wz_det3*cv_det3/cv_stea3

 Tz_det3 = 1/Wz_det3
 Tp_det3 = 1/Wp_det3
 taud_det3 = Tf
 taun_det3 = taud_det3*Wt2_det3/Wt1_det3


H_PID = HPD * HPI
 HdPID = H_PID*Hf_VR2
Ho_det3 = HdPID/(1+HdPID)
bode(Ho_det3) % => Wb = 213
figure
step(Ho_det3); title('Regulator PID')
figure
t = 0:0.016:10;
lsim(Ho_det3,t,t); title('Regulator PID')
%% CALCULUL REGULATOARELOR IN CAZUL REGLARII IN CASCADA


Kap = 20 + 0.25*n
tauap = 10^(-2)
K1 = 0.25 + 0.01*n
K2 = 5 + 0.1*n
Tm1 = 0.05 + 2*0.001*n
Tm2 = 0.5 + 0.01*n
KTomega = 0.1
TTomega = 0.01
TQ1 = 5
TTM = 5
KTM = 0.12 + 0.01 * n
TB = 60
K = 0.9
T = 10
omegasteaTC = 1
taum = 60/omegasteaTC
TG = 2
KG = 0.16
Kpg = 0.01
KceKv= 0.025
Tv = 4
Kc = 200 + 2*n
KtetaT = 0.8
KtetaZ = 0.3
TtetaZ = 120
Tc = 600 + 5*n
tauc = 0.15 * Tc
TtetaT = 100 + 2*n
taut = 0.05*TtetaT
 Ktetam = 0.16
 Ttetam = 4
 Ktetac = 0.1
 Ttetac = 16
 Mr = 0
 timp_mort = 1/(1+tauap)
 TMstelat = Tm2 + TTomega
 KMstelat = KTomega * K2


 KTM_stea = KTM/KTomega
 TTM_stea = TTM - TTomega

 K_stea = K*KG
 T_stea = T + timp_mort + TG

 Kf = KceKv*Kc*Ktetac
 Tf = Tv + Tc + Ttetac

 Ktetat_stea = KtetaT / Ktetac
 Ttetat_stea = TtetaT - Ttetac

 %% CALCULUL REGULATORULUI SISTEMULUI DE REGLARE A DEBITULUI (Qm)

 num = [Tm1*TMstelat Tm1+TMstelat 1];
 den = [2*tauap*Kap*K1*KMstelat 0]; 
 HR_omega = tf(num,den)
 HO_omega = tf(1,[2*(tauap^2) 2*tauap 1])
 Homegao = tf(1,[2*tauap 1])

H_Prod = tf(KTM_stea , [TTM_stea 1])

TTM_dublustea = TTM_stea + 2*tauap

num1 = [T_stea^2 T_stea];
den1 = [2*TTM_dublustea*T_stea*KTM_stea*K_stea 0];
HRQ = tf(num1, den1)


tau1 = Tm1
tau2 = TMstelat
VRomega = 1/(2*tauap*Kap*K1*KMstelat)

VRQ = T_stea/(2*TTM_dublustea*KTM_stea*K_stea)
tauiQ = T_stea

Mr = 0.1
Tsum = tauap

num2 = [KMstelat*2*Mr*T*Tsum KMstelat*2*Mr*Tsum];
den2 = [conv([TMstelat 1],[2*tauap*tauap 2*tauap 1])];
uomega = tf(num2,den2)
iq = tf([KTM_stea],[TTM_stea 1])

figure
t=0:0.1:10
step_response = step(uomega,t)
plot(t,step_response)

Hf1 = tf([Kap*K1],[tauap*Tm1 tauap+Tm1 1])
Hf2 = tf([KMstelat],[TMstelat 1])
Hf3 = tf([KTM_stea*K_stea],[TTM_stea*T_stea TTM_stea+T_stea 1])

Hf_reg = tf(12.8991,[conv([0.01 1 0],[12.994 17.98 1])])
figure
bode(Hf_reg)

W_prim = 0.0311
modul = 468.27

VR = 1/modul
Ti = 4/W_prim
 HReg = tf([VR*Ti Ti],[Ti 0])

 uomega1 = Hf2/(1+HReg*Hf1*Hf2*Hf3)

figure
step(0.1*uomega1)
G1 = feedback(uomega1, 0.1);
step(G1);

%% CALCULUL REGULATORULUI SISTEMULUI DE REGLARE A TEMPERATURII (Qm)

VRC = (0.9*tauc)/(Kf*Tf)
tauic = 3.3*tauc

num_temp = [2.3*tauc*VRC*Kf VRC*Kf]
den_temp = [3.3*(tauc^2)*Tf 3.3*tauc*(tauc+Tf) 3.3*tauc+2.3*tauc*VRC*Kf VRC*Kf]
HOC_temp = tf(num_temp , den_temp)
p1 = -0.0106
p2 = -0.0009
p3 = -0.0007

T1_prim = -0.0011
T2_prim = -0.0006
B = 1
K_prim = VRC*Kf/B
T0C = T1_prim+T2_prim-2.3*tauc

H0C_prim = tf(K_prim , [T0C 1])


[mag, w] = bode(H0C_prim)
VRMAX = max(mag)
w0 = w(find(abs(w) == min(abs(w))))

T0 = 2*pi/w0

Vrm = VRMAX*0.75
tauim = 0.6*T0
taudm = 0.1*T0

% Proportional (P) controller
P = tf([Vrm], [1]);

% Integral (I) controller
I = tf([1], [tauim, 0]);

% Derivative (D) controller
D = tf([taudm, 0], [1]);

% PID controller
PID = P + I + D

Hrm = P
Te = 0.1*(TtetaT+Ttetam)

H_num_om = Hrm*H0C_prim*tf([Ktetat_stea*Ttetam Ktetat_stea],1,'IOdelay', taut);

H1 = tf([Ttetat_stea*Ttetam Ttetat_stea+Ttetam],1)
H2 = tf([Hrm*H0C_prim*Ktetat_stea*Ktetam],1 ,'IOdelay', taut) 
H_den_om = H1+H2
Hom = H_num_om/H_den_om
Hom_discret = c2d(Hom, Te , 'zoh')

figure
step(Hom)
hold on
step(Hom_discret)
%% CALCULUL REGULATORULUI CU PREDICTIE


Kap = 20 + 0.25*n
tauap = 10^(-2)
K1 = 0.25 + 0.01*n
K2 = 5 + 0.1*n
Tm1 = 0.05 + 2*0.001*n
Tm2 = 0.5 + 0.01*n
KTomega = 0.1
TTomega = 0.01
TQ1 = 5
TTM = 5
KTM = 0.12 + 0.01 * n
TB = 60
K = 0.9
T = 10
omegasteaTC = 1
taum = 60/omegasteaTC
TG = 2
KG = 0.16
Kpg = 0.01
KceKv= 0.025
Tv = 4
Kc = 200 + 2*n
KtetaT = 0.8
KtetaZ = 0.3
TtetaZ = 120
Tc = 600 + 5*n
tauc = 0.15 * Tc
TtetaT = 100 + 2*n
taut = 0.05*TtetaT
 Ktetam = 0.16
 Ttetam = 4
 Ktetac = 0.1
 Ttetac = 16
 Mr = 0
 timp_mort = 1/(1+tauap)
 TMstelat = Tm2 + TTomega
 KMstelat = KTomega * K2


 KTM_stea = KTM/KTomega
 TTM_stea = TTM - TTomega

 K_stea = K*KG
 T_stea = T + timp_mort + TG

 Kf = KceKv*Kc*Ktetac
 Tf = Tv + Tc + Ttetac

 Ktetat_stea = KtetaT / Ktetac
 Ttetat_stea = TtetaT - Ttetac

 T1_stea = T + TG
num = [KTM_stea*K_stea]
den = [conv(conv([2*tauap 1],[TTM_stea 1]),[T1_stea 1])]
 Hf = tf(num, den, 'IOdelay', timp_mort)
Hf_prim = tf(num, den)


 Tmin = 2*tauap + TTM_stea
 Tr = 2*Tmin
 TB = 4*Tmin
 B = tf(1 , [4*Tmin*Tmin 4*Tmin 0])

 HR1 = B/Hf_prim

 VR = T1_stea/(4*KTM_stea*K_stea*Tmin)
 taui = T1_stea

 HR = tf([VR*taui VR],[taui 0])

 H0 = tf(1,[4*Tmin*Tmin 4*Tmin 1], 'IOdelay', timp_mort)

 % bode(H0)
 tr_prim = 6*Tmin
 tr = 6*Tmin + timp_mort

 % figure
 % bode(Hf_prim)

 w_secund = 0.235
 modul_Hf = 0.0222
 VR_barat = 1/modul_Hf
 taui_barat = 4/w_secund

 HR_barat = tf([VR_barat*taui_barat VR_barat],[taui_barat 0])

 figure;
 step(H0)

