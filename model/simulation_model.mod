$Global

$Prob
- PK/PD simulation script

$CMT  @annotated
ABS  : Absorption compartment
CENT : Central compartment (mg)
PERI : Peripheral (mg)
RESP : Reponse compartment

$PARAM @annotated
TVKA     : 0.459        : Absorption rate
vc       : 67.9         : Individual VC
k20      : 0.178        : Elimination
k23      : 0.368        : Cmt transferring 2 -> 3 
k32      : 0.0484       : Cmt transferring 3 -> 2
MICBL    : 1.00         : MIC_C1D4 baseline
EMAX     : 10           : TVEMAX
EC50     : 200          : TVEC50
KOUT     : 0.0314       : TVKOUT
HILL     : 2            : TVHILL

$MAIN
double ik20   = k20;
double ik23   = k23;
double ik32   = k32;
double KAVAR = TVKA;
double Emax  = EMAX;
double Ec50  = EC50;
double Kout  = KOUT;
double Hill  = HILL;

double Kin   = Kout * MICBL;

RESP_0 = MICBL;

$ODE
dxdt_ABS  = -KAVAR*ABS;
dxdt_CENT = KAVAR*ABS - k20*CENT - k23*CENT + k32*PERI;
dxdt_PERI = k23*CENT - k32*PERI;

double CP = CENT/vc;

double EFFECT = 1 + (Emax * pow(CP, Hill)) / (pow(Ec50, Hill) + pow(CP, Hill)); 
dxdt_RESP = Kin*EFFECT - Kout*RESP;

$TABLE
double DV_RESP = RESP;
double DV_PK   = log(CENT/vc);

$CAPTURE
DV_RESP DV_PK
