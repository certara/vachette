# -------------------------------------------------------------------------
#  Sponsor           : Merck
#  Project Number    : MRK-FTE-VRSV-640
#  Project Root Path : Local
# -------------------------------------------------------------------------
#  Program : vachette-models-v33.R
#  Author  : Jos Lommerse - Certara
#  Date    : 17 January 2022
#  Purpose : mrgsolve model codes
# -------------------------------------------------------------------------
#  Software : R version 4.1.2 (2021-11-01)
#  Platform : ThinkPad P15s Gen 2 Model 20W6008AUS
#  Environment : Windows-10, 11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz
# -------------------------------------------------------------------------

library(mrgsolve)

# ------ model code -------------

iv <- '

[ PARAM ] @annotated
TVCL   : 1    : Clearance (L/hr)
TVV    : 5    : Volume of distribution (L)
WT     : 70   : Weight (kg)

[ CMT ] CENT

[ MAIN ]
double CL   = TVCL*pow(WT/70,0.75)*exp(ECL);
double V    = TVV*exp(EV);

[ OMEGA ] @name IIV @annotated
ECL   : 0.0 : Eta-CL
EV    : 0.0 : Eta-V

[ SIGMA ] @name SGMA @annotated
PROP   : 0.0 : Prop-err
ADD    : 0.0 : Add-err

[ ODE ]
dxdt_CENT = -(CL/V)*CENT;    // i.v.

[ TABLE ]
double IPRED  = CENT/V;
double DV     = IPRED*(1+PROP)+ADD;

while(DV < 0) {
  simeps();
  DV = IPRED*(1+PROP)+ADD;
}

[ CAPTURE ] WT CP=IPRED CP_RUV=DV ECL EV PROP ADD;

'


indirect.response <- '

[ PARAM ] @annotated
TVCL   : 1    : Clearance (L/hr)
TVV    : 5    : Volume of distribution (L)
TVEC50 : 2    : EC50
TVKIN  : 2    : Synthesis rate
TVBSLN : 10   : Baseline
WT     : 70   : Weight (kg)

[ CMT ] CENT RESP

[ GLOBAL ]
#define INHIB (1-CENT/(CENT+EC50))

[ MAIN ]

double CL   = TVCL*exp(ECL);
double V    = TVV*exp(EV);
double EC50 = TVEC50*pow(WT/70,1.25)*exp(EEC50);
double KIN  = TVKIN*pow(WT/70,1.25)*exp(EKIN);
double BSLN = TVBSLN*exp(EBSLN);
double KOUT = KIN/BSLN;

RESP_0 = KIN/KOUT;

[ OMEGA ] @name IIV @annotated
ECL   : 0.0 : Eta-CL
EV    : 0.0 : Eta-V
EEC50 : 0.0 : Eta-EC50
EKIN  : 0.0 : Eta-KIN
EBSLN : 0.0 : Eta-BSLN

[ SIGMA ] @name SGMA @annotated
PROP   : 0.0 : Prop-err
ADD    : 0.0 : Add-err

[ ODE ]
dxdt_CENT = -(CL/V)*CENT;    // i.v.
dxdt_RESP = KIN*INHIB - KOUT*RESP;

[ TABLE ]
double IPRED    = RESP;
double RESP_RUV = IPRED*(1+PROP)+ADD;

while(RESP_RUV < 0) {
simeps();
RESP_RUV = IPRED*(1+PROP)+ADD;
}

[ CAPTURE ] CP = CENT/V WT RESP_RUV TVCL ECL;

'

oral1cmt <- '

[ PARAM ] @annotated
TVCL : 1    : Clearance (L/hr)
TVV  : 35   : Volume of distribution (L)
TVKA : 0.5  : Absorption rate constant (1/hr)
TVF  : 1    : Bioavailability
WT   : 70   : Weight (kg)

[ CMT ] GUT CENT

[ SET ] delta=1

[ MAIN ]

double CL = TVCL*pow(WT/70,0.75)*exp(ECL);
double V  = TVV*exp(EV);
double KA = TVKA*pow(WT/70,0.75)*exp(EKA);

F_GUT  = TVF*pow(WT/70,0.9);

[ OMEGA ] @name IIV @annotated
ECL   : 0.0 : Eta-CL
EV    : 0.0 : Eta-V
EKA   : 0.0 : Eta-KA

[ SIGMA ] @name SGMA @annotated
PROP   : 0.0 : Prop-err
ADD    : 0.0 : Add-err

[ ODE ]
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (CL/V)*CENT;

[ TABLE ]
double IPRED  = CENT/V;
double DV     = IPRED*(1+PROP)+ADD;

while(DV < 0) {
simeps();
DV = IPRED*(1+PROP)+ADD;
}

[ CAPTURE ] WT CP=IPRED CP_RUV=DV;

'

oral1cmt.wt.age <- '

[ PARAM ] @annotated
TVCL : 1    : Clearance (L/hr)
TVV  : 35   : Volume of distribution (L)
TVKA : 0.5  : Absorption rate constant (1/hr)
TVF  : 1    : Bioavailability
WT   : 70   : Weight (kg)
AGE  : 30   : Age (yrs)

[ CMT ] GUT CENT

[ SET ] delta=1

[ MAIN ]

double CL = TVCL*pow(AGE/30,0.75)*exp(ECL);
double V  = TVV;
double KA = TVKA*pow(WT/70,0.75)*exp(EKA);

F_GUT  = TVF*pow(WT/70,0.9);

[ OMEGA ] @name IIV @annotated
ECL   : 0.0 : Eta-CL
EV    : 0.0 : Eta-V
EKA   : 0.0 : Eta-KA

[ SIGMA ] @name SGMA @annotated
PROP   : 0.0 : Prop-err
ADD    : 0.0 : Add-err

[ ODE ]
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (CL/V)*CENT;

[ TABLE ]
double IPRED  = CENT/V;
double DV     = IPRED*(1+PROP)+ADD;

while(DV < 0) {
simeps();
DV = IPRED*(1+PROP)+ADD;
}

[ CAPTURE ] WT AGE CP=IPRED CP_RUV=DV;

'

pembro <- '
// -------------------------------------------------------------------------
//  Sponsor           : Pembrolizumab model for VACHETTE
//  Project Number    : MRK-FTE-VRSV-640
//  Project Root Path : Local
// -------------------------------------------------------------------------
//  Program : Pembrolizumab
//  Author  : Nele Mueller-Plock - Certara
//  Date    : 08 Mar 2022
//  Purpose : 1. script contains model code for published Pembrolizumab PK model (Ahmadi, CPT Pharmacometrics Syst. Pharmacol. (2017) 6, 49xxxxxxx57), to be sourced
// -------------------------------------------------------------------------
//  Software : Version 2021.09.1 - © 2009-2021 RStudio, PBC
//             R version 4.1.2 (2021-11-01)
//  Environment : Windows NT 10.0; Win64; x64
//

$CMT CENT PERI

$PARAM

TVCL = 0.219687,
TVV1 = 3.47684,
TVQ  = 0.794715,
TVV2 = 4.05996,
TVEXCLQ = 0.594596,
TVEXV1V2 = 0.489134,
RESERR = -0.272147

CLALB = -0.906561,
CLBSLD = 0.0871905,
CLEGFR = 0.134997,
CLSEX = -0.151971,
CLADIAGN = 0.144816,
CLBECOGN = -0.0738877,
CLIPIP = 0.139662,

V1ALB = -0.207954,
V1SEX = -0.134075,
V1IPIP = 0.0736625,

SEX = 1,
IPIP = 0,
ALB = 39.6,
EGFR = 88.47,
BSLD = 89.60,
BECOG = 1,
ADIAGN = 1,
WT = 76.8,             //mean weight as in footnote of Table 3 in publication

$GLOBAL

double MWT = 76.8;     //mean weight as in footnote of Table 3 in publication

$MAIN

F_CENT = 1;

//covariates
double CV1SEX = (SEX == 2 ? (1 + V1SEX) : 1);
double CV1IPIP = (IPIP == 1 ? (1 + V1IPIP) : 1);
double CV1ALB = pow((ALB/39.60),V1ALB);

double V1COV=CV1ALB*CV1IPIP*CV1SEX;

double CCLSEX = (SEX == 2 ? (1 + CLSEX) : 1);
double CCLIPIP = (IPIP == 1 ? (1 + CLIPIP) : 1);
double CCLEGFR = pow((EGFR/88.47),CLEGFR);
double CCLBSLD = pow((BSLD/89.60),CLBSLD);
double CCLBECOGN = (BECOG == 0 ? (1 + CLBECOGN) : 1);
double CCLALB = pow((ALB/39.60),CLALB);
double CCLADIAGN = (ADIAGN == 2 ? (1 + CLADIAGN) : 1);

double CLCOV=CCLADIAGN*CCLALB*CCLBECOGN*CCLBSLD*CCLEGFR*CCLIPIP*CCLSEX;

//typical parameters
double TCL = TVCL * pow((WT/MWT),TVEXCLQ) * CLCOV;
double TV1 = TVV1 * pow((WT/MWT),TVEXV1V2) * V1COV;
double TQ = TVQ * pow((WT/MWT),TVEXCLQ);
double TV2 = TVV2 * pow((WT/MWT),TVEXV1V2);

//individual parameters
double CLi = TCL * exp(etaCLQ);
double V1i = TV1 * exp(etaV);
double V2i = TV2 * exp(etaV);
double Qi = TQ * exp(etaCLQ);


double K10 = CLi/V1i;
double K12 = Qi/V1i;
double K21 = Qi/V2i;

$ODE

double CP = CENT/V1i;

dxdt_CENT  = -K12*CENT + K21*PERI - K10*CENT;
dxdt_PERI = K12*CENT - K21*PERI;

$OMEGA @annotated @block
etaCLQ: 0.1341 : ETA on clearance and Q
etaV: 0.0159127  0.0417008 : ETA on volumes

$SIGMA 1

$TABLE
double CONC   = CENT/V1i;
double CONC_RUV = exp(log(CENT/V1i) + EPS(1)*RESERR);


$CAPTURE CONC CONC_RUV SEX IPIP ALB EGFR BSLD BECOG ADIAGN
'



pembro2 <- '
// -------------------------------------------------------------------------
//  Sponsor           : Pembrolizumab model for VACHETTE
//  Project Number    : MRK-FTE-VRSV-640
//  Project Root Path : Local
// -------------------------------------------------------------------------
//  Program : Pembrolizumab
//  Author  : Nele Mueller-Plock - Certara
//  Date    : 08 Mar 2022
//  Purpose : 1. script contains model code for published Pembrolizumab PK model (Ahmadi, CPT Pharmacometrics Syst. Pharmacol. (2017) 6, 49????"57), to be sourced
// -------------------------------------------------------------------------
//  Software : Version 2021.09.1 - © 2009-2021 RStudio, PBC
//             R version 4.1.2 (2021-11-01)
//  Environment : Windows NT 10.0; Win64; x64
//

$CMT CENT PERI

$PARAM

TVCL = 0.219687,
TVV1 = 3.47684,
TVQ  = 0.794715,
TVV2 = 4.05996,
TVEXCLQ = 0.594596,
TVEXV1V2 = 0.489134,
//original: RESERR = -0.272147
RESERR = -0.272147/5

CLALB = -0.906561,
CLBSLD = 0.0871905,
CLEGFR = 0.134997,
CLSEX = -0.151971,
CLADIAGN = 0.144816,
CLBECOGN = -0.0738877,
CLIPIP = 0.139662,

V1ALB = -0.207954,
V1SEX = -0.134075,
V1IPIP = 0.0736625,

SEX = 1,
IPIP = 0,
ALB = 39.6,
EGFR = 88.47,
BSLD = 89.60,
BECOG = 1,
ADIAGN = 1,
WT = 76.8,             //mean weight as in footnote of Table 3 in publication

$GLOBAL

double MWT = 76.8;     //mean weight as in footnote of Table 3 in publication

$MAIN

F_CENT = 1;

//covariates
double CV1SEX = (SEX == 2 ? (1 + V1SEX) : 1);
double CV1IPIP = (IPIP == 1 ? (1 + V1IPIP) : 1);
double CV1ALB = pow((ALB/39.60),V1ALB);

double V1COV=CV1ALB*CV1IPIP*CV1SEX;

double CCLSEX = (SEX == 2 ? (1 + CLSEX) : 1);
double CCLIPIP = (IPIP == 1 ? (1 + CLIPIP) : 1);
double CCLEGFR = pow((EGFR/88.47),CLEGFR);
double CCLBSLD = pow((BSLD/89.60),CLBSLD);
double CCLBECOGN = (BECOG == 0 ? (1 + CLBECOGN) : 1);
double CCLALB = pow((ALB/39.60),CLALB);
double CCLADIAGN = (ADIAGN == 2 ? (1 + CLADIAGN) : 1);

double CLCOV=CCLADIAGN*CCLALB*CCLBECOGN*CCLBSLD*CCLEGFR*CCLIPIP*CCLSEX;

//typical parameters
double TCL = TVCL * pow((WT/MWT),TVEXCLQ) * CLCOV;
double TV1 = TVV1 * pow((WT/MWT),TVEXV1V2) * V1COV;
double TQ = TVQ * pow((WT/MWT),TVEXCLQ);
double TV2 = TVV2 * pow((WT/MWT),TVEXV1V2);

//individual parameters
double CLi = TCL * exp(etaCLQ);
double V1i = TV1 * exp(etaV);
double V2i = TV2 * exp(etaV);
double Qi = TQ * exp(etaCLQ);


double K10 = CLi/V1i;
double K12 = Qi/V1i;
double K21 = Qi/V2i;

$ODE

double CP = CENT/V1i;

dxdt_CENT  = -K12*CENT + K21*PERI - K10*CENT;
dxdt_PERI = K12*CENT - K21*PERI;

$OMEGA @annotated @block
//original: etaCLQ: 0.1341               : ETA on clearance and Q
//original: etaV:   0.0159127  0.0417008 : ETA on volumes
etaCLQ: 0.067         : ETA on clearance and Q
etaV:   0.008  0.0208 : ETA on volumes

$SIGMA 1

$TABLE
double CONC   = CENT/V1i;
double CONC_RUV = exp(log(CENT/V1i) + EPS(1)*RESERR);


$CAPTURE CONC CONC_RUV SEX IPIP ALB EGFR BSLD BECOG ADIAGN
'



pembro3 <- '
// -------------------------------------------------------------------------
//  Sponsor           : Pembrolizumab model for VACHETTE
//  Project Number    : MRK-FTE-VRSV-640
//  Project Root Path : Local
// -------------------------------------------------------------------------
//  Program : Pembrolizumab
//  Author  : Nele Mueller-Plock / Jos Lommerse  - Certara
//  Date    : 14 FEB 2023
//  Purpose : 1. script contains model code for published Pembrolizumab PK model (Ahmadi, CPT Pharmacometrics Syst. Pharmacol. (2017) 6, 49xxxxxxx57), to be sourced
// -------------------------------------------------------------------------
//  Software : Version 2021.09.1 - © 2009-2021 RStudio, PBC
//             R version 4.1.2 (2021-11-01)
//  Environment : Windows NT 10.0; Win64; x64
//

$CMT CENT PERI

$PARAM

TVCL = 0.219687,
TVV1 = 3.47684,
TVQ  = 0.794715,
TVV2 = 4.05996,
TVEXCLQ = 0.594596,
TVEXV1V2 = 0.489134,
RESERR = -0.272147

CLALB = -0.906561,
CLBSLD = 0.0871905,
CLEGFR = 0.134997,
CLSEX = -0.151971,
CLADIAGN = 0.144816,
CLBECOGN = -0.0738877,
CLIPIP = 0.139662,

V1ALB = -0.207954,
V1SEX = -0.134075,
V1IPIP = 0.0736625,

SEX = 1,
IPIP = 0,
ALB = 39.6,
EGFR = 88.47,
BSLD = 89.60,
BECOG = 1,
ADIAGN = 1,
WT = 76.8,             //mean weight as in footnote of Table 3 in publication

$GLOBAL

double MWT = 76.8;     //mean weight as in footnote of Table 3 in publication

$MAIN

F_CENT = 1;

//covariates
double CV1SEX = (SEX == 2 ? (1 + V1SEX) : 1);
double CV1IPIP = (IPIP == 1 ? (1 + V1IPIP) : 1);
double CV1ALB = pow((ALB/39.60),V1ALB);

double V1COV=CV1ALB*CV1IPIP*CV1SEX;

double CCLSEX = (SEX == 2 ? (1 + CLSEX) : 1);
double CCLIPIP = (IPIP == 1 ? (1 + CLIPIP) : 1);
double CCLEGFR = pow((EGFR/88.47),CLEGFR);
double CCLBSLD = pow((BSLD/89.60),CLBSLD);
double CCLBECOGN = (BECOG == 0 ? (1 + CLBECOGN) : 1);
double CCLALB = pow((ALB/39.60),CLALB);
double CCLADIAGN = (ADIAGN == 2 ? (1 + CLADIAGN) : 1);

double CLCOV=CCLADIAGN*CCLALB*CCLBECOGN*CCLBSLD*CCLEGFR*CCLIPIP*CCLSEX;

//typical parameters
double TCL = TVCL * pow((WT/MWT),TVEXCLQ) * CLCOV;
double TV1 = TVV1 * pow((WT/MWT),TVEXV1V2) * V1COV;
double TQ = TVQ * pow((WT/MWT),TVEXCLQ);
double TV2 = TVV2 * pow((WT/MWT),TVEXV1V2);

//individual parameters
double CLi = TCL * exp(etaCLQ);
double V1i = TV1 * exp(etaV);
double V2i = TV2 * exp(etaV);
double Qi = TQ * exp(etaCLQ);


double K10 = CLi/V1i;
double K12 = Qi/V1i;
double K21 = Qi/V2i;

$ODE

double CP = CENT/V1i;

dxdt_CENT  = -K12*CENT + K21*PERI - K10*CENT;
dxdt_PERI = K12*CENT - K21*PERI;

[ OMEGA ] @name IIV @annotated
etaCLQ   : 0.0 : etaCLQ
etaV     : 0.0 : etaV

// Prop error only, via RESERR
//[ SIGMA ] @name SGMA @annotated
//PROP   : 0.0 : Prop-err
//ADD    : 0.0 : Add-err

//$OMEGA @annotated @block
//etaCLQ: 0.1341 : ETA on clearance and Q
//etaV: 0.0159127  0.0417008 : ETA on volumes

$SIGMA 1

$TABLE
double CONC     = CENT/V1i;
double CONC_RUV = exp(log(CENT/V1i) + EPS(1)*RESERR);


$CAPTURE CONC CONC_RUV SEX IPIP ALB EGFR BSLD BECOG ADIAGN
'

oral2cmt <- '

[ PARAM ] @annotated
TVCL : 1    : Clearance (L/hr)
TVV2  : 20   : Central volume of distribution (L)
TVQ  : 20   : Distribution Clearance (L/h)
TVV3 : 20   : Peripheral volume (L)
TVKA : 4   : Absorption rate constant (1/hr)
TVF  : 1    : Bioavailability
WT   : 70   : Weight (kg)

[ CMT ] GUT CENT PERIPH

[ SET ] delta=1

[ MAIN ]

double CL = TVCL*pow(WT/70,0.75)*exp(ECL);
double V2  = TVV2*exp(EV);
double Q = TVQ*pow(WT/70,0.75)*exp(ECL);
double V3  = TVV3*exp(EV);
double KA = TVKA*pow(WT/70,0.75)*exp(EKA);

F_GUT  = TVF*pow(WT/70,0.9);

[ OMEGA ] @name IIV @annotated
ECL   : 0.0 : Eta-CL
EV    : 0.0 : Eta-V
EKA   : 0.0 : Eta-KA

[ SIGMA ] @name SGMA @annotated
PROP   : 0.0 : Prop-err
ADD    : 0.0 : Add-err

[ ODE ]
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (Q/V2)*CENT + (Q/V3)*PERIPH - (CL/V2)*CENT;
dxdt_PERIPH = (Q/V2)*CENT - (Q/V3)*PERIPH;

[ TABLE ]
double IPRED  = CENT/V2;
double DV     = IPRED*(1+PROP)+ADD;

while(DV < 0) {
simeps();
DV = IPRED*(1+PROP)+ADD;
}

[ CAPTURE ] WT CP=IPRED CP_RUV=DV;

'

oral1cmt.misfit <- '

[ PARAM ] @annotated
TVCL : 1.17    : Clearance (L/hr)
TVV  : 37.2   : Volume of distribution (L)
TVKA : 21.5  : Absorption rate constant (1/hr)
TVF  : 1    : Bioavailability
WT   : 70   : Weight (kg)

[ CMT ] GUT CENT

[ SET ] delta=1

[ MAIN ]

double CL = TVCL*pow(WT/70,0.75)*exp(ECL);
double V  = TVV*exp(EV);
double KA = TVKA;

F_GUT  = TVF*pow(WT/70,0.9);

[ OMEGA ] @name IIV @annotated
ECL   : 0.0824 : Eta-CL
EV    : 0.0973 : Eta-V
EKA   : 0.00001 : Eta-KA

[ SIGMA ] @name SGMA @annotated
PROP   : 0.0575 : Prop-err
ADD    : 0.0 : Add-err

[ ODE ]
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (CL/V)*CENT;

[ TABLE ]
double IPRED  = CENT/V;
double DV     = IPRED*(1+PROP)+ADD;

while(DV < 0) {
simeps();
DV = IPRED*(1+PROP)+ADD;
}

[ CAPTURE ] WT CP=IPRED CP_RUV=DV;

'
