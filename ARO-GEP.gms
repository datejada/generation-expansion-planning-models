*-------------------------------------------------------------------------
*        Universidad Pontificia Comillas de Madrid
*        Optimization Techniques
*        Diego Alejandro Tejada Arango
*-------------------------------------------------------------------------

$TITLE Two-Stage Adaptive Robust Optimization (ARO)-Generation Expansion Planning (GEP) 

* ========================================================================
* SETS DEFINITION
* ========================================================================
SETS
   p       "time periods (e.g., hours)      " / h01 *h24 /
   g       "generation          technologies" / wind, solar, ccgt, ocgt /
   r(g)    "subset of renewable technologies" / wind, solar/
   l       "iterations                      " / l01 * l10 /
   ll(l)   "iterations subset               "
;
* ========================================================================
* PARAMETERS AND SCALARS
* ========================================================================
SCALARS
pWeight        "weight of representative period [days]" /365/
pENSCost       "energy not supplied cost    [kEUR/MWh]" /0.180/
pEXCCost       "excess of energy    cost    [kEUR/MWh]" /0.180/
pSP_UncertBudg "subproblem: Uncertainty budget        " /0.5  /
pBdTol         "tolerance   for Bender's decomposition" / 1e-6/
pLB            "lower bound for Bender's decomposition" 
pUB            "upper bound for Bender's decomposition"
;
PARAMETERS
pDemand(p)  "demand per time period [MW]"
/
h01	950
h02	870
h03	814
h04	779
h05	758
h06	751
h07	779
h08	834
h09	902
h10	956
h11	1010
h12	1023
h13	1018
h14	1010
h15	980
h16	965
h17	963
h18	997
h19	1093
h20	1114
h21	1115
h22	1107
h23	1053
h24	1035
/
;
TABLE pGenInfo(g,*) "generation information"
*       kEUR/MWh   kEUR/MW/year   MW
         VarCost    InvCost     UnitCap
ocgt      0.070         25         100
ccgt      0.050         40         400
wind      0.001         70          50 
solar     0.000         50          10   
;
TABLE pMinProf(p,g) "mininum generation profile [p.u.]"
* low wind, low solar
        wind   solar
h01     0.11   0.00
h02     0.11   0.00
h03     0.11   0.00
h04     0.11   0.00
h05     0.10   0.00
h06     0.10   0.00
h07     0.10   0.00
h08     0.09   0.01
h09     0.09   0.12
h10     0.09   0.28
h11     0.09   0.42
h12     0.09   0.51
h13     0.10   0.53
h14     0.12   0.50
h15     0.14   0.40
h16     0.15   0.23
h17     0.16   0.05
h18     0.16   0.00
h19     0.16   0.00
h20     0.15   0.00
h21     0.14   0.00
h22     0.13   0.00
h23     0.12   0.00
h24     0.12   0.00
;
TABLE pMaxProf(p,g) "maximum generation profile [p.u.]"
* high wind, high solar
        wind   solar
h01     0.68   0.00
h02     0.69   0.00
h03     0.70   0.00
h04     0.71   0.00
h05     0.73   0.00
h06     0.74   0.02
h07     0.75   0.12
h08     0.76   0.30
h09     0.77   0.50
h10     0.78   0.66
h11     0.79   0.78
h12     0.80   0.83
h13     0.81   0.83
h14     0.81   0.78
h15     0.80   0.68
h16     0.79   0.53
h17     0.78   0.35
h18     0.77   0.17
h19     0.76   0.04
h20     0.75   0.00
h21     0.74   0.00
h22     0.74   0.00
h23     0.74   0.00
h24     0.74   0.00
;
* profiles for thermal units
pMinProf(p,g)$[NOT r(g)] = 0 ;
pMaxProf(p,g)$[NOT r(g)] = 1 ;

PARAMETERS
* input parameters
pVarCost(g) "variable production cost      [kEUR/MWh]    "
pInvCost(g) "annual investment cost per MW [kEUR/MW/year]"
pUnitCap(g) "capacity per unit             [MW]          "
pGenAvai(g) "generation availability       [0-1]"

* parameters for Bender's Decomposition
pSP_InstalUnits  (  g  ) "subproblem: number of installed generation units  [N]"
pSP_InstalUnits_L(  g,l) "subproblem: number of installed generation units at iteration l [N]"
pSP_BigM_MaxProd (p,g  ) "subproblem: aux Big-M for Complementary Slackness Condition of max production"
pSP_BigM_MinProd (p,g  ) "subproblem: aux Big-M for Complementary Slackness Condition of min production"
pSP_BigM_MaxENS  (p    ) "subproblem: aux Big-M for Complementary Slackness Condition of max ENS       "
pSP_BigM_MinENS  (p    ) "subproblem: aux Big-M for Complementary Slackness Condition of min ENS       "
pSP_BigM_MaxEXC  (p    ) "subproblem: aux Big-M for Complementary Slackness Condition of max EXC       "
pSP_BigM_MinEXC  (p    ) "subproblem: aux Big-M for Complementary Slackness Condition of min EXC       "
pMP_GenAvai      (  g,l) "Master Problem: subproblem optimal value at iteration l for uncertain parameter"
pBdBounds        (l,*  ) "lower and upper bound at each iteration"
;
* ========================================================================
* VARIABLES
* ========================================================================
INTEGER VARIABLES
   vMP_InstalUnits(g)   "number of installed generation units [N]"
;
BINARY VARIABLES
   vSP_CSC_MaxProd(p,g) "aux binary for Complementary Slackness Condition of max production"
   vSP_CSC_MinProd(p,g) "aux binary for Complementary Slackness Condition of min production"
   vSP_CSC_MaxENS (p  ) "aux binary for Complementary Slackness Condition of max ENS       "
   vSP_CSC_MinENS (p  ) "aux binary for Complementary Slackness Condition of min ENS       "
   vSP_CSC_MaxEXC (p  ) "aux binary for Complementary Slackness Condition of max EXC       "
   vSP_CSC_MinEXC (p  ) "aux binary for Complementary Slackness Condition of min EXC       "
;
POSITIVE VARIABLE
   vMP_Product    (p,g,l) "Master Problem: generation production        [MW]  "
   vMP_ENS        (p  ,l) "Master Problem: energy not supplied          [MW]  "
   vMP_EXC        (p  ,l) "Master Problem: excess of energy             [MW]  "
   vSP_Product    (p,g  ) "subproblem:     generation production        [MW]  "
   vSP_ENS        (p    ) "subproblem:     energy not supplied          [MW]  "
   vSP_EXC        (p    ) "subproblem:     excess of energy             [MW]  "
   vSP_GenAvai    (g    ) "subproblem: uncertain generation availability[p.u.]"
   vSP_DualMaxProd(p,g  ) "dual variable of max production constraint   [kEUR/MWh]"
   vSP_DualMinProd(p,g  ) "dual variable of min production constraint   [kEUR/MWh]"  
   vSP_DualMaxENS (p    ) "dual variable of max ENS        constraint   [kEUR/MWh]"
   vSP_DualMinENS (p    ) "dual variable of min ENS        constraint   [kEUR/MWh]"
   vSP_DualMaxEXC (p    ) "dual variable of max EXC        constraint   [kEUR/MWh]"
   vSP_DualMinEXC (p    ) "dual variable of min EXC        constraint   [kEUR/MWh]"
;
FREE VARIABLES
   vMP_ObjFun      "Master Problem: Investment + Recourse function [kEUR]"
   vMP_Theta       "recourse function                              [kEUR]"
   vSP_ObjFun      "Subproblem: operating cost                     [kEUR]"
   vSP_DualBal(p)  "dual variable of balance equation              [kEUR/MWh]"
;
* ========================================================================
* EQUATIONS AND MODEL DEFINITION
* ========================================================================
EQUATIONS
   eMP_ObjFun               "Master Problem: Investment + Recourse function [kEUR]"
   eMP_BendersCuts  (    l) "Master Problem: Bender's Cuts                  [kEUR]"
   eMP_Balance      (p  ,l) "Master Problem: power balance constriant         [MW]"
   eMP_MaxProd      (p,g,l) "Master Problem: generation production constraint [MW]"
   eSP_ObjFun               "Subproblem: operating cost                     [kEUR]"
   eSP_UncertSet            "Subproblem: uncertainty set constraint           [MW]"
   eSP_Balance      (p    ) "Subproblem: power balance constriant             [MW]"
   eSP_MaxProd      (p,g  ) "Subproblem: generation production constraint     [MW]"
   eSP_MinProd      (p,g  ) "Subproblem: generation production constraint     [MW]"
   eSP_MaxENS       (p    ) "Subproblem: maximum ENS constraint               [MW]"
   eSP_MinENS       (p    ) "Subproblem: minimum ENS constraint               [MW]"
   eSP_MaxEXC       (p    ) "Subproblem: maximum EXC constraint               [MW]"
   eSP_MinEXC       (p    ) "Subproblem: minimum EXC constraint               [MW]"
   eSP_DL_DProd     (p,g  ) "Subproblem: lagrangian derivative respect to generation production [kEUR/MWh]"
   eSP_DL_DENS      (p    ) "Subproblem: lagrangian derivative respect to ENS                   [kEUR/MWh]"
   eSP_DL_DEXC      (p    ) "Subproblem: lagrangian derivative respect to ENS                   [kEUR/MWh]"
   eSP_CSC_MaxProd_a(p,g  ) "Subproblem: Complementary Slackness Condition max production part a[kEUR/MWh]"
   eSP_CSC_MaxProd_b(p,g  ) "Subproblem: Complementary Slackness Condition max production part b[kEUR/MWh]"
   eSP_CSC_MinProd_a(p,g  ) "Subproblem: Complementary Slackness Condition min production part a[kEUR/MWh]"
   eSP_CSC_MinProd_b(p,g  ) "Subproblem: Complementary Slackness Condition min production part b[kEUR/MWh]"
   eSP_CSC_MaxENS_a (p    ) "Subproblem: Complementary Slackness Condition max ENS        part a[kEUR/MWh]"
   eSP_CSC_MaxENS_b (p    ) "Subproblem: Complementary Slackness Condition max ENS        part b[kEUR/MWh]"
   eSP_CSC_MinENS_a (p    ) "Subproblem: Complementary Slackness Condition min ENS        part a[kEUR/MWh]"
   eSP_CSC_MinENS_b (p    ) "Subproblem: Complementary Slackness Condition min ENS        part b[kEUR/MWh]"
   eSP_CSC_MaxEXC_a (p    ) "Subproblem: Complementary Slackness Condition max EXC        part a[kEUR/MWh]"
   eSP_CSC_MaxEXC_b (p    ) "Subproblem: Complementary Slackness Condition max EXC        part b[kEUR/MWh]"
   eSP_CSC_MinEXC_a (p    ) "Subproblem: Complementary Slackness Condition min EXC        part a[kEUR/MWh]"
   eSP_CSC_MinEXC_b (p    ) "Subproblem: Complementary Slackness Condition min EXC        part b[kEUR/MWh]"
;
   eMP_ObjFun..
    vMP_ObjFun =E=
        + SUM[g, pInvCost(g)*pUnitCap(g)*vMP_InstalUnits(g)]
        + vMP_Theta
;
   eMP_BendersCuts(ll)..
    vMP_Theta =G=
        pWeight*[
                 + SUM[(p,g),pVarCost(g)*vMP_Product(p,g,ll)]
                 + SUM[(p  ),pENSCost   *vMP_ENS    (p,  ll)]
                 + SUM[(p  ),pEXCCost   *vMP_EXC    (p,  ll)]
                 ]
;
   eMP_Balance(p,ll)..
    SUM[g,vMP_Product(p,g,ll)] + vMP_ENS(p,ll) =E= pDemand(p) + vMP_EXC(p,ll)
;
   eMP_MaxProd (p,g,ll)..
    vMP_Product(p,g,ll) =L= [pMinProf(p,g)+[pMaxProf(p,g)-pMinProf(p,g)]*pMP_GenAvai(g,ll)] *pUnitCap(g)*vMP_InstalUnits(g)
;
   eSP_ObjFun..
    vSP_ObjFun =E=
        pWeight*[
                 + SUM[(p,g),pVarCost(g)*vSP_Product(p,g)]
                 + SUM[(p  ),pENSCost   *vSP_ENS    (p  )]
                 + SUM[(p  ),pEXCCost   *vSP_EXC    (p  )]
                 ]
;
   eSP_UncertSet..
    +                  SUM[g,pGenAvai(g)-vSP_GenAvai(g)] =L=
    + pSP_UncertBudg * SUM[g,pGenAvai(g)               ]
;
   eSP_Balance(p)..
    SUM[g,vSP_Product(p,g)] + vSP_ENS(p) =E= pDemand(p) + vSP_EXC(p)
;
   eSP_MaxProd (p,g)..
    vSP_Product(p,g) =L= [pMinProf(p,g)+[pMaxProf(p,g)-pMinProf(p,g)]*vSP_GenAvai(g)] * pUnitCap(g) * pSP_InstalUnits(g)
;
   eSP_MinProd (p,g)..
   -vSP_Product(p,g) =L= 0
;
   eSP_MaxENS(p)..
    vSP_ENS  (p) =L= pDemand(p)
;
   eSP_MinENS(p)..
   -vSP_ENS  (p) =L= 0
;
   eSP_MaxEXC(p)..
    vSP_EXC  (p) =L= pDemand(p)
;
   eSP_MinEXC(p)..
   -vSP_EXC  (p) =L= 0
;
   eSP_DL_DProd(p,g)..
    pWeight*pVarCost(g) - vSP_DualBal(p) + vSP_DualMaxProd(p,g) - vSP_DualMinProd(p,g) =E= 0
;
   eSP_DL_DENS (p  )..
    pWeight*pENSCost    - vSP_DualBal(p) + vSP_DualMaxENS (p  ) - vSP_DualMinENS (p  ) =E= 0
;
   eSP_DL_DEXC (p  )..
    pWeight*pEXCCost    + vSP_DualBal(p) + vSP_DualMaxEXC (p  ) - vSP_DualMinEXC (p  ) =E= 0
;
   eSP_CSC_MaxProd_a(p,g)..
    vSP_DualMaxProd(p,g) =L= pSP_BigM_MaxProd(p,g) * vSP_CSC_MaxProd(p,g)
;
   eSP_CSC_MaxProd_b(p,g)..
    [pMinProf(p,g)+[pMaxProf(p,g)-pMinProf(p,g)]*vSP_GenAvai(g)] * pUnitCap(g) * pSP_InstalUnits(g) - vSP_Product(p,g) =L= pSP_BigM_MaxProd(p,g) *[1 - vSP_CSC_MaxProd(p,g)]
;
   eSP_CSC_MinProd_a(p,g)..
    vSP_DualMinProd (p,g) =L= pSP_BigM_MinProd(p,g) *     vSP_CSC_MinProd(p,g)
;
   eSP_CSC_MinProd_b(p,g)..
    vSP_Product     (p,g) =L= pSP_BigM_MinProd(p,g) *[1 - vSP_CSC_MinProd(p,g)]
;
   eSP_CSC_MaxENS_a (p  )..
    vSP_DualMaxENS  (p  ) =L= pSP_BigM_MaxENS (p  ) *     vSP_CSC_MaxENS (p  )
;
   eSP_CSC_MaxENS_b (p  )..
    pDemand(p)-vSP_ENS(p) =L= pSP_BigM_MaxENS (p  ) *[1 - vSP_CSC_MaxENS (p  )]
;
   eSP_CSC_MinENS_a (p  )..
    vSP_DualMinENS  (p  ) =L= pSP_BigM_MinENS (p  ) *     vSP_CSC_MinENS (p  )
;
   eSP_CSC_MinENS_b (p  )..
    vSP_ENS         (p  ) =L= pSP_BigM_MinENS (p  ) *[1 - vSP_CSC_MinENS (p  )]
;
   eSP_CSC_MaxEXC_a (p  )..
    vSP_DualMaxEXC  (p  ) =L= pSP_BigM_MaxEXC (p  ) *     vSP_CSC_MaxEXC (p  )
;
   eSP_CSC_MaxEXC_b (p  )..
   pDemand(p)-vSP_EXC(p) =L= pSP_BigM_MaxEXC (p  ) *[1 - vSP_CSC_MaxEXC (p  )]
;
   eSP_CSC_MinEXC_a (p  )..
    vSP_DualMinEXC  (p  ) =L= pSP_BigM_MinEXC (p  ) *     vSP_CSC_MinEXC (p  )
;
   eSP_CSC_MinEXC_b (p  )..
    vSP_EXC         (p  ) =L= pSP_BigM_MinEXC (p  ) *[1 - vSP_CSC_MinEXC (p  )]
;

MODEL Master     / eMP_ObjFun, eMP_BendersCuts, eMP_Balance, eMP_MaxProd / ;
MODEL Subproblem / all -  Master/ ;

* ========================================================================
* PARAMETER DEFINITION
* ========================================================================

pVarCost(g) = pGenInfo(g,'VarCost') ;
pInvCost(g) = pGenInfo(g,'InvCost') ;
pUnitCap(g) = pGenInfo(g,'UnitCap') ;
pGenAvai(g) = 1                     ;

pSP_BigM_MaxProd(p,g) = 5e3 ;
pSP_BigM_MinProd(p,g) = 5e3 ;
pSP_BigM_MaxEXC (p  ) = 5e3 ;
pSP_BigM_MinEXC (p  ) = 5e3 ;
pSP_BigM_MaxENS (p  ) = 5e3 ;
pSP_BigM_MinENS (p  ) = pDemand(p) ;

* ========================================================================
* CONSTRAINTS AS BOUNDS OF VARIABLES
* ========================================================================

*  Using the problems inputs
vMP_InstalUnits.UP(g)$[NOT r(g)] = CEIL[SMAX[p,pDemand(p)/[pUnitCap(g)              +1]]]; 
vMP_InstalUnits.UP(g)$[    r(g)] = CEIL[SMAX[p,pDemand(p)/[pUnitCap(g)*pMinProf(p,g)+1]]]; 

vMP_Product.UP(p,g,l) = pMaxProf(p,g) * pUnitCap(g) * vMP_InstalUnits.UP(g) ;
vMP_ENS.UP    (p,  l) = pDemand (p  ) ;

vSP_GenAvai.UP(g) = pGenAvai(g) ;
vSP_GenAvai.LO(g) = 0           ;

* ========================================================================
* OPTIONS AND INITIAL VALUES
* ========================================================================

* to allow CPLEX correctly detect rays in an infeasible problem
* only simplex method can be used and no preprocessing neither scaling options
* optimality and feasibility tolerances are very small to avoid primal degeneration

FILE COPT / cplex.opt /
;
PUT  COPT PUTCLOSE 'ScaInd -1' / 'LPMethod 1' / 'PreInd 0' / 'EpOpt 1e-9' / 'EpRHS 1e-9' / 'IIS 1' /;
;
Subproblem.OptFile = 1 ;
;
* parameters initialization
vMP_Theta.lo     = -9999;
pLB              = -INF ;
pUB              =  INF ;
pMP_GenAvai(g,l) =  0   ;

LL(l) = NO ;

* option to find the solution to optimality
OPTION optcr=0;
*OPTION limrow=100,limcol=100;
* ========================================================================
* BENDERS DECOMPOSITION
* ========================================================================

LOOP(l $[ABS(1-pLB/pUB) > pBdTol],

*  Solve master problem
   Solve Master using MIP minimizing vMP_ObjFun ;

*  Update the lower bound
   pLB = vMP_ObjFun.L ;
   
*  Fix optimal values of first-stage variables
   pSP_InstalUnits  (g  ) = vMP_InstalUnits.L(g) ;
   pSP_InstalUnits_L(g,l) = vMP_InstalUnits.L(g) ;

*  Solve master problem
   Solve Subproblem using MIP maximizing vSP_ObjFun ;

*  Update the upper bound
   pUB = MIN[pUB, vMP_ObjFun.L + vSP_ObjFun.L - vMP_Theta.L] ;

*  Update the iteration counter
   ll(l) = YES ;

*  Set optimal values from the subproblem's solution
   pMP_GenAvai(g,l) = vSP_GenAvai.L(g) ; 

*  reset bounds for the recourse function
   vMP_Theta.lo = -INF ;
   vMP_Theta.up =  INF ;

*  Store bounds for reporting
   pBdBounds(l,'lower') = pLB ;
   pBdBounds(l,'upper') = pUB ;
   
* for debugging
   if (Subproblem.ModelStat = 10, pLB = pUB)

   )
;

* gdx with all results
execute_unload 'TwoStage-ARO-GEP.gdx'
