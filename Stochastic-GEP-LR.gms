*-------------------------------------------------------------------------
*        Universidad Pontificia Comillas de Madrid
*        Optimization Techniques
*        Diego Alejandro Tejada Arango
*-------------------------------------------------------------------------

$TITLE Two-Stage Stochastic Generation Expansion Planning - Using Lagrangian Relaxation (LR)

* ========================================================================
* SETS DEFINITION
* ========================================================================
SETS
   p       "time periods (e.g., hours)      " /h01 *h24 /
   sc      "uncertainty scenarios           " /sc01,sc02,sc03/
   g       "generation          technologies" / wind, solar, ccgt, ocgt /
   r(g)    "subset of renewable technologies" / wind, solar/
   l       "iterations                      " / l001 * l100 /
   ll(l)   "iterations subset               "
* dinamic set (to be defined depending on the input data)
   sca(sc) "active scenarios"
;
* ========================================================================
* PARAMETERS AND SCALARS
* ========================================================================
SCALARS
pWeight  "weight of representative period [days]" /365/
pENSCost "energy not supplied cost    [kEUR/MWh]" /0.180/
;
PARAMETER
pScProb(sc) "scenario probability [p.u.]"
/sc01 0.2, sc02 0.5, sc03 0.3/

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
TABLE pRenProf(p,r,sc) "renewable profile [p.u.]"
* sc01 -> low wind, high solar; sc02 -> avg wind and solar; sc03 -> high wind and low solar
	wind.sc01	wind.sc02	wind.sc03	solar.sc01	solar.sc02	solar.sc03
h01	  0.11	      0.54	       0.68	       0.00	        0.00	        0.00
h02	  0.11	      0.54	       0.69	       0.00	        0.00	        0.00
h03	  0.11	      0.53	       0.70	       0.00	        0.00	        0.00
h04	  0.11	      0.52	       0.71	       0.00	        0.00	        0.00
h05	  0.10	      0.51	       0.73	       0.00	        0.00	        0.00
h06	  0.10	      0.50	       0.74	       0.02	        0.00	        0.00
h07	  0.10	      0.48	       0.75	       0.12	        0.01	        0.00
h08	  0.09	      0.47	       0.76	       0.30	        0.07	        0.01
h09	  0.09	      0.46	       0.77	       0.50	        0.20	        0.12
h10	  0.09	      0.45	       0.78	       0.66	        0.36	        0.28
h11	  0.09	      0.45	       0.79	       0.78	        0.50	        0.42
h12	  0.09	      0.45	       0.80	       0.83	        0.57	        0.51
h13	  0.10	      0.43	       0.81	       0.83	        0.59	        0.53
h14	  0.12	      0.41	       0.81	       0.78	        0.54	        0.50
h15	  0.14	      0.38	       0.80	       0.68	        0.44	        0.40
h16	  0.15	      0.35	       0.79	       0.53	        0.29	        0.23
h17	  0.16	      0.34	       0.78	       0.35	        0.13	        0.05
h18	  0.16	      0.35	       0.77	       0.17	        0.03	        0.00
h19	  0.16	      0.36	       0.76	       0.04	        0.00	        0.00
h20	  0.15	      0.38	       0.75	       0.00	        0.00	        0.00
h21	  0.14	      0.41	       0.74	       0.00	        0.00	        0.00
h22	  0.13	      0.43	       0.74	       0.00	        0.00	        0.00
h23	  0.12	      0.46	       0.74	       0.00	        0.00	        0.00
h24	  0.12	      0.48	       0.74	       0.00	        0.00	        0.00
;
* parameters for Benders decomposition
PARAMETERS
   pDifference                "difference in dual variables among iterations"
   pLRTol                     "LR relative tolerance      " / 1e-6 /                     
   pZ_Lower                   "lower bound                " / -INF /                              
   pZ_Upper                   "upper bound                " /  INF /
   pZ_Bounds_L     (l,*     ) "bounds as each iteration   "
   pInstalUnits_L  (l,  g,sc) "installed units       in iteration l"
   pProduct_L      (l,p,g,sc) "production            in iteration l"
   pENS_L          (l,p,  sc) "energy not supplied   in iteration l"
   pLambda_L       (l,  g,sc) "lagrangian multiplier in iteration l"
;
* ========================================================================
* VARIABLES
* ========================================================================
INTEGER VARIABLE
   vInstalUnits(g,sc) "number of installed generation units [N]"
;
POSITIVE VARIABLE
   vProduct(p,g,sc)   "generation production per scenario [MW]"
   vENS    (p,  sc)   "energy not supplied   per scenario [MW]"   
;
FREE VARIABLES
* complete problem variables
   vTotalCost    "Total Cost = Investment + Operation [kEUR]"
   vInvesCost    "Total investment Cost               [kEUR]"   
   vOperaCost    "Total operating  Cost               [kEUR]"
* variables for lagrangian relaxation
   vLambda(g,sc) "lagrangian multiplier                     "
   vLRObjFun     "lagrangian term in the objective function "
   vTheta        "recourse function                         "
;
* ========================================================================
* EQUATIONS AND MODEL DEFINITION
* ========================================================================
EQUATIONS
* Complete problem equations 
   eTotalCost         "Total Cost = Investment + Operation [kEUR]"
   eInvesCost         "Total investment Cost               [kEUR]"   
   eOperaCost         "Total operating  Cost               [kEUR]"   
   eBalance  (p,  sc) "power balance constriant            [MW]  "
   eRenProd  (p,g,sc) "renewable  production constriant    [MW]  "
   eMaxProd  (p,g,sc) "generation production constraint    [MW]  "
   eInvestNAC(  g,sc) "Investment Nonanticipativity constraints  "
* equations for lagrangian relaxation
   eLRObjFun          "lagrangian term in the objective function    "
   eLRCuts   (l     ) "lagrangian cuts at iteration l               "
;
eTotalCost.. vTotalCost =E= vInvesCost + vOperaCost + vLRObjFun
;
eInvesCost.. vInvesCost =E= SUM[(g,sca(sc)), pScProb(sc)*pGenInfo(g,'InvCost')*pGenInfo(g,'UnitCap')*vInstalUnits(g,sc)]
;
eOperaCost.. vOperaCost =E= pWeight * SUM[(p,g,sca(sc)),
                                            pScProb(sc)*[
                                             + pGenInfo(g,'VarCost')*vProduct(p,g,sc)
                                             + pENSCost             *vENS    (p,  sc)]]
;
eLRObjFun.. vLRObjFun =E= SUM[(g,sca(sc)),
                               vLambda.l(g,sc)*[vInstalUnits(g,sc)-vInstalUnits(g,sc++1)]]
;
eBalance(p,sca(sc))..
    SUM[g,vProduct(p,g,sc)] + vENS(p,sc) =E= pDemand(p) 
;
eRenProd(p,r,sca(sc))..
    vProduct(p,r,sc) =L= pRenProf(p,r,sc) * pGenInfo(r,'UnitCap')*vInstalUnits(r,sc)
;
eMaxProd(p,g,sca(sc))$[NOT r(g)]..
    vProduct(p,g,sc) =L=                    pGenInfo(g,'UnitCap')*vInstalUnits(g,sc)
;
* it is possible to reduce one equation using: $[ORD(sc)<CARD(sc)]
* however, it is easier to separate the problem for the LR using the cycling operator ++
eInvestNAC(g,sca(sc))..
    vInstalUnits(g,sc) =E= vInstalUnits(g,sc++1)
;
eLRCuts(ll)..
    vTheta
        =l= + SUM[(g,sca(sc)), pScProb(sc)*pGenInfo(g,'InvCost')*pGenInfo(g,'UnitCap')*pInstalUnits_L(ll,g,sc)] 
            + pWeight * SUM[(p,g,sca(sc)),
                                            pScProb(sc)*[
                                             + pGenInfo(g,'VarCost')*pProduct_L(ll,p,g,sc)
                                             + pENSCost             *pENS_L    (ll,p,  sc)]]

            + SUM[(g,sca(sc)),
                               vLambda(g,sc)*[pInstalUnits_L(ll,g,sc)-pInstalUnits_L(ll,g,sc++1)]]
;
MODEL Master_LR     / eLRCuts / ;
MODEL Subproblem_LR / eTotalCost, eInvesCost, eOperaCost, eLRObjFun, eBalance, eRenProd, eMaxProd            / ;
MODEL Complete      / eTotalCost, eInvesCost, eOperaCost,            eBalance, eRenProd, eMaxProd, eInvestNAC/ ;

* ========================================================================
* OPTIONS AND INITIAL VALUES
* ========================================================================

* to allow CPLEX correctly detect rays in an infeasible problem
* only simplex method can be used and no preprocessing neither scaling options
* optimality and feasibility tolerances are very small to avoid primal degeneration

FILE COPT / cplex.opt / ;
PUT  COPT putclose 'ScaInd -1' / 'LPMethod 1' / 'PreInd 0' / 'EpOpt 1e-9' / 'EpRHS 1e-9' / ;

Subproblem_LR.OptFile = 1 ; Complete.OptFile = 1 ;

* active  uncertainty scenarios with probability
sca(sc) $[pScProb(sc)] = YES ;

* parameters initialization
pDifference                 = INF;
LL              (l        ) = NO ;
pInstalUnits_L  (l,  g,sca) = 0  ;
pProduct_L      (l,p,g,sca) = 0  ;
pENS_L          (l,p,  sca) = 0  ;

* lambda is free variable because it comes from equility constraint
* Therefore, we need to put both upper and lower bounds
vLambda.up      (    g,sca) = +1e9;
vLambda.lo      (    g,sca) = 0;

*  it is important to put bounds on the subproblem's variables
*  to avoid an unbounded subproblem. If the variables
*  do not have initial bounds, one can impose maximum bounds
*  that makes sense to the problem. Unfixing the variables
*  is also necessary to let GAMS optimize againg. 
vInstalUnits.LO(g,sca) = 0 ;
*  Naive approach (all to a big-number 9999)
*vInstalUnits.UP(g,sca) = 9999 ; 
*  Using the problems inputs
   vInstalUnits.UP(g,sca) $[NOT pGenInfo(g,'UnitCap')]= 0 ;
   vInstalUnits.UP(g,sca) $[    pGenInfo(g,'UnitCap') AND NOT r(g)]= CEIL[SMAX[ p    ,pDemand(p)/ pGenInfo(g,'UnitCap')                    ]];
   vInstalUnits.UP(r,sca) $[    pGenInfo(r,'UnitCap')             ]= CEIL[SMAX[(p,sc),pDemand(p)/[pGenInfo(r,'UnitCap')*pRenProf(p,r,sc)+1]]]; 

* option to find the solution to optimality
OPTION optcr=0;

* ========================================================================
* Lagrangian Relaxation
* ========================================================================

LOOP (l $(pDifference > pLRTol),

*  solving master problem
   IF (ORD(l) = 1,
      vLambda.l(g,sca) = 0
   ELSE
      SOLVE Master_LR using LP maximizing vTheta
   ) ;
   pLambda_L(l,g,sca) = vLambda.l(g,sca) ;

*  solving subproblem
   solve Subproblem_LR using MIP minimizing vTotalCost ;

   pInstalUnits_L  (l,  g,sca) = vInstalUnits.l  (  g,sca) ;   
   pProduct_L      (l,p,g,sca) = vProduct.l      (p,g,sca) ;
   pENS_L          (l,p,  sca) = vENS.l          (p,  sca) ;

   pZ_Lower = vTotalCost.l ;
   pZ_Upper = vTheta.l     ;

*  save the bounds of the next iteration (to report results)
   pZ_Bounds_L(l+1,'lower') = pZ_Lower ;
   pZ_Bounds_L(l+1,'upper') = pZ_Upper ;

*  increasing the set of LR cuts
   ll(l) = YES ;

*  updating the multiplier difference
   pDifference $[ORD(l) > 1] = SUM[(g,sca), abs(pLambda_L(l,g,sca)-pLambda_L(l-1,g,sca))]
) ;

* optimal solution of the LR (linear combination of all iterations)
* this solution is infeasible in the complete problem
vInstalUnits.l(  g,sca) = sum[ll, eLRCuts.m(ll)*pInstalUnits_L(ll,  g,sca)];
vProduct.l    (p,g,sca) = sum[ll, eLRCuts.m(ll)*pProduct_L    (ll,p,g,sca)];
vENS.l        (p,  sca) = sum[ll, eLRCuts.m(ll)*pENS_L        (ll,p,  sca)];
vTotalCost.l            = sum[ll, eLRCuts.m(ll)*[
                                + SUM[(g,sca(sc)), pScProb(sc)*pGenInfo(g,'InvCost')*pGenInfo(g,'UnitCap')*pInstalUnits_L(ll,g,sc)] 
                                + pWeight * SUM[(p,g,sca(sc)),
                                            pScProb(sc)*[
                                             + pGenInfo(g,'VarCost')*pProduct_L(ll,p,g,sc)
                                             + pENSCost             *pENS_L    (ll,p,  sc)]]
                                                ]
                             ]
;

* result parameters
PARAMETERS
pInstalCap(g,sc) "installed capacity       [MW]     "
pScPrices (p,sc) "scenario  prices         [EUR/MWh]"
pEVPrices (p   ) "expected value of prices [EUR/MWh]"
;
pInstalCap(g,sc) = pGenInfo(g,'UnitCap')*vInstalUnits.L(g,sc)
;
pScPrices (p,sca(sc)) =                       eBalance.M(p,sc)  *1e3 / [pWeight * pScProb(sc)];
pEVPrices (p        ) = SUM[sc, pScProb(sc) * pScPrices (p,sc)]                               ;

* gdx with all results
execute_unload 'TwoStageStochGEP-LR.gdx'

$ontext
NOTE: LR optimal solution is not necessarily feasible in the complicating constraints.
 To solve this situation, it is necessary to introduce a penalty in the o.f. of the subproblem,
 a.k.a Augmented Lagrangian.
$offtext
;
*$stop
* ========================================================================
* COMPLETE MODEL SOLUTION FOR VALITATION
* ========================================================================

* for the complete problem, we don't use the LR term in the objective function
vLRObjFun.fx = 0;
* solve the problem
SOLVE Complete USING MIP MINIMIZING vTotalCost
;

* result parameters
pInstalCap(g,sc) = pGenInfo(g,'UnitCap')*vInstalUnits.L(g,sc)
;
pScPrices (p,sca(sc)) =                       eBalance.M(p,sc)  *1e3 / [pWeight * pScProb(sc)];
pEVPrices (p        ) = SUM[sc, pScProb(sc) * pScPrices (p,sc)]                               ;

* gdx with all results
execute_unload 'TwoStageStochGEP-Complete.gdx'
