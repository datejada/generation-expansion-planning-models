*-------------------------------------------------------------------------
*        Universidad Pontificia Comillas de Madrid
*        Optimization Techniques
*        Diego Alejandro Tejada Arango
*-------------------------------------------------------------------------

$TITLE Two-Stage Stochastic Generation Expansion Planning - Using Benders

* ========================================================================
* SETS DEFINITION
* ========================================================================
SETS
   p       "time periods (e.g., hours)      " /h01 *h24 /
   sc      "uncertainty scenarios           " /sc01,sc02,sc03/
   g       "generation          technologies" / wind, solar, ccgt, ocgt /
   r(g)    "subset of renewable technologies" / wind, solar/
   l       "iterations                      " / l01 * l50 /
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
   pBdTol                     "relative Benders' tolerance" / 1e-6 /                     
   pZ_Lower                   "lower bound                " / -INF /                              
   pZ_Upper                   "upper bound                " /  INF /
   pZ_Bounds_L     (l,*     ) "bounds as each iteration   "
   pInstalUnits_L  (l,  g   ) "first stage variables values               in iteration l"
   pDualRenProd_L  (l,p,g,sc) "dual variables of second stage constraints in iteration l"
   pDualMaxProd_L  (l,p,g,sc) "dual variables of second stage constraints in iteration l"
   pDelta          (l       ) "cut type (feasibility 0 optimality 1)      in iteration l"
   pZ2_L           (l       ) "subproblem objective function value        in iteration l"
;
* ========================================================================
* VARIABLES
* ========================================================================
INTEGER VARIABLE
   vInstalUnits(g)    "number of installed generation units [N]"
;
POSITIVE VARIABLE
   vProduct(p,g,sc) "generation production per scenario [MW]"
   vENS    (p,  sc) "energy not supplied   per scenario [MW]"   
;
FREE VARIABLES
* complete problem variables
   vTotalCost    "Total Cost = Investment + Operation [kEUR]"
   vInvesCost    "Total investment Cost               [kEUR]"   
   vOperaCost    "Total operating  Cost               [kEUR]"
* variables for Benders decomposition
   vZ1           "first  stage objective function           "
   vZ2           "second stage objective function           "
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
* equations for Benders decomposition
   eZ1                "first  stage     objective function       "
   eZ2                "second stage     objective function       "
   eBdCuts   (l     ) "Benders cuts at iteration l               "
;
eZ1..        vZ1        =E= vInvesCost + vTheta
;                       
eZ2..        vZ2        =E=            + vOperaCost
;
eTotalCost.. vTotalCost =E= vInvesCost + vOperaCost
;
eInvesCost.. vInvesCost =E= SUM[g, pGenInfo(g,'InvCost')*pGenInfo(g,'UnitCap')*vInstalUnits(g)]
;
eOperaCost.. vOperaCost =E= pWeight * SUM[(p,g,sca),
                                            pScProb(sca)*[
                                             + pGenInfo(g,'VarCost')*vProduct(p,g,sca)
                                             + pENSCost             *vENS    (p,  sca)]]
;
eBalance(p,sca(sc))..
    SUM[g,vProduct(p,g,sc)] + vENS(p,sc) =E= pDemand(p) 
;
eRenProd(p,r,sca(sc))..
    vProduct(p,r,sc) =L= pRenProf(p,r,sc) * pGenInfo(r,'UnitCap')*vInstalUnits(r)
;
eMaxProd(p,g,sca(sc))$[NOT r(g)]..
    vProduct(p,g,sc) =L=                    pGenInfo(g,'UnitCap')*vInstalUnits(g)
;

eBdCuts(ll)..
    pDelta(ll)* vTheta
        =G= + pZ2_L(ll) 
            - SUM[(p,r,sc)           , pDualRenProd_L(ll,p,r,sc) * pRenProf(p,r,sc) * pGenInfo(r,'UnitCap') * [pInstalUnits_L(ll,r) - vInstalUnits(r)]]
            - SUM[(p,g,sc)$[NOT r(g)], pDualMaxProd_L(ll,p,g,sc)                    * pGenInfo(g,'UnitCap') * [pInstalUnits_L(ll,g) - vInstalUnits(g)]]
;

MODEL Master_Bd     / eZ1       , eInvesCost, eBdCuts                                 / ;
MODEL Subproblem_Bd / eZ2                   , eOperaCost, eBalance, eRenProd, eMaxProd/ ;
MODEL Complete      / eTotalCost, eInvesCost, eOperaCost, eBalance, eRenProd, eMaxProd/ ;

* ========================================================================
* OPTIONS AND INITIAL VALUES
* ========================================================================

* to allow CPLEX correctly detect rays in an infeasible problem
* only simplex method can be used and no preprocessing neither scaling options
* optimality and feasibility tolerances are very small to avoid primal degeneration

FILE COPT / cplex.opt /
;
PUT  COPT PUTCLOSE 'ScaInd -1' / 'LPMethod 1' / 'PreInd 0' / 'EpOpt 1e-9' / 'EpRHS 1e-9' / ;
;
Subproblem_Bd.OptFile = 1 ;
;
* parameters initialization
vTheta.fx                  = 0 ;
LL              (l       ) = NO; 
pInstalUnits_L  (l,  g   ) = 0 ;
pDualRenProd_L  (l,p,g,sc) = 0 ;
pDualMaxProd_L  (l,p,g,sc) = 0 ;
pDelta          (l       ) = 0 ;
pZ2_L           (l       ) = 0 ;

* option to find the solution to optimality
OPTION optcr=0;

* active  uncertainty scenarios with probability
sca(sc) $[pScProb(sc)] = YES ;

* ========================================================================
* BENDERS DECOMPOSITION
* ========================================================================

LOOP(l $[ABS(1-pZ_Lower/pZ_Upper) > pBdTol],

*  solving master problem
   solve Master_Bd using MIP minimizing vZ1 ;

*  storing the master solution
   pInstalUnits_L(l,g) = vInstalUnits.L(g) ;

*  fixing first-stage variables and solving subproblem
   vInstalUnits.FX (g) = vInstalUnits.L(g) ;

*  solving subproblem
   solve Subproblem_Bd using RMIP minimizing vZ2 ;

*  storing parameters to build a new Benders' cut
   if (Subproblem_Bd.ModelStat = 4,
      pDelta(l) = 0 ;
      pZ2_L (l) = Subproblem_Bd.SumInfes ;
   else
*     updating lower and upper bound
      pZ_Lower =              vZ1.L ;
      pZ_Upper = MIN(pZ_Upper,vZ1.L - vTheta.L + vZ2.L) ;

      pDelta(l) =     1 ;
      pZ2_L (l) = vZ2.L ;
   ) ;

   vTheta.lo = -INF ;
   vTheta.up =  INF ;
   pDualRenProd_L (l,p,g,sc) = eRenProd.M(p,g,sc);
   pDualMaxProd_L (l,p,g,sc) = eMaxProd.M(p,g,sc);

*  it is important to put bounds on the master's variables
*  to avoid an unbounded master problem. If the variables
*  do not have initial bounds, one can impose maximum bounds
*  that makes sense to the problem. Unfixing the variables
*  is also necessary to let GAMS optimize againg. 
   vInstalUnits.LO(g)                             = 0 ;
*  Naive approach (all to a big-number 9999)
*   vInstalUnits.UP(g)                             = 9999 ; 
*  Using the problems inputs
   vInstalUnits.UP(g) $[NOT pGenInfo(g,'UnitCap')]= 0 ;
   vInstalUnits.UP(g) $[    pGenInfo(g,'UnitCap') AND NOT r(g)]= CEIL[SMAX[ p    ,pDemand(p)/ pGenInfo(g,'UnitCap')                    ]];
   vInstalUnits.UP(r) $[    pGenInfo(r,'UnitCap')             ]= CEIL[SMAX[(p,sc),pDemand(p)/[pGenInfo(r,'UnitCap')*pRenProf(p,r,sc)+1]]]; 

*  increase the set of Benders' cuts
   LL(l) = YES ;
  
*  save the bounds of the next iteration (to report results)
   pZ_Bounds_L(l+1,'lower') = pZ_Lower ;
   pZ_Bounds_L(l+1,'upper') = pZ_Upper ;
) ;

* result parameters
PARAMETERS
pInstalCap(g   ) "installed capacity       [MW]     "
pScPrices (p,sc) "scenario  prices         [EUR/MWh]"
pEVPrices (p   ) "expected value of prices [EUR/MWh]"
;
pInstalCap(g) = pGenInfo(g,'UnitCap')*vInstalUnits.L(g)
;
pScPrices (p,sca(sc)) =                       eBalance.M(p,sc)  *1e3 / [pWeight * pScProb(sc)];
pEVPrices (p        ) = SUM[sc, pScProb(sc) * pScPrices (p,sc)]                               ;

* gdx with all results
execute_unload 'TwoStageStochGEP-Benders.gdx'

*$stop
* ========================================================================
* COMPLETE MODEL SOLUTION FOR VALITATION
* ========================================================================

SOLVE Complete USING MIP MINIMIZING vTotalCost
;

* result parameters
pInstalCap(g) = pGenInfo(g,'UnitCap')*vInstalUnits.L(g)
;
pScPrices (p,sca(sc)) =                       eBalance.M(p,sc)  *1e3 / [pWeight * pScProb(sc)];
pEVPrices (p        ) = SUM[sc, pScProb(sc) * pScPrices (p,sc)]                               ;

* gdx with all results
execute_unload 'TwoStageStochGEP-Complete.gdx'
