*-------------------------------------------------------------------------
*        Universidad Pontificia Comillas de Madrid
*        Optimization Techniques
*        Diego Alejandro Tejada Arango
*-------------------------------------------------------------------------

$TITLE Worst Case Scenario -Generation Expansion Planning (GEP) 

* ========================================================================
* SETS DEFINITION
* ========================================================================
SETS
   p       "time periods (e.g., hours)      " /h01 *h24 /
   g       "generation          technologies" / wind, solar, ccgt, ocgt /
   r(g)    "subset of renewable technologies" / wind, solar/

;
* ========================================================================
* PARAMETERS AND SCALARS
* ========================================================================
SCALARS
pWeight  "weight of representative period [days]" /365/
pENSCost "energy not supplied cost    [kEUR/MWh]" /0.180/
;
PARAMETER
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
TABLE pProf(p,g) "generation profile [p.u.]"
* worst case -> low wind, low solar
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
* profile for thermal units is equal to 1
pProf(p,g)$[NOT r(g)] = 1
;
* ========================================================================
* VARIABLES
* ========================================================================
INTEGER VARIABLE
   vInstalUnits(g)    "number of installed generation units [N]"
;
POSITIVE VARIABLE
   vProduct(p,g) "generation production per scenario [MW]"
   vENS    (p  ) "energy not supplied   per scenario [MW]"   
;
FREE VARIABLES
   vTotalCost         "Total Cost = Investment + Operation [kEUR]"
   vInvesCost         "Total investment Cost               [kEUR]"   
   vOperaCost         "Total operating  Cost               [kEUR]"   
;
* ========================================================================
* EQUATIONS AND MODEL DEFINITION
* ========================================================================
EQUATIONS
   eTotalCost      "Total Cost = Investment + Operation [kEUR]"
   eInvesCost      "Total investment Cost               [kEUR]"   
   eOperaCost      "Total operating  Cost               [kEUR]"   
   eBalance  (p  ) "power balance constriant            [MW]  "
   eMaxProd  (p,g) "generation production constraint    [MW]  "
;
   eTotalCost.. vTotalCost =E= vInvesCost + vOperaCost
;
   eInvesCost.. vInvesCost =E= SUM[g, pGenInfo(g,'InvCost')*pGenInfo(g,'UnitCap')*vInstalUnits(g)]
;
   eOperaCost.. vOperaCost =E= pWeight * SUM[(p,g),[
                                                    + pGenInfo(g,'VarCost')*vProduct(p,g)
                                                    + pENSCost             *vENS    (p  )
                                                    ]
                                             ]
;
   eBalance(p).. SUM[g,vProduct(p,g)] + vENS(p) =E= pDemand(p) 
;
   eMaxProd(p,g)..
        + vProduct(p,g) =L= pProf(p,g) * pGenInfo(g,'UnitCap')*vInstalUnits(g)
;

MODEL Determ_WC_GEP / eTotalCost, eInvesCost, eOperaCost, eBalance, eMaxProd /
;

* ========================================================================
* DETERMINISTIC WORST CASE SCENARIO
* ========================================================================
* option to find the solution to optimality
OPTION optcr=0;

SOLVE Determ_WC_GEP USING MIP MINIMIZING vTotalCost
;

* result parameters
PARAMETERS
pInstalCap(g) "installed capacity       [MW]     "
pPrices   (p) "electricity  prices      [EUR/MWh]"
;

pInstalCap(g) = pGenInfo  (g,'UnitCap')*vInstalUnits.L(g) ;
pPrices   (p) = eBalance.M(p          )*1e3 / pWeight     ;

* gdx with all results
execute_unload 'WorstCaseScenario-GEP.gdx'