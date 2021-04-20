*-------------------------------------------------------------------------
*        Universidad Pontificia Comillas de Madrid
*        Optimization Techniques
*        Diego Alejandro Tejada Arango
*-------------------------------------------------------------------------

$TITLE Multi-Stage Stochastic Generation Expansion Planning

* NOTE: For the sake of simplicity, we assume that all monetary values are referred to the
* same point in time, and thus, it is not necessary to multiply by discount rates.

* ========================================================================
* SETS DEFINITION
* ========================================================================
SETS
   y       "years                           " /y01 *y05 /
   p       "time periods (e.g., hours)      " /h01 *h24 /
   sc      "            scenarios           " /sc01,sc02,sc03/
   g       "generation          technologies" / wind, solar, ccgt, ocgt /
   r(g)    "subset of renewable technologies" / wind, solar/
* dinamic set (to be defined depending on the input data)
   sca(sc) "active scenarios"
   ya (y ) "active years    "
;
ALIAS (y,yy)
* ========================================================================
* PARAMETERS AND SCALARS
* ========================================================================
SCALARS
pWeight  "weight of representative period [days]" /365/
pENSCost "energy not supplied cost    [kEUR/MWh]" /0.180/
;
PARAMETER
pCumDemIncr(y ) "cumulative yearly demand increase [p.u.]"
pOrder     (y ) "ordinal of the year                     "
pDemIncr   (y ) "     yearly demand increase       [p.u.]"
/
y01 0.00
y02 0.01
y03 0.01
y04 0.01
y05 0.01
/
pScProb    (sc) "scenario probability           [p.u.]"
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

* ========================================================================
* VARIABLES
* ========================================================================
INTEGER VARIABLE
   vInstalUnits (y,g) "           number of installed generation units [N]"
   vCumInstUnits(y,g) "cumulative number of installed generation units [N]"
;
POSITIVE VARIABLE
   vProduct(y,p,g,sc) "generation production per scenario [MW]"
   vENS    (y,p,  sc) "energy not supplied   per scenario [MW]"   
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
   eTotalCost           "Total Cost = Investment + Operation [kEUR]"
   eInvesCost           "Total investment Cost               [kEUR]"   
   eOperaCost           "Total operating  Cost               [kEUR]"
   eCumInUnit(y,  g   ) "cumulative Installed units          [N ]  "
   eBalance  (y,p,  sc) "power balance constriant            [MW]  "
   eRenProd  (y,p,g,sc) "renewable  production constriant    [MW]  "
   eMaxProd  (y,p,g,sc) "generation production constraint    [MW]  "
;
    eTotalCost.. vTotalCost =E= vInvesCost + vOperaCost
;
    eInvesCost.. vInvesCost =E= SUM[(ya(y),g), [CARD(y)-ORD(y)+1]*pGenInfo(g,'InvCost')*pGenInfo(g,'UnitCap')*vInstalUnits(y,g)]
;
    eOperaCost.. vOperaCost =E= pWeight * SUM[(ya(y),p,g,sca(sc)),
                                               pScProb(sc)*[
                                                + pGenInfo(g,'VarCost')*vProduct(y,p,g,sc)
                                                + pENSCost             *vENS    (y,p,  sc)]]
;
   eCumInUnit(ya(y),g)..
        vCumInstUnits(y,g) =E= SUM[yy $[pOrder(yy) <= ORD(y)], vInstalUnits(yy,g)]
;
   eBalance(ya(y),p,sca(sc))..
        SUM[g,vProduct(y,p,g,sc)] + vENS(y,p,sc) =E= pDemand(p) * pCumDemIncr(y)
;
   eRenProd(ya(y),p,r,sca(sc))..
        vProduct(y,p,r,sc) =L= pRenProf(p,r,sc) * pGenInfo(r,'UnitCap')*vCumInstUnits(y,r)
;
   eMaxProd(ya(y),p,g,sca(sc))$[not r(g)]..
        vProduct(y,p,g,sc) =L=                    pGenInfo(g,'UnitCap')*vCumInstUnits(y,g)
;

MODEL TwoStageStochGEP / all /
;
* ========================================================================
* MODEL SOLUTION AND RESULTS
* ========================================================================
* option to find the solution to optimality
OPTION optcr=0;

* active only uncertainty scenarios with probability
sca(sc) $[pScProb(sc)] = YES ;

* compute the cumulative yearly demand growth
ya         (y) $[ORD(y)=1 OR pDemIncr(y)] = YES                   ;
pOrder     (y) = ORD(y)                                           ;
pCumDemIncr(y) = PROD[yy $[pOrder(yy) <= ORD(y)], 1+pDemIncr(yy)] ;

SOLVE TwoStageStochGEP USING MIP MINIMIZING vTotalCost
;

* result parameters
PARAMETERS
pInstalCap(y,g   ) "installed capacity       [MW]     "
pScPrices (y,p,sc) "scenario  prices         [EUR/MWh]"
pEVPrices (y,p   ) "expected value of prices [EUR/MWh]"
;
pInstalCap(y,g) = pGenInfo(g,'UnitCap')*vInstalUnits.L(y,g)
;
pScPrices (y,p,sca(sc)) =                       eBalance.M(y,p,sc)  *1e3 / [pWeight * pScProb(sc)];
pEVPrices (y,p        ) = SUM[sc, pScProb(sc) * pScPrices (y,p,sc)]                               ;

* gdx with all results
execute_unload 'MultiStageStochGEP.gdx'