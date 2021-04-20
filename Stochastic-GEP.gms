*-------------------------------------------------------------------------
*        Universidad Pontificia Comillas de Madrid
*        Optimization Techniques
*        Diego Alejandro Tejada Arango
*-------------------------------------------------------------------------

$TITLE Two-Stage Stochastic Generation Expansion Planning

* ========================================================================
* SETS DEFINITION
* ========================================================================
SETS
   p       "time periods (e.g., hours)      " /h01 *h24 /
   sc      "            scenarios           " /sc01,sc02,sc03,sc_exp/
   scu(sc) "uncertainty scenarios           " /sc01,sc02,sc03/
   sce(sc) "expected value of scenarios     " /sc_exp/
   g       "generation          technologies" / wind, solar, ccgt, ocgt /
   r(g)    "subset of renewable technologies" / wind, solar/
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
   vTotalCost         "Total Cost = Investment + Operation [kEUR]"
   vInvesCost         "Total investment Cost               [kEUR]"   
   vOperaCost         "Total operating  Cost               [kEUR]"   
;
* ========================================================================
* EQUATIONS AND MODEL DEFINITION
* ========================================================================
EQUATIONS
   eTotalCost         "Total Cost = Investment + Operation [kEUR]"
   eInvesCost         "Total investment Cost               [kEUR]"   
   eOperaCost         "Total operating  Cost               [kEUR]"   
   eBalance  (p,  sc) "power balance constriant            [MW]  "
   eRenProd  (p,g,sc) "renewable  production constriant    [MW]  "
   eMaxProd  (p,g,sc) "generation production constraint    [MW]  "
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
   eMaxProd(p,g,sca(sc))$[not r(g)]..
        vProduct(p,g,sc) =L=                    pGenInfo(g,'UnitCap')*vInstalUnits(g)
;

MODEL TwoStageStochGEP / all /
;
* ========================================================================
* MODEL SOLUTION AND RESULTS
* ========================================================================
* option to find the solution to optimality
OPTION optcr=0;

* active only uncertainty scenarios with probability
sca(sc) $[scu(sc) AND pScProb(sc)] = YES ;

SOLVE TwoStageStochGEP USING MIP MINIMIZING vTotalCost
;

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
execute_unload 'TwoStageStochGEP.gdx'
*$stop
* ========================================================================
* STOCHASTIC MEASURES
* ========================================================================
;
PARAMETERS
pRP                 "Recourse Problem solution 'here and now'                               "
pWS                 "Wait and See solution                                                  "
pEVPI               "Expected Value of Perfect Information                                  "
pEV                 "Expect Value Problem result                                            "
pEEV                "Expected Result of Expect Value Problem                                "
pVSS                "Value of the Stochastic solution                                       "
pScProb_aux  (  sc) "auxiliary parameter of scenario probabilities                          "
pDetTotCost  (  sc) "total cost         of each deterministic scenario                      "
pDetENS      (  sc) "energy not supply  of each deterministic scenario                      "
pDetInstalCap(g,sc) "installed capacity of each deterministic scenario                      "
pEEVTotCost  (  sc) "total cost         of each deterministic scenario using the EV solution"
pEVVENS      (  sc) "energy not supply  of each deterministic scenario using the EV solution"
pEEVInstalCap(g,sc) "installed capacity of each deterministic scenario using the EV solution"
pEVInstalCap (g   ) "installed capacity of expected value problem                           "
;

* 1) save the solution from the stochastic model
pRP = vTotalCost.L ;

* 2) for the WS value we need the solution of each deterministic scenario
* 2.1) store the initial values of prob in a aux parameter and clean the initial values
pScProb_aux(sca(sc)) = pScProb(sc) ;
pScProb    (    sc ) = 0           ;
* 2.2)loop over the uncertanty scenarios solving once at the time
sca(sc) = NO ;
loop(scu,
    pScProb(scu) = 1   ;
    sca    (scu) = YES ;
    SOLVE TwoStageStochGEP USING MIP MINIMIZING vTotalCost ;
    pDetTotCost  (  scu) = vTotalCost.L ;
    pDetInstalCap(g,scu) = pGenInfo(g,'UnitCap')*vInstalUnits.L(g);
    pDetENS      (  scu) = SUM[p,vENS.L(p,scu)];
    pScProb      (  scu) = 0  ;
    sca          (  scu) = NO ;
);
* 2.3) we determine the wait and see solution
pWS = sum[scu, pScProb_aux(scu)*pDetTotCost(scu)]
;
* 3) we determine the EVPI
pEVPI = pRP - pWS 
;
* 4) determine the expected value of scenarios
pRenProf(p,r,sce) = SUM[scu,pScProb_aux(scu)*pRenProf(p,r,scu)]
;
* 5) solve the model for the expected value and determine EEV
pScProb(sce) = 1   ;
sca    (sce) = YES ;
SOLVE TwoStageStochGEP USING MIP MINIMIZING vTotalCost ;
pEV = vTotalCost.L ;
pDetInstalCap(g,sce) = pGenInfo(g,'UnitCap')*vInstalUnits.L(g)
;
* 6) save the solution of the expected value
pEVInstalCap(g) = vInstalUnits.L(g)
;
* 7) solve each deterministic scenario with this solution
sca(sc) = NO ;
loop(scu,
    pScProb(scu) = 1   ;
    sca    (scu) = YES ;
    vInstalUnits.FX(g) = pEVInstalCap(g) ;
    SOLVE TwoStageStochGEP USING MIP MINIMIZING vTotalCost ;
    pEEVTotCost  (  scu) = vTotalCost.L ;
    pEVVENS      (  scu) = SUM[p,vENS.L(p,scu)];
    pScProb      (  scu) = 0  ;
    sca          (  scu) = NO ;
);
* 8) we determine the EEV
pEEV = sum[scu, pScProb_aux(scu)*pEEVTotCost(scu)]
;
* 9) Determine VSS
pVSS = pEEV - pRP
;
* gdx with stochastic measures
execute_unload 'TwoStageStochGEP-Stoch-Measures.gdx'
pRP, pWS, pEVPI, pEV, pEEV, pVSS, pDetTotCost, pDetENS
pDetInstalCap, pEEVTotCost, pEVVENS, pEVInstalCap, pRenProf
; 
