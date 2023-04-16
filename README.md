# Generation Expansion Planning (GEP) models considering uncertainties on renewable energy resources (RES)

## Files description
The following files solve the GEP problem for three scenarios of wind and solar production using different approaches:

* **[Stochastic-GEP.gms](Stochastic-GEP.gms)**: Two-Stage Stochastic Generation Expansion Planning
* **[Stochastic-GEP-Benders.gms](Stochastic-GEP-Benders.gms)**: Two-Stage Stochastic Generation Expansion Planning - Using Benders
* **[Stochastic-GEP-Benders-Multicut.gms](Stochastic-GEP-Benders-Multicut.gms)**: Two-Stage Stochastic Generation Expansion Planning - Using Benders Multicut
* **[Stochastic-GEP-LR.gms](Stochastic-GEP-LR.gms)**: Two-Stage Stochastic Generation Expansion Planning - Using Lagrangian Relaxation (LR)
* **[Stochastic-GEP-Multistage.gms](Stochastic-GEP-Multistage.gms)**: Multi-Stage Stochastic Generation Expansion Planning
* **[SRO-GEP.gms](SRO-GEP.gms)**: Single-Stage Static Robust Optimization (SRO) or Scenario-based RO -Generation Expansion Planning (GEP)
* **[ARO-GEP.gms](ARO-GEP.gms)**: Two-Stage Adaptive Robust Optimization (ARO)-Generation Expansion Planning (GEP)
* **[WCS-GEP.gms](WCS-GEP.gms)**: Worst Case Scenario -Generation Expansion Planning (GEP)
* **[Bilevel-Centralized-GEPM.gms](Bilevel-Centralized-GEPM.gms)**: The Bilevel formulation Deterministic Single-Node Static Generation Expansion Planning Model (GEPM)

The models are developed in [GAMS](https://www.gams.com/) and solved with [CPLEX](https://www.ibm.com/analytics/cplex-optimizer), but you could use any other solver (e.g., [GUROBI](https://www.gurobi.com/), [Cbc](https://github.com/coin-or/Cbc)).

## GEP Formulation

### Indices
| **Name** | **Description**                   |
|----------|-----------------------------------|
| $p$      | time periods                      |
| $g$      | generation technologies           |
| $r(g)$   | subset of renewable techonologies |
| $sc$     | scenarios                         |

### Parameters
| **Name**   | **Domains** | **Description**                                             |
|------------|-------------|-------------------------------------------------------------|
| $pVOLL   $ |             | Value of Lost Load [\$/MWh]                                 |
| $pWeight $ |             | Representative period weight [hours]                        |
| $pInvCost$ | $g$         | Investment cost [\$/MW]                                     |
| $pVarCost$ | $g$         | Variable production cost [\$/MWh]                           |
| $pUnitCap$ | $g$         | Capacity per each invested unit [MW/unit]                   |
| $pRenProf$ | $r,p,sc$    | Renewable profile (e.g., load factor) [p.u.]                |
| $pDemand $ | $p$         | Demand [MW]                                                 |
| $pScProb $ | $sc$        | Scenario probability [p.u.]                                 |

### Variables
| **Name**    | **Domains** | **Description**              |
|-------------|-------------|------------------------------|
| $vTotCost $ |             | Total system cost [\$]       |
| $vInvCost $ |             | Total investment cost [\$]   |
| $vOpeCost $ |             | Total operating cost [\$]    |
| $vGenInv  $ | $g$         | Generation investment [1..N] |
| $vGenProd $ | $g,p,sc$    | Generation production [MW]   |
| $vLossLoad$ | $p,sc$      | Loss of load [MW]            |

### Equations
| **Name**                                    | **Domains** | **Description**                    |
|---------------------------------------------|-------------|------------------------------------|
| [eObjFun](#eobjfun)                         |             | Total system cost      [\$]        |
| [eInvCost](#einvcost)                       |             | Total investment cost      [\$]    |
| [eOpeCost](#eopecost)                       |             | Total operating cost      [\$]     |
| [eBalance](#ebalance)                       | $p,sc$      | Power system balance   [MWh]       |
| [eMaxProd](#emaxprod)                       | $g,p,sc$    | Maximum generation production [MW] |
| [eRenProd](#erenprod)                       | $r,p,sc$    | Maximum renewable production [MW]  |

#### *eObjFun*
$$
\displaystyle{\min{vTotCost = vInvCost + vOpeCost}}
$$

#### *eInvCost*
$$
vInvCost = \displaystyle \sum_{g}(pInvCost_{g} \cdot pUnitCap_{g} \cdot vGenInv_{g})
$$

#### *eOpeCost*
$$
vOpeCost = pWeight \cdot {\left(\displaystyle \sum_{sc}pScProb_{sc}\cdot{\left(\sum_{g,p}pVarCost_{g} \cdot vGenProd_{g,p,sc} + \sum_{p,sc}pVOLL \cdot vLossLoad_{p,sc}\right)}\right)}
$$

#### *eBalance*
$$
\displaystyle \sum_{g}vGenProd_{g,p,sc} + vLossLoad_{p,sc} = pDemand_{p} \quad \forall{p,sc} 
$$

#### *eMaxProd*
$$
vGenProd_{g,p,sc} \leq pUnitCap_{g} \cdot vGenInv_{g} \quad \forall{g,p,sc} 
$$

#### *eRenProd*
$$
vGenProd_{r,p,sc} \leq pRenProf_{r,p,sc} \cdot pUnitCap_{r} \cdot vGenInv_{r} \quad \forall{r,p,sc} 
$$

#### *Bounds*
$vGenProd_{g,p,sc}\geq 0 ~ \forall g, p, sc $

$vLossLoad_{p,sc}\geq 0 ~ \forall p, sc $

$vGenInv_{g} \in \mathbb{Z}^{+} ~ \forall g $

## References
The main references to model the optimization problems are:

[1] [Optimization Techniques by Andrés Ramos Galán](https://pascua.iit.comillas.edu/aramos/OT.htm)

[2] [A. J. Conejo, L. Baringo, S. J. Kazempour and A. S. Siddiqui, Investment in Electricity Generation and Transmission, Cham, Zug, Switzerland:Springer, 2016.](https://link.springer.com/book/10.1007/978-3-319-29501-5)

[3] [Sun X.A., Conejo A.J. (2021) Static Robust Optimization. In: Robust Optimization in Electric Energy Systems. International Series in Operations Research & Management Science, vol 313. Springer, Cham.]( https://doi.org/10.1007/978-3-030-85128-6_2)
