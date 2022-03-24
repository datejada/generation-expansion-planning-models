# Generation Expansion Planning (GEP) models considering uncertainties on renewable energy resources (RES)

The following files solve the GEP problem for three scenarios of wind and solar production using different approaches:

* **Stochastic-GEP.gms**: Two-Stage Stochastic Generation Expansion Planning
* **Stochastic-GEP-Benders.gms**: Two-Stage Stochastic Generation Expansion Planning - Using Benders
* **Stochastic-GEP-Benders-Multicut.gms**: Two-Stage Stochastic Generation Expansion Planning - Using Benders Multicut
* **Stochastic-GEP-LR.gms**: Two-Stage Stochastic Generation Expansion Planning - Using Lagrangian Relaxation (LR)
* **Stochastic-GEP-Multistage.gms**: Multi-Stage Stochastic Generation Expansion Planning
* **SRO-GEP.gms**: Single-Stage Static Robust Optimization (SRO) or Scenario-based RO -Generation Expansion Planning (GEP)
* **ARO-GEP.gms**: Two-Stage Adaptive Robust Optimization (ARO)-Generation Expansion Planning (GEP)
* **WCS-GEP.gms**: Worst Case Scenario -Generation Expansion Planning (GEP)
* **Bilevel-Centralized-GEPM.gms**: The Bilevel formulation Deterministic Single-Node Static Generation Expansion Planning Model (GEPM)

The models are developed in [GAMS](https://www.gams.com/) and solved with [CPLEX](https://www.ibm.com/analytics/cplex-optimizer), but you could use any other solver (e.g., [GUROBI](https://www.gurobi.com/), [Cbc](https://github.com/coin-or/Cbc)).

The main references to model the optimization problems are:

[1] [Optimization Techniques by Andrés Ramos Galán](https://pascua.iit.comillas.edu/aramos/OT.htm)

[2] [A. J. Conejo, L. Baringo, S. J. Kazempour and A. S. Siddiqui, Investment in Electricity Generation and Transmission, Cham, Zug, Switzerland:Springer, 2016.](https://link.springer.com/book/10.1007/978-3-319-29501-5)

[3] [Sun X.A., Conejo A.J. (2021) Static Robust Optimization. In: Robust Optimization in Electric Energy Systems. International Series in Operations Research & Management Science, vol 313. Springer, Cham.]( https://doi.org/10.1007/978-3-030-85128-6_2)
