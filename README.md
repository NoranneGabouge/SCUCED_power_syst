# Power Systems Optimization
## Security Constrained Unit Commitment Economic Dispatch problem solved with Benders decomposition

### SC-UCED problem 

The problem addressed in this project is the Unit Commitment
and Economic Dispatch, which is fundamental in the operation of electrical grids. The
setting of our problem is the following: before the beginning of an operation period, the
power system operator receives the demand forecast in each time period and price offers
of each generating unit in the system for the entire operation horizon. We consider here
an operation period of one day and time intervals of one hour. With this information, the
operator must solve the following two problems:
- The Unit Commitment problem (UC) consists in selecting which generating units
will be turned on and which ones will be left unused, in each time period. The
operator must strive here to minimize fixed operation cost while committing enough
generators to supply all demand, satisfy operational requirements such as spinning
reserves and minimum running times and satisfy system requirements, which amounts
to guaranteeing that the Economic Dispatch problem will be feasible.
- The Economic Dispatch problem (ED) consists in determining the output of each generating unit that has been selected in the UC. The constraints that must be satisfied
here are the physical constraints of the electric system: mainly verifying that trans-
mission line capacities will not be exceeded in the base case or in any of a select set
of contingencies. In this stage, the operator aims at minimizing the variable cost of
generating electricity.

Since we are not only considering the operation of the system in its normal state but also in
the case that one of a select set of contingencies occurs, we speak of a Security-Constrained
Unit Commitment and Economic Dispatch (SC-UCED). 

### Scope
This project considers the SC-UCED problem formulated as a Mixed Integer Program (MIP) and solves it using Benders decomposition.
We are limiting our approach to solving the SC-UCED problem for one day with hourly intervals (half-hourly intervals for some computational experiments), 
applying Benders' decomposition and relying on a simplified linear formulation of power flow (so-called DC Power Flow).
Costs are modeled as a piecewise linear convex function and we assume fast decoupled power flow without losses

 The original problem formulation is decomposed into three layers, each one of
them corresponding to one specific component of the problem (UC master problem,
ED subproblems, SC sub-subproblems). 

We conduct computational experiments and
compare the results to the performance of an off-the-shelf solver.

Data available in ```data``` folder obtained from IEEE test cases with different number of buses

Completed: Spring 2020  
Co-authors: Tomas Valencia Zuluaga, Tim Schmidtlein
