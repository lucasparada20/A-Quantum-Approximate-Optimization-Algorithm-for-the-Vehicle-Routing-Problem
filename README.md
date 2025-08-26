# A-Quantum-Approximate-Optimization-Algorithm-for-the-Vehicle-Routing-Problem

## Model Formulation

## QUBO Formulation for the VRP

We model the Vehicle Routing Problem (VRP) as a Quadratic Unconstrained Binary Optimization (QUBO) as follows. 

**(1) Objective (total travel cost)**

$\min \ \sum_{(i,j)\in A} c_{ij}\,x_{ij}$

where:
- \(x_{ij}\in\{0,1\}\) is 1 if arc \((i,j)\) is used,
- \(c_{ij}\) is the travel cost on arc \((i,j)\),
- \(A\) is the set of allowed arcs.

**(2) Outgoing arc constraint**  
$$
\sum_{j \in V} x_{ij} = 1 \quad \forall i \in N
$$

**(3) Incoming arc constraint**  
$$
\sum_{j \in V} x_{ji} = 1 \quad \forall i \in N
$$

**(4) Number of vehicles leaving the depot**  
$$
\sum_{j \in A(0)} x_{0j} = k
$$

**(5) Number of vehicles returning to the depot**  
$$
\sum_{j \in A(0)} x_{j0} = k
$$

Here $$A(0)$$ is the set of arcs of the depot, $$N$$ is the set of customer nodes and $$k$$ is the number of vehicles.

