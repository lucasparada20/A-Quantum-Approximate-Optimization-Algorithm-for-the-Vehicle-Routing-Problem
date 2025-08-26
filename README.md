# A-Quantum-Approximate-Optimization-Algorithm-for-the-Vehicle-Routing-Problem

## Model Formulation

We model the Vehicle Routing Problem (VRP) as a Quadratic Unconstrained Binary Optimization (QUBO). Let graph $G=(V,A)$, with $V=N\cup\{0\}$ be the set of nodes including the depot $\{0\}$, and $A(0)$ be the set of arcs that connect the depot. We are given a fleet of $k$ vehicles. The model follows:

**(1) Objective (total travel cost)**
$\min \ \sum_{(i,j)\in A} c_{ij}x_{ij}$

Where:
- $x_{ij}\in $\{0,1\}$ is 1 if arc $\(i,j\)$ is used.
- $c_{ij}$ is the travel cost on arc $\(i,j\)$.

**(2) Outgoing arc constraint**  
$\sum_{j \in V} x_{ij} = 1 \quad \forall i \in N$

**(3) Incoming arc constraint**  
$\sum_{j \in V} x_{ji} = 1 \quad \forall i \in N$

**(4) Number of vehicles leaving the depot**  
$\sum_{j \in A(0)} x_{0j} = k$

**(5) Number of vehicles returning to the depot**  
$\sum_{j \in A(0)} x_{j0} = k$

```math
\begin{aligned}
\min_{x \in \{0,1\}} \quad 
& \underbrace{\sum_{i \in V} \sum_{\substack{j \in V \\ j \neq i}} c_{ij}\, x_{ij}}_{\text{(1) travel cost}} \\[1em]
&+ A \underbrace{\sum_{i \in N} \left( \sum_{\substack{j \in V \\ j \neq i}} x_{ij} - 1 \right)^{2}}_{\text{(2) outgoing arc constraints}} \\[1em]
&+ A \underbrace{\sum_{i \in N} \left( \sum_{\substack{j \in V \\ j \neq i}} x_{ji} - 1 \right)^{2}}_{\text{(3) incoming arc constraints}} \\[1em]
&+ A \underbrace{\left( \sum_{j \in V \setminus \{0\}} x_{0j} - k \right)^{2}}_{\text{(4) depot out-degree}} \\[1em]
&+ A \underbrace{\left( \sum_{j \in V \setminus \{0\}} x_{j0} - k \right)^{2}}_{\text{(5) depot in-degree}}
\end{aligned}
```

## QUBO Formulation for the VRP


