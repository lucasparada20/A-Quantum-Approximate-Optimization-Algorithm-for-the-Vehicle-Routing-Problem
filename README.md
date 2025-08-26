# A-Quantum-Approximate-Optimization-Algorithm-for-the-Vehicle-Routing-Problem

## Model Formulation

We model the Vehicle Routing Problem (VRP) as a Quadratic Unconstrained Binary Optimization (QUBO) following the approach of these [authors](https://ieeexplore.ieee.org/document/9774961). Let graph $G=(V,A)$, with $V = N \cup \{0\}$ be the set of nodes including the depot $\{0\}$, and $A(0)$ be the set of arcs that connect the depot. We are given a fleet of $k$ vehicles. The model follows:

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

In practice, this is not a VRP as the subtour elimination constraints are lacking. The problem in (1)-(5) resembles an assignment problem, but it is a simple enough problem to embed in a QAOA for operations research practitioners transitioning into quantum information technologies.

## QUBO Formulation for the VRP

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

There is a C++ `main.cpp` file that will build and solve the QUBO formulation using the Cplex library. The code requires a hardcoded distance matrix and the parameter $k$. In the current implementation, there are two distance matrices hardcoded, so you can uncomment one to test it. To compile and run the code, you first need to give the path to your Cplex library in the `Makefile` here:

```bash
CPLEX_DIR = /some/path/to/Cplex
```
Next, compile, link, and run:

```bash
make
./vrp_qubo
```

## The Python Source Code

A modular implementation of QAOA is provided in the `src` directory. The modules are the following.

* `Qubo.py`: Builds the QUBO model for a given distance matrix and number of vehicles using the following methods. `build_vrp_qubo()` resembles the C++ implementation in `main.cpp`, and returns the quadratic term matrix $Q$, the linear term array $h$, and the constant term. Next, `qubo_to_hamiltonian()` builds a Pennylane Hamiltonian object from the QUBO matrices. Lastly, the function `compute_qubo_value()` is a helper method that computes the Ising cost of a given solution provided in a bitstring Numpy ndarray.
* `Qaoa.py`: Uses PennyLane's `default.qubit` simulator backend and SciPy optimizers to minimize the QUBO Hamiltonian. It also performs multiple optimization restarts to help avoid poor local minima.
  * **`qaoa_layer(gamma, beta, H)`** – Applies one QAOA layer by evolving under the problem Hamiltonian `H` for time `gamma`, followed by RX rotations with angle `2*beta` on all qubits.  
  * **`qaoa_circuit(params, H, n_wires)`** – Builds the full QAOA circuit with `p` alternating layers of cost and mixing unitaries, starting from the uniform superposition state \(|+\rangle^{\otimes n}\).  
  * **`run(qubo_lin, qubo_quad, const_term, var_map, p=1, n_shots=10000)`** – Converts a QUBO into a Hamiltonian, executes the QAOA circuit on PennyLane’s `default.qubit` simulator, and uses classical optimization with multiple restarts to minimize the expected cost, returning samples, the best energy, and the optimal parameters.  


