# Simulating magic circuits with Clifford-augmented MPS
The stabilizer formalism allows a subset of quantum circuits, the *Clifford circuits* made up of Hadamard, CNOT and phase gates, to be efficiently simulated on classical computers. However, universal computation is only possible by adding T gates to the set of simulated operations. The CAMPS (Clifford-augmented matrix product state) ansatz for quantum states joins the computational power of the stabilizer formalism and matrix product states, but in a naive implementation still struggles with simulating T gates due to an exponential increase in the reuqired bond dimension.

This project, a stepping stone to my Master's thesis, implements an analytical disentangling algorithm that keeps the MPS bond dimension to 1 for the first N simulated T gates, as described in the papers by Fux *et al.* (2025) and by Liu and Clark (2025).

## References
G. E. Fux, B. Béri, R. Fazio, and E. Tirrito, “Disentangling Magic States with Classically Simulable Quantum Circuits,” Phys. Rev. Lett., vol. 135, no. 26, Dec. 2025, doi: [10.1103/ggp1-byj1](https://doi.org/10.1103/ggp1-byj1). 

Z. Liu and B. K. Clark, “Classical simulability of Clifford+T circuits with Clifford-augmented matrix product states,” Aug. 26, 2025, arXiv:2412.17209. doi: [10.48550/arXiv.2412.17209](https://doi.org/10.48550/arXiv.2412.17209).
