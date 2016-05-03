# XYZ_plaquette_deco_SSE

Performs SSE for the XYZ model on the Kagome lattice

-lattice_so.f Generates lattices
-dloop.f90 performs a Stochastic Series Expansion (SSE) QMC algorithm in the global Sz basis designed 
to study the XYZ Hamiltonian at finite temperature using a 2+1-dimensional simulation cell. 
Within the SSE, the Hamiltonian was implemented with a triangular plaquette breakup, which helps 
ergodicity in the regime where Jz/J±± is large. Using this Hamiltonian breakup, the standard SSE-directed 
loop equations were modified to include sampling of off-diagonal operators of the type. The resulting 
algorithm is highly efficient, scaling linearly in the number of lattice sites V and the inverse 
temperature β. This scaling is modified to V2β in the cases where a full q-dependent structure factor 
measurement is required. The program was implemented in Fortran and verified by comparing results for 
small clusters with exact diagonalization data  
