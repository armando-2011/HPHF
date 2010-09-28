//High Performance Hartree Fock
//Based on / inspired by:
//Szabo and Ostlund, "Modern Quantum Chemistry"

#include <stdio.h>
#include<math.h>


int main() {
	
	// The SCF procedure is a 12 step process, because anything less is for alcoholics.
	//
	// 1:
	// Specify a set of nuclear coordinates (R_A), atomic numbers
	// (Z_A), number of electrons (N), and a basis set (phi_m)
	//
	// 2:
	// Calculate the required molecular integrals: S_mn, Hcore_mn, and (mn|ls)
	//
	// 3:
	// Diagonalize the overlap matrix S, and obtain a transforamtion matrix X using either
	// X = Sqrt(X) = U Sqrt(s) U_dagger
	// X = U Sqrt(s)
	// where s is the diagonalized form of S
	//
	// 4:
	// Obtain a guess for the Density Matrix P
	//
	// 5:
	// Calculate G from from P and (mn|ls)
	// G_mn = SUM_ls P_ls [ (mn|sl) - 0.5 (ml|sn) ]
	//
	// 6:
	// Add G to Hcore to get the Fock matrix F
	//
	// 7:
	// Calculate the transformed Fock matrix
	// F' = X_dagger F X
	//
	// 8:
	// Diagonalize F' to obtain C' and E
	//
	// 9:
	// Calcualte C
	// C = X C'
	//
	// 10:
	// Form a new density matrix P from C
	// P_mn = 2 SUM_a^(N/2) C_ma Cstar_na
	//
	// 11:
	// Check if the new density matrix is the same as the old density matrix (within a threshold).
	// If yes, move on. If not, go back to step 5 with the new P.
	//
	// 12:
	// After the proces is converged, use C, P, and F to calcualte
	// expectation values and other quantities of interest

	
	return 0;
}
