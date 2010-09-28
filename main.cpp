// High Performance Hartree Fock
// Based on / inspired by:
// Szabo and Ostlund, "Modern Quantum Chemistry"

// STATUS : Hardcoded for H2

// C includes
#include <cstdio>
#include <cmath>
// C++ includes
#include <iostream>
#include <string>
#include <vector>
//#include <Eigen/Eigen>

using namespace std;
//using namespace Eigen;

struct atom
{
	string element;
	double charge;
	double position[3];
	double Nelec;
	vector<double> coeff;
	vector<double> exponent;
};

int main() {
	
	// The SCF procedure is a 12 step process, because anything less is for alcoholics.
	//
	// 1:
	// Specify a set of nuclear coordinates (R_A), atomic numbers
	// (Z_A), number of electrons (N), and a basis set (phi_m)
	
	// I will hardcode in H2 to start, 
	int NAtoms;
	NAtoms = 2;
	double charge;
	charge = 0;

	atom atoms[NAtoms];

	atoms[0].element = "Hydrogen";
	atoms[0].charge = 1;
	atoms[0].Nelec = 1;
	atoms[0].position[0] = 0;
	atoms[0].position[1] = 0;
	atoms[0].position[2] = 0;

	atoms[1].element = "Hydrogen";
	atoms[1].charge = 1;
	atoms[1].Nelec = 1;
	atoms[1].position[0] = 0.75;
	atoms[1].position[1] = 0;
	atoms[1].position[2] = 0;

	// Robot Roll-Call
	for (int i = 0; i < NAtoms; i++) {
		cout << "Atom " <<  i+1 << ":\n"
		     << "Name: " <<  atoms[i].element << "\n"	
		     << "Charge: " << atoms[i].charge << "\n"
		     << "Position: " << atoms[i].position[0] << " " 
		     << atoms[i].position[1] << " " << atoms[i].position[2]
		     << "\n   ---   \n"
		     ;
	}
	
	double N_total;
	N_total = 0;
	for (int i = 0; i < NAtoms; i++){
		N_total = N_total + atoms[i].Nelec;
	}
	N_total = N_total + charge;
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
