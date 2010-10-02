// High Performance Hartree Fock
// Based on / inspired by:
// Szabo and Ostlund, "Modern Quantum Chemistry"

// STATUS : Hardcoded for H2

// C includes
#include <cstdio>
#include <cmath>
// C++ includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

const double pi = 3.1415926535898;


struct atom
{
	string element;
	double charge;
	VectorXd position;
	int Nelec;
	VectorXd coeff;    
	VectorXd exponent;
};

// Include our functions
#include "functions.h"

int main() {
	bool printinfo;
	printinfo = false;
	
	// The SCF procedure is a 12 step process, because anything less is for alcoholics.
	
	
	// 1:
	// Specify a set of nuclear coordinates (R_A), 
	// atomic numbers (Z_A),
	// number of electrons (N), 
	// and a basis set (phi_m)
	
	// I will hardcode in H2 to start, 
	int NAtoms;
	NAtoms = 2;
	double charge;
	charge = 0;

	atom atoms[NAtoms];

	atoms[0].element = "Hydrogen";
	atoms[0].charge = 1;
	atoms[0].Nelec = 1;
	atoms[0].position = VectorXd::Zero(3);
	atoms[0].position << 0, 0, 0;

	atoms[1].element = "Hydrogen";
	atoms[1].charge = 1;
	atoms[1].Nelec = 1;
	atoms[1].position = VectorXd::Zero(3);
	atoms[1].position << 0.75, 0, 0;

	// Robot Roll-Call
	if (printinfo){
	for (int i = 0; i < NAtoms; i++) {
		cout << "Atom " <<  i+1 << ":\n"
		     << "Name: " <<  atoms[i].element << "\n"	
		     << "Charge: " << atoms[i].charge << "\n"
		     << "Position: " << atoms[i].position.transpose() // << " " 
		     //<< atoms[i].position[1] << " " << atoms[i].position[2]
		     << "\n   ---   \n"
		     ;
	}
	}

	// Get the total number of electrons	
	int N_total;
	N_total = 0;
	for (int i = 0; i < NAtoms; i++){
		N_total = N_total + atoms[i].Nelec;
	}
	N_total = N_total + charge;

	// We are working with CLOSED SHELL HF at the moment
	if ( N_total % 2 != 0){
		cout << "You do not have an even number of electrons. \nAborting";
			return(1);
	}
	// Tell the user something
	if (printinfo){
	cout << "The Number of electrons is " << N_total << "\n";
	}

	// Lets Hard-Code in the STO-6G basis set for Hydrogen.
	// Taken from https://bse.pnl.gov/bse/portal
	for (int i = 0; i < NAtoms; i++){
		atoms[i].exponent = VectorXd::Zero(6);
		atoms[i].coeff = VectorXd::Zero(6);
		atoms[i].exponent << 35.52322122, 6.513143725, 1.822142904, 0.625955266, 0.243076747, 0.100112428;
		atoms[i].coeff << 0.00916359628, 0.04936149294, 0.16853830490, 0.37056279970, 0.41649152980, 0.13033408410;
	}
	
	// The STO basis sets from the EMSL have had their exponents scaled, but not their coefficients.
	// This has been hard-coded for 6 (j loop). Will need to be changed later.
	for (int i = 0; i < NAtoms; i++) {
		for (int j = 0; j < 6; j++){
			atoms[i].coeff[j] = atoms[i].coeff[j] * 
				pow((2 * atoms[i].exponent[j] / pi), 0.75);
		}
	}

	// Spit something out
	if (printinfo){
	cout << "Coefficients:\n" << atoms[0].coeff.transpose() << "\n";
	cout << "Exponents:\n" << atoms[0].exponent.transpose() << "\n";
	}
	if (printinfo){
		cout.precision(10);
		//I need to coeff. to higher precison for comparison
		for (int i=0; i<6; i++){
			cout << atoms[0].coeff[i] << endl;
		}
	}


	//
	// 2:
	// Calculate the required molecular integrals: S_mn, Hcore_mn, and (mn|ls)
	
	// The STO-6G basis set is a linear combination of 6 Gaussians,
	// but technically it is only one basis function. Therefore our
	// system here will have 2 total basis functions (i.e. a 2x2)
	int NBasis (2);


	// Set up our overlap Matrix
	MatrixXd S(NBasis,NBasis);

	for (int i = 0; i < NBasis; i++){
		for (int j = 0; j < NBasis; j++){
			S(i,j) = overlap(atoms[i].coeff, atoms[i].exponent, 
					 atoms[j].coeff, atoms[j].exponent, 
					 atoms[i].position[0], atoms[j].position[0] );
		}
	}

	// Use SVD do diagonalize S
	bool solved;
	VectorXd b(NBasis);
	VectorXd x(NBasis);
	solved =  S.svd().solve(b, &x); // Test to make sure it is solvable
	if (! solved){
		cout << "Cannot diagonalize overlap matrix!!!";
		return(0);
	}
	SVD<MatrixXd> lala(S);
	lala.solve(b, &x);
	cout << lala.matrixU() << endl;
	cout << lala.matrixV() << endl;
	cout << lala.matrixV()*S*lala.matrixV() << endl;
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
