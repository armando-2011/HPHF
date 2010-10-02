// Functions used

double overlap(VectorXd coeffA, VectorXd expA, 
	       VectorXd coeffB, VectorXd expB,
	       double centerA, double centerB){
	double SAB (0);
	// We are working with Gaussians*Gaussians, 
	// so if we know the exponenets and coefficents we have an
	// analytical expression for the integral.
	
	//Find out how many Gaussians are on centers A and B
	int NGA, NGB;
	NGA = coeffA.rows();
	NGB = coeffB.rows();

	double R_BA;
	R_BA = centerB - centerA;
	// Set some temp variables
	double d, f, g, h;
	for (int i = 0; i < NGA; i++){
		for (int j = 0; j < NGB; j++){
			d = coeffA[i]*coeffB[j];
			f = (expA[i] + expB[j]);
			g = -2*expB[j]*R_BA;
			h = -expB[j]*R_BA*R_BA;
			SAB = SAB + d*sqrt(pi/f)*exp((g*g/(4*f)) + h);
		}
	}

	// Results for A=B on mathematica = 0.485213642 8094333
	// Results from this algorithm    = 0.485213642 7543415
	return(SAB);
}


//////////////////////////////////////////////////////
// Calculate the f function FO only (s-type orbitals)
// Used in V() and TWOE()
double FO(double arg){
	double F;
	if (arg < 1e-6){
		F = 1 - arg/3;
	}
	else {
		F = sqrt(pi/arg)*erf(sqrt(arg))/2;
	}
	return F;
}

// Calculate the overlaps of un-normalized primitives
double S(double A, double B, double R_ab){
	double overlap;
	overlap = pow((pi/(A + B)),1.5) * exp(-A*B*R_ab/(A + B));
	return(overlap);
}

// Calculate the kinteic energy integrals for un-normalized primitives
double T(double A, double B, double R_ab){
	double kinetic;
	kinetic = A*B/(A+B)*(3 - 2*A*B*R_ab/(A+B))*pow(pi/(A+B),1.5);
	return kinetic;
}

// Calculate the un-normalized nuclear attraction integrals
double V(double A, double B, double R_ab, double Rcp2, double Zc){
	double N;
	N = 2*pi/(A+B)*FO((A+B)*Rcp2)*exp(-A*B*R_ab/(A+B));
	N = -N*Zc;
	return N;
}

// Two electron integrals
double TWOE(double A, double B, double C, double D,double Rab2, double Rcd2, double Rpq2){
	double integral;
	integral = 2*pow(pi,2.5) / 
		( (A+B) * (C+B)*sqrt(A+B+C+D)) *
		FO( (A+B) * (C+D) *Rpq2 / (A+B+C+D)) *
		exp( -A*B*Rab2/(A+B) - C*D*Rcd2/(C+D));
	return integral;
}
