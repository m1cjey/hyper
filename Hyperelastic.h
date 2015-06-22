#ifndef HYPERELASTIC
#define HYPERELASTIC

class hyperelastic
{
public:

	unsigned ID;
	int NEI[200];
	int NEI_w[200];
	int N;
	int N_w;
	double lambda;
	double n_q[DIMENSION];
	double half_p[DIMENSION];
	double stress[DIMENSION][DIMENSION];
	double P_wall[DIMENSION];
	double differential_p[DIMENSION];
	double p[DIMENSION];
	double ang_p[DIMENSION];
	double Ai[DIMENSION][DIMENSION];
	double inverse_Ai[DIMENSION][DIMENSION];
	double t_inverse_Ai[DIMENSION][DIMENSION];
	double Fi[DIMENSION][DIMENSION];
	double t_inverse_Fi[DIMENSION][DIMENSION];
	double J;
	double differential_Fi[DIMENSION][DIMENSION];
	double t_differential_Fi[DIMENSION][DIMENSION];
	double vis_force[DIMENSION];
	double old_F[DIMENSION];
	int stress0;
	int highest;	//à¯Ç¡í£ÇËééå±âêÕóp15/2/8




};

class hyperelastic2
{
public:
	unsigned ID;
	double DFDlambda;
	double wiin;
	double DgDq[DIMENSION];
	double newton_DgDq[DIMENSION];
	double aiin[DIMENSION];
	double n0ij[DIMENSION];
	double spl_f;
};

#endif