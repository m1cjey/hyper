#include "stdafx.h"	

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,bool repetation,double **F);
void renew_lambda(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t);
void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double **F);
void calc_transposed_inverse_matrix(double **M,bool transport,bool inversion);
double calc_det(double **M,int N);
double calc_det3(double **M);
void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER);
void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void calc_inverse_matrix_for_NR(int N, double *a);
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F);
void calc_F(vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1);
void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t,double **F);
void inverse(double **a,int N);
void ludcmp(double **a,int N,int *index,double *d);
void lubksb(double **a,int N,int *index,double b[]);
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,double **F);
void contact_judge_hyper(mpsconfig &CON,vector<mpselastic> &PART, vector<hyperelastic> &HYPER,int t);
void contact_judge_hyper2(mpsconfig CON, vector<mpselastic> &PART, vector<hyperelastic> &HYPER, int hyper_number, int t);
void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1, int t);
void transpose(double **M,double **N);
void output_newton_data1(double *fx, double *DfDx, double *n_rx, double *n_ry, double *n_rz, int hyper_number, int count, int t);
void output_newton_data2(double E, double *XX, int hyper_number, int count, int t);
void calc_gravity(mpsconfig CON,vector<hyperelastic> &HYPER,int hyper_number);
void calculation_vec_norm(vector<mpselastic> PART, vector<hyperelastic> &HYPER, int hyper_number,int particle_number,int t);
void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t);
void contact_judge(mpsconfig &CON, vector<mpselastic> PART,vector<hyperelastic> &HYPER,double max_h,int t);
void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void iccg2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);
void CG3D(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X);


void calc_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int t,double **F)
{	
	int h_num=HYPER.size();
	int p_num=PART.size();
	double Dt=CON.get_dt();
	cout<<"Hypercalculation starts."<<endl;
	double max_h=PART[p_num-h_num].r[A_Z];
	//calc_gravity(CON,HYPER,h_num);
	
//	contact_judge(CON,PART,HYPER,max_h,t);

	if(t==1)
	{
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=0;
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=PART[i].r[D];
		for(int k=p_num-h_num+1;k<p_num;k++)	if(max_h<PART[k].r[A_Z])	max_h=PART[k].r[A_Z];
		calc_constant(CON,PART,HYPER,HYPER1);
		calc_stress(CON,HYPER);
	}
	
	if(CON.get_flag_vis()==ON)	calc_vis_f(CON,PART,HYPER,HYPER1,0,t);

	if(t==1 || t%CON.get_interval()==0)
	{
		output_hyper_data(PART,HYPER,HYPER1,t);
		momentum_movie_AVS(CON,t,PART,HYPER,F);
		output_energy(CON,PART,HYPER,t);
	}


	newton_raphson(CON,PART,HYPER,HYPER1,t,F);

	calc_half_p(CON,PART,HYPER,HYPER1,0,F);

	calc_F(PART,HYPER,HYPER1);

	calc_stress(CON,HYPER);
	
	calc_differential_p(CON,PART,HYPER,HYPER1,F);

	renew_lambda(CON,HYPER,HYPER1,t);

	calc_half_p(CON,PART,HYPER,HYPER1,1,F);

	cout<<"Hypercalculation ends."<<endl;
}

void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
//	cout<<"初期値計算";

	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);	//考慮が必要かもしれない
	double mi=V*CON.get_hyper_density();
	int h_num=HYPER.size();
	int model=CON.get_model_number();
	cout<<"V"<<V<<endl;

	////初期運動量
	//曲げねじり
	/*if(model==21)
	{		
		int t=30,b=2;
		double max=0,min=0;

		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double Y=PART[i].q0[A_Y];
			double X=PART[i].q0[A_X];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=mi*(-t*part_p*part_p*part_p*Y+b*(3*part_p*part_p-1));
			HYPER[i].p[A_Y]=mi*t*part_p*part_p*part_p*X;
			HYPER[i].p[A_Z]=0;
		}
	}*/

	//曲げ
	if(model==21)
	{
		int b=2;
		double max=0,min=0;

		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=mi*b*(3*part_p*part_p-1);
		}
	}//*/

	//ねじり
	/*if(model==21)
	{
		int t=30;
		double max=0,min=0;

		for(int i=0;i<h_num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
		}
		double H=max-min+le;
		cout<<H;
		//double H=1.8;
		for(int i=0;i<h_num;i++)	
		{
			double Z=PART[i].q0[A_Z];
			double Y=PART[i].q0[A_Y];
			double X=PART[i].q0[A_X];
			double part_p=(Z/H)*2;
			HYPER[i].p[A_X]=-mi*t*part_p*part_p*part_p*Y;
			HYPER[i].p[A_Y]=mi*t*part_p*part_p*part_p*X;
			HYPER[i].p[A_Z]=0;
		}
	}//*/


	//回転
	if(model==22)
	{
		for(int i=0;i<h_num;i++)
		{
			PART[i].p[A_X]=mi*0.4*(PART[i].q0[A_Z]-PART[i].q0[A_Y]);
			PART[i].p[A_Y]=mi*0.4*(PART[i].q0[A_X]-PART[i].q0[A_Z]);
			PART[i].p[A_Z]=mi*0.4*(PART[i].q0[A_Y]-PART[i].q0[A_X]);
		}
	}

	////角運動量計算
	for(int i=0;i<h_num;i++)
	{
		HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
		HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
		HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
	}


	////近傍粒子の記憶とaiin,wiin,Aiの計算
	for(int i=0;i<h_num;i++)
	{
		int N=0;
		double dis=0;
		double wiin=0;
		double aiin[DIMENSION];
		for(int j=0;j<h_num;j++)
		{
			wiin=0;
			aiin[A_X]=PART[j].q0[A_X]-PART[i].q0[A_X];	aiin[A_Y]=PART[j].q0[A_Y]-PART[i].q0[A_Y];	aiin[A_Z]=PART[j].q0[A_Z]-PART[i].q0[A_Z];
			
			HYPER1[i*h_num+j].aiin[A_X]=aiin[A_X];
			HYPER1[i*h_num+j].aiin[A_Y]=aiin[A_Y];
			HYPER1[i*h_num+j].aiin[A_Z]=aiin[A_Z];

			dis=sqrt(aiin[A_X]*aiin[A_X]+aiin[A_Y]*aiin[A_Y]+aiin[A_Z]*aiin[A_Z]);
			if(dis<r && j!=i)
			{	
				wiin=kernel4(r,dis);
				HYPER[i].NEI[N]=j;
				N++;
			}
			else	wiin=0;

			HYPER[i].pnd+=wiin;
			HYPER1[i*h_num+j].wiin=wiin;
		}
		HYPER[i].N=N;
	}
	
	////Ai, Fi関連の計算
	double **p_Ai=new double *[DIMENSION];
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Ai[D]=new double [DIMENSION];
		p_Fi[D]=new double[DIMENSION];
	}
	for(int i=0;i<h_num;i++)
	{
		//Aiの計算
		double Ai[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		Ai[0][0]=0;	Ai[0][1]=0;	Ai[0][2]=0;
		Ai[1][0]=0;	Ai[1][1]=0;	Ai[1][2]=0;
		Ai[2][0]=0;	Ai[2][1]=0;	Ai[2][2]=0;

		int Ni=HYPER[i].N;
		double a[DIMENSION];
		for(int j=0;j<Ni;j++)	
		{
			int k=HYPER[i].NEI[j];
			double w=HYPER1[i*h_num+k].wiin;		
			double a[DIMENSION]={HYPER1[i*h_num+k].aiin[A_X],	a[A_Y]=HYPER1[i*h_num+k].aiin[A_Y], a[A_Z]=HYPER1[i*h_num+k].aiin[A_Z]};
			
			Ai[0][0]+=w*a[0]*a[0];		Ai[0][1]+=w*a[0]*a[1];		Ai[0][2]+=w*a[0]*a[2];
			Ai[1][0]+=w*a[1]*a[0];		Ai[1][1]+=w*a[1]*a[1];		Ai[1][2]+=w*a[1]*a[2];
			Ai[2][0]+=w*a[2]*a[0];		Ai[2][1]+=w*a[2]*a[1];		Ai[2][2]+=w*a[2]*a[2];
		}
		HYPER[i].Ai[0][0]=Ai[0][0];	HYPER[i].Ai[0][1]=Ai[0][1];	HYPER[i].Ai[0][2]=Ai[0][2];
		HYPER[i].Ai[1][0]=Ai[1][0];	HYPER[i].Ai[1][1]=Ai[1][1];	HYPER[i].Ai[1][2]=Ai[1][2];
		HYPER[i].Ai[2][0]=Ai[2][0];	HYPER[i].Ai[2][1]=Ai[2][1];	HYPER[i].Ai[2][2]=Ai[2][2];
		
		//inverse_Ai,t_inverse_Aiの計算
		p_Ai[0][0]=Ai[0][0];	p_Ai[0][1]=Ai[0][1];	p_Ai[0][2]=Ai[0][2];
		p_Ai[1][0]=Ai[1][0];	p_Ai[1][1]=Ai[1][1];	p_Ai[1][2]=Ai[1][2];
		p_Ai[2][0]=Ai[2][0];	p_Ai[2][1]=Ai[2][1];	p_Ai[2][2]=Ai[2][2];
		inverse(p_Ai,DIMENSION);
		HYPER[i].inverse_Ai[0][0]=p_Ai[0][0];		HYPER[i].inverse_Ai[0][1]=p_Ai[0][1];		HYPER[i].inverse_Ai[0][2]=p_Ai[0][2];
		HYPER[i].inverse_Ai[1][0]=p_Ai[1][0];		HYPER[i].inverse_Ai[1][1]=p_Ai[1][1];		HYPER[i].inverse_Ai[1][2]=p_Ai[1][2];
		HYPER[i].inverse_Ai[2][0]=p_Ai[2][0];		HYPER[i].inverse_Ai[2][1]=p_Ai[2][1];		HYPER[i].inverse_Ai[2][2]=p_Ai[2][2];

		HYPER[i].t_inverse_Ai[0][0]=p_Ai[0][0];		HYPER[i].t_inverse_Ai[0][1]=p_Ai[1][0];		HYPER[i].t_inverse_Ai[0][2]=p_Ai[2][0];
		HYPER[i].t_inverse_Ai[1][0]=p_Ai[0][1];		HYPER[i].t_inverse_Ai[1][1]=p_Ai[1][1];		HYPER[i].t_inverse_Ai[1][2]=p_Ai[2][1];
		HYPER[i].t_inverse_Ai[2][0]=p_Ai[0][2];		HYPER[i].t_inverse_Ai[2][1]=p_Ai[1][2];		HYPER[i].t_inverse_Ai[2][2]=p_Ai[2][2];
		
		//Fiの計算
		p_Fi[0][0]=Ai[0][0]*p_Ai[0][0]+Ai[0][1]*p_Ai[1][0]+Ai[0][2]*p_Ai[2][0];	p_Fi[0][1]=Ai[0][0]*p_Ai[0][1]+Ai[0][1]*p_Ai[1][1]+Ai[0][2]*p_Ai[2][1];	p_Fi[0][2]=Ai[0][0]*p_Ai[0][2]+Ai[0][1]*p_Ai[1][2]+Ai[0][2]*p_Ai[2][2];
		p_Fi[1][0]=Ai[1][0]*p_Ai[0][0]+Ai[1][1]*p_Ai[1][0]+Ai[1][2]*p_Ai[2][0];	p_Fi[1][1]=Ai[1][0]*p_Ai[0][1]+Ai[1][1]*p_Ai[1][1]+Ai[1][2]*p_Ai[2][1];	p_Fi[1][2]=Ai[1][0]*p_Ai[0][2]+Ai[1][1]*p_Ai[1][2]+Ai[1][2]*p_Ai[2][2];
		p_Fi[2][0]=Ai[2][0]*p_Ai[0][0]+Ai[2][1]*p_Ai[1][0]+Ai[2][2]*p_Ai[2][0];	p_Fi[2][1]=Ai[2][0]*p_Ai[0][1]+Ai[2][1]*p_Ai[1][1]+Ai[2][2]*p_Ai[2][1];	p_Fi[2][2]=Ai[2][0]*p_Ai[0][2]+Ai[2][1]*p_Ai[1][2]+Ai[2][2]*p_Ai[2][2];		

		HYPER[i].Fi[0][0]=p_Fi[0][0];	HYPER[i].Fi[0][1]=p_Fi[0][1];	HYPER[i].Fi[0][2]=p_Fi[0][2];	
		HYPER[i].Fi[1][0]=p_Fi[1][0];	HYPER[i].Fi[1][1]=p_Fi[1][1];	HYPER[i].Fi[1][2]=p_Fi[1][2];	
		HYPER[i].Fi[2][0]=p_Fi[2][0];	HYPER[i].Fi[2][1]=p_Fi[2][1];	HYPER[i].Fi[2][2]=p_Fi[2][2];	
		//Jの計算
		double J=calc_det3(p_Fi);
		HYPER[i].J=J;
//		cout<<"HYPER["<<i<<"].J="<<HYPER[i].J<<endl;

		//t_inverse_Fiの計算
		inverse(p_Fi,DIMENSION);
		HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];		HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];		HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];		HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];		HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];		HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];		HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];
	}
	for(int D=0;D<DIMENSION;D++)
	{
		delete[]	p_Ai[D];
		delete[]	p_Fi[D];
	}
	delete[]	p_Ai;
	delete[]	p_Fi;


	////n0ijの計算
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N;

		double p_n0ij[DIMENSION]={0,0,0};
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];			
			p_n0ij[A_X]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_X];
			p_n0ij[A_Y]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_Y];
			p_n0ij[A_Z]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[A_Z];

			HYPER1[i*h_num+k].n0ij[A_X]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_X][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_X][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_X][2]*HYPER1[i*h_num+k].aiin[2]);
			HYPER1[i*h_num+k].n0ij[A_Y]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_Y][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_Y][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_Y][2]*HYPER1[i*h_num+k].aiin[2]);
			HYPER1[i*h_num+k].n0ij[A_Z]=V*HYPER1[k*h_num+i].wiin*(HYPER[k].t_inverse_Ai[A_Z][0]*HYPER1[i*h_num+k].aiin[0]+HYPER[k].t_inverse_Ai[A_Z][1]*HYPER1[i*h_num+k].aiin[1]+HYPER[k].t_inverse_Ai[A_Z][2]*HYPER1[i*h_num+k].aiin[2]);
		}
		HYPER1[i*h_num+i].n0ij[A_X]=V*(HYPER[i].t_inverse_Ai[A_X][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_X][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_X][2]*p_n0ij[2]);
		HYPER1[i*h_num+i].n0ij[A_Y]=V*(HYPER[i].t_inverse_Ai[A_Y][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_Y][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_Y][2]*p_n0ij[2]);
		HYPER1[i*h_num+i].n0ij[A_Z]=V*(HYPER[i].t_inverse_Ai[A_Z][0]*p_n0ij[0]+HYPER[i].t_inverse_Ai[A_Z][1]*p_n0ij[1]+HYPER[i].t_inverse_Ai[A_Z][2]*p_n0ij[2]);
	}
		
	////DgDqの計算
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];
			HYPER1[k*h_num+i].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+k].n0ij[2]);		
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
	}	
//	cout<<"----------OK"<<endl;
 }


/////ニュートンラフソン法 
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t,double **F)
{
	/////fx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照
	/////DfDx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照

	int calc_type=1;//ニュートラフソンの反復方法 0:偏微分項の逆行列をそのまま求める　1:線形方程式を利用

	//OPENMPがオンであれば書き込まれる。threadsはpcで扱える最大スレッド数
	#ifdef _OPENMP
    printf("OpenMP : On, threads = %d\n", omp_get_max_threads());
	#endif

	//最大スレッド数で計算し続けると過負荷でCPUが非常に熱くなるため、並列化数を指定している
	//なお、最大スレッドが12のときは8～10ぐらいが目安
	omp_set_num_threads(8);
	//pn=2;//test,とりあえず2元でとけるかどうか確認 
	//////////////////　f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm

	int h_num=HYPER.size();
	double *fx=new double [h_num];//関数値。
	double *DfDx=new double [h_num*h_num];//関数の偏微分値。
	double *XX=new double [h_num];//現在の解。	
	double *XX_old=new double [h_num];//1ステップ前の解。
	double ep=1e-5;//収束判定
	double E=1;//現在の誤差
/*	double start=0;
	double end=0;
	double newton_t=0;*/
	int count=0;//反復回数
	double d;
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double Dt=CON.get_interval();
	double sum=0;
	double E_old=0;
	int dec_flag=ON;

	#pragma omp parallel for
	for(int i=0;i<h_num;i++)
	{
		XX[i]=1;
		XX_old[i]=0;
		fx[i]=0;
		for(int j=0;j<h_num;j++)	DfDx[i*h_num+j]=0;
	}


	//	for(int i=0; i<N; i++) XX[i]=1;///初期値を与える。とりあえず1で
	cout<<"NR法開始";
//	start=clock();
	while(E>ep)
	{
		count++;
		for(int i=0; i<h_num; i++)	XX_old[i]=XX[i];	//解を記憶

//		if(count==1)	for(int i=0;i<h_num;i++)	for(int j=0;j<h_num;j++)	for(int D=0;D<DIMENSION;D++)	HYPER1[i*h_num+j].newton_DgDq[D]=HYPER1[i*N+j].DgDq[D];

		calc_newton_function(CON,PART,HYPER,HYPER1,XX,fx,DfDx,h_num,count,t,F);


/*		//現在の関数値を求める
		if(count==1) cout<<fx[0]<<" "<<fx[1]<<endl;
		//現在の偏微分値を求める
		//calc_DfDx(XX)////現在の偏微分値を求める。超弾性体ならば、calc_DgDq()などで求められるはず
		DfDx[0*N+0]=2*XX[0];
		DfDx[0*N+1]=2*XX[1];
		DfDx[1*N+0]=2*XX[0]/9;
		DfDx[1*N+1]=2*XX[1];
		if(count==1) cout<<DfDx[0]<<" "<<DfDx[1]<<" "<<DfDx[2]<<" "<<DfDx[3]<<endl;*/

		///値の更新
		if(calc_type==0)//逆行列を利用 逆行列が求まりさえすれば速いはず
		{
			calc_inverse_matrix_for_NR(h_num,DfDx);

			for(int i=0; i<h_num; i++) 
			{
				d=0; //変化量
				for(int j=0; j<h_num; j++)	d+=DfDx[i*h_num+j]*fx[j];
				XX[i]-=d;
			}
		}
		else if(calc_type==1)//逆行列を用いない、安定するはずだが、遅くなるはず
		{	
			
			int *b_ind=new int [h_num*h_num];
			double *b_val=new double [h_num*h_num];
			int *ptr=new int [h_num+1];
			int all_ind_num=0;
			
			ptr[0]=0;
			for(int i=0;i<h_num;i++)
			{
				ptr[i+1]=0;
				for(int j=0;j<h_num;j++)
				{
					b_ind[i*h_num+j]=0;
					b_val[i*h_num+j]=0;
				}
			}

			
			for(int i=0;i<h_num;i++)
			{
				int ind_num=0;
				for(int j=0;j<h_num;j++)
				{
					if(DfDx[i*h_num+j]!=0)
					{
						b_ind[all_ind_num+ind_num]=j;	
						b_val[all_ind_num+ind_num]=DfDx[i*h_num+j];
						ind_num++;
						ptr[i+1]=ind_num+ptr[i];
					}			
				}
				//cout<<"ptr"<<i+1<<" "<<ptr[i+1]<<endl;
				all_ind_num+=ind_num;			
			}

			cout<<all_ind_num<<endl;
			int *ind=new int [all_ind_num];
			double *val=new double [all_ind_num];
			for(int i=0;i<all_ind_num;i++)
			{
				ind[i]=b_ind[i];
				val[i]=b_val[i];
			}
			//for(int i=0;i<h_num;i++)	cout<<"ptr"<<i<<"	"<<ptr[i]<<endl;
			CG3D(&CON,val,ind,ptr,h_num,fx,all_ind_num,XX);
			///iccg2(&CON,val,ind,ptr,h_num,fx,all_ind_num,XX);
			//BiCGStab2_method(&CON,val,ind,ptr,h_num,fx,all_ind_num,XX);
			//for(int i=0;i<h_num;i++)	cout<<"lambda"<<i<<"="<<XX[i]<<endl;

			delete[] b_ind;
			delete[] b_val;
			delete[] ind;
			delete[] val;
			delete[] ptr;

			/*gauss(DfDx,fx,h_num);
			for(int i=0;i<h_num;i++)	XX[i]-=fx[i];//*0.5*mi/(Dt*Dt)*V*fx[i];*/
		}

		//誤差の評価
		E_old=E;
		E=0;
		sum=0;
		for(int i=0; i<h_num; i++)
		{
			E+=fabs(XX[i]-XX_old[i]);
			sum+=fabs(XX[i]);
		}
		E/=sum;


		if(count==1 || count%200==0)
		{		
			/*
			cout<<"XX_old["<<i<<"]-d=X["<<i<<"]	";
			cout<<XX_old[i]<<" - ";
			cout<<d<<" = ";
			cout<<XX[i]<<endl;*/

			cout<<"反復回数	"<<count<<" E="<<E<<endl;
			//output_newton_data2(E,XX,h_num,count,t);

		}
		if(count>CON.get_nr())	break;
		else if(dec_flag==ON)	if(E_old-E<0)	break;	
	}
//	end=clock();
//	newton_t=(end-start)/CLOCKS_PER_SEC;

	cout<<"反復完了";
	#pragma omp parallel for
	for(int i=0;i<h_num;i++) HYPER[i].lambda=XX[i];

//	for(int i=0;i<N;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
	delete[]	fx;
	delete[]	DfDx;
	delete[]	XX;
	delete[]	XX_old;
	cout<<"---------- OK"<<endl;
}


void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t,double **F)
{
	clock_t t3=clock();

	int h_num=hyper_number;
	int flag_vis=CON.get_flag_vis();
	bool flag_FEM=CON.get_FEM_flag();
	int flag_G=CON.get_flag_G();
	int model_num=CON.get_model_number();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double density=CON.get_hyper_density();
	double *n_rx=new double[h_num];
	double *n_ry=new double[h_num];
	double *n_rz=new double[h_num];

	double **n_DgDq_x=new double *[h_num];
	double **n_DgDq_y=new double *[h_num];
	double **n_DgDq_z=new double *[h_num];
	for(int i=0;i<h_num;i++)
	{
		n_rx[i]=0;
		n_ry[i]=0;
		n_rz[i]=0;
		n_DgDq_x[i]=new double [h_num];
		n_DgDq_y[i]=new double [h_num];
		n_DgDq_z[i]=new double [h_num];
	}
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			n_DgDq_x[i][j]=0;
			n_DgDq_y[i][j]=0;
			n_DgDq_z[i][j]=0;
		}
	}

	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double [DIMENSION];

	////位置座標の更新	
	for(int i=0;i<h_num;i++)
	{
/*		if(model_num==30||model_num==23)
		{
			if(PART[i].q0[A_Z]!=0)
			{
				//half_pの計算
				double p_half_p[DIMENSION]={0,0,0};

				int Ni=HYPER[i].N;
				for(int j=0;j<Ni;j++)
				{	
					int k=HYPER[i].NEI[j];
					p_half_p[A_X]+=(HYPER[k].stress[0][0]-lambda[k])*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[0][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[0][2]*HYPER1[k*h_num+i].DgDq[2];
					p_half_p[A_Y]+=HYPER[k].stress[1][0]*HYPER1[k*h_num+i].DgDq[0]+(HYPER[k].stress[1][1]-lambda[k])*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[1][2]*HYPER1[k*h_num+i].DgDq[2];
					p_half_p[A_Z]+=HYPER[k].stress[2][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[2][1]*HYPER1[k*h_num+i].DgDq[1]+(HYPER[k].stress[2][2]-lambda[k])*HYPER1[k*h_num+i].DgDq[2];
				}//jに関するfor文の終わり
				p_half_p[A_X]+=(HYPER[i].stress[0][0]-lambda[i])*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
				p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-lambda[i])*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
				p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-lambda[i])*HYPER1[i*h_num+i].DgDq[2];
		
				//重力の影響
				if(flag_G==ON)	p_half_p[A_Z]-=9.8*mi;
				//粘性項の影響
				if(flag_vis==ON)
				{
					p_half_p[A_X]+=density*HYPER[i].vis_force[A_X];
					p_half_p[A_Y]+=density*HYPER[i].vis_force[A_Y];
					p_half_p[A_Z]+=density*HYPER[i].vis_force[A_Z];
				}
				//磁場の考慮
				if(flag_FEM==ON)
				{
					p_half_p[A_X]+=F[A_X][i]*V;//density;
					p_half_p[A_Y]+=F[A_Y][i]*V;//density;
					p_half_p[A_Z]+=F[A_Z][i]*V;//density;
				}
				//位置座標の計算
				n_rx[i]=PART[i].r[A_X]+Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;
				n_ry[i]=PART[i].r[A_Y]+Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;
				n_rz[i]=PART[i].r[A_Z]+Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
			}
			else
			{
				n_rx[i]=PART[i].q0[A_X];
				n_ry[i]=PART[i].q0[A_Y];
				n_rz[i]=PART[i].q0[A_Z];
			}
		}
		else*/
		{
			//half_pの計算
			double p_half_p[DIMENSION]={0,0,0};

			int Ni=HYPER[i].N;
			for(int j=0;j<Ni;j++)
			{	
				int k=HYPER[i].NEI[j];
				p_half_p[A_X]+=(HYPER[k].stress[0][0]-lambda[k])*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[0][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[0][2]*HYPER1[k*h_num+i].DgDq[2];
				p_half_p[A_Y]+=HYPER[k].stress[1][0]*HYPER1[k*h_num+i].DgDq[0]+(HYPER[k].stress[1][1]-lambda[k])*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[1][2]*HYPER1[k*h_num+i].DgDq[2];
				p_half_p[A_Z]+=HYPER[k].stress[2][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[2][1]*HYPER1[k*h_num+i].DgDq[1]+(HYPER[k].stress[2][2]-lambda[k])*HYPER1[k*h_num+i].DgDq[2];
			}//jに関するfor文の終わり
			p_half_p[A_X]+=(HYPER[i].stress[0][0]-lambda[i])*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
			p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-lambda[i])*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
			p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-lambda[i])*HYPER1[i*h_num+i].DgDq[2];
		
			//重力の影響
			if(flag_G==ON)	p_half_p[A_Z]-=9.8*mi;
			//粘性項の影響
			if(flag_vis==ON)
			{
				p_half_p[A_X]+=density*HYPER[i].vis_force[A_X];
				p_half_p[A_Y]+=density*HYPER[i].vis_force[A_Y];
				p_half_p[A_Z]+=density*HYPER[i].vis_force[A_Z];
			}
			//磁場の考慮
			if(flag_FEM==ON)
			{
				p_half_p[A_X]+=F[A_X][i]*V;//density;
				p_half_p[A_Y]+=F[A_Y][i]*V;//density;
				p_half_p[A_Z]+=F[A_Z][i]*V;//density;
			}
			//位置座標の計算
			n_rx[i]=PART[i].r[A_X]+Dt*(HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X])/mi;
			n_ry[i]=PART[i].r[A_Y]+Dt*(HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y])/mi;
			n_rz[i]=PART[i].r[A_Z]+Dt*(HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z])/mi;
		}
	}


	////DgDqとfxの更新
	for(int i=0;i<h_num;i++)
	{
		//Fiの計算
		int Ni=HYPER[i].N;
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};	

		for(int in=0;in<Ni;in++)
		{
			int inn=HYPER[i].NEI[in];
			double w=HYPER1[i*h_num+inn].wiin;

			fi[0][0]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[0][1]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[0][2]+=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[A_Z];
			fi[1][0]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[1][1]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[1][2]+=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[A_Z];
			fi[2][0]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_X];
			fi[2][1]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_Y];
			fi[2][2]+=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[A_Z];
		}

		p_Fi[0][0]=fi[0][0]*HYPER[i].inverse_Ai[0][0]+fi[0][1]*HYPER[i].inverse_Ai[1][0]+fi[0][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[0][1]=fi[0][0]*HYPER[i].inverse_Ai[0][1]+fi[0][1]*HYPER[i].inverse_Ai[1][1]+fi[0][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[0][2]=fi[0][0]*HYPER[i].inverse_Ai[0][2]+fi[0][1]*HYPER[i].inverse_Ai[1][2]+fi[0][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[1][0]=fi[1][0]*HYPER[i].inverse_Ai[0][0]+fi[1][1]*HYPER[i].inverse_Ai[1][0]+fi[1][2]*HYPER[i].inverse_Ai[2][0];	
		p_Fi[1][1]=fi[1][0]*HYPER[i].inverse_Ai[0][1]+fi[1][1]*HYPER[i].inverse_Ai[1][1]+fi[1][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[1][2]=fi[1][0]*HYPER[i].inverse_Ai[0][2]+fi[1][1]*HYPER[i].inverse_Ai[1][2]+fi[1][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[2][0]=fi[2][0]*HYPER[i].inverse_Ai[0][0]+fi[2][1]*HYPER[i].inverse_Ai[1][0]+fi[2][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[2][1]=fi[2][0]*HYPER[i].inverse_Ai[0][1]+fi[2][1]*HYPER[i].inverse_Ai[1][1]+fi[2][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[2][2]=fi[2][0]*HYPER[i].inverse_Ai[0][2]+fi[2][1]*HYPER[i].inverse_Ai[1][2]+fi[2][2]*HYPER[i].inverse_Ai[2][2];

		//Jの計算
		double J=calc_det3(p_Fi);

		//fxの計算
		fx[i]=V*(1-J);//1-J;//

		//t_inverse_Fiの計算
		inverse(p_Fi,DIMENSION);

		//DgDqの計算
		for(int j=0;j<Ni;j++)
		{
			int k=HYPER[i].NEI[j];
			n_DgDq_x[i][k]=J*(p_Fi[0][0]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_y[i][k]=J*(p_Fi[0][1]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[k*h_num+i].n0ij[2]);
			n_DgDq_z[i][k]=J*(p_Fi[0][2]*HYPER1[k*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[k*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[k*h_num+i].n0ij[2]);
		}
		n_DgDq_x[i][i]=J*(p_Fi[0][0]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][0]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][0]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_y[i][i]=J*(p_Fi[0][1]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][1]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][1]*HYPER1[i*h_num+i].n0ij[2]);
		n_DgDq_z[i][i]=J*(p_Fi[0][2]*HYPER1[i*h_num+i].n0ij[0]+p_Fi[1][2]*HYPER1[i*h_num+i].n0ij[1]+p_Fi[2][2]*HYPER1[i*h_num+i].n0ij[2]);
	}

	////DfDxの更新
/*	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			double DFDlambda=0;
			for(int k=0;k<h_num;k++)	DFDlambda+=n_DgDq_x[i][k]*HYPER1[j*h_num+k].DgDq[A_X]+n_DgDq_y[i][k]*HYPER1[j*h_num+k].DgDq[A_Y]+n_DgDq_z[i][k]*HYPER1[j*h_num+k].DgDq[A_Z];
			if(DFDlambda!=0)	DfDx[i*h_num+j]=-Dt*Dt*0.5/mi*DFDlambda;//-DFDlambda;//
		}
	}//*/
	
	
	for(int k=0;k<h_num;k++)
	{
		double DFDlambda=0;
		int Nk=HYPER[k].N;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				DfDx[in*h_num+jn]-=Dt*Dt*0.5/mi*(n_DgDq_x[in][k]*HYPER1[jn*h_num+k].DgDq[A_X]+n_DgDq_y[in][k]*HYPER1[jn*h_num+k].DgDq[A_Y]+n_DgDq_z[in][k]*HYPER1[jn*h_num+k].DgDq[A_Z]);
			}
			DfDx[k*h_num+in]-=Dt*Dt*0.5/mi*(n_DgDq_x[k][k]*HYPER1[in*h_num+k].DgDq[A_X]+n_DgDq_y[k][k]*HYPER1[in*h_num+k].DgDq[A_Y]+n_DgDq_z[k][k]*HYPER1[in*h_num+k].DgDq[A_Z]);
			DfDx[in*h_num+k]-=Dt*Dt*0.5/mi*(n_DgDq_x[in][k]*HYPER1[k*h_num+k].DgDq[A_X]+n_DgDq_y[in][k]*HYPER1[k*h_num+k].DgDq[A_Y]+n_DgDq_z[in][k]*HYPER1[k*h_num+k].DgDq[A_Z]);
		}
		DfDx[k*h_num+k]-=Dt*Dt*0.5/mi*(n_DgDq_x[k][k]*HYPER1[k*h_num+k].DgDq[A_X]+n_DgDq_y[k][k]*HYPER1[k*h_num+k].DgDq[A_Y]+n_DgDq_z[k][k]*HYPER1[k*h_num+k].DgDq[A_Z]);		
	}//*/

	////出力
//	if(count%200==0 && count>CON.get_nr()/2)
//	if(t==1||t%CON.get_interval()==0)	if(count%200==0||count==1)	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);
//	if(t==1||t%(10*CON.get_interval())==0)	if(count%200==0||count==1)	output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);

	ofstream t_loge("time_log_newton.dat", ios::app);
	clock_t t4=clock();
	t_loge<<"step="<<t<<", count="<<count<<", time="<<1000*(t4-t3)/CLOCKS_PER_SEC<<"[e-3sec]"<<endl;
	t_loge.close();


	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;

	for(int i=0;i<h_num;i++)
	{
		delete[]	n_DgDq_x[i];
		delete[]	n_DgDq_y[i];
		delete[]	n_DgDq_z[i];
	}
	delete[]	n_DgDq_x;
	delete[]	n_DgDq_y;
	delete[]	n_DgDq_z;

	delete[]	n_rx;
	delete[]	n_ry;
	delete[]	n_rz;
}

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,bool repetation,double **F)
{
//	if(repetation==0)	cout<<"仮の運動量＆位置座標計算";
//	else	cout<<"運動量計算";

	int h_num=HYPER.size();
	double Dt=CON.get_dt();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	int flag_vis=CON.get_flag_vis();
	bool flag_FEM=CON.get_FEM_flag();
	int flag_G=CON.get_flag_G();
	double density=CON.get_hyper_density();
	int model_num=CON.get_model_number();

	for(int i=0;i<h_num;i++)
	{
		double p_half_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{		
			int k=HYPER[i].NEI[j];
			p_half_p[A_X]+=(HYPER[k].stress[0][0]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[0][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[0][2]*HYPER1[k*h_num+i].DgDq[2];
			p_half_p[A_Y]+=HYPER[k].stress[1][0]*HYPER1[k*h_num+i].DgDq[0]+(HYPER[k].stress[1][1]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[1][2]*HYPER1[k*h_num+i].DgDq[2];
			p_half_p[A_Z]+=HYPER[k].stress[2][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[2][1]*HYPER1[k*h_num+i].DgDq[1]+(HYPER[k].stress[2][2]-HYPER[k].lambda)*HYPER1[k*h_num+i].DgDq[2];
		}//jに関するfor文の終わり	
		p_half_p[A_X]+=(HYPER[i].stress[0][0]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[0][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[0][2]*HYPER1[i*h_num+i].DgDq[2];
		p_half_p[A_Y]+=HYPER[i].stress[1][0]*HYPER1[i*h_num+i].DgDq[0]+(HYPER[i].stress[1][1]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[1][2]*HYPER1[i*h_num+i].DgDq[2];
		p_half_p[A_Z]+=HYPER[i].stress[2][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[2][1]*HYPER1[i*h_num+i].DgDq[1]+(HYPER[i].stress[2][2]-HYPER[i].lambda)*HYPER1[i*h_num+i].DgDq[2];

		//重力項
		if(flag_G==ON)	p_half_p[A_Z]-=9.8*mi;
		//粘性項
		if(flag_vis==ON)
		{
			p_half_p[A_X]+=density*HYPER[i].vis_force[A_X];
			p_half_p[A_Y]+=density*HYPER[i].vis_force[A_Y];
			p_half_p[A_Z]+=density*HYPER[i].vis_force[A_Z];
		}
		//磁力項
		if(flag_FEM==ON)
		{
			p_half_p[A_X]+=F[A_X][i]*V;//density;
			p_half_p[A_Y]+=F[A_Y][i]*V;//density;
			p_half_p[A_Z]+=F[A_Z][i]*V;//density;
		/*	if(i==0){
			cout<<" 1="<<mi;
			cout<<" 2="<<F[A_X][i]*mi;
			cout<<" 3="<<F[A_X][i]/density<<endl;
			cout<<" 4="<<HYPER[i].vis_force[A_X];
			cout<<" 5="<<9.8*mi<<endl;}
		}*/
		}
		if(repetation==0)
		{

			//half_pの更新
			HYPER[i].half_p[A_X]=HYPER[i].p[A_X]+Dt*0.5*p_half_p[A_X];
			HYPER[i].half_p[A_Y]=HYPER[i].p[A_Y]+Dt*0.5*p_half_p[A_Y];
			HYPER[i].half_p[A_Z]=HYPER[i].p[A_Z]+Dt*0.5*p_half_p[A_Z];//
			//位置座標の更新
			/*if(model_num==30||model_num==23)
			{
				if(PART[i].q0[A_Z]!=0)
				{
					PART[i].r[A_X]+=Dt*HYPER[i].half_p[A_X]/mi;
					PART[i].r[A_Y]+=Dt*HYPER[i].half_p[A_Y]/mi;
					PART[i].r[A_Z]+=Dt*HYPER[i].half_p[A_Z]/mi;					
				}
				else
				{
					PART[i].r[A_X]=PART[i].q0[A_X];
					PART[i].r[A_Y]=PART[i].q0[A_Y];
					PART[i].r[A_Z]=PART[i].q0[A_Z];
				}
			}
			else*/
			{
				PART[i].r[A_X]+=Dt*HYPER[i].half_p[A_X]/mi;
				PART[i].r[A_Y]+=Dt*HYPER[i].half_p[A_Y]/mi;
				PART[i].r[A_Z]+=Dt*HYPER[i].half_p[A_Z]/mi;
			}
		}
		else
		{
			//運動量の更新
			HYPER[i].p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*p_half_p[A_X];
			HYPER[i].p[A_Y]=HYPER[i].half_p[A_Y]+Dt*0.5*p_half_p[A_Y];
			HYPER[i].p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*p_half_p[A_Z];////
			//速度の更新
			PART[i].u[A_X]=HYPER[i].half_p[A_X]/mi;
			PART[i].u[A_Y]=HYPER[i].half_p[A_Y]/mi;
			PART[i].u[A_Z]=HYPER[i].half_p[A_Z]/mi;
			//角運動量の更新
			HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
			HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
			HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
		}
	}//iに関するfor文の終わり
	
//	cout<<"----------OK"<<endl;
}

void calc_F(vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1)
{
//	cout<<"Fi計算";
	////Fiの更新
	int h_num=HYPER.size();

	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];

	for(int i=0;i<h_num;i++)
	{
		double fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
		//Fiの計算

		int Ni=HYPER[i].N;	

		for(int in=0;in<Ni;in++)
		{
			int inn=HYPER[i].NEI[in];
			double w=HYPER1[i*h_num+inn].wiin;
			double a[DIMENSION]={HYPER1[i*h_num+inn].aiin[A_X],	HYPER1[i*h_num+inn].aiin[A_Y],	HYPER1[i*h_num+inn].aiin[A_Z]};
			
			fi[0][0]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_X];	fi[0][1]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_Y];	fi[0][2]+=w*(PART[inn].r[A_X]-PART[i].r[A_X])*a[A_Z];
			fi[1][0]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_X];	fi[1][1]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_Y];	fi[1][2]+=w*(PART[inn].r[A_Y]-PART[i].r[A_Y])*a[A_Z];
			fi[2][0]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_X];	fi[2][1]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_Y];	fi[2][2]+=w*(PART[inn].r[A_Z]-PART[i].r[A_Z])*a[A_Z];
		}

		p_Fi[0][0]=fi[0][0]*HYPER[i].inverse_Ai[0][0]+fi[0][1]*HYPER[i].inverse_Ai[1][0]+fi[0][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[0][1]=fi[0][0]*HYPER[i].inverse_Ai[0][1]+fi[0][1]*HYPER[i].inverse_Ai[1][1]+fi[0][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[0][2]=fi[0][0]*HYPER[i].inverse_Ai[0][2]+fi[0][1]*HYPER[i].inverse_Ai[1][2]+fi[0][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[1][0]=fi[1][0]*HYPER[i].inverse_Ai[0][0]+fi[1][1]*HYPER[i].inverse_Ai[1][0]+fi[1][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[1][1]=fi[1][0]*HYPER[i].inverse_Ai[0][1]+fi[1][1]*HYPER[i].inverse_Ai[1][1]+fi[1][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[1][2]=fi[1][0]*HYPER[i].inverse_Ai[0][2]+fi[1][1]*HYPER[i].inverse_Ai[1][2]+fi[1][2]*HYPER[i].inverse_Ai[2][2];
		p_Fi[2][0]=fi[2][0]*HYPER[i].inverse_Ai[0][0]+fi[2][1]*HYPER[i].inverse_Ai[1][0]+fi[2][2]*HYPER[i].inverse_Ai[2][0];
		p_Fi[2][1]=fi[2][0]*HYPER[i].inverse_Ai[0][1]+fi[2][1]*HYPER[i].inverse_Ai[1][1]+fi[2][2]*HYPER[i].inverse_Ai[2][1];
		p_Fi[2][2]=fi[2][0]*HYPER[i].inverse_Ai[0][2]+fi[2][1]*HYPER[i].inverse_Ai[1][2]+fi[2][2]*HYPER[i].inverse_Ai[2][2];
		
		HYPER[i].Fi[0][0]=p_Fi[0][0];	HYPER[i].Fi[0][1]=p_Fi[0][1];	HYPER[i].Fi[0][2]=p_Fi[0][2];	
		HYPER[i].Fi[1][0]=p_Fi[1][0];	HYPER[i].Fi[1][1]=p_Fi[1][1];	HYPER[i].Fi[1][2]=p_Fi[1][2];	
		HYPER[i].Fi[2][0]=p_Fi[2][0];	HYPER[i].Fi[2][1]=p_Fi[2][1];	HYPER[i].Fi[2][2]=p_Fi[2][2];	
		//Jの計算
		double J=calc_det3(p_Fi);
	//	for(int i=0;i<h_num;i++)	cout<<"J["<<i<<"]="<<J<<endl;
		HYPER[i].J=J;

		//t_inverse_Fiの計算
		inverse(p_Fi,DIMENSION);
		HYPER[i].t_inverse_Fi[0][0]=p_Fi[0][0];	HYPER[i].t_inverse_Fi[0][1]=p_Fi[1][0];	HYPER[i].t_inverse_Fi[0][2]=p_Fi[2][0];
		HYPER[i].t_inverse_Fi[1][0]=p_Fi[0][1];	HYPER[i].t_inverse_Fi[1][1]=p_Fi[1][1];	HYPER[i].t_inverse_Fi[1][2]=p_Fi[2][1];
		HYPER[i].t_inverse_Fi[2][0]=p_Fi[0][2];	HYPER[i].t_inverse_Fi[2][1]=p_Fi[1][2];	HYPER[i].t_inverse_Fi[2][2]=p_Fi[2][2];
	}
	
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;

	//calculation of DgDq
	for(int i=0;i<h_num;i++)
	{
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{			
			int k=HYPER[i].NEI[j];
			HYPER1[k*h_num+i].DgDq[A_X]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Y]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+k].n0ij[2]);
			HYPER1[k*h_num+i].DgDq[A_Z]=HYPER[k].J*(HYPER[k].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+k].n0ij[0]+HYPER[k].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+k].n0ij[1]+HYPER[k].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+k].n0ij[2]);
			//cout<<"i"<<i<<"j"<<k<<"	"<<HYPER1[k*h_num+i].DgDq[A_X]<<","<<HYPER1[k*h_num+i].DgDq[A_Y]<<","<<HYPER1[k*h_num+i].DgDq[A_Z]<<endl;
		}
		HYPER1[i*h_num+i].DgDq[A_X]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_X][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_X][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_X][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Y]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Y][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Y][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Y][2]*HYPER1[i*h_num+i].n0ij[2]);
		HYPER1[i*h_num+i].DgDq[A_Z]=HYPER[i].J*(HYPER[i].t_inverse_Fi[A_Z][0]*HYPER1[i*h_num+i].n0ij[0]+HYPER[i].t_inverse_Fi[A_Z][1]*HYPER1[i*h_num+i].n0ij[1]+HYPER[i].t_inverse_Fi[A_Z][2]*HYPER1[i*h_num+i].n0ij[2]);
		//cout<<"i"<<i<<"j"<<i<<"	"<<HYPER1[i*h_num+i].DgDq[A_X]<<","<<HYPER1[i*h_num+i].DgDq[A_Y]<<","<<HYPER1[i*h_num+i].DgDq[A_Z]<<endl;
	}
//	cout<<"----------OK"<<endl;

}

void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER)
{
//	cout<<"応力計算";
	int h_num=HYPER.size();

	double **d_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	d_Fi[D]=new double [DIMENSION];

	double c10=CON.get_c10();
	double c01=CON.get_c01();
	
	double b[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double bb[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};

	for(int j=0;j<h_num;j++)
	{	
		double J=HYPER[j].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[j].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[j].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[j].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[j].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[j].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[j].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[j].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[j].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[j].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[j].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[j].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[j].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[j].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[j].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[j].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[j].Fi[2][2];
		}

		b[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[0][1]*d_Fi[0][1]+d_Fi[0][2]*d_Fi[0][2];
		b[0][1]=d_Fi[0][0]*d_Fi[1][0]+d_Fi[0][1]*d_Fi[1][1]+d_Fi[0][2]*d_Fi[1][2];
		b[0][2]=d_Fi[0][0]*d_Fi[2][0]+d_Fi[0][1]*d_Fi[2][1]+d_Fi[0][2]*d_Fi[2][2];
		b[1][0]=d_Fi[1][0]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[0][2];
		b[1][1]=d_Fi[1][0]*d_Fi[1][0]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[1][2]*d_Fi[1][2];
		b[1][2]=d_Fi[1][0]*d_Fi[2][0]+d_Fi[1][1]*d_Fi[2][1]+d_Fi[1][2]*d_Fi[2][2];
		b[2][0]=d_Fi[2][0]*d_Fi[0][0]+d_Fi[2][1]*d_Fi[0][1]+d_Fi[2][2]*d_Fi[0][2];
		b[2][1]=d_Fi[2][0]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[1][2];
		b[2][2]=d_Fi[2][0]*d_Fi[2][0]+d_Fi[2][1]*d_Fi[2][1]+d_Fi[2][2]*d_Fi[2][2];

		bb[0][0]=b[0][0]*b[0][0]+b[0][1]*b[1][0]+b[0][2]*b[2][0];
		bb[0][1]=b[0][0]*b[0][1]+b[0][1]*b[1][1]+b[0][2]*b[2][1];
		bb[0][2]=b[0][0]*b[0][2]+b[0][1]*b[1][2]+b[0][2]*b[2][2];
		bb[1][0]=b[1][0]*b[0][0]+b[1][1]*b[1][0]+b[1][2]*b[2][0];
		bb[1][1]=b[1][0]*b[0][1]+b[1][1]*b[1][1]+b[1][2]*b[2][1];
		bb[1][2]=b[1][0]*b[0][2]+b[1][1]*b[1][2]+b[1][2]*b[2][2];
		bb[2][0]=b[2][0]*b[0][0]+b[2][1]*b[1][0]+b[2][2]*b[2][0];
		bb[2][1]=b[2][0]*b[0][1]+b[2][1]*b[1][1]+b[2][2]*b[2][1];
		bb[2][2]=b[2][0]*b[0][2]+b[2][1]*b[1][2]+b[2][2]*b[2][2];

		double trace_b=b[0][0]+b[1][1]+b[2][2];
		double trace_bb=bb[0][0]+bb[1][1]+bb[2][2];
	
		HYPER[j].stress[0][0]=2/J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
		HYPER[j].stress[0][1]=2/J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
		HYPER[j].stress[0][2]=2/J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
		HYPER[j].stress[1][0]=2/J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
		HYPER[j].stress[1][1]=2/J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
		HYPER[j].stress[1][2]=2/J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
		HYPER[j].stress[2][0]=2/J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
		HYPER[j].stress[2][1]=2/J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
		HYPER[j].stress[2][2]=2/J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));
		
		/*cout<<"stress="<<endl;
		cout<<HYPER[j].stress[0][0]<<"	"<<HYPER[j].stress[0][1]<<"	"<<HYPER[j].stress[0][2]<<endl;
		cout<<HYPER[j].stress[1][0]<<"	"<<HYPER[j].stress[1][1]<<"	"<<HYPER[j].stress[1][2]<<endl;
		cout<<HYPER[j].stress[2][0]<<"	"<<HYPER[j].stress[2][1]<<"	"<<HYPER[j].stress[2][2]<<endl;*/
	}

	for(int D=0;D>DIMENSION;D++)	delete[]	d_Fi[D];
	delete[]	d_Fi;

//	cout<<"----------OK"<<endl;
}

void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,double **F)
{
//	cout<<"運動量微分値計算";

	int h_num=HYPER.size();
	int flag_vis=CON.get_flag_vis();
	int flag_G=CON.get_flag_G();
	int flag_FEM=CON.get_FEM_flag();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double density=CON.get_hyper_density();
	for(int i=0;i<h_num;i++)
	{
		double p_differential_p[DIMENSION]={0,0,0};
		int Ni=HYPER[i].N;
		for(int j=0;j<Ni;j++)
		{				
			int k=HYPER[i].NEI[j];
			p_differential_p[A_X]+=HYPER[k].stress[A_X][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_X][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_X][2]*HYPER1[k*h_num+i].DgDq[2];
			p_differential_p[A_Y]+=HYPER[k].stress[A_Y][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_Y][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_Y][2]*HYPER1[k*h_num+i].DgDq[2];	
			p_differential_p[A_Z]+=HYPER[k].stress[A_Z][0]*HYPER1[k*h_num+i].DgDq[0]+HYPER[k].stress[A_Z][1]*HYPER1[k*h_num+i].DgDq[1]+HYPER[k].stress[A_Z][2]*HYPER1[k*h_num+i].DgDq[2];	
		}
		HYPER[i].differential_p[A_X]=HYPER[i].half_p[A_X]+Dt*0.5*(p_differential_p[A_X]+HYPER[i].stress[A_X][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_X][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_X][2]*HYPER1[i*h_num+i].DgDq[2]);
		HYPER[i].differential_p[A_Y]=HYPER[i].half_p[A_X]+Dt*0.5*(p_differential_p[A_X]+HYPER[i].stress[A_Y][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Y][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Y][2]*HYPER1[i*h_num+i].DgDq[2]);	
		HYPER[i].differential_p[A_Z]=HYPER[i].half_p[A_Z]+Dt*0.5*(p_differential_p[A_Z]+HYPER[i].stress[A_Z][0]*HYPER1[i*h_num+i].DgDq[0]+HYPER[i].stress[A_Z][1]*HYPER1[i*h_num+i].DgDq[1]+HYPER[i].stress[A_Z][2]*HYPER1[i*h_num+i].DgDq[2]);////p_differential_p[A_Z];////
		
		//重力影響
		if(flag_G==ON)	HYPER[i].differential_p[A_Z]-=Dt*0.5*9.8*mi;
		if(flag_vis==ON)
		{
			HYPER[i].differential_p[A_X]+=Dt*0.5*HYPER[i].vis_force[A_X]*density;
			HYPER[i].differential_p[A_Y]+=Dt*0.5*HYPER[i].vis_force[A_Y]*density;
			HYPER[i].differential_p[A_Z]+=Dt*0.5*HYPER[i].vis_force[A_Z]*density;
		}
		if(flag_FEM==ON)
		{
			HYPER[i].differential_p[A_X]+=Dt*0.5*F[A_X][i]*V;//density;
			HYPER[i].differential_p[A_Y]+=Dt*0.5*F[A_Y][i]*V;//density;
			HYPER[i].differential_p[A_Z]+=Dt*0.5*F[A_Z][i]*V;//density;
		}
	}
//	cout<<"----------OK"<<endl;

}

void renew_lambda(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int t)
{
	clock_t t3=clock();

//	cout<<"Lambda計算";

	int h_num=HYPER.size();
	double le=CON.get_distancebp();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mk=V*CON.get_hyper_density();

	double *N_Left=new double[h_num*h_num];
	double *N_Right=new double[h_num];
	for(int i=0;i<h_num;i++)
	{
		N_Right[i]=0;
		for(int j=0;j<h_num;j++)	N_Left[j*h_num+i]=0;
	}

	
/*	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			double N_left=0;
			for(int k=0;k<h_num;k++)
			{
				N_left+=HYPER1[i*h_num+k].DgDq[0]*HYPER1[j*h_num+k].DgDq[0]+HYPER1[i*h_num+k].DgDq[1]*HYPER1[j*h_num+k].DgDq[1]+HYPER1[i*h_num+k].DgDq[2]*HYPER1[j*h_num+k].DgDq[2];
				//cout<<"i"<<i<<"j"<<j<<"k"<<k<<"	"<<HYPER1[i*h_num+k].DgDq[0]<<"	"<<HYPER1[j*h_num+k].DgDq[0]<<"	"<<HYPER1[i*h_num+k].DgDq[1]<<"	"<<HYPER1[j*h_num+k].DgDq[1]<<"	"<<HYPER1[i*h_num+k].DgDq[2]<<"	"<<HYPER1[j*h_num+k].DgDq[2]<<endl;
				//cout<<"i"<<i<<"j"<<j<<"k"<<k<<"	"<<HYPER1[i*h_num+k].DgDq[0]*HYPER1[j*h_num+k].DgDq[0]+HYPER1[i*h_num+k].DgDq[1]*HYPER1[j*h_num+k].DgDq[1]+HYPER1[i*h_num+k].DgDq[2]*HYPER1[j*h_num+k].DgDq[2]<<endl;
			}
			N_Left[i*h_num+j]=Dt*0.5*N_left;//Dt/2/mk*N_left;
		}//jに関するfor文の終わり
	}//iに関するfor文の終わり
	for(int i=0;i<h_num;i++)
	{
		double N_right=0;
		for(int j=0;j<h_num;j++)
		{
			N_right+=HYPER[j].differential_p[0]*HYPER1[i*h_num+j].DgDq[0]+HYPER[j].differential_p[1]*HYPER1[i*h_num+j].DgDq[1]+HYPER[j].differential_p[2]*HYPER1[i*h_num+j].DgDq[2];
			//cout<<"i"<<i<<"j"<<j<<"	"<<HYPER[j].differential_p[0]*HYPER1[i*h_num+j].DgDq[0]+HYPER[j].differential_p[1]*HYPER1[i*h_num+j].DgDq[1]+HYPER[j].differential_p[2]*HYPER1[i*h_num+j].DgDq[2]<<endl;
		}
		N_Right[i]=N_right;//1/mk*N_right;	
	}//*/
	
	for(int k=0;k<h_num;k++)
	{
		int Nk=HYPER[k].N;
		double N_right=0;
		for(int i=0;i<Nk;i++)
		{
			int in=HYPER[k].NEI[i];
			for(int j=0;j<Nk;j++)
			{
				int jn=HYPER[k].NEI[j];
				N_Left[in*h_num+jn]+=Dt*0.5*(HYPER1[in*h_num+k].DgDq[0]*HYPER1[jn*h_num+k].DgDq[0]+HYPER1[in*h_num+k].DgDq[1]*HYPER1[jn*h_num+k].DgDq[1]+HYPER1[in*h_num+k].DgDq[2]*HYPER1[jn*h_num+k].DgDq[2]);
			}
			N_Left[in*h_num+k]+=Dt*0.5*(HYPER1[in*h_num+k].DgDq[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER1[in*h_num+k].DgDq[1]*HYPER1[k*h_num+k].DgDq[1]+HYPER1[in*h_num+k].DgDq[2]*HYPER1[k*h_num+k].DgDq[2]);
			N_Left[k*h_num+in]+=Dt*0.5*(HYPER1[k*h_num+k].DgDq[0]*HYPER1[in*h_num+k].DgDq[0]+HYPER1[k*h_num+k].DgDq[1]*HYPER1[in*h_num+k].DgDq[1]+HYPER1[k*h_num+k].DgDq[2]*HYPER1[in*h_num+k].DgDq[2]);
			N_right+=HYPER[in].differential_p[0]*HYPER1[k*h_num+in].DgDq[0]+HYPER[in].differential_p[1]*HYPER1[k*h_num+in].DgDq[1]+HYPER[in].differential_p[2]*HYPER1[k*h_num+in].DgDq[2];
		}//jに関するfor文の終わり
		N_Right[k]=N_right+HYPER[k].differential_p[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER[k].differential_p[1]*HYPER1[k*h_num+k].DgDq[1]+HYPER[k].differential_p[2]*HYPER1[k*h_num+k].DgDq[2];//1/mk*N_right;	
		N_Left[k*h_num+k]+=Dt*0.5*(HYPER1[k*h_num+k].DgDq[0]*HYPER1[k*h_num+k].DgDq[0]+HYPER1[k*h_num+k].DgDq[1]*HYPER1[k*h_num+k].DgDq[1]+HYPER1[k*h_num+k].DgDq[2]*HYPER1[k*h_num+k].DgDq[2]);
	}//iに関するfor文の終わり*/

/*
	for(int i=0;i<h_num;i++)	cout<<"N_Right["<<i<<"]="<<N_Right[i]<<endl;
	cout<<endl;
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)	cout<<"N_Left["<<i<<"]["<<j<<"]="<<N_Left[i*h_num+j]<<endl;
		cout<<endl;
	}*/

	//lambdaを求める

	/*ofstream fl("N_Left.csv");
	ofstream fr("N_Reft.csv");

	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			fl<<N_Left[i*h_num+j]<<",";
		}
		fl<<endl;
		fr<<N_Right[i]<<endl;
	}
	fl.close();
	fr.close();//*/
		
	/*gauss(N_Left,N_Right,h_num);
	for(int i=0;i<h_num;i++)	HYPER[i].lambda=N_Right[i];//*/

	
	int *b_ind=new int [h_num*h_num];
	double *b_val=new double [h_num*h_num];
	int *ptr=new int [h_num+1];

	int all_ind_num=0;
	for(int i=0;i<h_num;i++)
	{
		ptr[i+1]=0;
		for(int j=0;j<h_num;j++)
		{
			b_ind[i*h_num+j]=0;
			b_val[i*h_num+j]=0;
		}
	}

	ptr[0]=0;
	for(int i=0;i<h_num;i++)
	{
		int ind_num=0;
		for(int j=0;j<h_num;j++)
		{
			if(N_Left[i*h_num+j]!=0)
			{
				b_ind[all_ind_num+ind_num]=j;	
				b_val[all_ind_num+ind_num]=N_Left[i*h_num+j];
				ind_num++;
				ptr[i+1]=ind_num+ptr[i];
			}			
		}
		//cout<<"ptr"<<i+1<<" "<<ptr[i+1]<<endl;
		all_ind_num+=ind_num;			
	}


	int *ind=new int [all_ind_num];
	double *val=new double [all_ind_num];
	for(int i=0;i<all_ind_num;i++)
	{
		ind[i]=b_ind[i];
		val[i]=b_val[i];
	}

	double *X=new double [h_num];
	//BiCGStab2_method(&CON,val,ind,ptr,h_num,N_Right,all_ind_num,X);
	//iccg2(&CON,val,ind,ptr,h_num,N_Right,all_ind_num,X);
	CG3D(&CON,val,ind,ptr,h_num,N_Right,all_ind_num,X);

	delete[] ind;
	delete[] val;
	delete[] ptr;

	for(int i=0;i<h_num;i++)	HYPER[i].lambda=X[i];

	delete[] X;//*/

//	for(int i=0;i<h_num;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;

	ofstream t_loge("time_log_gauss.dat", ios::app);
	clock_t t4=clock();
	t_loge<<"step="<<t<<", time="<<1000*(t4-t3)/CLOCKS_PER_SEC<<"[e-3sec]"<<endl;
	t_loge.close();

	delete [] N_Left;
	delete [] N_Right;

//	cout<<"----------OK"<<endl;
}

//detを求める関数　※自己流のため自信な
double calc_det(double **M,int N)
{
	double det=0;
	double det_plus;
	double det_minus;
	int	a=0;

	for(int j=0;j<N;j++)
	{
		det_plus=1;
		for(int i=0;i<N;i++)	det_plus*=M[i%N][(j+i)%N];
		det+=det_plus;

		det_minus=1;
		for(int i=0;i<N;i++)
		{
			if(j-i<0)
			{
				a=j-i;
				a+=N;
				det_minus*=M[i%N][a%N];
			}
			else
			{
				det_minus*=M[i%N][(j-i)%N];
			}
		}
		det-=det_minus;
	}
	return det;
}

double calc_det3(double **M)
{
	double det=0;
	det+=M[0][0]*M[1][1]*M[2][2]+M[0][1]*M[1][2]*M[2][0]+M[0][2]*M[1][0]*M[2][1];
	det-=M[0][2]*M[1][1]*M[2][0]+M[0][1]*M[1][0]*M[2][2]+M[0][0]*M[1][2]*M[2][1];
	return det;
}


//逆行列を求める関数 function.hの圧力計算用関数に類似した名前の関数があったのですみわけ
void calc_inverse_matrix_for_NR(int N, double *a)
{
	//N=未知数
	double buf=0;

	double *inv_a=new  double[N*N];					//逆行列格納
	
	//単位行列を作る
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			if(i==j) inv_a[j*N+i]=1;
			else		inv_a[j*N+i]=0;
		}
	}

	//掃き出し法
	for(int i=0;i<N;i++)
	{
		if(a[i*N+i]<DBL_EPSILON) cout<<"a[i][i]=0"<<endl;
		buf=1/a[i*N+i];
		for(int j=0;j<N;j++)
		{
			a[j*N+i]*=buf;
			inv_a[j*N+i]*=buf;
		}
		for(int j=0;j<N;j++)
		{
			if(i!=j)
			{
				buf=a[i*N+j];//a[j][i];
				for(int k=0;k<N;k++)
				{
					a[k*N+j]-=a[k*N+i]*buf;//a[j][k]-=a[i][k]*buf;
					inv_a[k*N+j]-=inv_a[k*N+i]*buf;
				}
			}
		}
	}

	for(int i=0;i<N;i++)	for(int j=0;j<N;j++) a[i*N+j]=inv_a[i*N+j];

	delete [] inv_a;
}



void inverse(double **a,int N)
{
	double d=0;
	double *col=new double [N];
	double **y=new double *[N];
	for(int i=0;i<N;i++)	y[i]=new double [N];
	int *index=new int [N];

	//int p_revision=0,m_revision=0;
	int fig=0,x=0;
	int	flag=OFF;

	for(int D=0;D<N;D++)
	{
		col[D]=0;
		index[D]=0;
		for(int D2=0;D2<N;D2++)	y[D][D2]=0;
	}

	for(int D=0;D<N;D++)	for(int D2=0;D2<N;D2++)	if(fabs(a[D][D2])<DBL_EPSILON)	a[D][D2]=0;

	for(int D=0;D<N;D++)
	{
		for(int D2=0;D2<N;D2++)
		{
			if(fabs(a[D][D2])>pow(10.0,5.0))
			{
				flag=ON;
				x=log10(fabs(a[D][D2]));
				if(x>fig)	fig=x;
			}
		}
	}	
	
	if(flag==ON)	for(int D=0;D<N;D++) for(int D2=0;D2<N;D2++) a[D][D2]/=pow(10.0,fig);

	ludcmp(a,N,index,&d);
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<N;i++)	col[i]=0.0;
		col[j]=1.0;
		lubksb(a,N,index,col);
		for(int i=0;i<N;i++)	y[i][j]=col[i];
	}

	for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	a[i][j]=y[i][j];
	if(flag==ON)	for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	a[i][j]/=pow(10.0,fig);

	for(int i=0;i<N;i++)	delete[]	y[i];
	delete[]	y;
	delete[]	col;
	delete[]	index;
}




void ludcmp(double **a,int N,int *index,double *d)
{
	int imax=0;
	int count=0;
	double max=0,temp=0,sum=0,dum=0;
	*d=1.0;
	vector<float> buf(N);

	for(int i=0;i<N;i++)	buf[i]=0;

	for(int i=0;i<N;i++)
	{
		max=0.0;
		for(int j=0;j<N;j++)
		{
			temp=fabs(a[i][j]);
			if(temp>max)	max=temp;
		}
		if(max==0.0&&count==0)
		{
			count++;
			cout<<"特異行列である"<<endl;
		}
		buf[i]=1.0/max;
	}

	for(int j=0;j<N;j++)
	{
		for(int i=0;i<j;i++)
		{
			sum=a[i][j];
			for(int k=0;k<i;k++)	sum-=a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		max=0.0;
		for(int i=j;i<N;i++)
		{
			sum=a[i][j];
			for(int k=0;k<j;k++)	sum-=a[i][k]*a[k][j];
			a[i][j]=sum;
			dum=buf[i]*fabs(sum);
			if(dum>=max)
			{
				max=dum;
				imax=i;
			}
		}

		if(j!=imax)
		{
			for(int k=0;k<N;k++)
			{
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d=-(*d);
			buf[imax]=buf[j];
		}
		index[j]=imax;
		if(fabs(a[j][j])<DBL_EPSILON)
		{
			cout<<"error:a["<<j<<"]["<<j<<"]=0"<<endl;
			a[j][j]=0;
		}
		if(j!=N-1)
		{
			dum=1.0/(a[j][j]);
			for(int i=j+1;i<N;i++)a[i][j]*=dum;
		}
	}
	buf.clear();
}

void lubksb(double **a,int N,int *index,double b[])
{
	int ip=0,ii=0,count=0;
	double sum=0;

	for(int i=0;i<N;i++)
	{
		ip=index[i];
		sum=b[ip];
		b[ip]=b[i];
		if(count)	for(int j=ii;j<i;j++)	sum-=a[i][j]*b[j];
		else if(sum)
		{
			count++;
			ii=i;
		}
		b[i]=sum;
	}

	for(int i=N-1;i>=0;i--)
	{
		sum=b[i];
		for(int j=i+1;j<N;j++)	sum-=a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}



//圧力など粒子の持つ情報をコンター図で表示する関数
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,double **F)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	int h_num=HYPER.size();
	double TIME=CON.get_step()*CON.get_dt();
	double le=CON.get_distancebp();
	int STEP=CON.get_step()/CON.get_interval()+1;		//出力する総ステップ数
	int step;
	int n=3;	//出力次元
	double mi=get_volume(&CON)*CON.get_hyper_density();

	if(CON.get_interval()==1)	step=t/CON.get_interval();
	else step=t/CON.get_interval()+1;
	

	if(t==1) 
	{
		ofstream fp("momentum.inp", ios::trunc);			
		fp<<STEP<<endl;//総ステップ数
		fp<<"data_geom"<<endl;
		fp.close();
	}

	//mainファイル書き込み
	ofstream fp("momentum.inp",ios :: app);
	fp<<"step"<<step<<" TIME="<<TIME<<endl;

	//fp<<step<<endl;
	
	//fp<<"data_geom"<<endl;
	//fp<<"step1"<<endl;
	//fp<<"step"<<t/CON->get_interval()+1<<" TIME="<<TIME<<endl;
	
	if(n==3)
	{
		fp<<h_num<<" "<<h_num<<endl;	//節点数と要素数出力
	
		//節点番号とその座標の出力 
		for(int i=0;i<h_num;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	
		//要素番号と要素形状の種類、そして要素を構成する節点番号出力
		for(int i=0;i<h_num;i++)	fp<<i<<"  0 pt "<<i<<endl;

		//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
		fp<<"10 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
		fp<<"10 1 1 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
		//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
		fp<<"p_x,"<<endl;
		fp<<"p_y,"<<endl;
		fp<<"p_z,"<<endl;
		fp<<"lambda,"<<endl;
		fp<<"F_x,"<<endl;
		fp<<"F_y,"<<endl;
		fp<<"F_z,"<<endl;
		fp<<"v_x,"<<endl;
		fp<<"v_y,"<<endl;
		fp<<"v_z,"<<endl;
		//fp<<"P,N/m^2"<<endl;
		//fp<<"value1,??"<<endl;

		//各節点の情報値入力
		for(int i=0;i<h_num;i++)
		{
			fp<<i<<" "<<HYPER[i].p[A_X]<<" "<<HYPER[i].p[A_Y]<<" "<<HYPER[i].p[A_Z]<<" "<<HYPER[i].lambda<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<" "<<HYPER[i].vis_force[A_X]<<" "<<HYPER[i].vis_force[A_Y]<<" "<<HYPER[i].vis_force[A_Z]<<endl;
			//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
		}
	}
	else if(n==2)
	{
		int num=0;
		int *Nxz=new int [h_num];
		fp<<h_num<<" "<<h_num<<endl;	//節点数と要素数出力
	
		//節点番号とその座標の出力 
		for(int i=0;i<h_num;i++)
		{
			if(PART[i].q0[A_Y]<1/2*le&&PART[i].q0[A_Y]>-1/2*le)
			{
				fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
				Nxz[num]=i;
				num++;
			}
		}
		//要素番号と要素形状の種類、そして要素を構成する節点番号出力
		for(int i=0;i<num;i++)	fp<<i<<"  0 pt "<<i<<endl;

		//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
//		fp<<"10 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
//		fp<<"10 1 1 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
		//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
		fp<<"7 0"<<endl;
		fp<<"7 1 1 1 1 1 1 1"<<endl;
		fp<<"v_x,"<<endl;
		fp<<"v_y,"<<endl;
		fp<<"v_z,"<<endl;
		fp<<"lambda,"<<endl;
		fp<<"F_x,"<<endl;
		fp<<"F_y,"<<endl;
		fp<<"F_z,"<<endl;
/*		fp<<"v_x,"<<endl;
		fp<<"v_y,"<<endl;
		fp<<"v_z,"<<endl;*/
		//fp<<"P,N/m^2"<<endl;
		//fp<<"value1,??"<<endl;

		//各節点の情報値入力
		for(int i=0;i<num;i++)
		{
			int j=Nxz[i];
			fp<<j<<" "<<1/mi*HYPER[j].p[A_X]<<" "<<1/mi*HYPER[j].p[A_Y]<<" "<<1/mi*HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<endl;
			//fp<<j<<" "<<HYPER[j].p[A_X]<<" "<<HYPER[j].p[A_Y]<<" "<<HYPER[j].p[A_Z]<<" "<<HYPER[j].lambda<<" "<<F[A_X][j]<<" "<<F[A_Y][j]<<" "<<F[A_Z][j]<<" "<<HYPER[j].vjs_force[A_X]<<" "<<HYPER[j].vjs_force[A_Y]<<" "<<HYPER[j].vjs_force[A_Z]<<endl;
			//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
		}
		delete[] Nxz;
	}
	fp.close();
}


void contact_judge_hyper(mpsconfig CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,int t)
{
	cout<<"壁影響計算_勾配モデル";

	//アルゴリズム
	// 0. i周辺の粒子数密度が増加した場合，影響半径内にある粒子を探索し，以下を行う
	// 1. 「接触の可能性がある粒子」（(PART[j].PND>PART[j].PND0)が真？）を調べる
	// 2. 圧力を置換
	// 3. 初期配置の粒子と重複しないように接触の可能性がある粒子との間で力を計算する
	int dim=3;
	double r=CON.get_h_dis();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double Dt=CON.get_dt();
	int p_num=PART.size();
	int h_num=HYPER.size();
	int w_num=p_num-h_num;

	double **w=new double *[h_num];
	double **dis=new double *[h_num];
	double **NEI_w=new double*[h_num];

	for(int i=0;i<h_num;i++)
	{
		w[i]=new double [w_num];
		dis[i]=new double [w_num];
		NEI_w[i]=new double[w_num];
	}

	if(t!=1)	for(int i=h_num;i<w_num+h_num;i++)
	{
		PART[i].r[0]+=PART[i].u[0]*Dt;
		PART[i].r[1]+=PART[i].u[1]*Dt;
		PART[i].r[2]+=PART[i].u[2]*Dt;
	}
		
//		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		for(int i=0;i<h_num;i++)
		{
			int N_w=0;
			double pnd=0;

			for(int j=0;j<w_num;j++)
			{
				int k=j+h_num;
			//初期配置の粒子とは普通に圧力勾配を計算（内力計算）
			//else if(PART[i].contact==true)
				double dis_temp0=0;
				double dis_temp1=0;
				double dis_temp=0;
				for(int D=0;D<DIMENSION;D++)
				{
					dis_temp0=PART[k].r[D]-PART[i].r[D];
					dis_temp1+=dis_temp0*dis_temp0;
				}
				dis_temp=sqrt(dis_temp1);
				if(dis_temp<r)
				{
					dis[i][N_w]=dis_temp;
					w[i][N_w]=kernel(r,dis_temp);
					pnd+=w[i][N_w];
					NEI_w[i][N_w]=k;
					N_w++;				
				}
				//現在位置での周辺粒子数を取得			
			}
			if(N_w>0)
			{
				for(int j=0;j<h_num;j++)
				{
					if(j!=i)
					{					
						double qiin[DIMENSION];
						for(int D=0;D<DIMENSION;D++)	qiin[D]=PART[j].r[D]-PART[i].r[D];
						double dis0=sqrt(qiin[A_X]*qiin[A_X]+qiin[A_Y]*qiin[A_Y]+qiin[A_Z]*qiin[A_Z]);
						double wiin=kernel(r,dis0);
						pnd+=wiin;
					}
				}
				HYPER[i].pnd=pnd;

				double gra_accel_i[DIMENSION];
				for(int D=0;D<DIMENSION;D++)	gra_accel_i[D]=0;

				for(int j=0;j<N_w;j++)
				{
					int k=NEI_w[i][j];
					double rij[DIMENSION];
					double accel_i[DIMENSION];
					for(int D=0;D<DIMENSION;D++)
					{
						rij[D]=0;
						accel_i[D]=0;

						rij[D]=PART[k].r[D]-PART[i].r[D];
						accel_i[D]*=rij[D]*w[i][j]/(dis[i][j]*dis[i][j]);
						gra_accel_i[D]+=accel_i[D];
					}
				}
				for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]-=Dt*mi*3/pnd*gra_accel_i[D];
			}
		}

	for(int i=0;i<h_num;i++)
	{
		delete[] w[i];
		delete[] dis[i];
		delete[] NEI_w[i];
	}
	delete[] w;
	delete[] dis;
	delete[]NEI_w;
	cout<<"----------OK"<<endl;
}

void contact_judge_hyper2(mpsconfig CON, vector<mpselastic> &PART, vector<hyperelastic> &HYPER, int hyper_number, int t)
{
	//アルゴリズム
	// 0. i周辺の粒子数密度が増加した場合，影響半径内にある粒子を探索し，以下を行う
	// 1. 「接触の可能性がある粒子」（(PART[j].PND>PART[j].PND0)が真？）を調べる
	// 2. 圧力を置換
	// 3. 初期配置の粒子と重複しないように接触の可能性がある粒子との間で力を計算する

	cout<<"壁影響計算_距離関数";

	int h_num=hyper_number;
	int p_num=PART.size();
	double r=CON.get_h_dis();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double Dt=CON.get_dt();
	
	double y_min=(atan(-10.0)+PI*0.5)/PI;
	double y_max=(atan(10.0)+PI*0.5)/PI;

	//壁が平面の場合
	/*
	if(t==1)	calculation_vec_norm(PART,HYPER,h_num,p_num,t);
	else if(t!=1||PART[h_num].u[A_Z]!=0)
	{
		for(int i=h_num;i<p_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].r[D]+=PART[i].u[D]*Dt;
		calculation_vec_norm(PART,HYPER,h_num,p_num,t);
	}*/

	stringstream s;
	s<<"./Wall/dp"<<t<<".dat";
	ofstream dp(s.str());

	ofstream s_pnd("pnd.csv", ios::app);
	if(t==1)
	{
		s_pnd<<"t"<<",";
		for(int i=0;i<h_num;i++)	s_pnd<<i<<",";
	}
	s_pnd<<t<<",";

	double vec_norm[DIMENSION]={0,0,1};
	//for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[h_num].vec_norm[D]=vec_norm[D];

	for(int i=0;i<h_num;i++)
	{
		double pnd=0;
		double max_dp[DIMENSION]={-2*HYPER[i].p[0]*vec_norm[0],	-2*HYPER[i].p[1]*vec_norm[1],	-2*HYPER[i].p[2]*vec_norm[2]};
		double dp_i[DIMENSION]={0,	0,	0};

		for(int j=h_num;j<p_num;j++)
		{
			double vec_dis[DIMENSION]={PART[j].r[0]-PART[i].r[0],	PART[j].r[1]-PART[i].r[1],	PART[j].r[2]-PART[i].r[2]};
			double dis=0;
			double	in_vec=vec_dis[0]*vec_dis[0]+vec_dis[1]*vec_dis[1]+vec_dis[2]*vec_dis[2];
			dis=sqrt(in_vec);
			if(dis<r)
			{
				double w=kernel(r,dis);
				pnd+=w;
				
				//高さ方向の距離計算
				double h=vec_norm[0]*vec_dis[0]+vec_norm[1]*vec_dis[1]+vec_norm[2]*vec_dis[2];
				double avs_h=fabs(h);
				double x=10-avs_h/r*20;
				double y=((atan(x)+PI*0.5)/PI-y_min)/(y_max-y_min);
				double dp_j[DIMENSION]={max_dp[0]*y,	max_dp[1]*y,	max_dp[2]*y};
				
				dp_i[A_X]+=dp_j[A_X]*w;
				dp_i[A_Y]+=dp_j[A_Y]*w;
				dp_i[A_Z]+=dp_j[A_Z]*w;
			}			
		}
		if(pnd!=0 && HYPER[i].p[A_Z]<=0)
		{
			HYPER[i].p[A_X]+=dp_i[A_X]/pnd;
			HYPER[i].p[A_Y]+=dp_i[A_Y]/pnd;
			HYPER[i].p[A_Z]+=dp_i[A_Z]/pnd;
		}

		s_pnd<<pnd<<",";
		dp<<i<<","<<dp_i[A_X]<<","<<dp_i[A_Y]<<","<<dp_i[A_Z]<<","<<endl;
	}
	s_pnd<<endl;
	s_pnd.close();
	dp.close();
	cout<<"----------OK"<<endl;
}

void calculation_vec_norm(vector<mpselastic> PART, vector<hyperelastic> &HYPER,int hyper_number,int particle_number,int t)
{
	int h_num=hyper_number;
	int p_num=particle_number;

	//法線ベクトル計算	壁が平面であること前提
	int id_norm[3];
	int count_max=0;
	double maxZ=0;
	double maxY=0;
	double maxX=0;
	for(int i=h_num;i<p_num;i++)
	{			
		if(PART[i].surface==ON)
		{
			if(maxZ<PART[i].r[A_Z])
			{
				maxZ=PART[i].r[A_Z];
				maxY=PART[i].r[A_Y];
				maxX=PART[i].r[A_X];
			}
			else if(maxZ=PART[i].r[A_Z])
			{
				if(maxY!=PART[i].r[A_Y]&&maxX!=PART[i].r[A_X])
				{
					maxZ=PART[i].r[A_Z];
					maxY=PART[i].r[A_Y];
					maxX=PART[i].r[A_X];
					id_norm[count_max]=i;
					cout<<"粒子番号"<<id_norm[count_max]<<"Z座標"<<maxZ<<endl;
					count_max++;
				}
			}
		}
		if(count_max==3)	break;
	}

	double *r=new double [3*3];
	double *vec=new double [DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		vec[D]=0;
		for(int D2=0;D2<DIMENSION;D2++)	r[D*DIMENSION+D2]=0;
	}

	double r_n[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{					
		r[0*DIMENSION+D]=PART[id_norm[0]].r[D];
		r[1*DIMENSION+D]=PART[id_norm[1]].r[D];
		r[2*DIMENSION+D]=PART[id_norm[2]].r[D];
		r_n[D]=(r[0*DIMENSION+D]+r[1*DIMENSION+D]+r[2*DIMENSION+D])/3;
		cout<<r_n[D];
	}
	cout<<endl;
	for(int nn=0;nn<3;nn++)	for(int D=0;D<DIMENSION;D++)	r[nn*DIMENSION+D]-=r_n[D];
	for(int nn=0;nn<3;nn++)	for(int D=0;D<DIMENSION;D++)	cout<<r[nn*DIMENSION+D]<<endl;;

	gauss(r,vec,DIMENSION);
	cout<<"x"<<vec[A_X]<<"y"<<vec[A_Y]<<"z"<<vec[A_Z]<<endl;
	double in_vec=vec[A_X]*vec[A_X]+vec[A_Y]*vec[A_Y]+vec[A_Z]*vec[A_Z];
	double ab_vec=sqrt(in_vec);
	for(int D=0;D<DIMENSION;D++)	HYPER[h_num].vec_norm[D]=vec[D]/ab_vec;
	cout<<"x"<<HYPER[h_num].vec_norm[A_X]<<"y"<<HYPER[h_num].vec_norm[A_Y]<<"z"<<HYPER[h_num].vec_norm[A_Z]<<endl;
	
	delete[]	r;
	delete[]	vec;
}

void contact_judge(mpsconfig &CON, vector<mpselastic> PART,vector<hyperelastic> &HYPER,double max_h,int t)
{
	int h_num=HYPER.size();
	double r=CON.get_h_dis();
	double Dt=CON.get_interval();

	for(int i=0;i<h_num;i++)
	{
		double qik_z=PART[i].r[A_Z]-max_h;
		if(qik_z<r)	HYPER[i].p[A_Z]-=HYPER[i].p[A_Z]/qik_z*Dt;
	}	
}


void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,int t)
{
	int h_num=HYPER.size();


	/*if(t==1)
	{
		////計算した各定数の出力
		ofstream ai("Ai.csv");
		ofstream inai("inverse_Ai.csv");
		ofstream aiin("aiin.csv");
		ofstream n0("noij.csv");
		ofstream wiin("wiin.csv");

		aiin<<"i"<<","<<"j"<<","<<"x"<<","<<"y"<<","<<"z"<<endl;
		n0<<"i"<<","<<"j"<<","<<"x"<<","<<"y"<<","<<"z"<<endl;

		for(int j=0;j<h_num;j++)	wiin<<j<<",";
		wiin<<endl;

		for(int i=0;i<h_num;i++)
		{
			inai<<i<<",";
			inai<<HYPER[i].inverse_Ai[0][A_X]<<","<<HYPER[i].inverse_Ai[0][A_Y]<<","<<HYPER[i].inverse_Ai[0][A_Z]<<","<<endl;		
			inai<<","<<HYPER[i].inverse_Ai[1][A_X]<<","<<HYPER[i].inverse_Ai[1][A_Y]<<","<<HYPER[i].inverse_Ai[1][A_Z]<<","<<endl;	
			inai<<","<<HYPER[i].inverse_Ai[2][A_X]<<","<<HYPER[i].inverse_Ai[2][A_Y]<<","<<HYPER[i].inverse_Ai[2][A_Z]<<","<<endl;	

			ai<<i<<",";
			ai<<HYPER[i].Ai[0][A_X]<<","<<HYPER[i].Ai[0][A_Y]<<","<<HYPER[i].Ai[0][A_Z]<<","<<endl;
			ai<<","<<HYPER[i].Ai[1][A_X]<<","<<HYPER[i].Ai[1][A_Y]<<","<<HYPER[i].Ai[1][A_Z]<<","<<endl;
			ai<<","<<HYPER[i].Ai[2][A_X]<<","<<HYPER[i].Ai[2][A_Y]<<","<<HYPER[i].Ai[2][A_Z]<<","<<endl;
			
			for(int j=0;j<h_num;j++)
			{
				aiin<<i<<","<<j<<","<<HYPER1[i*h_num+j].aiin[A_X]<<","<<HYPER1[i*h_num+j].aiin[A_Y]<<","<<HYPER1[i*h_num+j].aiin[A_Z]<<","<<endl;
				n0<<i<<","<<j<<","<<HYPER1[i*h_num+j].n0ij[A_X]<<","<<HYPER1[i*h_num+j].n0ij[A_Y]<<","<<HYPER1[i*h_num+j].n0ij[A_Z]<<","<<endl;
				wiin<<i<<","<<HYPER1[i*h_num+j].wiin<<",";
			}
			wiin<<endl;
		}
		ai.close();
		inai.close();
		n0.close();
		wiin.close();
	}*/

/*	stringstream ss_dgdq;
	ss_dgdq<<"./DgDq/DgDq"<<t<<".csv";
	ofstream dg(ss_dgdq.str());
	ofstream d_p("d_P.csv", ios::app);
	ofstream h_p("h_P.csv", ios::app);*/
/*	ofstream stress("stress.csv", ios::app);*/
/*	ofstream J("J.csv", ios::app);
	ofstream ti_Fi("ti_Fi.csv", ios::app);
	ofstream Fi("Fi.csv", ios::app);*/
	ofstream lam("lambda.csv", ios::app);
	lam<<HYPER[182].lambda<<","<<HYPER[911].lambda<<","<<HYPER[1640].lambda<<","<<HYPER[2369].lambda<<endl;
	lam.close();

	ofstream j("J.csv", ios::app);
	j<<HYPER[182].J<<","<<HYPER[911].J<<","<<HYPER[1640].J<<","<<HYPER[2369].J<<endl;
	j.close();


	ofstream p("P.csv", ios::app);
	p<<sqrt(HYPER[182].p[A_X]*HYPER[182].p[A_X]+HYPER[182].p[A_Z]*HYPER[182].p[A_Z])<<",";
	p<<sqrt(HYPER[911].p[A_X]*HYPER[911].p[A_X]+HYPER[911].p[A_Z]*HYPER[911].p[A_Z])<<",";
	p<<sqrt(HYPER[1640].p[A_X]*HYPER[1640].p[A_X]+HYPER[1640].p[A_Z]*HYPER[1640].p[A_Z])<<",";	
	p<<sqrt(HYPER[2369].p[A_X]*HYPER[2369].p[A_X]+HYPER[2369].p[A_Z]*HYPER[2369].p[A_Z])<<endl;
	p.close();

	/*
	if(t==1)
	{
		p<<"t"<<",";		
		d_p<<"t"<<",";
		h_p<<"t"<<",";/
		lam<<"t"<<",";
		stress<<"t"<<",";
		ti_Fi<<"t"<<",";
		Fi<<"t"<<",";
		J<<"t"<<",";//

		for(int i=0;i<h_num;i++)
		{
			p<<i<<","<<","<<",";
			lam<<i<<",";
			d_p<<i<<","<<","<<",";
			h_p<<i<<","<<","<<",";
			stress<<i<<","<<","<<",";
			ti_Fi<<i<<","<<","<<",";
			Fi<<i<<","<<","<<",";
			J<<i<<",";
		}
		p<<endl;
		lam<<endl;
		d_p<<endl;
		h_p<<endl;
		stress<<endl;
		ti_Fi<<endl;
		Fi<<endl;
		J<<endl;
	}

	p<<t<<",";
	lam<<t<<",";
	d_p<<t<<",";
	h_p<<t<<",";
	stress<<t<<",";
	ti_Fi<<t<<",";
	Fi<<t<<",";
	J<<t<<",";
	for(int i=0;i<h_num;i++)
	{
		
		p<<HYPER[i].p[A_X]<<","<<HYPER[i].p[A_Y]<<","<<HYPER[i].p[A_Z]<<",";
		lam<<HYPER[i].lambda<<",";
/*		d_p<<HYPER[i].differential_p[A_X]<<","<<HYPER[i].differential_p[A_Y]<<","<<HYPER[i].differential_p[A_Z]<<","<<endl;
		h_p<<HYPER[i].half_p[A_X]<<","<<HYPER[i].half_p[A_Y]<<","<<HYPER[i].half_p[A_Z]<<","<<endl;

		stress<<HYPER[i].stress[A_X][A_X]<<","<<HYPER[i].stress[A_X][A_Y]<<","<<HYPER[i].stress[A_X][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_X][A_X]<<","<<HYPER[i].t_inverse_Fi[A_X][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_X][A_Z]<<",";	
		Fi<<HYPER[i].Fi[0][0]<<","<<HYPER[i].Fi[0][1]<<","<<HYPER[i].Fi[0][2]<<",";
		J<<HYPER[i].J<<",";
		for(int j=0;j<h_num;j++)	dg<<i<<","<<j<<","<<HYPER1[j*h_num+i].DgDq[A_X]<<","<<HYPER1[j*h_num+i].DgDq[A_Y]<<","<<HYPER1[j*h_num+i].DgDq[A_Z]<<","<<endl;
		
	}
	stress<<endl<<",";
	ti_Fi<<endl<<",";
	Fi<<endl<<",";
	for(int i=0;i<h_num;i++)
	{
		stress<<HYPER[i].stress[A_Y][A_X]<<","<<HYPER[i].stress[A_Y][A_Y]<<","<<HYPER[i].stress[A_Y][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_Y][A_X]<<","<<HYPER[i].t_inverse_Fi[A_Y][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_Y][A_Z]<<",";
		Fi<<HYPER[i].Fi[1][0]<<","<<HYPER[i].Fi[1][1]<<","<<HYPER[i].Fi[1][2]<<",";
	}
	stress<<endl<<",";
	ti_Fi<<endl<<",";
	Fi<<endl<<",";
	for(int i=0;i<h_num;i++)
	{
		stress<<HYPER[i].stress[A_Z][A_X]<<","<<HYPER[i].stress[A_Z][A_Y]<<","<<HYPER[i].stress[A_Z][A_Z]<<",";
		ti_Fi<<HYPER[i].t_inverse_Fi[A_Z][A_X]<<","<<HYPER[i].t_inverse_Fi[A_Z][A_Y]<<","<<HYPER[i].t_inverse_Fi[A_Z][A_Z]<<",";
		Fi<<HYPER[i].Fi[2][0]<<","<<HYPER[i].Fi[2][1]<<","<<HYPER[i].Fi[2][2]<<",";
	}

	p<<endl;
	lam<<endl;
	d_p<<endl;
	h_p<<endl;
	stress<<endl;
	ti_Fi<<endl;
	Fi<<endl;
	J<<endl;

	dg.close();
	stress.close();
	ti_Fi.close();
	Fi.close();
	d_p.close();
	h_p.close();
	J.close();
	p.close();
	lam.close();*/
}

void output_newton_data1(double *fx, double *DfDx, double *n_rx, double *n_ry, double *n_rz,int hyper_number,int count, int t)
{
	int h_num=hyper_number;

	stringstream ss_r;
	ss_r<<"./Newton_raphson/position"<<t<<".xlsx";
	stringstream ss_Df;
	ss_Df<<"./Newton_raphson/DfDx "<<t<<".xlsx";
	stringstream ss_fx;
	ss_fx<<"./Newton_raphson/fx"<<t<<".xlsx";
		
	if(count==1)
	{
		ofstream init0(ss_r.str(), ios::trunc);
		ofstream init1(ss_Df.str(), ios::trunc);
		ofstream init3(ss_fx.str(), ios::trunc);
	
		init0.close();
		init1.close();
		init3.close();
	}

	ofstream r(ss_r.str(), ios::app);
	ofstream sfx(ss_fx.str(), ios::app);

	ofstream Df(ss_Df.str(), ios::app);

	if(count==1)
	{
		r<<"count"<<",";
		sfx<<"count"<<",";
		for(int i=0;i<h_num;i++)
		{
			r<<i<<","<<","<<",";
			sfx<<i<<",";
		}
		r<<endl;
		sfx<<endl;
	}
	r<<count<<",";
	sfx<<count<<",";

	Df<<"count"<<","<<count<<endl;
	for(int i=0;i<h_num;i++)
	{
		r<<n_rx[i]<<","<<n_ry[i]<<","<<n_rz[i]<<",";
		sfx<<fx[i]<<",";
		for(int j=0;j<h_num;j++) Df<<DfDx[i*h_num+j]<<",";
		Df<<endl;

	}
	r<<endl;
	sfx<<endl;

	r.close();
	sfx.close();

	Df.close();
}

void output_newton_data2(double E, double *XX, int hyper_number, int count, int t)
{
	int h_num=hyper_number;
	stringstream ss_E;
	ss_E<<"./Newton_raphson/E"<<t<<".xlsx";
	
	stringstream ss_lam;
	ss_lam<<"./Newton_raphson/lambda"<<t<<".xlsx";
		
	if(count==1)
	{
		ofstream init0(ss_E.str(), ios::trunc);
		ofstream init1(ss_lam.str(), ios::trunc);
	
		init0.close();
		init1.close();
	}

	ofstream e(ss_E.str(), ios::app);
	ofstream lam(ss_lam.str(), ios::app);

	if(count==1)
	{
		e<<"反復回数"<<","<<"E"<<endl;
		lam<<"反復回数"<<","<<"lambda"<<endl;
		for(int i=0;i<h_num;i++)	lam<<","<<i;
		lam<<endl;
	}
	
	e<<count<<","<<E<<endl;

	lam<<count;
	for(int i=0;i<h_num;i++)	lam<<","<<XX[i];
	lam<<endl;

	e.close();
	lam.close();
}

void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER,int t)
{
//	cout<<"弾性ポテンシャル計算";
	int h_num=HYPER.size();
	int p_num=PART.size();
	double d_Fi[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	vector<double>	W;
	W.reserve(h_num);
	for(int i=0;i<h_num;i++)	W.emplace_back(0);
	double dC[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	double dC2[DIMENSION][DIMENSION]={{0,0,0},{0,0,0},{0,0,0}};
	/*
	for(int j=0;j<4;j++)
	{
		if(j==0)	i=182;
		else if(j==1) i=911;
		else if(j==2) i=1640;
		else if(j==3) i=2369;

		double J=HYPER[i].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[i].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[i].Fi[2][2];
		}

		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

		dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
		dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
		dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
		dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
		dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
		dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
		dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
		dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
		dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
		W[i]=c10*(Ic-3)+c01*(IIc-3);
	}*/

	
	for(int i=0;i<h_num;i++)
	{
		double J=HYPER[i].J;
		if(J<0){
			d_Fi[0][0]=-1/pow(-J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=-1/pow(-J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=-1/pow(-J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=-1/pow(-J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=-1/pow(-J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=-1/pow(-J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=-1/pow(-J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=-1/pow(-J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=-1/pow(-J,1/3)*HYPER[i].Fi[2][2];
		}
		else
		{
			d_Fi[0][0]=1/pow(J,1/3)*HYPER[i].Fi[0][0];	d_Fi[0][1]=1/pow(J,1/3)*HYPER[i].Fi[0][1];	d_Fi[0][2]=1/pow(J,1/3)*HYPER[i].Fi[0][2];
			d_Fi[1][0]=1/pow(J,1/3)*HYPER[i].Fi[1][0];	d_Fi[1][1]=1/pow(J,1/3)*HYPER[i].Fi[1][1];	d_Fi[1][2]=1/pow(J,1/3)*HYPER[i].Fi[1][2];
			d_Fi[2][0]=1/pow(J,1/3)*HYPER[i].Fi[2][0];	d_Fi[2][1]=1/pow(J,1/3)*HYPER[i].Fi[2][1];	d_Fi[2][2]=1/pow(J,1/3)*HYPER[i].Fi[2][2];
		}

		dC[0][0]=d_Fi[0][0]*d_Fi[0][0]+d_Fi[1][0]*d_Fi[1][0]+d_Fi[2][0]*d_Fi[2][0];
		dC[0][1]=d_Fi[0][0]*d_Fi[0][1]+d_Fi[1][0]*d_Fi[1][1]+d_Fi[2][0]*d_Fi[2][1];
		dC[0][2]=d_Fi[0][0]*d_Fi[0][2]+d_Fi[1][0]*d_Fi[1][2]+d_Fi[2][0]*d_Fi[2][2];
		dC[1][0]=d_Fi[0][1]*d_Fi[0][0]+d_Fi[1][1]*d_Fi[1][0]+d_Fi[2][1]*d_Fi[2][0];
		dC[1][1]=d_Fi[0][1]*d_Fi[0][1]+d_Fi[1][1]*d_Fi[1][1]+d_Fi[2][1]*d_Fi[2][1];
		dC[1][2]=d_Fi[0][1]*d_Fi[0][2]+d_Fi[1][1]*d_Fi[1][2]+d_Fi[2][1]*d_Fi[2][2];
		dC[2][0]=d_Fi[0][2]*d_Fi[0][0]+d_Fi[1][2]*d_Fi[1][0]+d_Fi[2][2]*d_Fi[2][0];
		dC[2][1]=d_Fi[0][2]*d_Fi[0][1]+d_Fi[1][2]*d_Fi[1][1]+d_Fi[2][2]*d_Fi[2][1];
		dC[2][2]=d_Fi[0][2]*d_Fi[0][2]+d_Fi[1][2]*d_Fi[1][2]+d_Fi[2][2]*d_Fi[2][2];

		dC2[0][0]=dC[A_X][0]*dC[0][A_X]+dC[A_X][1]*dC[1][A_X]+dC[A_X][2]*dC[2][A_X];
		dC2[0][1]=dC[A_X][0]*dC[0][A_Y]+dC[A_X][1]*dC[1][A_Y]+dC[A_X][2]*dC[2][A_Y];
		dC2[0][2]=dC[A_X][0]*dC[0][A_Z]+dC[A_X][1]*dC[1][A_Z]+dC[A_X][2]*dC[2][A_Z];
		dC2[1][0]=dC[A_Y][0]*dC[0][A_X]+dC[A_Y][1]*dC[1][A_X]+dC[A_Y][2]*dC[2][A_X];
		dC2[1][1]=dC[A_Y][0]*dC[0][A_Y]+dC[A_Y][1]*dC[1][A_Y]+dC[A_Y][2]*dC[2][A_Y];
		dC2[1][2]=dC[A_Y][0]*dC[0][A_Z]+dC[A_Y][1]*dC[1][A_Z]+dC[A_Y][2]*dC[2][A_Z];
		dC2[2][0]=dC[A_Z][0]*dC[0][A_X]+dC[A_Z][1]*dC[1][A_X]+dC[A_Z][2]*dC[2][A_X];
		dC2[2][1]=dC[A_Z][0]*dC[0][A_Y]+dC[A_Z][1]*dC[1][A_Y]+dC[A_Z][2]*dC[2][A_Y];
		dC2[2][2]=dC[A_Z][0]*dC[0][A_Z]+dC[A_Z][1]*dC[1][A_Z]+dC[A_Z][2]*dC[2][A_Z];

		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		double Ic=trace_dC;
		double IIc=0.50*(trace_dC*trace_dC-trace_dC2);
		W[i]=c10*(Ic-3)+c01*(IIc-3);
	}//*/
//	cout<<"----------OK"<<endl;

	ofstream e("E.csv", ios::app);
	ofstream e_T("E_T.csv", ios::app);
	ofstream e_g("E_g.csv", ios::app);
	ofstream e_W("E_W.csv", ios::app);
	ofstream e_lam("E_lam.csv", ios::app);

	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
/*	double sum_e_T=0;
	double sum_e_g=0;
	double sum_e_lam=0;
	double sum_e_W=0;
	double sum_e=0;
	double vv=0;
	double energy=0;

	for(int j=0;j<4;j++)
	{
		if(j==0)	i=182;
		else if(j==1) i=911;
		else if(j==2) i=1640;
		else if(j==3) i=2369;

		vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		//energy=0.5/mi*vv+mi*9.8*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*9.8*PART[i].r[A_Z]<<",";
		e_W<<W[i]*V<<",";
		e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";
	}*/
	
	if(t==1)
	{
		e<<"t"<<",";
		e_T<<"t"<<",";
		e_g<<"t"<<",";
		e_W<<"t"<<",";
		e_lam<<"t"<<",";
		for(int i=0;i<h_num;i++)
		{
			e<<i<<",";
			e_T<<i<<",";
			e_g<<i<<",";
			e_W<<i<<",";
			e_lam<<i<<",";
		}
		e<<endl;
		e_T<<endl;
		e_g<<endl;
		e_W<<endl;
		e_lam<<endl;
	}

	e<<t<<",";
	e_T<<t<<",";
	e_g<<t<<",";
	e_W<<t<<",";
	e_lam<<t<<",";
	double sum_e_T=0;
	double sum_e_g=0;
	double sum_e_lam=0;
	double sum_e_W=0;
	double sum_e=0;
	double vv=0;
	double energy=0;

	for(int i=0;i<h_num;i++)
	{
		vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		//energy=0.5/mi*vv+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		energy=0.5/mi*vv+mi*9.8*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_T+=0.5/mi*vv;
		sum_e_g+=mi*9.8*PART[i].r[A_Z];
		sum_e_lam+=HYPER[i].lambda*(1-HYPER[i].J)*V;
		sum_e_W+=W[i]*V;
		sum_e+=energy;
		e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*9.8*PART[i].r[A_Z]<<",";
		e_W<<W[i]*V<<",";
		e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";
	}
	e<<sum_e<<endl;
	e_T<<sum_e_T<<endl;
	e_g<<sum_e_g<<endl;
	e_W<<sum_e_W<<endl;
	e_lam<<sum_e_lam<<endl;//*/

	e.close();
	e_T.close();
	e_g.close();
	e_W.close();
	e_lam.close();
}

void calc_gravity(mpsconfig CON,vector<hyperelastic> &HYPER,int hyper_number)
{
	int h_num=hyper_number;
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	//重力降下
	double Dt=CON.get_dt();
	for(int i=0;i<h_num;i++)	HYPER[i].p[A_Z]-=9.8*Dt;
}

void transpose(double **M, double **N)
{
	N[0][0]=M[0][0];
	N[1][1]=M[1][1];
	N[2][2]=M[2][2];
	N[0][1]=M[1][0];
	N[0][2]=M[2][0];
	N[1][0]=M[0][1];
	N[1][2]=M[2][1];
	N[2][0]=M[0][2];
	N[2][1]=M[1][2];
}


void BiCGStab2_method(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :非ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解
	
	int count=0;
	double pk=1;
	double E=1;//誤差
	double alp,beta,rr,w,ita;

	double *r=new double [pn];	//残差
	double *P=new double [pn];	//探索ベクトル
	double *Pj=new double [pn];
	double *rj=new double [pn];
	double *AP=new double [pn];
	double *AtPj=new double [pn];
	double *e=new double [pn];
	double *Ae=new double [pn];
	double *y=new double [pn];
	double *u=new double [pn];
	double *W=new double [pn];
	double *Z=new double [pn];
	double *e_pre=new double [pn];
	beta=0;
	w=0;
	ita=0;

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0;
		r[n]=B[n];
	}

	for(int n=0;n<pn;n++)
	{
		P[n]=r[n];//初期化
		rj[n]=r[n];
		Pj[n]=rj[n];
		AP[n]=0;
		y[n]=0;
		u[n]=0;
		W[n]=0;
		Z[n]=0;
		e[n]=0;
		e_pre[n]=0;//1ステップ前のe[]
	}
	double rr0=0;
	for(int n=0;n<pn;n++) rr0+=r[n]*r[n];
	cout<<"rr0="<<rr0<<endl;
	 cout<<"BiCGstab2法スタート  -----未知数="<<pn<<"  ---";
	 unsigned int time=GetTickCount();
	double ep=CON->get_FEMCGep();//収束判定

	while(E>ep)// EP=CON->get_CGep();//収束判定(convergence test)
	{
		count++;

		for(int n=0;n<pn;n++)
		{
			//P[n]=r[n]+beta*(P[n]-w*AP[n]);
			P[n]=r[n]+beta*(P[n]-u[n]);
		} 

		////pk(rとrjの内積)を求める
		pk=0;
		for(int n=0;n<pn;n++) pk+=r[n]*rj[n];

		//////////////alpを求める
		for(int n=0;n<pn;n++)
		{      
			AP[n]=0;
			//cout<<"ptrn"<<ptr[n]<<"ptrn+1"<<ptr[n+1]<<endl;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			//for(int m=0;m<pn;m++) AP[n]+=A[n][m]*P[m];
		}
		double APrj=0;
		for(int n=0;n<pn;n++)  APrj+=rj[n]*AP[n];
		alp=pk/APrj;
		//cout<<"alp="<<alp<<" APrj="<<APrj<<endl;
		//////////////////////

		for(int n=0;n<pn;n++) y[n]=e[n]-r[n]-alp*W[n]+alp*AP[n];

		for(int n=0;n<pn;n++)
		{
			e_pre[n]=e[n];//変更前の値を記憶
			e[n]=r[n]-alp*AP[n];
		}

		for(int n=0;n<pn;n++)
		{
			Ae[n]=0;
			//for(int m=0;m<pn;m++) Ae[n]+=A[n][m]*e[m];
			for(int m=ptr[n];m<ptr[n+1];m++) Ae[n]+=val[m]*e[ind[m]];
		}

		if(count%2!=0)
		{
			double e_dot_ae  = 0;
			double ae_dot_ae = 0;
			for(int n=0;n<pn;n++)
			{
				e_dot_ae+=e[n]*Ae[n];
				ae_dot_ae+=Ae[n]*Ae[n];
			}
			w = e_dot_ae / ae_dot_ae;
			ita=0;
		}
		else
		{
			double e_dot_ae  = 0;
			double y_dot_y=0;
			double y_dot_e=0;
			double Ae_dot_y=0;
			double Ae_dot_Ae=0;
			for(int n=0;n<pn;n++) e_dot_ae+=e[n]*Ae[n];
			for(int n=0;n<pn;n++) y_dot_y+=y[n]*y[n];
			for(int n=0;n<pn;n++) y_dot_e+=y[n]*e[n];
			for(int n=0;n<pn;n++) Ae_dot_y+=Ae[n]*y[n];
			for(int n=0;n<pn;n++) Ae_dot_Ae+=Ae[n]*Ae[n];

			double co=Ae_dot_Ae*y_dot_y-Ae_dot_y*Ae_dot_y;

			w=y_dot_y*e_dot_ae-y_dot_e*Ae_dot_y;
			w/=co;

			ita=Ae_dot_Ae*y_dot_e-Ae_dot_y*e_dot_ae;
			ita/=co;
		}
		

		for(int n=0;n<pn;n++) 
		{
			u[n]=w*AP[n]+ita*(e_pre[n]-r[n]+beta*u[n]);
			Z[n]=w*r[n]+ita*Z[n]-alp*u[n];
			X[n]+=alp*P[n]+Z[n];
			r[n]=e[n]-ita*y[n]-w*Ae[n];
		}

		///beta
		beta=0;
		for(int n=0;n<pn;n++) beta+=r[n]*rj[n];
		beta/=w*APrj;
		
		//cout<<"beta="<<beta<<endl;
		//////////////////////

		//W[n]
		for(int n=0;n<pn;n++) W[n]=Ae[n]+beta*AP[n];
		
		//////////////////誤差
		rr=0;
		for(int n=0;n<pn;n++) rr+=r[n]*r[n];
		E=rr/rr0;
		//E=sqrt(rr);
		if(count==1||count%10==0)	cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////
	}
	
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;

	delete [] r;
	delete [] Pj;
	delete [] rj;
	delete [] AP;
	delete [] AtPj;
	delete [] P;
	delete [] e;
	delete [] Ae;
	delete [] y;
	delete [] u;
	delete [] W;
	delete [] Z;
	delete [] e_pre;
}


void iccg2(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	double accel=0.87;//CON->get_CGaccl();//加速ファクタ
	
	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *val2=new double [num2];
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];
	double *r=new double [pn];
	double *P=new double [pn];

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
		r[k]=0;
		P[k]=0;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;//加速ﾌｧｸﾀ
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];//ゼロ要素の値
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg法
	double alp,beta;
	double rLDLt_r;
	double E=1;//誤差
	double *AP = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *D1 = new double [pn];//D行列
	
	/////不完全コレスキｰ分解
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
	        int i=ind2[m];//列番号
	        if(i==0)
			{
				val2[m]=val2[m];
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
		    
			}
			if(i>0 && i<k)
			{
				double sum=0;
				
				for(int j=ptr2[k];j<m;j++)
				{	
					for(int J=ptr2[i];J<ptr2[i+1];J++)//i行目のなかから列の一致するものを探している。少し手間か？
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				double sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
				else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
				D1[k]=1/val2[m];
				//if(val2[m]>0) cout<<"EE"<<endl;
            }
	    }
	}    
	///不完全コレスキー分解完了/////////*/

	///列を基準にした配列に値を代入
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=val2[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 
		else
		{
		    double sum=0;
		    /////////        
		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y[ind2[m]];//式（3.78）
		    int m=ptr2[i+1]-1;
		    y[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]がもとまった。
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y[i]-D1[i]*sum;	
	}
	/////////////////*/
	
	for(int n=0;n<pn;n++){ P[n]=LDLt_r[n];
	cout<<"P"<<n<<"	"<<P[n]<<endl;
	}

	cout<<"ICCG法:未知数="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=CON->get_FEMCGep();//収束判定
	rLDLt_r=0;
	for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//最初のrLDLt_rだけここで求める
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alpを求める
		double PAP=0;
		#pragma omp parallel for reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		cout<<"alp="<<alp<<endl;
		//////////////////////
		E=0;
		#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++)
		{
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P 更新後の場所
			r[n]-=alp*AP[n];// r=r-alp*AP       更新後の残差
			E+=r[n]*r[n];						//更新後の誤差
		}
		E=sqrt(E);
		cout<<"E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/val2[0]; //式（3.77） 新
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//対角成分は除くからptr[i+1]-1
			    {
			        sum+=val2[m]*y[ind2[m]];//式（3.78）
			    }
			    int m=ptr2[i+1]-1;
			    y[i]=(r[i]-sum)/val2[m];
			}
		}////y[i]がもとまった。
	
		/////////LDLt_r[i]を求める
		for(int i=pn-1;i>=0;i--)
		{
		    double sum=0;
			for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
			
		    LDLt_r[i]=y[i]-D1[i]*sum;	
		}
		/////////////////*/
	
		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		beta=beta*rLDLt_r;
		/////////////////*/
		
		///////////////////// P=r+beta*P
		for(int n=0;n<pn;n++) P[n]=LDLt_r[n]+beta*P[n];//iccg
		if(count==1||count%10==0)	cout<<"E="<<E<<" count="<<count<<endl;
	}
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";
		
	delete [] AP;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;
	delete [] r;
	delete [] P;
	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	//*count2=count;//反復回数を格納して返す
}

void CG3D(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	cout<<"CG法:未知数="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=CON->get_FEMCGep();//収束判定

	double alp,beta;
	double E=1;//誤差
	double *AP = new double [pn];
	double *r=new double [pn];
	double *P=new double [pn];

	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alpを求める
		double PAP=0;
		double Pr=0;
		#pragma omp parallel for reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
			Pr+=P[n]*r[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=Pr/PAP;
		cout<<"alp="<<alp<<endl;
		//////////////////////
		E=0;
		#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++)
		{
			X[n]+=alp*P[n];// X(k+1)=X(k)+alp*P 更新後の場所
			r[n]-=alp*AP[n];// r=r-alp*AP       更新後の残差
			E+=r[n]*r[n];						//更新後の誤差
		}
		E=sqrt(E);
		if(count==1||count%10==0)	cout<<"E="<<E<<" count="<<count<<endl;
		////////////////////////
		
		///////////////////////beta
		double rAP=0;
		for(int n=0;n<pn;n++)	rAP+=r[n]*AP[n];
		beta=-rAP/PAP;
		cout<<"beta="<<beta<<endl;
		for(int n=0;n<pn;n++)	P[n]=r[n]+beta*P[n];

	}
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";
			
	delete [] AP;
	delete [] r;
	delete [] P;

}

hyperelastic::hyperelastic()
{
	for(int i=0;i<200;i++)
	{
		NEI[i]=0;
	}
	N=0;
	lambda=1;
	J=0;
	pnd=0;
	for(int D=0;D<DIMENSION;D++)
	{
		vec_norm[D]=0;
		half_p[D]=0;
		differential_p[D]=0;
		p[D]=0;	
		ang_p[D]=0;
		vis_force[D]=0;
		for(int D2=0;D2<DIMENSION;D2++)
		{
			stress[D][D2]=0;
			Ai[D][D2]=0;
			t_inverse_Ai[D][D2]=0;
			t_inverse_Fi[D][D2]=0;
			Fi[D][D2]=0;
		}
	}
}

hyperelastic2::hyperelastic2()
{
	wiin=0;
	spl_f=0;
	for(int D=0;D<DIMENSION;D++)
	{
		DgDq[D]=0;
		aiin[D]=0;
		n0ij[D]=0;
	}
}

