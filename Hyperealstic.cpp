#include "stdafx.h"		


void calc_hyper(mpsconfig &CON,vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int hyper_number,int t)
{	
	int N=hyper_number;
//	
	double Dt=CON.get_dt()*CON.get_interval();

	if(t==1)
	{
//		system("mkdir Momentum");
		system("mkdir d_Momentum");
		system("mkdir h_Momentum");
		system("mkdir Newton_raphson");
		system("mkdir Position");
		system("mkdir System");
		system("mkdir Ai");
		system("mkdir Fi");
		system("mkdir Lambda");
		system("mkdir renew_Lambda");
		system("mkdir Hy_stress");
	}

/*	stringstream sp1;
	sp1<<"./Momentum/momentum"<<t<<".dat";
	string filename=sp1.str();
	ofstream pf1(filename);*/

	stringstream ss;
	ss<<"./d_Momentum/differential_momentum"<<t<<".dat";
	string fl=ss.str();
	ofstream sdm(fl);

	stringstream ss2;
	ss2<<"./h_Momentum/half_momentum"<<t<<".dat";
	string fl2=ss2.str();
	ofstream shm(fl2);

	stringstream ss3;
	ss3<<"./Position/position"<<t<<".dat";
	string fl3=ss3.str();
	ofstream po(fl3);

	stringstream ss4;
	ss4<<"./Hy_stress/Stress"<<"t"<<t<<".csv";
	cout<<PART.size();
	string fl4=ss4.str();
	ofstream stress(fl4);

	stringstream ss5;
	ss5<<"./Fi/Fi "<<"t"<<t<<".csv";
	string fl5=ss5.str();
	ofstream fi(fl5);

	stringstream ss6;
	ss6<<"./Fi/det_Fi "<<"t"<<t<<".dat";
	string fl6=ss6.str();
	ofstream dfi(fl6);

	stringstream ss7;
	ss7<<"./renew_Lambda/renew_lambda "<<"t"<<t<<".dat";
	string fl7=ss7.str();
	ofstream lam(fl7);

	if(t==1)
	{
		for(int i=0;i<N;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=0;
		for(int i=0;i<N;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=PART[i].r[D];
		
		if(CON.get_ttest()==ON)	//引っ張り試験用15/2/8
		{
			double highest=0;
			for(int i=0;i<N;i++)	if(highest<PART[i].r[A_Z])	highest=PART[i].r[A_Z];
			for(int i=0;i<N;i++)
			{
				if(PART[i].r[A_Z]==highest)	HYPER[i].highest=ON;
				else
				{
					HYPER[i].highest=OFF;
				}
			}
		}

		calc_constant(CON,PART,HYPER,HYPER1,N);

/*		pf1<<"Momentum"<<t<<endl;
		for(int i=0;i<N;i++)	pf1<<i<<"	"<<HYPER[i].p[A_X]<<"	"<<HYPER[i].p[A_Y]<<"	"<<HYPER[i].p[A_Z]<<endl;
		pf1<<endl;*/

		calc_stress(CON,PART,HYPER,HYPER1,N);

		po<<"Position"<<t<<endl;		
		stress<<"Stress"<<t<<endl;
		dfi<<"det_Fi"<<t<<endl;
		fi<<"Fi"<<t<<endl;	
		for(int i=0;i<N;i++)
		{
			po<<i<<"	";
			for(int D=0;D<DIMENSION-1;D++)	po<<PART[i].r[D]<<"	";
			po<<PART[i].r[A_Z]<<endl;
	
			if(HYPER[i].stress0==ON)	stress<<i<<","<<0<<endl;
			else
			{				
				for(int D=0;D<DIMENSION;D++)
				{
					stress<<i<<",";
					for(int D2=0;D2<DIMENSION;D2++)	stress<<HYPER[i].stress[D][D2]<<",";
					stress<<endl;
				}
				stress<<endl;
			}

			dfi<<i<<"	"<<HYPER[i].J<<endl;
			fi<<i<<endl;
			for(int D=0;D<DIMENSION;D++)
			{
					for(int D2=0;D2<DIMENSION;D2++)	fi<<HYPER[i].Fi[D][D2]<<",";
					fi<<endl;
			}
		}
	}

	if(t==1 || t%CON.get_interval()==0)	momentum_movie_AVS(CON,t,PART,HYPER,N);

	if(t!=1)
	{
/*		for(int i=0;i<N;i++)
		{
			cout<<"PART["<<i<<"].r= ";
			for(int D=0;D<DIMENSION;D++)	cout<<PART[i].r[D]<<"	";
			cout<<endl;
		}
		cout<<endl;*/

		po<<"Position"<<t<<endl;		
		stress<<"Stress"<<t<<endl;
		dfi<<"det_Fi"<<t<<endl;
		fi<<"Fi"<<t<<endl;	
		for(int i=0;i<N;i++)
		{
			po<<i<<"	";
			for(int D=0;D<DIMENSION-1;D++)	po<<PART[i].r[D]<<"	";
			po<<PART[i].r[A_Z]<<endl;
			if(HYPER[i].stress0==ON)	stress<<i<<","<<0<<endl;
			else
			{	
				for(int D=0;D<DIMENSION;D++)	
				{
					stress<<i<<",";
					for(int D2=0;D2<DIMENSION;D2++)	stress<<HYPER[i].stress[D][D2]<<",";
					stress<<endl;
				}
				stress<<endl;
			}
			dfi<<i<<"	"<<HYPER[i].J<<endl;
			fi<<i<<endl;
			for(int D=0;D<DIMENSION;D++)
			{
					for(int D2=0;D2<DIMENSION;D2++)	fi<<HYPER[i].Fi[D][D2]<<",";
					fi<<endl;
			}
			fi<<endl;
		}
	}
	newton_raphson(CON,PART,HYPER,HYPER1,N,t);
	cout<<"Newton_raphson is ended."<<endl;

	calc_half_p(CON,PART,HYPER,HYPER1,N,0,t);

	for(int i=0;i<N;i++)
	{
		shm<<"half_p"<<i<<"	";
		for(int D=0;D<DIMENSION;D++)	shm<<HYPER[i].half_p[D]<<"	";
		shm<<endl;
	}
	
/*	rf<<"Renewed position"<<endl;
	for(int i=0;i<N;i++)
	{
		rf<<"PART["<<i<<"].r=";
		for(int D=0;D<DIMENSION;D++)	rf<<PART[i].r[D]<<"	";
		rf<<endl;
	}
	rf<<endl;*/

	calc_F(PART,HYPER,HYPER1,N,t);

	calc_stress(CON,PART,HYPER,HYPER1,N);

/*	for(int i=0;i<N;i++)
	{
			cout<<"HYPER["<<i<<"].stress="<<endl;
			for(int D=0;D<DIMENSION;D++)	
			{
				for(int D2=0;D2<DIMENSION;D2++)	cout<<HYPER[i].stress[D][D2]<<" ";
				cout<<endl;
			}
			cout<<endl;
	}
	cout<<endl;*/

	double p_DgDq[3];
	for(int j=0;j<N;j++)	for(int i=0;i<N;i++)	for(int D=0;D<DIMENSION;D++)	HYPER1[j*N+i].DgDq[D]=0;
	//calculation of DgDq
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<N;i++)
		{
			for(int D=0;D<DIMENSION;D++)
			{
				p_DgDq[D]=0;
				for(int D2=0;D2<DIMENSION;D2++)	p_DgDq[D]+=HYPER[j].t_inverse_Fi[D][D2]*HYPER1[i*N+j].n0ij[D2];
			}
/*			cout<<"p_DgDq["<<j<<"]["<<i<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<p_DgDq[D]<<" ";
			cout<<endl;*/
			for(int D=0;D<DIMENSION;D++)	HYPER1[j*N+i].DgDq[D]=HYPER[j].J*p_DgDq[D];
/*			cout<<"J["<<j<<"]="<<HYPER[j].J<<endl;
			cout<<"DgDq["<<j<<"*N+"<<i<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[j*N+i].DgDq[D]<<" ";
			cout<<endl;*/
		}
//		cout<<endl;
	}

/*	for(int j=0;j<N;j++)
	{
		for(int i=0;i<N;i++)
		{
			cout<<"HYPER1["<<j<<"]["<<i<<"].DgDq=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[j*N+i].DgDq[D]<<" ";
			cout<<endl;
		}
		cout<<endl;
	}
	cout<<endl;*/

	calc_differential_p(CON,PART,HYPER,HYPER1,N);
	cout<<"Differential_p is calculated."<<endl;

	sdm<<"Differential_p"<<t<<endl;
	for(int i=0;i<N;i++)	sdm<<i<<"	"<<HYPER[i].differential_p[A_X]<<"	"<<HYPER[i].differential_p[A_Y]<<"	"<<HYPER[i].differential_p[A_Z]<<endl;


	renew_lambda(CON,PART,HYPER,HYPER1,N);
	
	lam<<"renew_lambda"<<t<<endl;
	for(int i=0;i<N;i++)	lam<<i<<"	"<<HYPER[i].lambda<<endl;

	cout<<"Lambda is renewed."<<endl;

	calc_half_p(CON,PART,HYPER,HYPER1,N,1,t);

/*	pf1<<"Renewed momentum"<<endl;
	for(int i=0;i<N;i++)	pf1<<"HYPER["<<i<<"].p="<<HYPER[i].p[A_X]<<"	"<<HYPER[i].p[A_Y]<<"	"<<HYPER[i].p[A_Z]<<endl;
	for(int i=0;i<N;i++)	cout<<"HYPER["<<i<<"].p="<<HYPER[i].p[A_X]<<"	"<<HYPER[i].p[A_Y]<<"	"<<HYPER[i].p[A_Z]<<endl;*/

	cout<<"Hypercalculation is ended."<<endl;

	stress.close();
	fi.close();
	dfi.close();
	lam.close();
	sdm.close();
	shm.close();
//	pf1.close();
	po.close();
}

void calc_constant(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number)
{
	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt()*CON.get_interval();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	int num=particle_number;

	if(CON.get_flag_vis()==OFF&&CON.get_FEM_flag()==OFF)	for(int i=0;i<num;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]=0;

	//曲げねじり
	if(CON.get_model_number()==21)
	{
		int t=30,b=2;
		double H=0,Z=0,X=0,Y=0,max=0,min=0;

		for(int i=0;i<num;i++)
		{
			if(max<PART[i].q0[A_Z])	max=PART[i].q0[A_Z];
			if(min>PART[i].q0[A_Z])	min=PART[i].q0[A_Z];
			H=max-min+le;
		}
		double part_p=0;

		for(int i=0;i<num;i++)	
		{
			Z=PART[i].q0[A_Z];
			Y=PART[i].q0[A_Y];
			X=PART[i].q0[A_X];
			part_p=(Z/H)*2;
			HYPER[i].p[A_X]=mi*(-t*part_p*part_p*part_p*Y+b*(3*part_p*part_p-1));
			HYPER[i].p[A_Y]=t*mi*part_p*part_p*part_p*X;
			HYPER[i].p[A_Z]=0;
		}
	}
	
	//回転
	if(CON.get_model_number()==22)
	{
		for(int i=0;i<num;i++)
		{
			PART[i].p[A_X]=mi*0.4*(PART[i].q0[A_Z]-PART[i].q0[A_Y]);
			PART[i].p[A_Y]=mi*0.4*(PART[i].q0[A_X]-PART[i].q0[A_Z]);
			PART[i].p[A_Z]=mi*0.4*(PART[i].q0[A_Y]-PART[i].q0[A_X]);
		}
	}

	//重力降下
	if(CON.get_model_number()==23)
	{
		for(int i=0;i<num;i++)
		{		
			HYPER[i].p[A_Z]=-9.8*mi*Dt;
//			HYPER[i].p[A_Z]=0;
			HYPER[i].p[A_X]=0;
			HYPER[i].p[A_Y]=0;
		}
	}

	//角運動量計算
	for(int i=0;i<particle_number;i++)
	{
		HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
		HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
		HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
	}

	stringstream ss;
	ss<<"./Ai/Ai.csv";
	string fl=ss.str();
	ofstream ai(fl);

	stringstream ss2;
	ss2<<"./Ai/inverse_Ai.csv";
	string fl2=ss2.str();
	ofstream inai(fl2);

	ofstream aii("aiin.csv");
	ofstream n0("noij.csv");
	ofstream wi("wiin.csv");
	
	ofstream dg("initial_DgDq.csv");

	////近傍粒子の記憶とaiin,wiin,Aiの計算
	double dis=0;
	int N=0;

	for(int i=0;i<num;i++)
	{
		HYPER[i].N=0;
		for(int j=0;j<200;j++)	HYPER[i].NEI[j]=0;
	}

	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			HYPER1[i*num+j].wiin=0;
			for(int D=0;D<DIMENSION;D++)	HYPER1[i*num+j].aiin[D]=0;
		}
	}
	
	for(int i=0;i<num;i++)
	{
		N=0;
		for(int j=0;j<num;j++)
		{
			for(int D=0;D<DIMENSION;D++)
			{
				HYPER1[i*num+j].aiin[D]=PART[j].q0[D]-PART[i].q0[D];
				aii<<HYPER1[i*num+j].aiin[D]<<",";
			}
			aii<<",";
			dis=sqrt(HYPER1[i*num+j].aiin[A_X]*HYPER1[i*num+j].aiin[A_X]+HYPER1[i*num+j].aiin[A_Y]*HYPER1[i*num+j].aiin[A_Y]+HYPER1[i*num+j].aiin[A_Z]*HYPER1[i*num+j].aiin[A_Z]);
			if(dis<r && j!=i)
			{				
//				HYPER1[i*num+j].wiin=kernel4(r,dis);	//一時消去15/2/10
				HYPER[i].NEI[N]=j;
				N++;
			}
			else
			{
				HYPER1[i*num+j].wiin=0;
			}
			wi<<HYPER1[i*num+j].wiin<<",";			
		}
		aii<<endl<<endl;
		wi<<endl;
		HYPER[i].N=N;
	}
	aii.close();
	wi.close();
	
	//Aiの計算
	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				HYPER[i].Ai[D][D2]=0;
				for(int j=0;j<num;j++)	HYPER[i].Ai[D][D2]+=HYPER1[i*num+j].wiin*HYPER1[i*num+j].aiin[D]*HYPER1[i*num+j].aiin[D2];
			}
		}
	}

	//inverse_Ai,t_inverse_Aiの計算
	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				HYPER[i].inverse_Ai[D][D2]=0;
				HYPER[i].t_inverse_Ai[D][D2]=0;
			}
		}
	}

	double **p_Ai=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Ai[D]=new double [DIMENSION];
	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Ai[D][D2]=0;

	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Ai[D][D2]=HYPER[i].Ai[D][D2];
		inverse(p_Ai,DIMENSION);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].inverse_Ai[D][D2]=p_Ai[D][D2];
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)p_Ai[D][D2]=p_Ai[D2][D];
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].t_inverse_Ai[D][D2]=p_Ai[D][D2];
	}
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Ai[D];
	delete[]	p_Ai;

	////Fiの計算
	//初期化
	for(int i=0;i<num;i++)
	{
		HYPER[i].J=0;
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)
		{
			HYPER[i].differential_Fi[D][D2]=0;
			HYPER[i].t_inverse_Fi[D][D2]=0;
		}
	}

	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=0;

	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				HYPER[i].Fi[D][D2]=0;
				for(int D3=0;D3<DIMENSION;D3++)	HYPER[i].Fi[D][D2]+=HYPER[i].Ai[D][D3]*HYPER[i].inverse_Ai[D3][D2];
			}
		}
	}
	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=HYPER[i].Fi[D][D2];
		HYPER[i].J=calc_det(p_Fi,DIMENSION);
//		cout<<"HYPER["<<i<<"].J="<<HYPER[i].J<<endl;
		if(HYPER[i].J>=0)	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].differential_Fi[D][D2]=1/pow(HYPER[i].J,1.0/3.0)*HYPER[i].Fi[D][D2];
		else
		{
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].differential_Fi[D][D2]=-1/pow(-HYPER[i].J,1.0/3.0)*HYPER[i].Fi[D][D2];
		}			
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=p_Fi[D2][D];
		inverse(p_Fi,DIMENSION);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].t_inverse_Fi[D][D2]=p_Fi[D][D2];
	}
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;

	//n0ijの計算
	double p_n0ij[DIMENSION];
	double p_n0ij_2[DIMENSION];
	double p_n0ij_3[DIMENSION];

	for(int i=0;i<num;i++)	for(int j=0;j<num;j++)	for(int D=0;D<DIMENSION;D++)	HYPER1[i*num+j].n0ij[D]=0;

	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			if(j!=i)
			{
	//			cout<<"n0["<<i<<"]["<<j<<"]=";
				for(int D=0;D<DIMENSION;D++)
				{
					p_n0ij[D]=0;
					for(int D2=0;D2<DIMENSION;D2++)	p_n0ij[D]+=HYPER[j].t_inverse_Ai[D][D2]*HYPER1[i*num+j].aiin[D2];
				}
				for(int D=0;D<DIMENSION;D++)	HYPER1[i*num+j].n0ij[D]=V*HYPER1[j*num+i].wiin*p_n0ij[D];
	//			for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[i*num+j].n0ij[D]<<" ";
			}
			else
			{
//				cout<<"n0["<<i<<"]["<<j<<"]=";
				for(int D=0;D<DIMENSION;D++)
				{
					p_n0ij_2[D]=0;
					for(int k=0;k<num;k++)		p_n0ij_2[D]+=HYPER1[k*num+i].wiin*HYPER1[i*num+k].aiin[D];
				}
//				for(int D=0;D<DIMENSION;D++)	cout<<"p_n0_2["<<i<<"]["<<j<<"]["<<D<<"]="<<p_n0ij_2[D]<<endl;
				for(int D=0;D<DIMENSION;D++)
				{
					p_n0ij_3[D]=0;
					for(int D2=0;D2<DIMENSION;D2++)	p_n0ij_3[D]+=HYPER[i].t_inverse_Ai[D][D2]*p_n0ij_2[D2];
				}
/*				for(int D=0;D<DIMENSION;D++)
				{
					for(int D2=0;D2<DIMENSION;D2++)	cout<<"t_i_Ai["<<D<<"]["<<D2<<"]="<<HYPER[i].t_inverse_Ai[D][D2]<<" ";
					cout<<endl;
				}
				for(int D=0;D<DIMENSION;D++)	cout<<"p_n0["<<i<<"]["<<j<<"]["<<D<<"]_3="<<p_n0ij_3[D]<<endl;*/
				for(int D=0;D<DIMENSION;D++)	HYPER1[i*num+j].n0ij[D]=V*p_n0ij_3[D];
//				for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[i*num+j].n0ij[D]<<" ";
			}
//			cout<<endl;
		}
//		cout<<endl;
	}
//	cout<<endl;

	
	//DgDqの計算
	double p_DgDq[3];
	for(int j=0;j<num;j++)	for(int i=0;i<num;i++)	for(int D=0;D<DIMENSION;D++)	HYPER1[j*num+i].DgDq[D]=0;
	for(int j=0;j<num;j++)
	{
		for(int i=0;i<num;i++)
		{
			for(int D=0;D<DIMENSION;D++)
			{
				p_DgDq[D]=0;
				for(int D2=0;D2<DIMENSION;D2++)	p_DgDq[D]+=HYPER[j].t_inverse_Fi[D][D2]*HYPER1[i*num+j].n0ij[D2];
			}
			for(int D=0;D<DIMENSION;D++)	HYPER1[j*num+i].DgDq[D]=HYPER[j].J*p_DgDq[D];			
		}
	}

	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				inai<<HYPER[i].inverse_Ai[D][D2]<<",";		
				ai<<HYPER[i].Ai[D][D2]<<",";
			}
			ai<<endl;
			inai<<endl;
		}
		ai<<endl<<endl;
		inai<<endl<<endl;
		
		for(int j=0;j<num;j++)
		{
//			cout<<"n0ij["<<i<<"]["<<j<<"]=";
			for(int D=0;D<DIMENSION;D++)	n0<<HYPER1[i*num+j].n0ij[D]<<",";
			n0<<",";
	//		for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[i*num+j].n0ij[D]<<" ";
//			cout<<endl;
			dg<<HYPER1[i*num+j].DgDq[A_X]<<","<<HYPER1[i*num+j].DgDq[A_Y]<<","<<HYPER1[i*num+j].DgDq[A_Z]<<","<<",";
		}
		dg<<endl<<endl;
		n0<<endl<<endl;
//		cout<<endl;
	}
	
	ai.close();
	inai.close();
	n0.close();
	dg.close();
}


/////ニュートンラフソン法 
void newton_raphson(mpsconfig &CON,vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,int t)
{
	/////fx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照
	/////DfDx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照

	int calc_type=0;//ニュートラフソンの反復方法 0:偏微分項の逆行列をそのまま求める　1:線形方程式を利用

	//pn=2;//test,とりあえず2元でとけるかどうか確認 
	//////////////////　f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm

	int N=particle_number;
	calc_type=1;
	double *fx=new double [N];//関数値。
	double *DfDx=new double [N*N];//関数の偏微分値。
	double *XX=new double [N];//現在の解。	
	double *XX_old=new double [N];//1ステップ前の解。
	double ep=1e-4;//収束判定
	double E=1;//現在の誤差
/*	double start=0;
	double end=0;
	double newton_t=0;*/
	int count=0;//反復回数
	double d;
	double V=get_volume(&CON);

	for(int i=0;i<N;i++)
	{
		if(t==1)	HYPER[i].lambda=1;
		XX[i]=0;
		XX_old[i]=0;
		fx[i]=0;
		for(int j=0;j<N;j++)	DfDx[i*N+j]=0;
	}

	stringstream ss;
	ss<<"./Newton_raphson/newton_raphson"<<t<<".dat";
	string filename=ss.str();
	ofstream newton(filename);

	//	for(int i=0; i<N; i++) XX[i]=1;///初期値を与える。とりあえず1で
//	cout<<"NR法開始----";
//	start=clock();
	while(E>ep)
	{
		count++;

/*		stringstream ss1;
		ss1<<"./Newton_raphson/fx "<<"t"<<t<<" count"<<count<<".dat";
		filename=ss1.str();
		ofstream fs1(filename);*/

/*		stringstream ss2;
		ss2<<"./Newton_raphson/DfDx "<<"t"<<t<<" count"<<count<<".csv";
		filename=ss2.str();
		ofstream fDfDx(filename);*/

/*		stringstream ss3;
		ss3<<"./Newton_raphson/inverse_DfDx "<<"t"<<t<<" count"<<count<<".csv";
		filename=ss3.str();
		ofstream i_fDfDx(filename);*/

		stringstream ss4;
		ss4<<"./Newton_raphson/lambda "<<"t"<<t<<" count"<<count<<".dat";
		string filename2=ss4.str();
		ofstream lam(filename2);

//		if(count==1)	for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	for(int D=0;D<DIMENSION;D++)	HYPER1[i*N+j].newton_DgDq[D]=HYPER1[i*N+j].DgDq[D];

			calc_newton_q(CON,PART,HYPER,HYPER1,N,count,t);
			//calculation of DgDq
			double	p_newton_DgDq[3];
			for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	for(int D=0;D<DIMENSION;D++)	HYPER1[j*N+i].newton_DgDq[D]=0;
			for(int i=0;i<N;i++)
			{
				for(int j=0;j<N;j++)
				{
					for(int D=0;D<DIMENSION;D++)
					{
						p_newton_DgDq[D]=0;
						for(int D2=0;D2<DIMENSION;D2++)	p_newton_DgDq[D]+=HYPER[j].t_inverse_Fi[D][D2]*HYPER1[i*N+j].n0ij[D2];
					}
					for(int D=0;D<DIMENSION;D++)	HYPER1[j*N+i].newton_DgDq[D]=HYPER[j].J*p_newton_DgDq[D];
				}
			}

		calc_lambda(CON,PART,HYPER,HYPER1,N,t,count);

		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				//現在の偏微分値を求める
				DfDx[i*N+j]=HYPER1[i*N+j].DFDlambda;
//				fDfDx<<DfDx[i*N+j]<<",";
			}
//			fDfDx<<endl;
			//現在の関数値を求める
			fx[i]=V*(1-HYPER[i].J);
//			cout<<"fx["<<i<<"]="<<fx[i]<<endl;
//			fs1<<fx[i]<<endl;
		}

		for(int i=0; i<N; i++)
		{
			XX[i]=HYPER[i].lambda;
			XX_old[i]=XX[i];	//解を記憶
		}

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
			calc_inverse_matrix_for_NR(N,DfDx);
/*			for(int i=0;i<N;i++)
			{
				for(int j=0;j<N;j++)
				{
					cout<<"DfDx["<<i<<"]["<<j<<"]="<<DfDx[i*N+j]<<endl;
					i_fDfDx<<DfDx[i*N+j]<<",";
				}
				cout<<endl;
				i_fDfDx<<endl;
			}*/

//			if(count==1) cout<<DfDx[0]<<" "<<DfDx[1]<<" "<<DfDx[2]<<" "<<DfDx[3]<<endl;
			for(int i=0; i<N; i++) 
			{
				cout<<"XX_old["<<i<<"]-d=X["<<i<<"]	";
				cout<<XX_old[i]<<" - ";
				d=0; //変化量
				for(int j=0; j<N; j++)	d+=DfDx[i*N+j]*fx[j];
				XX[i]-=d;
				cout<<d<<" = ";
				cout<<XX[i]<<endl;
//				if(count==1) cout<<d<<endl;
			}
		}
		else if(calc_type==1)//逆行列を用いない、安定するはずだが、遅くなるはず
		{
			gauss(DfDx,fx,N);
			if(count%100==1)
			{
				stringstream ss5;
				ss5<<"./Newton_raphson/transition "<<"t"<<t<<" count"<<count<<".dat";
				string filename3=ss5.str();
				ofstream tr(filename3);

				tr<<"Transition "<<"t"<<t<<"count"<<count<<endl;
				for(int i=0; i<N; i++) 
				{
					tr<<"XX_old["<<i<<"]-r=X["<<i<<"]	";
					tr<<XX_old[i]<<" - ";
					XX[i]-=fx[i];
					tr<<fx[i]<<" = ";
					tr<<XX[i]<<endl;
	//				if(count==1) cout<<d<<endl;
				}
				tr.close();
			}
			else
			{
				for(int i=0;i<N;i++)	XX[i]-=fx[i];
			}
		}

//		for(int i=0;i<N;i++)	cout<<"XX["<<i<<"]="<<XX[i]<<endl;
		for(int i=0;i<N;i++)	HYPER[i].lambda=XX[i];
//		for(int i=0;i<N;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
		lam<<"lambda"<<"t"<<t<<"count"<<count<<endl;
		for(int i=0;i<N;i++)	lam<<i<<"	"<<HYPER[i].lambda<<endl;

		//誤差の評価	
		E=0;
		for(int i=0; i<N; i++) E+=fabs(XX[i]-XX_old[i]);
		double sum;
		sum=0;
		for(int i=0;i<N;i++)	sum+=fabs(XX[i]);
		E/=sum;

		newton<<"反復回数"<<count<<"E"<<" "<<E<<endl;
		if(count%100==1)	cout<<"反復回数="<<count<<" E="<<E<<endl;

		lam.close();
//		fs1.close();
//		fDfDx.close();
//		i_fDfDx.close();
		if(count>CON.get_nr())	break;	//15/2/8

	}
//	end=clock();


//	newton_t=(end-start)/CLOCKS_PER_SEC;
	for(int i=0;i<N;i++)	newton<<"反復完了"<<endl;
//	cout<<"完了 反復回数="<<count<<"X=("<<XX[0]<<","<<XX[1]<<")"<<endl;

	delete[]	fx;
	delete[]	DfDx;
	delete[]	XX;
	delete[]	XX_old;

}

void calc_lambda(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,int t,int count)
{
	double Dt=CON.get_dt()*CON.get_interval();
	double V=get_volume(&CON);
	double mk=V*CON.get_hyper_density();
	double DFDlambda;

/*	stringstream	s2;
	s2<<"./Newton_raphson/newton_DgDq "<<"t"<<t<<" count"<<count<<".csv";
	string filename=s2.str();
	ofstream	fs2(filename);*/

	for(int i=0;i<particle_number;i++)
	{
		for(int j=0;j<particle_number;j++)
		{
			HYPER1[j*particle_number+i].DFDlambda=0;
			DFDlambda=0;
			for(int k=0;k<particle_number;k++)
			{		
/*				for(int D=0;D<DIMENSION;D++)	fs2<<HYPER1[i*particle_number+j].newton_DgDq[D]<<",";
				fs2<<",";*/
				for(int D=0;D<DIMENSION;D++)	DFDlambda+=HYPER1[i*particle_number+k].newton_DgDq[D]*HYPER1[j*particle_number+k].DgDq[D];
			}
			HYPER1[i*particle_number+j].DFDlambda=-Dt*Dt/2/mk*DFDlambda;
//			cout<<"DFDlmabda["<<i<<"]["<<j<<"]="<<HYPER1[i*particle_number+j].DFDlambda<<endl;
		}
//		fs2<<endl;
	}
//	fs2.close();

/*	stringstream	s;
	s<<"./Newton_raphson/DFDlambda "<<"t"<<t<<" count"<<count<<".csv";
	filename=s.str();
	ofstream	fs(filename);
			
	for(int i=0;i<particle_number;i++)
	{
		for(int j=0;j<particle_number;j++)	fs<<HYPER1[i*particle_number+j].DFDlambda<<",";
		fs<<endl;
	}
	fs.close();*/
	
}


void calc_newton_q(mpsconfig &CON,vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,int count,int t)
{
	double Dt=CON.get_dt()*CON.get_interval();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();

	//half_pの計算
	double p_half_p[3][3];
	double p_half_p2[3];
	double p_half_p3[3];
	double n_half_p[3];

	for(int i=0;i<particle_number;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[i].n_q[D]=0;
	for(int D=0;D<DIMENSION;D++)	n_half_p[D]=0;

/*	stringstream s5;
	s5<<"./Newton_raphson/partial_half_p "<<"t"<<t<<" count"<<count<<".dat";
	filename=s5.str();
	ofstream fss5(filename);

	stringstream s6;
	s6<<"./Newton_raphson/half_p "<<"t"<<t<<" count"<<count<<".dat";
	filename=s6.str();
	ofstream fss6(filename);*/
	
	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	p_half_p3[D]=0;
		for(int j=0;j<particle_number;j++)
		{		
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_half_p[D][D2]=0;
			if(HYPER[j].stress0==ON)	for(int D=0;D<DIMENSION;D++)	p_half_p[D][D]=-HYPER[j].lambda;
			else
			{
				for(int D=0;D<DIMENSION;D++)
				{
					for(int D2=0;D2<DIMENSION;D2++)
					{
						p_half_p[D][D2]=HYPER[j].stress[D][D2];
						if(D2==D)	p_half_p[D2][D2]-=HYPER[j].lambda;
					}
				}
			}
			for(int D=0;D<DIMENSION;D++)
			{	
				p_half_p2[D]=0;
				for(int D2=0;D2<DIMENSION;D2++)	p_half_p2[D]+=p_half_p[D][D2]*HYPER1[j*particle_number+i].DgDq[D2];
			}
			for(int D=0;D<DIMENSION;D++)	p_half_p3[D]+=p_half_p2[D];
		}//jに関するfor文の終わり			
//		cout<<"partial_half["<<i<<"]="<<Dt/2*p_half_p3[A_X]<<" "<<Dt/2*p_half_p3[A_Y]<<" "<<Dt/2*p_half_p3[A_Z]<<endl;


		////位置座標の更新
		for(int D=0;D<DIMENSION;D++)	HYPER[i].n_q[D]=PART[i].r[D]+Dt*n_half_p[D]/mi;
		if(CON.get_ttest()==ON&&HYPER[i].highest==ON)	HYPER[i].n_q[A_Z]=PART[i].q0[A_Z];	//引っ張り試験解析用15/2/8
/*		fss5<<"partial_half_p["<<i<<"]=";
		for(int D=0;D<DIMENSION;D++)fss5<<p_half_p3[D]<<"	";
		fss5<<endl;
		fss6<<"HYPER["<<i<<"].n_half_p="<<n_half_p[A_X]<<"	"<<n_half_p[A_Y]<<"	"<<n_half_p[A_Z]<<endl;*/
	}
/*	fss5.close();
	fss6.close();*/


	if(count%100==1)
	{
		stringstream s4;
		s4<<"./Newton_raphson/position "<<"t"<<t<<" count"<<count<<".dat";
		string filename=s4.str();
		ofstream fss4(filename);

		fss4<<"PART.n_q"<<"t"<<t<<"count"<<count<<endl;
		for(int i=0;i<particle_number;i++)	fss4<<i<<"	"<<HYPER[i].n_q[A_X]<<"	"<<HYPER[i].n_q[A_Y]<<"	"<<HYPER[i].n_q[A_Z]<<endl;		
		fss4.close();
	}

	//Fの計算
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double [DIMENSION];
	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=0;

/*	stringstream s;
	s<<"./Fi/Fi "<<"t"<<t<<" count"<<count<<".csv";
	filename=s.str();
	ofstream fss(filename);


	stringstream s3;
	s3<<"./Fi/t_inverse_Fi "<<"t"<<t<<" count"<<count<<".csv";
	filename=s3.str();
	ofstream fss3(filename);*/

	////qiinとFiの更新
	//初期化
	double fi[3][3];
	double fi2[3][3];
	double fi0[3][3];

	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi0[D][D2]=0;

	for(int i=0;i<particle_number;i++)
	{
		HYPER[i].J=0;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				HYPER[i].t_inverse_Fi[D][D2]=0;
				HYPER[i].Fi[D][D2]=0;			
			}
		}
	}

	//計算
	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]=0;
		for(int in=0;in<HYPER[i].N;in++)
		{
/*			cout<<"w["<<i<<"]["<<HYPER[i].NEI[in]<<"]="<<HYPER1[i*particle_number+HYPER[i].NEI[in]].wiin<<endl;
			cout<<"q=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER[HYPER[i].NEI[in]].n_q[D]-HYPER[i].n_q[D]<<" ";
			cout<<endl;
			cout<<"a["<<i<<"]["<<HYPER[i].NEI[in]<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[i*particle_number+HYPER[i].NEI[in]].aiin[D]<<" ";
			cout<<endl;*/
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi0[D][D2]=HYPER1[i*particle_number+HYPER[i].NEI[in]].wiin*(HYPER[HYPER[i].NEI[in]].n_q[D]-HYPER[i].n_q[D])*HYPER1[i*particle_number+HYPER[i].NEI[in]].aiin[D2];
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]+=fi0[D][D2];
				/*			cout<<"qiin=";
				for(int D=0;D<DIMENSION;D++)	cout<<PART[HYPER[i].NEI[in]].r[D]-PART[i].r[D]<<" ";
				cout<<endl;
				cout<<"aiin=";
				for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[i*particle_number+HYPER[i].NEI[in]].aiin[D]<<" ";
				cout<<endl;
				cout<<endl;*/
		}
//				cout<<endl;	

/*		cout<<"fi["<<i<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				cout<<fi[D][D2]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;

		cout<<"inverse_A["<<i<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				cout<<HYPER[i].inverse_Ai[D][D2]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;*/

		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				fi2[D][D2]=0;
				for(int D3=0;D3<DIMENSION;D3++)	fi2[D][D2]+=fi[D][D3]*HYPER[i].inverse_Ai[D3][D2];
			}
		}
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].Fi[D][D2]=fi2[D][D2];
	}

	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=HYPER[i].Fi[D][D2];
		HYPER[i].J=calc_det(p_Fi,DIMENSION);
/*		cout<<"J["<<i<<"]="<<HYPER[i].J<<endl;

		cout<<"Fi["<<i<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)	cout<<p_Fi[D*DIMENSION+D2]<<" ";
			cout<<endl;
		}*/

		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=p_Fi[D2][D];
		inverse(p_Fi,DIMENSION);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].t_inverse_Fi[D][D2]=p_Fi[D][D2];		
/*	cout<<"t_inverse_Fi["<<i<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)	cout<<p_Fi[D*DIMENSION+D2]<<" ";
			cout<<endl;
		}
		cout<<endl;*/
	}

/*	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				fss<<HYPER[i].Fi[D][D2]<<",";
				fss3<<HYPER[i].t_inverse_Fi[D][D2]<<",";
			}
			fss<<endl;
			fss3<<endl;
		}
		fss<<endl<<endl;
		fss3<<endl<<endl;
	}
	fss.close();
	fss3.close();*/
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;
}

void calc_half_p(mpsconfig &CON,vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,bool repetation,int t)
{
	double Dt=CON.get_dt()*CON.get_interval();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();

	double p_half_p[3][3];
	double p_half_p2[3];
	double p_half_p3[3];
//	int flag_wall=0;	//運動量の与え方を壁を通り抜けた粒子のみから全体に改良　15/2/5
	int wall=0;//試験作12/2/5

	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	p_half_p3[D]=0;
		for(int j=0;j<particle_number;j++)
		{				
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_half_p[D][D2]=0;
			if(HYPER[j].stress0==ON)	for(int D=0;D<DIMENSION;D++)	p_half_p[D][D]=-HYPER[j].lambda;
			else
			{
				for(int D=0;D<DIMENSION;D++)
				{
					for(int D2=0;D2<DIMENSION;D2++)
					{
						p_half_p[D2][D2]=HYPER[j].stress[D2][D2];
						if(D2==D)	p_half_p[D2][D2]-=HYPER[j].lambda;
					}
				}
			}			
			for(int D=0;D<DIMENSION;D++)
			{	
				p_half_p2[D]=0;
				for(int D2=0;D2<DIMENSION;D2++)	p_half_p2[D]+=p_half_p[D][D2]*HYPER1[j*particle_number+i].DgDq[D2];
			}
			for(int D=0;D<DIMENSION;D++)	p_half_p3[D]+=p_half_p2[D];
		}//jに関するfor文の終わり			
//		cout<<"partial_half["<<i<<"]="<<Dt/2*p_half_p3[A_X]<<" "<<Dt/2*p_half_p3[A_Y]<<" "<<Dt/2*p_half_p3[A_Z]<<endl;

		if(repetation==0)
		{
			for(int D=0;D<DIMENSION;D++)	HYPER[i].half_p[D]=HYPER[i].p[D]+Dt/2*p_half_p3[D];
			////位置座標の更新
			for(int D=0;D<DIMENSION;D++)	PART[i].r[D]+=Dt*HYPER[i].half_p[D]/mi;
			if(CON.get_ttest()==ON&&HYPER[i].highest==ON)	PART[i].r[A_Z]=PART[i].q0[A_Z];	//引っ張り試験解析用15/2/8
		}
		else
		{
			for(int D=0;D<DIMENSION;D++)
			{
				if(CON.get_flag_vis()==OFF&&CON.get_FEM_flag()==OFF)	HYPER[i].p[D]=0;
				HYPER[i].ang_p[D]=0;
			}

			for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]=HYPER[i].half_p[D]+Dt/2*p_half_p3[D];

			HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
			HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
			HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];

			if((PART[i].r[A_Z]<CON.get_r_z_wall()*le)&&(CON.get_flag_wall()==ON))//試験作15/2/5
			{
				HYPER[wall].Nw=i;
				wall++;
			}	
//			for(int D=0;D<DIMENSION;D++)cout<<"p["<<i<<"]="<<HYPER[i].half_p[D]<<"+"<<Dt<<"/2*"<<p_half_p3[D]<<" = "<<HYPER[i].p[D]<<endl;
		}
	}//iに関するfor文の終わり
	
	if(repetation!=0&&CON.get_flag_wall()==ON&&wall>0)	//試験的に壁モデルを作成15/2/5
	{
		for(int i=0;i<particle_number;i++)	
		{
			int count=0;
			for(int j=0;j<wall;j++)	
			{
				int k=HYPER[j].Nw;
				double dr_x=PART[i].r[A_X]-PART[k].r[A_X];
				double dr_y=PART[i].r[A_Y]-PART[k].r[A_Y];
				if(dr_x<DBL_EPSILON&&dr_y<DBL_EPSILON)	count++;
			}
			if(count>0)	HYPER[i].p[A_Z]-=HYPER[i].p[A_Z];
		}
	}

//	if(flag_wall>0)	for(int i=0;i<particle_number;i++)	HYPER[i].p[A_Z]-=HYPER[i].p[A_Z];	//運動量の与え方を壁を通り抜けた粒子のみから全体に改良　15/2/5
/*	if(repetation==0)
	{
		for(int i=0;i<particle_number;i++)
		{
			cout<<"half_p["<<i<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER[i].half_p[D]<<" ";
			cout<<endl;
		}
		cout<<endl;

		for(int i=0;i<particle_number;i++)
		{
			cout<<"PART["<<i<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<PART[i].r[D]<<" ";
			cout<<endl;
		}
		cout<<endl;
	}*/
}

void calc_F(vector<mpselastic>&PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number,int t)
{
	////qiinとFiの更新
	//初期化
	double fi[3][3];
	double fi2[3][3];
	double fi0[3][3];

	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi0[D][D2]=0;

	for(int i=0;i<particle_number;i++)
	{
		HYPER[i].J=0;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				HYPER[i].t_inverse_Fi[D][D2]=0;
				HYPER[i].differential_Fi[D][D2]=0;
				HYPER[i].Fi[D][D2]=0;			
			}
		}
	}

	//計算
	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]=0;
		for(int in=0;in<HYPER[i].N;in++)
		{
/*			cout<<"w["<<i<<"]["<<HYPER[i].NEI[in]<<"]="<<HYPER1[i*particle_number+HYPER[i].NEI[in]].wiin<<endl;
			cout<<"q=";
			for(int D=0;D<DIMENSION;D++)	cout<<PART[HYPER[i].NEI[in]].r[D]-PART[i].r[D]<<" ";
			cout<<endl;
			cout<<"a["<<i<<"]["<<HYPER[i].NEI[in]<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[i*particle_number+HYPER[i].NEI[in]].aiin[D]<<" ";
			cout<<endl;*/
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi0[D][D2]=HYPER1[i*particle_number+HYPER[i].NEI[in]].wiin*(PART[HYPER[i].NEI[in]].r[D]-PART[i].r[D])*HYPER1[i*particle_number+HYPER[i].NEI[in]].aiin[D2];
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]+=fi0[D][D2];
				/*			cout<<"qiin=";
				for(int D=0;D<DIMENSION;D++)	cout<<PART[HYPER[i].NEI[in]].r[D]-PART[i].r[D]<<" ";
				cout<<endl;
				cout<<"aiin=";
				for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[i*particle_number+HYPER[i].NEI[in]].aiin[D]<<" ";
				cout<<endl;
				cout<<endl;*/
		}
//				cout<<endl;	
/*		cout<<"fi["<<i<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				cout<<fi[D][D2]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;

		cout<<"inverse_A["<<i<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				cout<<HYPER[i].inverse_Ai[D][D2]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;*/

		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				fi2[D][D2]=0;
				for(int D3=0;D3<DIMENSION;D3++)	fi2[D][D2]+=fi[D][D3]*HYPER[i].inverse_Ai[D3][D2];
			}
		}
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].Fi[D][D2]=fi2[D][D2];
	}

	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=0;

	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=HYPER[i].Fi[D][D2];
		HYPER[i].J=calc_det(p_Fi,DIMENSION);
/*		cout<<"Fi["<<i<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)	cout<<p_Fi[D*DIMENSION+D2]<<" ";
			cout<<endl;
		}*/
		if(HYPER[i].J>=0)	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].differential_Fi[D][D2]=1/pow(HYPER[i].J,1.0/3.0)*p_Fi[D][D2];
		else
		{
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].differential_Fi[D][D2]=-1/pow(-HYPER[i].J,1.0/3.0)*p_Fi[D][D2];
		}
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=p_Fi[D2][D];
				
		inverse(p_Fi,DIMENSION);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].t_inverse_Fi[D][D2]=p_Fi[D][D2];		
/*		cout<<"t_inverse_Fi["<<i<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)	cout<<p_Fi[D*DIMENSION+D2]<<" ";
			cout<<endl;
		}
		cout<<endl;*/
	}
	for(int D=0;D<DIMENSION;D++)	delete[]	p_Fi[D];
	delete[]	p_Fi;

	stringstream s;
	s<<"./Fi/Fi "<<"t"<<t<<".csv";
	string filename=s.str();
	ofstream fss(filename);

	stringstream s2;
	s2<<"./Fi/differential_Fi "<<"t"<<t<<".csv";
	string filename2=s2.str();
	ofstream fss2(filename2);

	stringstream s3;
	s3<<"./Fi/t_inverse_Fi "<<"t"<<t<<".csv";
	string filename3=s3.str();
	ofstream fss3(filename3);

	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				fss<<HYPER[i].Fi[D][D2]<<",";
				fss2<<HYPER[i].differential_Fi[D][D2]<<",";
				fss3<<HYPER[i].t_inverse_Fi[D][D2]<<",";
			}
			fss<<endl;
			fss2<<endl;
			fss3<<endl;
		}
		fss<<endl<<endl;
		fss2<<endl<<endl;
		fss3<<endl<<endl;
	}
	fss.close();
	fss2.close();
	fss3.close();
	
//	for(int i=0;i<particle_number;i++)	cout<<"J["<<i<<"]="<<HYPER[i].J<<endl;

}

void calc_stress(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>HYPER1,int particle_number)
{
	double c10=CON.get_c10();
	double c01=CON.get_c01();

	for(int i=0;i<particle_number;i++)
	{
		HYPER[i].stress0=OFF;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				HYPER[i].stress[D][D2]=0;
				HYPER[i].t_differential_Fi[D][D2]=0;
			}
		}
	}

	double b[3][3];
	double bb[3][3];
	double trace_b,trace_bb;
	
	for(int j=0;j<particle_number;j++)	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[j].t_differential_Fi[D][D2]=HYPER[j].differential_Fi[D2][D];

	for(int j=0;j<particle_number;j++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				b[D][D2]=0;				
				for(int D3=0;D3<DIMENSION;D3++)	b[D][D2]+=HYPER[j].differential_Fi[D][D3]*HYPER[j].t_differential_Fi[D3][D2];
			}
		}
/*		cout<<"b["<<j<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				cout<<b[D][D2]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
		cout<<"t_differential_F["<<j<<"]="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				cout<<HYPER[j].differential_Fi[D2][D]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
		cout<<"differential_F="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				cout<<HYPER[j].differential_Fi[D][D2]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;*/
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				bb[D][D2]=0;				
				for(int D3=0;D3<DIMENSION;D3++)	bb[D][D2]+=b[D][D3]*b[D3][D2];
			}
		}
/*		cout<<"bb="<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				cout<<bb[D][D2]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;*/

		trace_b=0;
		trace_bb=0;
		for(int D=0;D<DIMENSION;D++)
		{
			trace_b+=b[D][D];
			trace_bb+=bb[D][D];
		}

//		cout<<"trace_b, trace_bb="<<trace_b<<" "<<trace_bb<<endl;
		
		double e=DBL_EPSILON;
		int c_b=0,c_b2=0;
		int c_bb=0,c_bb2=0;
		double b2[3][3];
		double bb2[3][3];

		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				b2[D][D2]=0;
				bb2[D][D2]=0;
			}
		}


		for(int D=0;D<DIMENSION;D++)
		{
			if(fabs(b[D][D]-1.0/3.0*trace_b)<e)
			{
				b2[D][D]=0.0;
				c_b2++;
			}
			else
			{
				b2[D][D]=b[D][D]-1.0/3.0*trace_b;
			}

			if(fabs(bb[D][D]-1.0/3.0*trace_bb)<e)
			{
				bb2[D][D]=0.0;
				c_bb2++;
			}
			else
			{
				bb2[D][D]=bb[D][D]-1.0/3.0*trace_bb;
			}
		}

		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				if(D2!=D)
				{
					if(fabs(b[D][D2])<e)	
					{
						b2[D][D2]=0.0;
						c_b++;
					}
					else
					{
						b2[D][D2]=b[D][D2];
					}
					if(fabs(bb[D][D2])<e)
					{
						bb2[D][D2]=0.0;
						c_bb++;
					}
					else
					{
						bb2[D][D2]=bb[D][D2];
					}
				}
			}
		}

		if((c_b+c_b2+c_bb+c_bb2)==18)	HYPER[j].stress0=ON;
		else if(((c_b+c_b2)==9)&&((c_bb+c_bb2)!=9))	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[j].stress[D][D2]=-c01*2/HYPER[j].J*bb2[D][D2];
		else if(((c_b+c_b2)!=9)&&((c_bb+c_bb2)==9))	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[j].stress[D][D2]=2/HYPER[j].J*(c10+c01*trace_b)*b2[D][D2];
		else
		{
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[j].stress[D][D2]=2/HYPER[j].J*((c10+c01*trace_b)*b2[D][D2]-c01*bb2[D][D2]);
		}
	}
}

void calc_differential_p(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number)
{
	double Dt=CON.get_dt()*CON.get_interval();

	double p_differential_p[3];
	double p_differential_p2[3];
	for(int i=0;i<particle_number;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[i].differential_p[D]=0;

	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	p_differential_p2[D]=0;
		for(int j=0;j<particle_number;j++)
		{							
	//			cout<<"partial_differential_p["<<j<<"]=0"<<endl;
			if(HYPER[j].stress0==OFF)
			{
				for(int D=0;D<DIMENSION;D++)
				{	
					p_differential_p[D]=0;
					for(int D2=0;D2<DIMENSION;D2++)	p_differential_p[D]+=HYPER[j].stress[D][D2]*HYPER1[j*particle_number+i].DgDq[D2];
//					cout<<p_differential_p[D]<<"	";
				}
//				cout<<endl;
/*				cout<<"DgDq";
				for(int D=0;D<DIMENSION;D++)	cout<<HYPER1[j*particle_number+i].DgDq[D];
				cout<<endl;

				cout<<"stress"<<endl;
				for(int D=0;D<DIMENSION;D++)
				{
					for(int D2=0;D2<DIMENSION;D2++)	cout<<HYPER[j].stress[D][D2]<<"	";
					cout<<endl;
				}
				cout<<endl;

				cout<<"partial_differential_p["<<j<<"]=";
				for(int D=0;D<DIMENSION;D++)cout<<partial_differential_p[D]<<" ";
				cout<<endl;*/
				for(int D=0;D<DIMENSION;D++)	p_differential_p2[D]+=p_differential_p[D];				
			}
/*			for(int D=0;D<DIMENSION;D++)	cout<<Dt/2*p_differential_p2[D]<<" ";
			cout<<endl;*/
		}		
		for(int D=0;D<DIMENSION;D++)	HYPER[i].differential_p[D]=HYPER[i].half_p[D]+Dt/2*p_differential_p2[D];
	}
}

void renew_lambda(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int particle_number)
{
	for(int i=0;i<particle_number;i++)	HYPER[i].lambda=0;

	double le=CON.get_distancebp();
	double Dt=CON.get_dt()*CON.get_interval();
	double V=get_volume(&CON);
	double mk=V*CON.get_hyper_density();

	double N_left;
	double N_right;
	double N_left2;
	double N_right2;
	double *N_Left=new double[particle_number*particle_number];
	double *N_Right=new double[particle_number];

	for(int i=0;i<particle_number;i++)
	{
		N_Right[i]=0;
		for(int j=0;j<particle_number;j++)	N_Left[j*particle_number+i]=0;
	}

	for(int i=0;i<particle_number;i++)
	{
		for(int j=0;j<particle_number;j++)
		{
			N_left2=0;
			for(int k=0;k<particle_number;k++)
			{
				N_left=0;
				for(int D=0;D<DIMENSION;D++)	N_left+=HYPER1[i*particle_number+k].DgDq[D]*HYPER1[j*particle_number+k].DgDq[D];
				N_left2+=N_left;
			}
			N_Left[i*particle_number+j]=Dt/2/mk*N_left2;
//			cout<<"N_Left["<<i<<"]["<<j<<"]= "<<Dt<<"/2/"<<mk<<"*"<<N_left2<<" = "<<N_Left[i*particle_number+j]<<endl;
		}//jに関するfor文の終わり
	}//iに関するfor文の終わり

	for(int i=0;i<particle_number;i++)
	{
		N_right2=0;
		for(int k=0;k<particle_number;k++)
		{
			N_right=0;
			for(int D=0;D<DIMENSION;D++)	N_right+=HYPER[k].differential_p[D]*HYPER1[i*particle_number+k].DgDq[D];
			N_right2+=N_right;
		}
		N_Right[i]=1/mk*N_right2;	
//		cout<<"N_Right["<<i<<"]= "<<"1/"<<mk<<"*"<<N_right2<<" = "<<N_Right[i]<<endl;
	}
/*
	for(int i=0;i<particle_number;i++)	cout<<"N_Right["<<i<<"]="<<N_Right[i]<<endl;
	cout<<endl;
	for(int i=0;i<particle_number;i++)
	{
		for(int j=0;j<particle_number;j++)	cout<<"N_Left["<<i<<"]["<<j<<"]="<<N_Left[i*particle_number+j]<<endl;
		cout<<endl;
	}*/

	//lambdaを求める
	gauss(N_Left,N_Right,particle_number);

	for(int i=0;i<particle_number;i++)	HYPER[i].lambda=N_Right[i];
	cout<<"lamda is renewed."<<endl;
//	for(int i=0;i<particle_number;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;

	delete [] N_Left;
	delete [] N_Right;
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
		if(max==0.0)
		{
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
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> &PART,vector<hyperelastic>&HYPER,int particle_number)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	
	double TIME=CON.get_step()*CON.get_dt();
	double le=CON.get_distancebp();
	int STEP=CON.get_step()/CON.get_interval()+1;		//出力する総ステップ数
	int n=particle_number;

		if(t==1) 
	{
		ofstream fp("momentum.inp");			
		fp<<CON.get_step()/CON.get_interval()+1<<endl;//総ステップ数
		fp<<"data_geom"<<endl;
		fp.close();
	}

	//mainファイル書き込み
	ofstream fp("momentum.inp",ios :: app);
	fp<<"step"<<t/CON.get_interval()+1<<" TIME="<<TIME<<endl;

	//fp<<step<<endl;
	
	//fp<<"data_geom"<<endl;
	//fp<<"step1"<<endl;
	//fp<<"step"<<t/CON->get_interval()+1<<" TIME="<<TIME<<endl;
	fp<<particle_number<<" "<<particle_number<<endl;	//節点数と要素数出力
	
	//節点番号とその座標の出力 
	for(int i=0;i<particle_number;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=0;i<particle_number;i++)
	{
		fp<<i<<"  0 pt "<<i<<endl;
	}

	//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
	fp<<"8 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
	fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	fp<<"p_x,"<<endl;
	fp<<"p_y,"<<endl;
	fp<<"p_z,"<<endl;
	fp<<"lambda,"<<endl;
	fp<<"J,"<<endl;
	fp<<"ap_x,"<<endl;
	fp<<"ap_y,"<<endl;
	fp<<"ap_z,"<<endl;

	//fp<<"P,N/m^2"<<endl;
	//fp<<"value1,??"<<endl;

	//各節点の情報値入力
	for(int i=0;i<particle_number;i++)
	{
		fp<<i<<" "<<HYPER[i].p[A_X]<<" "<<HYPER[i].p[A_Y]<<" "<<HYPER[i].p[A_Z]<<" "<<HYPER[i].lambda<<" "<<HYPER[i].J<<" "<<HYPER[i].ang_p[A_X]<<" "<<HYPER[i].ang_p[A_Y]<<" "<<HYPER[i].ang_p[A_Z]<<" "<<endl;
		//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
	}

	cout<<"OK"<<endl;
	fp.close();

}

















