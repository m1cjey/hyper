#include "stdafx.h"

void calc_vis_f(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>&HYPER1,int hyper_number,int t)
{
	int num=hyper_number;
	double d_nator[3];
	for(int D=0;D<DIMENSION;D++)	d_nator[D]=0;

	double n_nator=0;
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double vis=CON.get_h_viscousity();

	if(t==1)	system("mkdir Viscousity");

	stringstream s;
	s<<"./Viscousity/vis_force"<<t<<".dat";
	string filename=s.str();
	ofstream fs(filename);

	if(CON.get_FEM_flag()==OFF&&t==1)	for(int i=0;i<num;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]=0;	//t==1‚ð‘I‘ð€‚É‰Á‚¦‚½15/2/9

	calc_spl_f(CON,PART,HYPER1,hyper_number);

	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)	HYPER[i].vis_force[D]=0;
		for(int j=0;j<num;j++)
		{
			for(int D=0;D<DIMENSION;D++)	d_nator[D]=6/mi*vis*(HYPER[j].p[D]-HYPER[i].p[D])*HYPER1[i*num+j].spl_f;	//Ž®’ù³‚·‚é‚à‚Ð‚Ç‚­‚È‚Á‚½15/2/3
			for(int D=0;D<DIMENSION;D++)	HYPER[i].vis_force[D]+=d_nator[D];
		}

		n_nator=0;
		for(int j=0;j<num;j++)
		{
			n_nator+=((PART[j].r[A_X]-PART[i].r[A_X])*(PART[j].r[A_X]-PART[i].r[A_X])+(PART[j].r[A_Y]-PART[i].r[A_Y])*(PART[j].r[A_Y]-PART[i].r[A_Y])+(PART[j].r[A_Z]-PART[i].r[A_Z])*(PART[j].r[A_Z]-PART[i].r[A_Z]))*HYPER1[i*num+j].spl_f;
		}
		for(int D=0;D<DIMENSION;D++)	HYPER[i].vis_force[D]/=n_nator;
	}

	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			HYPER[i].p[D]+=HYPER[i].vis_force[D];
			fs<<HYPER[i].vis_force[D]<<" ";
		}
		fs<<endl;
	}
	fs<<endl;
	fs.close();
}

void calc_spl_f(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic2>&HYPER1,int hyper_number)
{
	int num=hyper_number;

	double h=CON.get_h_dis();
	double r=0,s=0,theta=0;

	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			theta=0;	//ŠÖ”C³15/2/8
			HYPER1[i*num+j].spl_f=0;
			r=sqrt(pow((PART[j].r[A_X]-PART[i].r[A_X]),2.0)+pow((PART[j].r[A_Y]-PART[i].r[A_Y]),2.0)+pow((PART[j].r[A_Z]-PART[i].r[A_Z]),2.0));
			s=r/(h/2);
			if(0<=s&&s<1)	theta=1/PI*(1-3/2*pow(s,2.0)+3/4*pow(s,3.0));
			if(1<=s&&s<2)	theta=1/PI/4*pow((2-s),3.0);
			HYPER1[i*num+j].spl_f=1/pow(h/2,3)*theta;	//ŠÖ”C³15/2/8
		}
	}
}
