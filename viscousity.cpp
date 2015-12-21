#include "stdafx.h"

void calc_vis_f(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic>&HYPER,vector<hyperelastic2>HYPER1,int rigid_number,int t)
{
	cout<<"”S«€ŒvŽZ";
	double v=CON.get_vis();
	double r=CON.get_h_dis();
	int h_num=HYPER.size();
	int r_num=rigid_number;
	int d=3;
	double mi=CON.get_hyper_density()*get_volume(&CON);
	double lambda=calclambda(CON);
	ofstream fv("viscousity.csv",ios::app);
	if(t==1)
	{
		fv<<",";
		for(int i=0;i<h_num-r_num;i++)
		{
			fv<<i<<","<<","<<","<<","<<","<<","<<",";
		}
		fv<<endl;
	}
	fv<<t<<",";

	for(int i=0;i<h_num-r_num;i++)
	{
		int Nh=0;
		double p_vis[3]={0,0,0};
		double n0=HYPER[i].pnd;
		double pnd=0;
		for(int j=0;j<h_num-r_num;j++)
		{
			if(j!=i)
			{
				double qij[3]={PART[j].r[A_X]-PART[i].r[A_X], PART[j].r[A_Y]-PART[i].r[A_Y], PART[j].r[A_Z]-PART[i].r[A_Z]};
				double dis=sqrt(qij[A_X]*qij[A_X]+qij[A_Y]*qij[A_Y]+qij[A_Z]*qij[A_Z]);
				if(dis<r)
				{
					double w=kernel4(r,dis);
					p_vis[A_X]+=(HYPER[j].p[A_X]-HYPER[i].p[A_X])*w;
					p_vis[A_Y]+=(HYPER[j].p[A_Y]-HYPER[i].p[A_Y])*w;
					p_vis[A_Z]+=(HYPER[j].p[A_Z]-HYPER[i].p[A_Z])*w;
				}
			}
		}
		HYPER[i].vis_force[A_X]=1/mi*v*2*d/lambda/n0*p_vis[A_X];
		HYPER[i].vis_force[A_Y]=1/mi*v*2*d/lambda/n0*p_vis[A_Y];
		HYPER[i].vis_force[A_Z]=1/mi*v*2*d/lambda/n0*p_vis[A_Z];
		/*HYPER[i].vis_force[A_X]=2*v*HYPER[i].pnd*p_vis_d[A_X]/p_vis_n;
		HYPER[i].vis_force[A_Y]=2*v*HYPER[i].pnd*p_vis_d[A_Y]/p_vis_n;
		HYPER[i].vis_force[A_Z]=2*v*HYPER[i].pnd*p_vis_d[A_Z]/p_vis_n;*/
		fv<<HYPER[i].vis_force[A_X]<<","<<HYPER[i].vis_force[A_Y]<<","<<HYPER[i].vis_force[A_Z]<<",";
	}

	fv<<endl;
	fv.close();

	cout<<"---------OK"<<endl;

	/*
	int h_num=hyper_number;
	double d_nator[3];
	for(int D=0;D<DIMENSION;D++)	d_nator[D]=0;

	double n_nator=0;
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double vis=CON.get_h_vis();

	if(t==1)	system("mkdir Viscousity");

	stringstream s;
	s<<"./Viscousity/vis_force"<<t<<".dat";
	string filename=s.str();
	ofstream fs(filename);

	if(CON.get_FEM_flag()==OFF&&t==1)	for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]=0;	//t==1‚ð‘I‘ð€‚É‰Á‚¦‚½15/2/9

	calc_spl_f(CON,PART,HYPER1,h_num);

	for(int i=0;i<h_num;i++)
	{
		for(int D=0;D<DIMENSION;D++)	HYPER[i].vis_force[D]=0;
		for(int j=0;j<h_num;j++)
		{
			for(int D=0;D<DIMENSION;D++)	d_nator[D]=6/mi*vis*(HYPER[j].p[D]-HYPER[i].p[D])*HYPER1[i*h_num+j].spl_f;	//Ž®’ù³‚·‚é‚à‚Ð‚Ç‚­‚È‚Á‚½15/2/3
			for(int D=0;D<DIMENSION;D++)	HYPER[i].vis_force[D]+=d_nator[D];
		}

		n_nator=0;
		for(int j=0;j<h_num;j++)
		{
			n_nator+=((PART[j].r[A_X]-PART[i].r[A_X])*(PART[j].r[A_X]-PART[i].r[A_X])+(PART[j].r[A_Y]-PART[i].r[A_Y])*(PART[j].r[A_Y]-PART[i].r[A_Y])+(PART[j].r[A_Z]-PART[i].r[A_Z])*(PART[j].r[A_Z]-PART[i].r[A_Z]))*HYPER1[i*h_num+j].spl_f;
		}
		for(int D=0;D<DIMENSION;D++)	HYPER[i].vis_force[D]/=n_nator;
	}

	for(int i=0;i<h_num;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			HYPER[i].p[D]+=HYPER[i].vis_force[D];
			fs<<HYPER[i].vis_force[D]<<" ";
		}
		fs<<endl;
	}
	fs<<endl;
	fs.close();*/
}

void calc_spl_f(mpsconfig &CON,vector<mpselastic>PART,vector<hyperelastic2>&HYPER1,int h_num)
{
	double h=CON.get_h_dis();
	double r=0,s=0,theta=0;

	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			theta=0;	//ŠÖ”C³15/2/8
			HYPER1[i*h_num+j].spl_f=0;
			r=sqrt(pow((PART[j].r[A_X]-PART[i].r[A_X]),2.0)+pow((PART[j].r[A_Y]-PART[i].r[A_Y]),2.0)+pow((PART[j].r[A_Z]-PART[i].r[A_Z]),2.0));
			s=r/(h/2);
			if(0<=s&&s<1)	theta=1/PI*(1-3/2*pow(s,2.0)+3/4*pow(s,3.0));
			if(1<=s&&s<2)	theta=1/PI/4*pow((2-s),3.0);
			HYPER1[i*h_num+j].spl_f=1/pow(h/2,3)*theta;	//ŠÖ”C³15/2/8
		}
	}
}


