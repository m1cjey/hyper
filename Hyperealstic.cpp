#include "stdafx.h"		

void calc_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int hyper_number,int t)
{	
	int N=hyper_number;

	if(t==1)
	{
		for(int i=0;i<N;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=0;
		for(int i=0;i<N;i++)for(int D=0;D<DIMENSION;D++)PART[i].q0[D]=PART[i].r[D];
		calc_constant(CON,PART,HYPER,HYPER1,N);
	}
	contact_judge_hyper(CON,PART,HYPER,t);
	cout<<"wall effect is calculated."<<endl;

	calc_stress(CON,HYPER,N);
	cout<<"stress is calculated."<<endl;

	newton_raphson(CON,PART,HYPER,HYPER1,N,t);
	cout<<"Newton_raphson is ended."<<endl;

	calc_half_p(CON,PART,HYPER,HYPER1,N,0,t);
	cout<<"half_p is calculated."<<endl;

	calc_F(PART,HYPER,HYPER1,N,t);
	cout<<"calc_F is calculated."<<endl;

	calc_stress(CON,HYPER,N);
	cout<<"calc_stress is renewed."<<endl;

	//calculation of DgDq
	double p_DgDq[3];	
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<N;i++)
		{
			for(int D=0;D<DIMENSION;D++)
			{
				HYPER1[j*N+i].DgDq[D]=0;
				p_DgDq[D]=0;
				for(int D2=0;D2<DIMENSION;D2++)	p_DgDq[D]+=HYPER[j].t_inverse_Fi[D][D2]*HYPER1[i*N+j].n0ij[D2];
			}
			for(int D=0;D<DIMENSION;D++)	HYPER1[j*N+i].DgDq[D]=HYPER[j].J*p_DgDq[D];
		}
	}

	calc_differential_p(CON,HYPER,HYPER1,N);
	cout<<"Differential_p is calculated."<<endl;

	renew_lambda(CON,HYPER,HYPER1,N);
	cout<<"Lambda is renewed."<<endl;

	calc_half_p(CON,PART,HYPER,HYPER1,N,1,t);
	cout<<"p is calculated."<<endl;

	if(t==1 || t%CON.get_interval()==0)
	{
		output_hyper_data(PART,HYPER,t);
		momentum_movie_AVS(CON,t,PART,HYPER,N);
	}

	cout<<"Hypercalculation is ended."<<endl;

}

void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int particle_number)
{
	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	int num=particle_number;

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

	ofstream ai("Ai.csv");
	ofstream inai("inverse_Ai.csv");

	ofstream aii("aiin.csv");
	ofstream n0("noij.csv");
	ofstream wi("wiin.csv");
	
	ofstream dg("initial_DgDq.csv");

	////近傍粒子の記憶とaiin,wiin,Aiの計算
	double dis=0;
	int N=0;

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
				HYPER1[i*num+j].wiin=kernel4(r,dis);	//一時消去15/2/10
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
	double a[DIMENSION];
	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)	
		{
			double w=HYPER1[i*num+j].wiin;		
			for(int D=0;D<DIMENSION;D++)	a[D]=HYPER1[i*num+j].aiin[D];
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].Ai[D][D2]+=w*a[D]*a[D2];
		}
	}

	//inverse_Ai,t_inverse_Aiの計算
	double **p_Ai=new double *[DIMENSION];
	double **M=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Ai[D]=new double [DIMENSION];
		M[D]=new double [DIMENSION];
	}
	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Ai[D][D2]=HYPER[i].Ai[D][D2];
		inverse(p_Ai,DIMENSION);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].inverse_Ai[D][D2]=p_Ai[D][D2];
		transpose(p_Ai,M);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].t_inverse_Ai[D][D2]=M[D][D2];
	}

	////Fiの計算
	double **p_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)	p_Fi[D]=new double[DIMENSION];

	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				p_Fi[D][D2]=0;
				for(int D3=0;D3<DIMENSION;D3++)		p_Fi[D][D2]+=HYPER[i].Ai[D][D3]*HYPER[i].inverse_Ai[D3][D2];
			}
		}		
		double J=calc_det(p_Fi,DIMENSION);
//		cout<<"HYPER["<<i<<"].J="<<HYPER[i].J<<endl;
		if(J>=0)	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].differential_Fi[D][D2]=1/pow(J,1.0/3.0)*p_Fi[D][D2];
		else
		{
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].differential_Fi[D][D2]=-1/pow(-J,1.0/3.0)*p_Fi[D][D2];
		}
		transpose(p_Fi,M);
		inverse(M,DIMENSION);

		HYPER[i].J=J;
		for(int D=0;D<DIMENSION;D++)for(int D2=0;D2<DIMENSION;D2++)HYPER[i].t_inverse_Fi[D][D2]=M[D][D2];
	}

	for(int D=0;D<DIMENSION;D++)
	{
		delete[]p_Ai[D];
		delete[]	p_Fi[D];
		delete[] M[D];
	}
	delete[]p_Ai;
	delete[]	p_Fi;
	delete[]M;

	//n0ijの計算
	double p_n0ij[DIMENSION];
	double p_n0ij_2[DIMENSION];
	double p_n0ij_3[DIMENSION];

	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			if(j!=i)
			{
				for(int D=0;D<DIMENSION;D++)
				{
					p_n0ij[D]=0;
					for(int D2=0;D2<DIMENSION;D2++)	p_n0ij[D]+=HYPER[j].t_inverse_Ai[D][D2]*HYPER1[i*num+j].aiin[D2];
				}
				for(int D=0;D<DIMENSION;D++)	HYPER1[i*num+j].n0ij[D]=V*HYPER1[j*num+i].wiin*p_n0ij[D];
			}
			else
			{
				for(int D=0;D<DIMENSION;D++)
				{
					p_n0ij_2[D]=0;
					for(int k=0;k<num;k++)		p_n0ij_2[D]+=HYPER1[k*num+i].wiin*HYPER1[i*num+k].aiin[D];
				}
				for(int D=0;D<DIMENSION;D++)
				{
					p_n0ij_3[D]=0;
					for(int D2=0;D2<DIMENSION;D2++)	p_n0ij_3[D]+=HYPER[i].t_inverse_Ai[D][D2]*p_n0ij_2[D2];
				}
				for(int D=0;D<DIMENSION;D++)	HYPER1[i*num+j].n0ij[D]=V*p_n0ij_3[D];
			}
		}
	}

	
	//DgDqの計算
	double p_DgDq[3];
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
			for(int D=0;D<DIMENSION;D++)	n0<<HYPER1[i*num+j].n0ij[D]<<",";
			n0<<",";
			dg<<HYPER1[i*num+j].DgDq[A_X]<<","<<HYPER1[i*num+j].DgDq[A_Y]<<","<<HYPER1[i*num+j].DgDq[A_Z]<<","<<",";
		}
		dg<<endl<<endl;
		n0<<endl<<endl;
	}
	
	ai.close();
	inai.close();
	n0.close();
	dg.close();
}


/////ニュートンラフソン法 
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int particle_number,int t)
{
	/////fx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照
	/////DfDx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照

	int calc_type=1;//ニュートラフソンの反復方法 0:偏微分項の逆行列をそのまま求める　1:線形方程式を利用

	//pn=2;//test,とりあえず2元でとけるかどうか確認 
	//////////////////　f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm



	int N=particle_number;
	double *fx=new double [N];//関数値。
	double *DfDx=new double [N*N];//関数の偏微分値。
	double *XX=new double [N];//現在の解。	
	double *XX_old=new double [N];//1ステップ前の解。
	double ep=1e-3;//収束判定
	double E=1;//現在の誤差
/*	double start=0;
	double end=0;
	double newton_t=0;*/
	int count=0;//反復回数
	double d;
	double V=get_volume(&CON);
	int dec_flag=OFF;

	for(int i=0;i<N;i++)
	{
		XX[i]=1;
		XX_old[i]=0;
		fx[i]=0;
		for(int j=0;j<N;j++)	DfDx[i*N+j]=0;
	}

	stringstream ss;
	ss<<"./Newton_raphson/newton_raphson"<<t<<".dat";
	ofstream newton(ss.str());

	//	for(int i=0; i<N; i++) XX[i]=1;///初期値を与える。とりあえず1で
	cout<<"NR法開始----";
//	start=clock();
	while(E>ep)
	{
		count++;
		for(int i=0; i<N; i++)	XX_old[i]=XX[i];	//解を記憶

//		if(count==1)	for(int i=0;i<N;i++)	for(int j=0;j<N;j++)	for(int D=0;D<DIMENSION;D++)	HYPER1[i*N+j].newton_DgDq[D]=HYPER1[i*N+j].DgDq[D];

		calc_newton_function(CON,PART,HYPER,HYPER1,XX,fx,DfDx,N,count,t);


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

			for(int i=0; i<N; i++) 
			{
				d=0; //変化量
				for(int j=0; j<N; j++)	d+=DfDx[i*N+j]*fx[j];
				XX[i]-=d;
			}
		}
		else if(calc_type==1)//逆行列を用いない、安定するはずだが、遅くなるはず
		{
			gauss(DfDx,fx,N);
			for(int i=0;i<N;i++)	XX[i]-=fx[i];
		}

		//誤差の評価
		double E_old=E;
		double sum=0;
		for(int i=0; i<N; i++)
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

			cout<<"反復回数="<<count<<" E="<<E<<endl;
			newton<<"反復回数"<<count<<"E"<<" "<<E<<endl;

//			if(count>CON.get_nr()/2)
			{
				stringstream ss4;
				ss4<<"./Newton_raphson/lambda "<<"t"<<t<<" count"<<count<<".dat";
				ofstream lam(ss4.str());

				lam<<"lambda"<<"t"<<t<<"count"<<count<<endl;

				for(int i=0; i<N; i++) lam<<i<<"	"<<XX[i]<<endl;	
				lam.close();
			}
		}
		if(count>CON.get_nr())	break;
	}
//	end=clock();
//	newton_t=(end-start)/CLOCKS_PER_SEC;

	cout<<"反復完了"<<endl;
	newton<<"反復完了"<<endl;

	for(int i=0;i<N;i++) HYPER[i].lambda=XX[i];
//	for(int i=0;i<N;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
	newton.close();
	delete[]	fx;
	delete[]	DfDx;
	delete[]	XX;
	delete[]	XX_old;
}

void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int num,int count,int t)
{

	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();

	//half_pの計算
	double p_half_p[3][3];
	double p_half_p2[3];
	double p_half_p3[3];
	double n_half_p[3];

	double *n_rx=new double[num];
	double *n_ry=new double[num];
	double *n_rz=new double[num];
	for(int i=0;i<num;i++)
	{
		n_rx[i]=0;
		n_ry[i]=0;
		n_rz[i]=0;
	}

	double **n_DgDq_x=new double *[num];
	double **n_DgDq_y=new double *[num];
	double **n_DgDq_z=new double *[num];
	for(int i=0;i<num;i++)
	{
		n_DgDq_x[i]=new double[num];
		n_DgDq_y[i]=new double [num];
		n_DgDq_z[i]=new double [num];
	}

	//Fの計算
	double **p_Fi=new double *[DIMENSION];
	double **M=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Fi[D]=new double [DIMENSION];
		M[D]=new double [DIMENSION];
	}
	
	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{		
			for(int D=0;D<DIMENSION;D++)
			{	
				for(int D2=0;D2<DIMENSION;D2++)
				{
					p_half_p[D][D2]=0;
					p_half_p[D][D2]=HYPER[j].stress[D][D2];
				}
				p_half_p[D][D]-=lambda[j];
			}
			for(int D=0;D<DIMENSION;D++)
			{	
				p_half_p2[D]=0;
				for(int D2=0;D2<DIMENSION;D2++)	p_half_p2[D]+=p_half_p[D][D2]*HYPER1[j*num+i].DgDq[D2];
			}
			for(int D=0;D<DIMENSION;D++)	p_half_p3[D]+=p_half_p2[D];
		}//jに関するfor文の終わり			
//		cout<<"partial_half["<<i<<"]="<<Dt/2*p_half_p3[A_X]<<" "<<Dt/2*p_half_p3[A_Y]<<" "<<Dt/2*p_half_p3[A_Z]<<endl;
		for(int D=0;D<DIMENSION;D++)	n_half_p[D]=HYPER[i].p[D]+Dt/2*p_half_p3[D];
		////位置座標の更新
		n_rx[i]=PART[i].r[A_X]+Dt*n_half_p[A_X]/mi;
		n_ry[i]=PART[i].r[A_Y]+Dt*n_half_p[A_Y]/mi;
		n_rz[i]=PART[i].r[A_Z]+Dt*n_half_p[A_Z]/mi;
	}

	for(int i=0;i<num;i++)
	{
	////qiinとFiの更新
		//初期化
		double fi[3][3];
		double fi0[3][3];
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)	
			{
				fi0[D][D2]=0;
				p_Fi[D][D2]=0;
			}
		}
		//DgDqの計算

		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]=0;
		int in_num=HYPER[i].N;

		for(int in=0;in<in_num;in++)
		{
			int inn=HYPER[i].NEI[in];
			double w=HYPER1[i*num+inn].wiin;
			for(int D=0;D<DIMENSION;D++)
			{
				fi0[A_X][D]=w*(n_rx[inn]-n_rx[i])*HYPER1[i*num+inn].aiin[D];
				fi0[A_Y][D]=w*(n_ry[inn]-n_ry[i])*HYPER1[i*num+inn].aiin[D];
				fi0[A_Z][D]=w*(n_rz[inn]-n_rz[i])*HYPER1[i*num+inn].aiin[D];
			}
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]+=fi0[D][D2];
		}

		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				p_Fi[D][D2]=0;
				for(int D3=0;D3<DIMENSION;D3++)	p_Fi[D][D2]+=fi[D][D3]*HYPER[i].inverse_Ai[D3][D2];
			}
		}
		double J=calc_det(p_Fi,DIMENSION);
		fx[i]=V*(1-J);
		transpose(p_Fi,M);
		inverse(M,DIMENSION);

		//calculation of DgDq
		for(int j=0;j<num;j++)
		{
			n_DgDq_x[i][j]=0;
			n_DgDq_y[i][j]=0;
			n_DgDq_z[i][j]=0;

			double	p_n_DgDq[3];
			for(int D=0;D<DIMENSION;D++)
			{
				p_n_DgDq[D]=0;
				for(int D2=0;D2<DIMENSION;D2++)	p_n_DgDq[D]+=M[D][D2]*HYPER1[j*num+i].n0ij[D2];
			}
			n_DgDq_x[i][j]=J*p_n_DgDq[A_X];
			n_DgDq_y[i][j]=J*p_n_DgDq[A_Y];
			n_DgDq_z[i][j]=J*p_n_DgDq[A_Z];
		}
	}

	for(int i=0;i<num;i++)
	{
		for(int j=0;j<num;j++)
		{
			double DFDlambda=0;
			DfDx[i*num+j]=0;
			for(int k=0;k<num;k++)
			{		
				DFDlambda+=n_DgDq_x[i][k]*HYPER1[j*num+k].DgDq[A_X];
				DFDlambda+=n_DgDq_y[i][k]*HYPER1[j*num+k].DgDq[A_Y];
				DFDlambda+=n_DgDq_z[i][k]*HYPER1[j*num+k].DgDq[A_Z];
			}
			DfDx[i*num+j]=-Dt*Dt/2/mi*DFDlambda;
		}
	}



//	if(count%200==0 && count>CON.get_nr()/2)
	if(count%200==0)
	{
		stringstream s;
		s<<"./Newton_raphson/position "<<"t"<<t<<" count"<<count<<".dat";
		ofstream fs(s.str());

		stringstream s2;
		s2<<"./Newton_raphson/fx "<<"t"<<t<<" count"<<count<<".dat";
		ofstream fs2(s2.str());
		
		stringstream s3;
		s3<<"./Newton_raphson/DfDx "<<"t"<<t<<" count"<<count<<".dat";
		ofstream fs3(s3.str());

		for(int i=0;i<num;i++)
		{
			fs<<i<<"	"<<n_rx[i]<<"	"<<n_ry[i]<<"	"<<n_rz[i]<<endl;		
			fs2<<i<<"	"<<fx[i]<<endl;
			for(int j=0;j<num;j++)
			{
				fs3<<DfDx[i*num+j]<<",";
			}
			fs3<<endl;
		}
		fs.close();
		fs2.close();
		fs3.close();
	}

	for(int D=0;D<DIMENSION;D++)
	{
		delete[]	p_Fi[D];
		delete[]M[D];
	}
	delete[]	p_Fi;
	delete[]M;

	for(int i=0;i<num;i++)
	{
		delete[]n_DgDq_x[i];
		delete[]n_DgDq_y[i];
		delete[]n_DgDq_z[i];
	}
	delete[]n_DgDq_x;
	delete[]n_DgDq_y;
	delete[]n_DgDq_z;

	delete[]n_rx;
	delete[]n_ry;
	delete[]n_rz;
}

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int particle_number,bool repetation,int t)
{
	double Dt=CON.get_dt();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();

	double p_half_p[3][3];
	double p_half_p2[3];
	double p_half_p3[3];
//	int flag_wall=0;	//運動量の与え方を壁を通り抜けた粒子のみから全体に改良　15/2/5

	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	p_half_p3[D]=0;
		for(int j=0;j<particle_number;j++)
		{				
			for(int D=0;D<DIMENSION;D++)
			{
				for(int D2=0;D2<DIMENSION;D2++)
				{
					p_half_p[D][D2]=0;
					p_half_p[D][D2]=HYPER[j].stress[D][D2];
				}
				p_half_p[D][D]-=HYPER[j].lambda;
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
		}
		else
		{
			for(int D=0;D<DIMENSION;D++)	HYPER[i].old_p[D]=HYPER[i].p[D];
			for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]=HYPER[i].half_p[D]+Dt/2*p_half_p3[D];

			HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
			HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
			HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];

//			for(int D=0;D<DIMENSION;D++)cout<<"p["<<i<<"]="<<HYPER[i].half_p[D]<<"+"<<Dt<<"/2*"<<p_half_p3[D]<<" = "<<HYPER[i].p[D]<<endl;
		}
	}//iに関するfor文の終わり
	/*
/*	if(repetation==0)
	{
		for(int i=0;i<particle_number;i++)
		{
			cout<<"half_p["<<i<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER[i].half_p[D]<<" ";
			cout<<endl;
		}
		cout<<endl;
		*/
}

void calc_F(vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int particle_number,int t)
{
	////qiinとFiの更新
	//初期化
	int num=particle_number;
	double fi[3][3];
	double fi0[3][3];

	double **p_Fi=new double *[DIMENSION];
	double **M=new double *[DIMENSION];

	for(int D=0;D<DIMENSION;D++)
	{
		p_Fi[D]=new double[DIMENSION];
		M[D]=new double[DIMENSION];
	}

	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		fi0[D][D2]=0;

	//計算
	for(int i=0;i<num;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]=0;
		int in_num=HYPER[i].N;
		for(int in=0;in<in_num;in++)
		{
			int inn=HYPER[i].NEI[in];
			double w=HYPER1[i*num+inn].wiin;
			double a[3];
			for(int D=0;D<DIMENSION;D++) a[D]=HYPER1[i*num+inn].aiin[D];
			for(int D=0;D<DIMENSION;D++)for(int D2=0;D2<DIMENSION;D2++)	fi0[D][D2]=w*(PART[inn].r[D]-PART[i].r[D])*a[D2];
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		fi[D][D2]+=fi0[D][D2];
		}
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				p_Fi[D][D2]=0;
				for(int D3=0;D3<DIMENSION;D3++)		p_Fi[D][D2]+=fi[D][D3]*HYPER[i].inverse_Ai[D3][D2];
			}
		}

		double J=calc_det(p_Fi,DIMENSION);
		HYPER[i].J=J;

		if(J>=0)	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		HYPER[i].differential_Fi[D][D2]=1/pow(J,1.0/3.0)*p_Fi[D][D2];
		else
		{
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		HYPER[i].differential_Fi[D][D2]=-1/pow(-J,1.0/3.0)*p_Fi[D][D2];
		}				
		inverse(p_Fi,DIMENSION);
		transpose(p_Fi,M);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		HYPER[i].t_inverse_Fi[D][D2]=M[D][D2];		
	}

	for(int D=0;D<DIMENSION;D++)
	{
		delete[]	p_Fi[D];
		delete[]M[D];
	}
	delete[]	p_Fi;
	delete[]M;

//	for(int i=0;i<num;i++)	cout<<"J["<<i<<"]="<<HYPER[i].J<<endl;

}

void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER,int particle_number)
{
	double c10=CON.get_c10();
	double c01=CON.get_c01();
	double b[3][3];
	double bb[3][3];
	double trace_b,trace_bb;
	
	for(int j=0;j<particle_number;j++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[j].t_differential_Fi[D][D2]=HYPER[j].differential_Fi[D2][D];

		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				b[D][D2]=0;				
				for(int D3=0;D3<DIMENSION;D3++)	b[D][D2]+=HYPER[j].differential_Fi[D][D3]*HYPER[j].t_differential_Fi[D3][D2];
			}
		}
		for(int D=0;D<DIMENSION;D++)
		{
			for(int D2=0;D2<DIMENSION;D2++)
			{
				bb[D][D2]=0;				
				for(int D3=0;D3<DIMENSION;D3++)	bb[D][D2]+=b[D][D3]*b[D3][D2];
			}
		}
		trace_b=0;
		trace_bb=0;
		for(int D=0;D<DIMENSION;D++)
		{
			trace_b+=b[D][D];
			trace_bb+=bb[D][D];
		}
		for(int D=0;D<DIMENSION;D++)
		{
			HYPER[j].stress[D][D]=2/HYPER[j].J*((c10+c01*trace_b)*(b[D][D]-1.0/3.0*trace_b)-c01*(bb[D][D]-1.0/3.0*trace_bb));
			for(int D2=0;D2<DIMENSION;D2++)	
			{
				if(D2!=D)	HYPER[j].stress[D][D2]=2/HYPER[j].J*((c10+c01*trace_b)*b[D][D2]-c01*bb[D][D2]);
			}
		}
	}
}

void calc_differential_p(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int particle_number)
{
	double Dt=CON.get_dt();

	double p_differential_p[3];
	double p_differential_p2[3];

	for(int i=0;i<particle_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)	p_differential_p2[D]=0;
		for(int j=0;j<particle_number;j++)
		{							
			for(int D=0;D<DIMENSION;D++)
			{	
				p_differential_p[D]=0;
				for(int D2=0;D2<DIMENSION;D2++)	p_differential_p[D]+=HYPER[j].stress[D][D2]*HYPER1[j*particle_number+i].DgDq[D2];
			}
			for(int D=0;D<DIMENSION;D++)	p_differential_p2[D]+=p_differential_p[D];				
		}		
		for(int D=0;D<DIMENSION;D++)	HYPER[i].differential_p[D]=HYPER[i].half_p[D]+Dt/2*p_differential_p2[D];
	}
}

void renew_lambda(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int particle_number)
{
	double le=CON.get_distancebp();
	double Dt=CON.get_dt();
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
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,int particle_number)
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
	fp<<n<<" "<<n<<endl;	//節点数と要素数出力
	
	//節点番号とその座標の出力 
	for(int i=0;i<n;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=0;i<n;i++)
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
	for(int i=0;i<n;i++)
	{
		fp<<i<<" "<<HYPER[i].p[A_X]<<" "<<HYPER[i].p[A_Y]<<" "<<HYPER[i].p[A_Z]<<" "<<HYPER[i].lambda<<" "<<HYPER[i].J<<" "<<HYPER[i].ang_p[A_X]<<" "<<HYPER[i].ang_p[A_Y]<<" "<<HYPER[i].ang_p[A_Z]<<" "<<endl;
		//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
	}

	cout<<"OK"<<endl;
	fp.close();

}


void contact_judge_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,int t)
{
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
		NEI_w[i]=new double[200];
	}

	if(t!=1)	for(int i=h_num;i<w_num+h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].r[D]+=PART[i].u[D]*Dt;
		
//		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		for(int i=0;i<h_num;i++)
		{
			int N_w=0;
			double pnd=0;

			for(int j=0;j<w_num;j++)
			{
				int k=j+h_num;
				dis[i][j]=0;
				w[i][j]=0;
			//初期配置の粒子とは普通に圧力勾配を計算（内力計算）
			//else if(PART[i].contact==true)
				double dis_temp0=0;
				double dis_temp1=0;
				double dis_temp=0;
				for(int D=0;D<DIMENSION;D++)
				{
					dis_temp0=PART[i].r[D]-PART[k].r[D];
					dis_temp1+=dis_temp0*dis_temp0;
				}
				dis_temp=sqrt(dis_temp1);
				if(dis_temp<r)
				{
					dis[i][j]=dis_temp;
					w[i][j]=kernel(r,dis_temp);
					NEI_w[i][N_w]=j;
					pnd+=w[i][j];
					N_w++;				
				}
				//現在位置での周辺粒子数を取得
				
			}
			if(N_w>0)
			{
				double gra_accel_i[DIMENSION];
				for(int D=0;D<DIMENSION;D++)	gra_accel_i[D]=0;

				for(int k=0;k<N_w;k++)
				{
					int j=NEI_w[i][k];
					double rij[DIMENSION];
					double accel_i[DIMENSION];
					for(int D=0;D<DIMENSION;D++)
					{
						rij[D]=0;
						accel_i[D]=0;

						rij[D]=PART[j].r[D]-PART[i].r[D];
						accel_i[D]=1/mi/Dt*(HYPER[i].old_p[D]-HYPER[i].p[D]);
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
}

void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,int t)
{
	int N=HYPER.size();

	stringstream s;
	s<<"./d_Momentum/differential_momentum"<<t<<".dat";
	ofstream sdm(s.str());

	stringstream s2;
	s2<<"./Position/position"<<t<<".dat";
	ofstream po(s2.str());

	stringstream s3;
	s3<<"./Hy_stress/Stress"<<"t"<<t<<".csv";
	ofstream stress(s3.str());

	stringstream s4;
	s4<<"./Fi/det_Fi "<<"t"<<t<<".dat";
	ofstream Fi_j(s4.str());
	
	stringstream s5;
	s5<<"./Fi/differential_Fi "<<"t"<<t<<".csv";
	ofstream d_Fi(s5.str());

	stringstream s6;
	s6<<"./Fi/t_inverse_Fi "<<"t"<<t<<".csv";
	ofstream i_t_Fi(s6.str());

	stringstream s7;
	s7<<"./renew_Lambda/renew_lambda "<<"t"<<t<<".dat";
	ofstream lam(s7.str());


	sdm<<"Differential_p"<<t<<endl;
	po<<"Position"<<t<<endl;		
	stress<<"Stress"<<t<<endl;
	Fi_j<<"det_Fi"<<t<<endl;
	lam<<""<<t<<endl;

	for(int i=0;i<N;i++)
	{

		sdm<<i<<"	"<<HYPER[i].differential_p[A_X]<<"	"<<HYPER[i].differential_p[A_Y]<<"	"<<HYPER[i].differential_p[A_Z]<<endl;
		po<<i<<"	";
		Fi_j<<i<<"	"<<HYPER[i].J<<endl;
		lam<<i<<"	"<<HYPER[i].lambda<<endl;
		for(int D=0;D<DIMENSION;D++)
		{
				po<<PART[i].r[D]<<"	";
				stress<<i<<",";
				for(int D2=0;D2<DIMENSION;D2++)
				{
					stress<<HYPER[i].stress[D][D2]<<",";
					d_Fi<<HYPER[i].differential_Fi[D][D2]<<",";
					i_t_Fi<<HYPER[i].t_inverse_Fi[D][D2]<<",";
				}
				stress<<endl;
				d_Fi<<endl;
				i_t_Fi<<endl;

		}
		po<<PART[i].r[A_Z]<<endl;
		stress<<endl;
		d_Fi<<endl<<endl;
		i_t_Fi<<endl<<endl;
	}
	stress.close();
	Fi_j.close();
	po.close();

	lam.close();
	sdm.close();
	d_Fi.close();
	i_t_Fi.close();
}

void transpose(double **M, double **N)
{
	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		N[D][D2]=0;
	for(int D=0;D<DIMENSION;D++)	N[D][D]=M[D][D];
	N[0][1]=M[1][0];
	N[0][2]=M[2][0];
	N[1][0]=M[0][1];
	N[1][2]=M[2][1];
	N[2][0]=M[0][2];
	N[2][1]=M[1][2];
}

hyperelastic::hyperelastic()
{
	for(int i=0;i<200;i++)
	{
		NEI[i]=0;
		NEI_w[i]=0;
	}
	N=0;
	N_w=0;
	lambda=1;
	J=0;
	for(int D=0;D<DIMENSION;D++)
	{
		half_p[D]=0;
		differential_p[D]=0;
		old_p[D]=0;
		p[D]=0;
		ang_p[D]=0;
		vis_force[D]=0;
		for(int D2=0;D2<DIMENSION;D2++)
		{
			stress[D][D2]=0;
			Ai[D][D2]=0;
			t_inverse_Ai[D][D2]=0;
			t_inverse_Fi[D][D2]=0;
			differential_Fi[D][D2]=0;
			t_differential_Fi[D][D2]=0;
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


