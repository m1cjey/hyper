#include "stdafx.h"		

void calc_hyper(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int hyper_number,int t)
{	
	int h_num=hyper_number;
	cout<<"Hypercalculation starts."<<endl;

	//calc_gravity(CON,HYPER,h_num);
	
	if(t==1)
	{
		for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	PART[i].q0[D]=0;
		for(int i=0;i<h_num;i++)for(int D=0;D<DIMENSION;D++)PART[i].q0[D]=PART[i].r[D];
		calc_constant(CON,PART,HYPER,HYPER1,h_num);
		calc_stress(CON,HYPER,h_num);
	}

	contact_judge_hyper2(CON,PART,HYPER,h_num,t);

	newton_raphson(CON,PART,HYPER,HYPER1,h_num,t);

	calc_half_p(CON,PART,HYPER,HYPER1,h_num,0,t);

	calc_F(PART,HYPER,HYPER1,h_num,t);

	calc_stress(CON,HYPER,h_num);

	calc_differential_p(CON,HYPER,HYPER1,h_num);

	renew_lambda(CON,HYPER,HYPER1,h_num);

	calc_half_p(CON,PART,HYPER,HYPER1,h_num,1,t);

	if(t==1 || t%CON.get_interval()==0)
	{
		output_hyper_data(PART,HYPER,HYPER1,t);
		momentum_movie_AVS(CON,t,PART,HYPER,h_num);
		output_energy(CON,PART,HYPER,h_num,t);
	}

	cout<<"Hypercalculation ends."<<endl;
}

void calc_constant(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> &HYPER1,int hyper_number)
{
	cout<<"初期値計算";

	double le=CON.get_distancebp();
	double r=CON.get_h_dis();
	double Dt=CON.get_dt();
	double V=get_volume(&CON);	//考慮が必要かもしれない
	double mi=V*CON.get_hyper_density();
	int h_num=hyper_number;
	int model=CON.get_model_number();

	////初期運動量
	//曲げねじり
	if(model==21)
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
			HYPER[i].p[A_X]=-t*part_p*part_p*part_p*Y+b*(3*part_p*part_p-1);
			HYPER[i].p[A_Y]=t*part_p*part_p*part_p*X;
			HYPER[i].p[A_Z]=0;
			for(int D=0;D<DIMENSION;D++) HYPER[i].p[D]*=mi;
		}
	}
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
		for(int j=0;j<h_num;j++)
		{
			double wiin=0;
			double aiin[DIMENSION];

			for(int D=0;D<DIMENSION;D++)
			{
				aiin[D]=PART[j].q0[D]-PART[i].q0[D];
				HYPER1[i*h_num+j].aiin[D]=aiin[D];
			}
			double dis=sqrt(aiin[A_X]*aiin[A_X]+aiin[A_Y]*aiin[A_Y]+aiin[A_Z]*aiin[A_Z]);
			if(dis<r && j!=i)
			{	
				wiin=kernel4(r,dis);	//一時消去15/2/10
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
	double a[DIMENSION];
	double **p_Ai=new double *[DIMENSION];
	double **p_Fi=new double *[DIMENSION];
	double **M=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Ai[D]=new double [DIMENSION];
		p_Fi[D]=new double[DIMENSION];
		M[D]=new double [DIMENSION];
	}
	
	for(int i=0;i<h_num;i++)
	{
		//Aiの計算
		for(int j=0;j<h_num;j++)	
		{
			double w=HYPER1[i*h_num+j].wiin;		
			for(int D=0;D<DIMENSION;D++)	a[D]=HYPER1[i*h_num+j].aiin[D];
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].Ai[D][D2]+=w*a[D]*a[D2];
		}

		//inverse_Ai,t_inverse_Aiの計算
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Ai[D][D2]=HYPER[i].Ai[D][D2];
		inverse(p_Ai,DIMENSION);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].inverse_Ai[D][D2]=p_Ai[D][D2];
		transpose(p_Ai,M);
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].t_inverse_Ai[D][D2]=M[D][D2];
		
		//Fiの計算
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=HYPER[i].Ai[D][0]*HYPER[i].inverse_Ai[0][D2]+HYPER[i].Ai[D][1]*HYPER[i].inverse_Ai[1][D2]+HYPER[i].Ai[D][2]*HYPER[i].inverse_Ai[2][D2];
			

		//Jの計算
		double J=calc_det(p_Fi,DIMENSION);
		HYPER[i].J=J;
//		cout<<"HYPER["<<i<<"].J="<<HYPER[i].J<<endl;

		//d_Fiの計算
		if(J>=0)	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].differential_Fi[D][D2]=1/pow(J,1.0/3.0)*p_Fi[D][D2];
		else
		{
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	HYPER[i].differential_Fi[D][D2]=-1/pow(-J,1.0/3.0)*p_Fi[D][D2];
		}
		//t_inverse_Fiの計算
		transpose(p_Fi,M);
		inverse(M,DIMENSION);
		for(int D=0;D<DIMENSION;D++)for(int D2=0;D2<DIMENSION;D2++)HYPER[i].t_inverse_Fi[D][D2]=M[D][D2];
	}
	for(int D=0;D<DIMENSION;D++)
	{
		delete[]	p_Ai[D];
		delete[]	p_Fi[D];
		delete[]	M[D];
	}
	delete[]	p_Ai;
	delete[]	p_Fi;
	delete[]	M;


	////n0ijの計算
	double p_n0ij[DIMENSION];
	double p_n0ij_2[DIMENSION];
	double p_n0ij_3[DIMENSION];
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			if(j!=i)
			{
				for(int D=0;D<DIMENSION;D++)	p_n0ij[D]=HYPER[j].t_inverse_Ai[D][0]*HYPER1[i*h_num+j].aiin[0]+HYPER[j].t_inverse_Ai[D][1]*HYPER1[i*h_num+j].aiin[1]+HYPER[j].t_inverse_Ai[D][2]*HYPER1[i*h_num+j].aiin[2];
				for(int D=0;D<DIMENSION;D++)	HYPER1[i*h_num+j].n0ij[D]=V*HYPER1[j*h_num+i].wiin*p_n0ij[D];
			}
			else
			{
				for(int D=0;D<DIMENSION;D++)
				{
					p_n0ij_2[D]=0;
					for(int k=0;k<h_num;k++)		p_n0ij_2[D]+=HYPER1[k*h_num+i].wiin*HYPER1[i*h_num+k].aiin[D];
				}
				for(int D=0;D<DIMENSION;D++)	p_n0ij_3[D]=HYPER[i].t_inverse_Ai[D][0]*p_n0ij_2[0]+HYPER[i].t_inverse_Ai[D][1]*p_n0ij_2[1]+HYPER[i].t_inverse_Ai[D][2]*p_n0ij_2[2];
				for(int D=0;D<DIMENSION;D++)	HYPER1[i*h_num+j].n0ij[D]=V*p_n0ij_3[D];
			}
		}
	}
		
	////DgDqの計算
	double p_DgDq[3];
	for(int j=0;j<h_num;j++)
	{
		for(int i=0;i<h_num;i++)
		{
			for(int D=0;D<DIMENSION;D++)
			{
				p_DgDq[D]=HYPER[j].t_inverse_Fi[D][0]*HYPER1[i*h_num+j].n0ij[0]+HYPER[j].t_inverse_Fi[D][1]*HYPER1[i*h_num+j].n0ij[1]+HYPER[j].t_inverse_Fi[D][2]*HYPER1[i*h_num+j].n0ij[2];
				HYPER1[j*h_num+i].DgDq[D]=HYPER[j].J*p_DgDq[D];			
			}
		}
	}

	cout<<"----------OK"<<endl;
}


/////ニュートンラフソン法 
void newton_raphson(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int hyper_number,int t)
{
	/////fx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照
	/////DfDx(N*N行列の各成分)は1次元配列で格納、(i,j)成分なら[j*N+i]で参照

	int calc_type=1;//ニュートラフソンの反復方法 0:偏微分項の逆行列をそのまま求める　1:線形方程式を利用

	//pn=2;//test,とりあえず2元でとけるかどうか確認 
	//////////////////　f1(x1,x2) = x1^2 + x2^2 -5 = 0 f2(x1,x2) = x1^2/9+ x2^2 -1 = 0  http://homepage1.nifty.com/gfk/excel_newton_ren.htm

	int h_num=hyper_number;
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
	int dec_flag=ON;

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

		calc_newton_function(CON,PART,HYPER,HYPER1,XX,fx,DfDx,h_num,count,t);


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
			gauss(DfDx,fx,h_num);
			for(int i=0;i<h_num;i++)	XX[i]-=fx[i];
		}

		//誤差の評価
		double E_old=E;
		E=0;
		double sum=0;
		for(int i=0; i<h_num; i++)
		{
			E+=fabs(XX[i]-XX_old[i]);
			sum+=fabs(XX[i]);
		}
		E/=sum;


		if(count==1 || count%CON.get_interval()==0)
		{		
			/*
			cout<<"XX_old["<<i<<"]-d=X["<<i<<"]	";
			cout<<XX_old[i]<<" - ";
			cout<<d<<" = ";
			cout<<XX[i]<<endl;*/

			cout<<"反復回数	"<<count<<" E="<<E<<endl;
			output_newton_data2(E,XX,h_num,count,t);

		}
		if(count>CON.get_nr())	break;
		else if(dec_flag==ON)	if(E_old-E<0)	break;	
	}
//	end=clock();
//	newton_t=(end-start)/CLOCKS_PER_SEC;

	cout<<"反復完了";

	for(int i=0;i<h_num;i++) HYPER[i].lambda=XX[i];
//	for(int i=0;i<N;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;
	delete[]	fx;
	delete[]	DfDx;
	delete[]	XX;
	delete[]	XX_old;
	cout<<"---------- OK"<<endl;
}


void calc_newton_function(mpsconfig &CON,vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,double *lambda,double *fx,double *DfDx,int hyper_number,int count,int t)
{
	int h_num=hyper_number;
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	int vec_g[DIMENSION];
	vec_g[A_X]=0,vec_g[A_Y]=0,vec_g[A_Z]=1;

	////位置座標の更新
	double p_half_p[3][3];
	double p_half_p2[3];
	double p_half_p3[3];
	double n_half_p[3];

	double *n_rx=new double[h_num];
	double *n_ry=new double[h_num];
	double *n_rz=new double[h_num];
	for(int i=0;i<h_num;i++)
	{
		n_rx[i]=0;
		n_ry[i]=0;
		n_rz[i]=0;
	}
	
	for(int i=0;i<h_num;i++)
	{
		//half_pの計算
		for(int D=0;D<DIMENSION;D++)	p_half_p3[D]=0;
		for(int j=0;j<h_num;j++)
		{	
			p_half_p[0][0]=HYPER[j].stress[0][0]-lambda[j];
			p_half_p[1][1]=HYPER[j].stress[1][1]-lambda[j];
			p_half_p[2][2]=HYPER[j].stress[2][2]-lambda[j];
			p_half_p[0][1]=HYPER[j].stress[0][1];
			p_half_p[0][2]=HYPER[j].stress[0][2];
			p_half_p[1][0]=HYPER[j].stress[1][0];
			p_half_p[1][2]=HYPER[j].stress[1][2];
			p_half_p[2][0]=HYPER[j].stress[2][0];
			p_half_p[2][1]=HYPER[j].stress[2][1];

			for(int D=0;D<DIMENSION;D++)	p_half_p2[D]=p_half_p[D][0]*HYPER1[j*h_num+i].DgDq[0]+p_half_p[D][1]*HYPER1[j*h_num+i].DgDq[1]+p_half_p[D][2]*HYPER1[j*h_num+i].DgDq[2];
			for(int D=0;D<DIMENSION;D++)	p_half_p3[D]+=p_half_p2[D];
		}//jに関するfor文の終わり
//		cout<<"partial_half["<<i<<"]="<<Dt/2*p_half_p3[A_X]<<" "<<Dt/2*p_half_p3[A_Y]<<" "<<Dt/2*p_half_p3[A_Z]<<endl;
		for(int D=0;D<DIMENSION;D++)	n_half_p[D]=HYPER[i].p[D]+Dt/2*(p_half_p3[D]-9.8*mi*vec_g[D]);
		
		//位置座標の計算
		n_rx[i]=PART[i].r[A_X]+Dt*n_half_p[A_X]/mi;
		n_ry[i]=PART[i].r[A_Y]+Dt*n_half_p[A_Y]/mi;
		n_rz[i]=PART[i].r[A_Z]+Dt*n_half_p[A_Z]/mi;
	}



	////DgDqとfxの更新
	double fi0[3][3];
	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi0[D][D2]=0;
	
	double **p_Fi=new double *[DIMENSION];
	double **M=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Fi[D]=new double [DIMENSION];
		M[D]=new double [DIMENSION];
	}

	double **n_DgDq_x=new double *[h_num];
	double **n_DgDq_y=new double *[h_num];
	double **n_DgDq_z=new double *[h_num];
	for(int i=0;i<h_num;i++)
	{
		n_DgDq_x[i]=new double[h_num];
		n_DgDq_y[i]=new double [h_num];
		n_DgDq_z[i]=new double [h_num];
	}

	for(int i=0;i<h_num;i++)
	{
		//Fiの計算
		int in_h_num=HYPER[i].N;
		double fi[3][3];
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]=0;
		
		for(int in=0;in<in_h_num;in++)
		{
			int inn=HYPER[i].NEI[in];
			double w=HYPER1[i*h_num+inn].wiin;
			for(int D=0;D<DIMENSION;D++)
			{
				fi0[A_X][D]=w*(n_rx[inn]-n_rx[i])*HYPER1[i*h_num+inn].aiin[D];
				fi0[A_Y][D]=w*(n_ry[inn]-n_ry[i])*HYPER1[i*h_num+inn].aiin[D];
				fi0[A_Z][D]=w*(n_rz[inn]-n_rz[i])*HYPER1[i*h_num+inn].aiin[D];
			}
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]+=fi0[D][D2];
		}
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=fi[D][0]*HYPER[i].inverse_Ai[0][D2]+fi[D][1]*HYPER[i].inverse_Ai[1][D2]+fi[D][2]*HYPER[i].inverse_Ai[2][D2];
		
		//Jの計算
		double J=calc_det(p_Fi,DIMENSION);
		
		//fxの計算
		fx[i]=V*(1-J);

		//t_inverse_Fiの計算
		transpose(p_Fi,M);
		inverse(M,DIMENSION);

		//DgDqの計算
		for(int j=0;j<h_num;j++)
		{
			n_DgDq_x[i][j]=0;
			n_DgDq_y[i][j]=0;
			n_DgDq_z[i][j]=0;

			double	p_n_DgDq[3];
			for(int D=0;D<DIMENSION;D++)	p_n_DgDq[D]=M[D][0]*HYPER1[j*h_num+i].n0ij[0]+M[D][1]*HYPER1[j*h_num+i].n0ij[1]+M[D][2]*HYPER1[j*h_num+i].n0ij[2];
			n_DgDq_x[i][j]=J*p_n_DgDq[A_X];
			n_DgDq_y[i][j]=J*p_n_DgDq[A_Y];
			n_DgDq_z[i][j]=J*p_n_DgDq[A_Z];
		}
	}

	////DfDxの更新
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			double DFDlambda=0;
			DfDx[i*h_num+j]=0;
			for(int k=0;k<h_num;k++)
			{		
				DFDlambda+=n_DgDq_x[i][k]*HYPER1[j*h_num+k].DgDq[A_X];
				DFDlambda+=n_DgDq_y[i][k]*HYPER1[j*h_num+k].DgDq[A_Y];
				DFDlambda+=n_DgDq_z[i][k]*HYPER1[j*h_num+k].DgDq[A_Z];
			}
			DfDx[i*h_num+j]=-Dt*Dt/2/mi*DFDlambda;
		}
	}


	////出力
//	if(count%200==0 && count>CON.get_nr()/2)
	if(count%CON.get_interval()==0||count==1)		output_newton_data1(fx,DfDx,n_rx,n_ry,n_rz,h_num,count,t);

	for(int D=0;D<DIMENSION;D++)
	{
		delete[]	p_Fi[D];
		delete[]	M[D];
	}
	delete[]	p_Fi;
	delete[]	M;

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

void calc_half_p(mpsconfig &CON,vector<mpselastic> &PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int hyper_number,bool repetation,int t)
{
	if(repetation==0)	cout<<"仮の運動量＆位置座標計算";
	else	cout<<"運動量計算";

	int h_num=hyper_number;
	double Dt=CON.get_dt();
	double le=CON.get_distancebp();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	int vec_g[DIMENSION];
	vec_g[A_X]=0,vec_g[A_Y]=0,vec_g[A_Z]=1;

	double p_half_p[3][3];
	double p_half_p2[3];
	double p_half_p3[3];

	for(int i=0;i<h_num;i++)
	{
		p_half_p3[0]=0,p_half_p3[1]=0,p_half_p3[2]=0;
		for(int j=0;j<h_num;j++)
		{				
			p_half_p[0][0]=HYPER[j].stress[0][0]-HYPER[j].lambda;
			p_half_p[1][1]=HYPER[j].stress[1][1]-HYPER[j].lambda;
			p_half_p[2][2]=HYPER[j].stress[2][2]-HYPER[j].lambda;
			p_half_p[0][1]=HYPER[j].stress[0][1];
			p_half_p[0][2]=HYPER[j].stress[0][2];
			p_half_p[1][0]=HYPER[j].stress[1][0];
			p_half_p[1][2]=HYPER[j].stress[1][2];
			p_half_p[2][0]=HYPER[j].stress[2][0];
			p_half_p[2][1]=HYPER[j].stress[2][1];

			for(int D=0;D<DIMENSION;D++)	p_half_p2[D]=p_half_p[D][0]*HYPER1[j*h_num+i].DgDq[0]+p_half_p[D][1]*HYPER1[j*h_num+i].DgDq[1]+p_half_p[D][2]*HYPER1[j*h_num+i].DgDq[2];
			for(int D=0;D<DIMENSION;D++)	p_half_p3[D]+=p_half_p2[D];
		}//jに関するfor文の終わり	

//		cout<<"partial_half["<<i<<"]="<<Dt/2*p_half_p3[A_X]<<" "<<Dt/2*p_half_p3[A_Y]<<" "<<Dt/2*p_half_p3[A_Z]<<endl;

		if(repetation==0)
		{
			//half_pの更新
			for(int D=0;D<DIMENSION;D++)	HYPER[i].half_p[D]=HYPER[i].p[D]+Dt/2*(p_half_p3[D]-mi*9.8*vec_g[D]);
			//位置座標の更新
			for(int D=0;D<DIMENSION;D++)	PART[i].r[D]+=Dt*HYPER[i].half_p[D]/mi;
		}

		else
		{
			HYPER[i].old_p[A_X]=HYPER[i].p[A_X],	HYPER[i].old_p[A_Y]=HYPER[i].p[A_Y],	HYPER[i].old_p[A_Z]=HYPER[i].p[A_Z];
			//運動量の更新
			for(int D=0;D<DIMENSION;D++)
			{
				HYPER[i].p[D]=HYPER[i].half_p[D]+Dt/2*(p_half_p3[D]-9.8*mi*vec_g[D]);
				PART[i].u[D]=HYPER[i].half_p[D]/mi;
			}

			HYPER[i].ang_p[A_X]=PART[i].r[A_Y]*HYPER[i].p[A_Z]-PART[i].r[A_Z]*HYPER[i].p[A_Y];
			HYPER[i].ang_p[A_Y]=PART[i].r[A_Z]*HYPER[i].p[A_X]-PART[i].r[A_X]*HYPER[i].p[A_Z];
			HYPER[i].ang_p[A_Z]=PART[i].r[A_X]*HYPER[i].p[A_Y]-PART[i].r[A_Y]*HYPER[i].p[A_X];
//			for(int D=0;D<DIMENSION;D++)cout<<"p["<<i<<"]="<<HYPER[i].half_p[D]<<"+"<<Dt<<"/2*"<<p_half_p3[D]<<" = "<<HYPER[i].p[D]<<endl;
		}
	}//iに関するfor文の終わり
	
/*	if(repetation==0)
	{
		for(int i=0;i<h_num;i++)
		{
			cout<<"half_p["<<i<<"]=";
			for(int D=0;D<DIMENSION;D++)	cout<<HYPER[i].half_p[D]<<" ";
			cout<<endl;
		}
		cout<<endl;
		*/
	cout<<"----------OK"<<endl;
}

void calc_F(vector<mpselastic> PART,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int hyper_number,int t)
{
	cout<<"Fi計算";
	////Fiの更新
	int h_num=hyper_number;
	double fi0[3][3];

	double **p_Fi=new double *[DIMENSION];
	double **M=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		p_Fi[D]=new double[DIMENSION];
		M[D]=new double[DIMENSION];
	}
	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		fi0[D][D2]=0;

	for(int i=0;i<h_num;i++)
	{
		//Fiの計算
		double fi[3][3];
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]=0;
		int in_h_num=HYPER[i].N;

		for(int in=0;in<in_h_num;in++)
		{
			int inn=HYPER[i].NEI[in];
			double w=HYPER1[i*h_num+inn].wiin;
			double a[3];
			for(int D=0;D<DIMENSION;D++) a[D]=HYPER1[i*h_num+inn].aiin[D];
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi0[D][D2]=w*(PART[inn].r[D]-PART[i].r[D])*a[D2];
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	fi[D][D2]+=fi0[D][D2];
		}

		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	p_Fi[D][D2]=fi[D][0]*HYPER[i].inverse_Ai[0][D2]+fi[D][1]*HYPER[i].inverse_Ai[1][D2]+fi[D][2]*HYPER[i].inverse_Ai[2][D2];


		//Jの計算
		double J=calc_det(p_Fi,DIMENSION);
	//	for(int i=0;i<h_num;i++)	cout<<"J["<<i<<"]="<<J<<endl;
		HYPER[i].J=J;

		//d_Fiの計算
		if(J>=0)	for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		HYPER[i].differential_Fi[D][D2]=1/pow(J,1.0/3.0)*p_Fi[D][D2];
		else
		{
			for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)		HYPER[i].differential_Fi[D][D2]=-1/pow(-J,1.0/3.0)*p_Fi[D][D2];
		}				

		//t_inverse_Fiの計算
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

	//calculation of DgDq
	double p_DgDq[3];	
	for(int j=0;j<h_num;j++)
	{
		for(int i=0;i<h_num;i++)
		{
			for(int D=0;D<DIMENSION;D++)
			{
				HYPER1[j*h_num+i].DgDq[D]=0;
				p_DgDq[D]=HYPER[j].t_inverse_Fi[D][0]*HYPER1[i*h_num+j].n0ij[0]+HYPER[j].t_inverse_Fi[D][1]*HYPER1[i*h_num+j].n0ij[1]+HYPER[j].t_inverse_Fi[D][2]*HYPER1[i*h_num+j].n0ij[2];
				HYPER1[j*h_num+i].DgDq[D]=HYPER[j].J*p_DgDq[D];
			}
		}
	}
	cout<<"----------OK"<<endl;

}

void calc_stress(mpsconfig &CON,vector<hyperelastic> &HYPER,int hyper_number)
{
	cout<<"応力計算";
	int h_num=hyper_number;

	double **d_Fi=new double *[DIMENSION];
	double **t_d_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		d_Fi[D]=new double	[DIMENSION];
		t_d_Fi[D]=new double [DIMENSION];
	}

	double c10=CON.get_c10();
	double c01=CON.get_c01();
	
	for(int j=0;j<h_num;j++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	d_Fi[D][D2]=HYPER[j].differential_Fi[D][D2];
		transpose(d_Fi,t_d_Fi);
		
		double b[3][3];
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	b[D][D2]=d_Fi[D][0]*t_d_Fi[0][D2]+d_Fi[D][1]*t_d_Fi[1][D2]+d_Fi[D][2]*t_d_Fi[2][D2];

		double bb[3][3];
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	bb[D][D2]=b[D][0]*b[0][D2]+b[D][1]*b[1][D2]+b[D][2]*b[2][D2];

		double trace_b=b[0][0]+b[1][1]+b[2][2];
		double trace_bb=bb[0][0]+bb[1][1]+bb[2][2];
	
		HYPER[j].stress[0][0]=2/HYPER[j].J*((c10+c01*trace_b)*(b[0][0]-1.0/3.0*trace_b)-c01*(bb[0][0]-1.0/3.0*trace_bb));
		HYPER[j].stress[1][1]=2/HYPER[j].J*((c10+c01*trace_b)*(b[1][1]-1.0/3.0*trace_b)-c01*(bb[1][1]-1.0/3.0*trace_bb));
		HYPER[j].stress[2][2]=2/HYPER[j].J*((c10+c01*trace_b)*(b[2][2]-1.0/3.0*trace_b)-c01*(bb[2][2]-1.0/3.0*trace_bb));
		HYPER[j].stress[0][1]=2/HYPER[j].J*((c10+c01*trace_b)*b[0][1]-c01*bb[0][1]);
		HYPER[j].stress[0][2]=2/HYPER[j].J*((c10+c01*trace_b)*b[0][2]-c01*bb[0][2]);
		HYPER[j].stress[1][0]=2/HYPER[j].J*((c10+c01*trace_b)*b[1][0]-c01*bb[1][0]);
		HYPER[j].stress[1][2]=2/HYPER[j].J*((c10+c01*trace_b)*b[1][2]-c01*bb[1][2]);
		HYPER[j].stress[2][0]=2/HYPER[j].J*((c10+c01*trace_b)*b[2][0]-c01*bb[2][0]);
		HYPER[j].stress[2][1]=2/HYPER[j].J*((c10+c01*trace_b)*b[2][1]-c01*bb[2][1]);
	}

	for(int D=0;D>DIMENSION;D++)
	{
		delete[]	d_Fi[D];
		delete[]	t_d_Fi[D];
	}
	delete[]	d_Fi;
	delete[]	t_d_Fi;

	cout<<"----------OK"<<endl;
}

void calc_differential_p(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int hyper_number)
{
	cout<<"運動量微分値計算";

	int h_num=hyper_number;
	int vec_g[DIMENSION];
	vec_g[A_X]=0,vec_g[A_Y]=0,vec_g[A_Z]=1;
	double Dt=CON.get_dt();
	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();
	double p_differential_p[3];
	double p_differential_p2[3];

	for(int i=0;i<h_num;i++)
	{
		for(int D=0;D<DIMENSION;D++)	p_differential_p2[D]=0;
		for(int j=0;j<h_num;j++)
		{							
			for(int D=0;D<DIMENSION;D++)	p_differential_p[D]=HYPER[j].stress[D][0]*HYPER1[j*h_num+i].DgDq[0]+HYPER[j].stress[D][1]*HYPER1[j*h_num+i].DgDq[1]+HYPER[j].stress[D][2]*HYPER1[j*h_num+i].DgDq[2];
			for(int D=0;D<DIMENSION;D++)	p_differential_p2[D]+=p_differential_p[D];				
		}		
		for(int D=0;D<DIMENSION;D++)	HYPER[i].differential_p[D]=HYPER[i].half_p[D]+Dt/2*(p_differential_p2[D]-9.8*mi*vec_g[D]);
	}
	cout<<"----------OK"<<endl;
}

void renew_lambda(mpsconfig &CON,vector<hyperelastic> &HYPER,vector<hyperelastic2> HYPER1,int hyper_number)
{
	cout<<"Lambda計算";

	int h_num=hyper_number;
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

	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)
		{
			double N_left2=0;
			for(int k=0;k<h_num;k++)
			{
				double N_left=HYPER1[i*h_num+k].DgDq[0]*HYPER1[j*h_num+k].DgDq[0]+HYPER1[i*h_num+k].DgDq[1]*HYPER1[j*h_num+k].DgDq[1]+HYPER1[i*h_num+k].DgDq[2]*HYPER1[j*h_num+k].DgDq[2];
				N_left2+=N_left;
			}
			N_Left[i*h_num+j]=Dt/2/mk*N_left2;
//			cout<<"N_Left["<<i<<"]["<<j<<"]= "<<Dt<<"/2/"<<mk<<"*"<<N_left2<<" = "<<N_Left[i*h_num+j]<<endl;
		}//jに関するfor文の終わり
	}//iに関するfor文の終わり

	for(int i=0;i<h_num;i++)
	{
		double N_right2=0;
		for(int k=0;k<h_num;k++)
		{
			double N_right=HYPER[k].differential_p[0]*HYPER1[i*h_num+k].DgDq[0]+HYPER[k].differential_p[1]*HYPER1[i*h_num+k].DgDq[1]+HYPER[k].differential_p[2]*HYPER1[i*h_num+k].DgDq[2];
			N_right2+=N_right;
		}
		N_Right[i]=1/mk*N_right2;	
//		cout<<"N_Right["<<i<<"]= "<<"1/"<<mk<<"*"<<N_right2<<" = "<<N_Right[i]<<endl;
	}
/*
	for(int i=0;i<h_num;i++)	cout<<"N_Right["<<i<<"]="<<N_Right[i]<<endl;
	cout<<endl;
	for(int i=0;i<h_num;i++)
	{
		for(int j=0;j<h_num;j++)	cout<<"N_Left["<<i<<"]["<<j<<"]="<<N_Left[i*h_num+j]<<endl;
		cout<<endl;
	}*/

	//lambdaを求める
	gauss(N_Left,N_Right,h_num);

	for(int i=0;i<h_num;i++)	HYPER[i].lambda=N_Right[i];
//	for(int i=0;i<h_num;i++)	cout<<"lambda["<<i<<"]="<<HYPER[i].lambda<<endl;

	delete [] N_Left;
	delete [] N_Right;

	cout<<"----------OK"<<endl;
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
void momentum_movie_AVS(mpsconfig &CON,int t,vector<mpselastic> PART,vector<hyperelastic> HYPER,int hyper_number)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	int h_num=hyper_number;
	double TIME=CON.get_step()*CON.get_dt();
	double le=CON.get_distancebp();
	int STEP=CON.get_step()/CON.get_interval()+1;		//出力する総ステップ数
	int step;

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
	fp<<h_num<<" "<<h_num<<endl;	//節点数と要素数出力
	
	//節点番号とその座標の出力 
	for(int i=0;i<h_num;i++) fp<<i<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	
	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=0;i<h_num;i++)
	{
		fp<<i<<"  0 pt "<<i<<endl;
	}

	//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
	fp<<"7 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
	fp<<"7 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	//fp<<"8 1 1 1 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	fp<<"p_x,"<<endl;
	fp<<"p_y,"<<endl;
	fp<<"p_z,"<<endl;
	fp<<"lambda,"<<endl;
	fp<<"ap_x,"<<endl;
	fp<<"ap_y,"<<endl;
	fp<<"ap_z,"<<endl;
	//fp<<"P,N/m^2"<<endl;
	//fp<<"value1,??"<<endl;

	//各節点の情報値入力
	for(int i=0;i<h_num;i++)
	{
		fp<<i<<" "<<HYPER[i].p[A_X]<<" "<<HYPER[i].p[A_Y]<<" "<<HYPER[i].p[A_Z]<<" "<<HYPER[i].lambda<<" "<<HYPER[i].ang_p[A_X]<<" "<<HYPER[i].ang_p[A_Y]<<" "<<HYPER[i].ang_p[A_Z]<<endl;
		//fp<<i<<" "<<NODE[i].depth<<" "<<NODE[i].L/le<<" "<<NODE[i].potential<<" "<<NODE[i].Fs<<" "<<endl;
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
	
	double y_min=(atan(-10.0)+PI/2)/PI;
	double y_max=(atan(10.0)+PI/2)/PI;

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
	double vec_norm[DIMENSION];
	vec_norm[A_X]=0;
	vec_norm[A_Y]=0;
	vec_norm[A_Z]=1;
	//for(int i=0;i<h_num;i++)	for(int D=0;D<DIMENSION;D++)	HYPER[h_num].vec_norm[D]=vec_norm[D];

	for(int i=0;i<h_num;i++)
	{
		double pnd=0;
		double max_dp[DIMENSION];
		for(int D=0;D<DIMENSION;D++)	max_dp[D]=-2*HYPER[i].p[D]*vec_norm[D];
		
		double dp_i[DIMENSION];
		dp_i[0]=0,	dp_i[1]=0,	dp_i[2]=0;

		for(int j=h_num;j<p_num;j++)
		{
			double vec_dis[DIMENSION];
			double dis=0;
			for(int D=0;D<DIMENSION;D++)	vec_dis[D]=PART[j].r[D]-PART[i].r[D];		
			double	in_vec=vec_dis[0]*vec_dis[0]+vec_dis[1]*vec_dis[1]+vec_dis[2]*vec_dis[2];
			dis=sqrt(in_vec);
			if(dis<r)
			{
				double w=kernel(r,dis);
				pnd+=w;
				
				//高さ方向の距離計算
				double dp_j[DIMENSION];
				double h=vec_norm[0]*vec_dis[0]+vec_norm[1]*vec_dis[1]+vec_norm[2]*vec_dis[2];
				double avs_h=fabs(h);
				double x=avs_h/r*20-10;
				double y=((atan(x)+PI/2)/PI-y_min)/(y_max-y_min);
				for(int D=0;D<DIMENSION;D++)
				{
					dp_j[D]=max_dp[D]*y;
					dp_i[D]+=dp_j[D]*w;
				}
			}			
		}
		if(pnd!=0 && HYPER[i].p[A_Z]<=0)	for(int D=0;D<DIMENSION;D++)	HYPER[i].p[D]+=dp_i[D]/pnd;
		s_pnd<<pnd<<",";
		dp<<i<<",";
		for(int D=0;D<DIMENSION;D++)	dp<<dp_i[D]<<",";
		dp<<endl;
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


void output_hyper_data(vector<mpselastic> PART,vector<hyperelastic> HYPER,vector<hyperelastic2> HYPER1,int t)
{
	int h_num=HYPER.size();

	if(t==1)
	{
		////計算した各定数の出力
		ofstream ai("Ai.csv");
		ofstream inai("inverse_Ai.csv");
		ofstream aiin("aiin.csv");
		ofstream n0("noij.csv");
		ofstream dg("initial_DgDq.csv");
		ofstream wiin("wiin.csv");

		aiin<<"粒子番号"<<","<<"周辺粒子番号"<<","<<"X変位"<<","<<"Y変位"<<","<<"Z変位"<<endl;
		n0<<"粒子番号"<<","<<"周辺粒子番号"<<","<<"X変位"<<","<<"Y変位"<<","<<"Z変位"<<endl;
		dg<<"粒子番号"<<","<<"周辺粒子番号"<<","<<"X変位"<<","<<"Y変位"<<","<<"Z変位"<<endl;

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

			aiin<<i<<",";
			n0<<i<<",";
			dg<<i<<",";
			wiin<<i<<",";	
			for(int j=0;j<h_num;j++)
			{
				aiin<<j<<",";
				n0<<j<<",";
				dg<<j<<",";
				aiin<<HYPER1[i*h_num+j].aiin[A_X]<<","<<HYPER1[i*h_num+j].aiin[A_Y]<<","<<HYPER1[i*h_num+j].aiin[A_Z]<<","<<endl;
				n0<<HYPER1[i*h_num+j].n0ij[A_X]<<","<<HYPER1[i*h_num+j].n0ij[A_Y]<<","<<HYPER1[i*h_num+j].n0ij[A_Z]<<","<<endl;
				dg<<HYPER1[i*h_num+j].DgDq[A_X]<<","<<HYPER1[i*h_num+j].DgDq[A_Y]<<","<<HYPER1[i*h_num+j].DgDq[A_Z]<<","<<endl;
				wiin<<HYPER1[i*h_num+j].wiin<<",";
			}
			wiin<<endl;
		}
		ai.close();
		inai.close();
		n0.close();
		dg.close();
		wiin.close();
	}



	ofstream d_p_X("d_P_X.csv", ios::app);
	ofstream d_p_Y("d_P_Y.csv", ios::app);
	ofstream d_p_Z("d_P_Z.csv", ios::app);

	ofstream h_p_X("h_P_X.csv", ios::app);
	ofstream h_p_Y("h_P_Y.csv", ios::app);
	ofstream h_p_Z("h_P_Z.csv", ios::app);

	ofstream lam("lambda.csv", ios::app);

	stringstream s2;
	s2<<"./Position/position"<<t<<".dat";
	ofstream po(s2.str());

	stringstream s3;
	s3<<"./Hy_stress/Stress"<<"t"<<t<<".csv";
	ofstream stress(s3.str());
	
	stringstream s5;
	s5<<"./Fi/differential_Fi "<<"t"<<t<<".csv";
	ofstream d_Fi(s5.str());

	stringstream s6;
	s6<<"./Fi/t_inverse_Fi "<<"t"<<t<<".csv";
	ofstream i_t_Fi(s6.str());

	ofstream p_X("P_X.csv", ios::app);
	ofstream p_Y("P_Y.csv", ios::app);
	ofstream p_Z("P_Z.csv", ios::app);
	ofstream p_c("P_corner.csv", ios::app);

	po<<"Position"<<t<<endl;		
	stress<<"Stress"<<t<<endl;

	double px=0,py=0,pz=0;

	if(t==1)
	{
		p_X<<"t"<<",";
		p_Y<<"t"<<",";
		p_Z<<"t"<<",";
		p_c<<"t"<<",";
		d_p_X<<"t"<<",";
		d_p_Y<<"t"<<",";
		d_p_Z<<"t"<<",";
		h_p_X<<"t"<<",";
		h_p_Y<<"t"<<",";
		h_p_Z<<"t"<<",";
		lam<<"t"<<",";

		p_c<<0<<","<<","<<","<<4<<","<<","<<","<<36<<","<<","<<","<<40<<endl;
		p_c<<","<<HYPER[0].p[A_X]<<","<<HYPER[0].p[A_Y]<<","<<HYPER[0].p[A_Z]<<",";
		p_c<<HYPER[4].p[A_X]<<","<<HYPER[4].p[A_Y]<<","<<HYPER[4].p[A_Z]<<",";
		p_c<<HYPER[36].p[A_X]<<","<<HYPER[36].p[A_Y]<<","<<HYPER[36].p[A_Z]<<",";
		p_c<<HYPER[40].p[A_X]<<","<<HYPER[40].p[A_Y]<<","<<HYPER[40].p[A_Z]<<",";
		
		for(int i=0;i<h_num;i++)
		{
			p_X<<i<<",";
			p_Y<<i<<",";
			p_Z<<i<<",";
			d_p_X<<i<<",";
			d_p_Y<<i<<",";
			d_p_Z<<i<<",";
			h_p_X<<i<<",";
			h_p_Y<<i<<",";
			h_p_Z<<i<<",";
			lam<<i<<",";
		}
		p_X<<endl;
		p_Y<<endl;
		p_Z<<endl;
		d_p_X<<endl;
		d_p_Y<<endl;
		d_p_Z<<endl;
		h_p_X<<endl;
		h_p_Y<<endl;
		h_p_Z<<endl;
		lam<<endl;
	}
	else
	{
		p_c<<t<<","<<HYPER[0].p[A_X]<<","<<HYPER[0].p[A_Y]<<","<<HYPER[0].p[A_Z]<<",";
		p_c<<HYPER[4].p[A_X]<<","<<HYPER[4].p[A_Y]<<","<<HYPER[4].p[A_Z]<<",";
		p_c<<HYPER[36].p[A_X]<<","<<HYPER[36].p[A_Y]<<","<<HYPER[36].p[A_Z]<<",";
		p_c<<HYPER[40].p[A_X]<<","<<HYPER[40].p[A_Y]<<","<<HYPER[40].p[A_Z]<<","<<endl;
	}

	p_X<<t<<",";
	p_Y<<t<<",";
	p_Z<<t<<",";
	d_p_X<<t<<",";
	d_p_Y<<t<<",";
	d_p_Z<<t<<",";
	h_p_X<<t<<",";
	h_p_Y<<t<<",";
	h_p_Z<<t<<",";
	lam<<t<<",";

	for(int i=0;i<h_num;i++)
	{
		
		px+=HYPER[i].p[A_X];
		py+=HYPER[i].p[A_Y];
		pz+=HYPER[i].p[A_Z];

		p_X<<HYPER[i].p[A_X]<<",";
		p_Y<<HYPER[i].p[A_Y]<<",";
		p_Z<<HYPER[i].p[A_Z]<<",";
		d_p_X<<HYPER[i].differential_p[A_X]<<",";
		d_p_Y<<HYPER[i].differential_p[A_Y]<<",";
		d_p_Z<<HYPER[i].differential_p[A_Z]<<",";
		h_p_X<<HYPER[i].half_p[A_X]<<",";
		h_p_Y<<HYPER[i].half_p[A_Y]<<",";
		h_p_Z<<HYPER[i].half_p[A_Z]<<",";

		po<<i<<"	";
		lam<<i<<","<<HYPER[i].lambda<<endl;
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
		po<<endl;
		stress<<endl;
		d_Fi<<endl<<endl;
		i_t_Fi<<endl<<endl;
	}

	p_X<<endl;
	p_Y<<endl;
	p_Z<<endl;
	d_p_X<<endl;
	d_p_Y<<endl;
	d_p_Z<<endl;
	h_p_X<<endl;
	h_p_Y<<endl;
	h_p_Z<<endl;
	lam<<endl;

	stress.close();
	po.close();

	lam.close();
	d_Fi.close();
	i_t_Fi.close();

	p_X.close();
	p_Y.close();
	p_Z.close();
	p_c.close();
	d_p_X.close();
	d_p_Y.close();
	d_p_Z.close();
	h_p_X.close();
	h_p_Y.close();
	h_p_Z.close();

	ofstream sum_p("sum_P.csv", ios::app);
	if(t==1)	sum_p<<"step"<<","<<"X"<<","<<"Y"<<","<<"Z"<<endl;
	sum_p<<t<<","<<px<<","<<py<<","<<pz<<endl;
	sum_p.close();
}

void output_newton_data1(double *fx, double *DfDx, double *n_rx, double *n_ry, double *n_rz,int hyper_number,int count, int t)
{
	int h_num=hyper_number;

	stringstream ss_r_X;
	ss_r_X<<"./Newton_raphson/r_X"<<t<<".dat";

	stringstream ss_r_Y;
	ss_r_Y<<"./Newton_raphson/r_Y"<<t<<".dat";

	stringstream ss_r_Z;
	ss_r_Z<<"./Newton_raphson/r_Z"<<t<<".dat";

	stringstream ss_fx;
	ss_fx<<"./Newton_raphson/fx "<<t<<".dat";
		
	if(count==1)
	{
		ofstream init0(ss_r_X.str(), ios::trunc);
		ofstream init1(ss_r_Y.str(), ios::trunc);
		ofstream init2(ss_r_Z.str(), ios::trunc);
		ofstream init3(ss_fx.str(), ios::trunc);
	
		init0.close();
		init1.close();
		init2.close();
		init3.close();
	}

	ofstream r_X(ss_r_X.str(), ios::app);
	ofstream r_Y(ss_r_Y.str(), ios::app);
	ofstream r_Z(ss_r_Z.str(), ios::app);
	ofstream sfx(ss_fx.str(), ios::app);

	stringstream ss_Df;
	ss_Df<<"./Newton_raphson/DfDx "<<"t"<<t<<" count"<<count<<".dat";
	ofstream Df(ss_Df.str());

	if(count==1)
	{
		r_X<<"反復回数"<<","<<"r_X"<<endl;
		r_Y<<"反復回数"<<","<<"r_Y"<<endl;
		r_Z<<"反復回数"<<","<<"r_Z"<<endl;
		sfx<<"反復回数"<<","<<"fx"<<endl;

		for(int i=0;i<h_num;i++)
		{
			r_X<<","<<i;
			r_Y<<","<<i;
			r_Z<<","<<i;
			sfx<<","<<i;
		}
		r_X<<endl;
		r_Y<<endl;
		r_Z<<endl;
		sfx<<endl;
	}
	
	r_X<<count;
	r_Y<<count;
	r_Z<<count;
	sfx<<count;

	Df<<"反復回数"<<","<<count<<endl;
	for(int i=0;i<h_num;i++)
	{
		r_X<<","<<n_rx[i];
		r_Y<<","<<n_ry[i];
		r_Z<<","<<n_rz[i];
		sfx<<","<<fx[i];
		for(int j=0;j<h_num;j++) Df<<DfDx[i*h_num+j]<<",";
		Df<<endl;

	}
	r_X<<endl;
	r_Y<<endl;
	r_Z<<endl;
	sfx<<endl;

	r_X.close();
	r_Y.close();
	r_Z.close();
	sfx.close();

	Df.close();
}
void output_newton_data2(double E, double *XX, int hyper_number, int count, int t)
{
	int h_num=hyper_number;
	stringstream ss_E;
	ss_E<<"./Newton_raphson/E"<<t<<".dat";
	
	stringstream ss_lam;
	ss_lam<<"./Newton_raphson/lambda"<<t<<".dat";
		
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

void output_energy(mpsconfig CON, vector<mpselastic> PART, vector<hyperelastic> HYPER, int hyper_number, int t)
{
	cout<<"弾性ポテンシャル計算";
	int h_num=hyper_number;
	int p_num=PART.size();

	double **d_Fi=new double *[DIMENSION];
	double **t_d_Fi=new double *[DIMENSION];
	for(int D=0;D<DIMENSION;D++)
	{
		d_Fi[D]=new double	[DIMENSION];
		t_d_Fi[D]=new double [DIMENSION];
	}

	double c10=CON.get_c10();
	double c01=CON.get_c01();
	double Ic=0, IIc=0;
	vector<double>	W;
	for(int i=0;i<h_num;i++)	W.push_back(0);

	for(int i=0;i<h_num;i++)
	{
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	d_Fi[D][D2]=HYPER[i].differential_Fi[D][D2];
		transpose(d_Fi,t_d_Fi);
		
		double dC[3][3];
		double dC2[3][3];
		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	dC[D][D2]=t_d_Fi[D][0]*d_Fi[0][D2]+t_d_Fi[D][1]*d_Fi[1][D2]+t_d_Fi[D][2]*d_Fi[2][D2];

		for(int D=0;D<DIMENSION;D++)	for(int D2=0;D2<DIMENSION;D2++)	dC2[D][D2]=dC[D][0]*dC[0][D2]+dC[D][1]*dC[1][D2]+dC[D][2]*dC[2][D2];

		double trace_dC=dC[0][0]+dC[1][1]+dC[2][2];
		double trace_dC2=dC2[0][0]+dC2[1][1]+dC2[2][2];

		Ic=trace_dC;
		IIc=0.5*(trace_dC*trace_dC-trace_dC2);
		W[i]=c10*(Ic-3)+c01*(IIc-3);
	}

	for(int D=0;D>DIMENSION;D++)
	{
		delete[]	d_Fi[D];
		delete[]	t_d_Fi[D];
	}
	delete[]	d_Fi;
	delete[]	t_d_Fi;

	cout<<"----------OK"<<endl;

	ofstream e("E.csv", ios::app);
	ofstream e_T("E_T.csv", ios::app);
	ofstream e_g("E_g.csv", ios::app);
	ofstream e_W("E_W.csv", ios::app);
	ofstream e_lam("E_lam.csv", ios::app);

	double V=get_volume(&CON);
	double mi=V*CON.get_hyper_density();

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
	for(int i=0;i<h_num;i++)
	{
		double vv=HYPER[i].p[0]*HYPER[i].p[0]+HYPER[i].p[1]*HYPER[i].p[1]+HYPER[i].p[2]*HYPER[i].p[2];
		double	energy=0.5/mi*vv+mi*9.8*PART[i].r[A_Z]+W[i]*V+HYPER[i].lambda*(1-HYPER[i].J)*V;
		e<<energy<<",";
		e_T<<0.5/mi*vv<<",";
		e_g<<mi*9.8*PART[i].r[A_Z]<<",";
		e_W<<W[i]*V<<",";
		e_lam<<HYPER[i].lambda*(1-HYPER[i].J)*V<<",";
	}
	e<<endl;
	e_T<<endl;
	e_g<<endl;
	e_W<<endl;
	e_lam<<endl;

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


