#include"stdafx.h"	
using namespace std;

//表面判定関数
//***********ハッシュ探索すれば早くなる？？？
//弾性体計算用にカスタマイズ（2012-11-05）
void freeon(elastic &ELAST, vector<mpselastic> &PART, double n0_4,int *INDEX,int **MESH, double *mindis, int t)
{
    //説明:freeon関数では各粒子の粒子数密度を求めるとともに、得られた密度から表面判定を行う。また、最低粒子間距離もついでにもとめている
	//freeon2より遅いが、そのぶん並列化が容易。
	double re=ELAST.get_re()*ELAST.get_distancebp();	//get_re()の値をそのまま使ってはダメ！！
	double mass=ELAST.get_particle_mass();
//	double volume=4*3.141592653*pow(re, 3)/3;
	int SIZE=ELAST.get_X_mesh()*ELAST.get_Y_mesh();
	
	//密度を求め、表面判定を行う。表面ならP=0にする
	//omp_set_num_threads(8);//スレッド数指定

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		PART[i].PND=0;
		PART[i].PND2=0;
		PART[i].N=0;
		PART[i].N2=0;
		PART[i].reset_current_distancebps();
		PART[i].reset_current_neighboursID();
		PART[i].reset_current_neighbours_position();
		PART[i].distance.clear();

		////粒子数密度測定
		double pnd=0.0;	//全粒子を考慮にいれた粒子数密度
		double pnd2=0;	//弾性体計算用粒子数密度
		double pnd4=0;	//表面判定用
		double sum=0;
		int N=0;		//周辺粒子数
		int N0=0;
			
		//空間に固有のID(INDEX, MESH)と粒子に固有のID(PART[i].index, PART[i].NEI[N])を照合して探索する
		for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
		{
			if(PART[i].index>=SIZE && PART[i].index<ELAST.get_number_of_mesh()-SIZE)
			{       
				for(int J=-1*ELAST.get_X_mesh();J<=ELAST.get_X_mesh();J+=ELAST.get_X_mesh())
				{
					for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
					{
						for(int L=0;L<INDEX[I+J+K];L++)
						{       
							int j=MESH[I+J+K][L];	//(I+J+K)番目の格子に含まれるL番目の粒子のIDを取得する
							if(j!=i)				//粒子IDが現在着目しているもの(i)と異なれば距離を計算する
							{ 
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
						
								if(dis<re)	
								{
									pnd+=kernel(re, dis);	//計算した距離が影響半径内であれば重み関数を計算する
									PART[i].NEI[N]=j;		//iの影響半径に存在する周辺粒子としてNEI[N]にIDを登録する。NEI[N]に格納されるIDはどのような順になる？
									PART[i].distance.insert(pair<int, double>(j, dis)); //周辺粒子との距離
									PART[i].set_current_neighboursID(j);	//最も重要
									PART[i].set_current_neighbours_position(X, Y, Z);
									PART[i].set_current_distancebps(dis);	//これいらない・・・
									N++;
									
								}
							}
						}
					}
				}	
			}
		}
		PART[i].PND=pnd;	//iの影響半径内に含まれる粒子数密度を取得
		PART[i].N=N;		//iの影響半径内に含まれる周辺粒子数を取得

//			PART[i].set_density(N*mass/volume);

//			cout<<"mass: "<<mass<<", N: "<<N<<", volume: "<<volume<<endl;
//			cout<<"density_i: "<<PART[i].get_density()<<", compensated density: "<<ELAST.get_density()<<endl;

	}//全粒子の表面判定が終了

//最小粒子間距離をもとめる

	if(ELAST.get_symp_flag()==ON)
	{
		double min=ELAST.get_distancebp();//最小粒子間距離
		int type1=0;
		int type2=0;
		double X1,Y1,Z1;

		//差し当たってこれいらない・・・コメントアウト
		//全ての粒子を総当りで検索している・・・bsearchを使って効率化する
		for(int i=0;i<PART.size();i++)
		{
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				if(dis<min)
				{
					type1=PART[i].type;
					type2=PART[j].type;
					min=dis;
					X1=PART[i].r[A_X];
					Y1=PART[i].r[A_Y];
					Z1=PART[i].r[A_Z];
				}
			}
		}//最小粒子間距離がもとまった
		cout<<"mindis="<<min<<" between "<<type1<<" & "<<type2<<endl;
		//cout<<"(x,y)="<<X1<<" , "<<Y1<<endl;
		*mindis=min;
	}
}

void freeon(mpsconfig &CON, vector<mpselastic> &PART, double n0_4,int *INDEX,int **MESH,double *mindis, int t)
{
    //説明:freeon関数では各粒子の粒子数密度を求めるとともに、得られた密度から表面判定を行う。また、最低粒子間距離もついでにもとめている
	//freeon2より遅いが、そのぶん並列化が容易。
	double re=CON.get_re()*CON.get_distancebp();	//get_re()の値をそのまま使ってはダメ！！
	double mass=CON.get_particle_mass();
//	double volume=4*3.141592653*pow(re, 3)/3;
	int SIZE=CON.get_X_mesh()*CON.get_Y_mesh();
	
	//密度を求め、表面判定を行う。表面ならP=0にする
	//omp_set_num_threads(8);//スレッド数指定

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		PART[i].PND=0;
		PART[i].PND2=0;
		PART[i].N=0;
		PART[i].N2=0;
		PART[i].reset_current_distancebps();
		PART[i].reset_current_neighboursID();
		PART[i].reset_current_neighbours_position();
		PART[i].distance.clear();

		////粒子数密度測定
		double pnd=0;	//全粒子を考慮にいれた粒子数密度
		double pnd2=0;	//弾性体計算用粒子数密度
		double pnd4=0;	//表面判定用
		double sum=0;
		int N=0;		//周辺粒子数
			
		//空間に固有のID(INDEX, MESH)と粒子に固有のID(PART[i].index, PART[i].NEI[N])を照合して探索する
		for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
		{
			if(PART[i].index>=SIZE && PART[i].index<CON.get_number_of_mesh()-SIZE)
			{       
				for(int J=-1*CON.get_X_mesh();J<=CON.get_X_mesh();J+=CON.get_X_mesh())
				{
					for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
					{
						for(int L=0;L<INDEX[I+J+K];L++)
						{       
							int j=MESH[I+J+K][L];	//(I+J+K)番目の格子に含まれるL番目の粒子のIDを取得する
							if(j!=i)				//粒子IDが現在着目しているもの(i)と異なれば距離を計算する
							{ 
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
						
								if(dis<re)	
								{
									if(PART[j].type!=WALL){	//壁粒子は初期周辺粒子に加えない
									pnd+=kernel(re,dis);	//計算した距離が影響半径内であれば重み関数を計算する
									PART[i].NEI[N]=j;		//iの影響半径に存在する周辺粒子としてNEI[N]にIDを登録する。NEI[N]に格納されるIDはどのような順になる？
									PART[i].distance.insert(pair<int, double>(j, dis)); //周辺粒子との距離
									PART[i].set_current_neighboursID(j);	//最も重要
									PART[i].set_current_neighbours_position(X, Y, Z);
									PART[i].set_current_distancebps(dis);	//これいらない・・・
									N++;
									}
									else {
									pnd+=kernel(re,dis);	//計算した距離が影響半径内であれば重み関数を計算する
									PART[i].NEI[N]=j;		//iの影響半径に存在する周辺粒子としてNEI[N]にIDを登録する。NEI[N]に格納されるIDはどのような順になる？
									PART[i].distance.insert(pair<int, double>(j, dis)); //周辺粒子との距離
									PART[i].set_current_neighboursID(j);	//最も重要
									PART[i].set_current_neighbours_position(X, Y, Z);
									PART[i].set_current_distancebps(dis);	//これいらない・・・
									N++;		//壁粒子は初期周辺粒子に加えない
									}//*/
								}
							}
						}
					}
				}	
			}
		}
		PART[i].PND=pnd;	//iの影響半径内に含まれる粒子数密度を取得
		PART[i].N=N;		//iの影響半径内に含まれる周辺粒子数を取得


	}//全粒子の表面判定が終了

//最小粒子間距離をもとめる

	double min=CON.get_distancebp();//最小粒子間距離
	int type1=0,type2=0;
	double X1,Y1,Z1;

	//差し当たってこれいらない・・・コメントアウト
	//全ての粒子を総当りで検索している・・・bsearchを使って効率化する
	for(int i=0;i<PART.size();i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<min)
			{
				type1=PART[i].type;
				type2=PART[j].type;
				min=dis;
				X1=PART[i].r[A_X];
				Y1=PART[i].r[A_Y];
				Z1=PART[i].r[A_Z];
			}
		}
	}//最小粒子間距離がもとまった
	cout<<"mindis="<<min<<" between "<<type1<<" & "<<type2<<endl;
	//cout<<"(x,y)="<<X1<<" , "<<Y1<<endl;
	*mindis=min;
	
}

///表面判定関数ver.2
void freeon2(mpsconfig &CON,vector<mpselastic> &PART,int particle_number,double n0_4,int *INDEX,int **MESH,double *mindis,int fluid_number,int out)
{
	///freeon2関数の説明:flag1[i]の導入により高速化。ただし、1CPUなら早いが、ﾏﾙﾁCPUによる並列化は厳しい
	double le=CON.get_distancebp();
	int SIZE=CON.get_X_mesh()*CON.get_Y_mesh();
	double d=2;
	if(CON.get_dimension()==3) d=3;

	if(CON.get_T_field()==ON && CON.get_insulate()==1)
	{
		out=particle_number;	//非断熱の温度場を計算する際は、OUTWALLも周辺粒子を格納しておく必要がある。よってout=particle_numberとすることでこれに対処
	}

	double *PND4=new double [out];//表面判定用粒子数密度
	int *flag1=new int [out];		//検査フラグ。0:未検査　1:検査済み
	///初期化
	#pragma omp parallel for
	for(int i=0;i<out;i++)
	{
		PART[i].PND=0;//初期化
	    PART[i].PND2=0;
	    PART[i].N=0;
	    PART[i].N2=0;
	    PART[i].N3=0;

		PND4[i]=0;
		flag1[i]=0;
	}
	
	//密度を求め、表面判定を行う。表面ならP=0にする
	for(int i=0;i<out;i++)//OUTWALL以外のCFD粒子。//OUTWALLの粒子数密度などはいらない
	{    
		////粒子数密度測定
		if(PART[i].index>=SIZE && PART[i].index<CON.get_number_of_mesh()-SIZE)
		{       
			for(int I=PART[i].index-1;I<=PART[i].index+1;I++)
			{       
				for(int J=-1*CON.get_X_mesh();J<=CON.get_X_mesh();J+=CON.get_X_mesh())
				{
					for(int K=-1*SIZE;K<=SIZE;K+=SIZE)
					{
						for(int L=0;L<INDEX[I+J+K];L++)
						{       
							int j=MESH[I+J+K][L];
							if(j<out)
							{
								if(flag1[j]==0 && j!=i)//まだ検査してないなら
								{
									double X=PART[j].r[A_X]-PART[i].r[A_X];
									double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
									double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
									double dis=sqrt(X*X+Y*Y+Z*Z);
							
									if(dis<=CON.get_re()*le)//勾配・発散
									{       
										double r=CON.get_re()*le;
										double w=kernel(r,dis);
										PART[i].PND+=w;
										PART[j].PND+=w;
										PART[i].NEI[PART[i].N]=j;
										PART[j].NEI[PART[j].N]=i;
										PART[i].N++;
										PART[j].N++;
									}
									if(dis<=CON.get_re2()*le)
									{       
										double r=CON.get_re2()*le;
										double w=kernel(r,dis);
							
										PART[i].PND2+=w;
										PART[j].PND2+=w;
										
										PART[i].NEI2[PART[i].N2]=j;
										PART[j].NEI2[PART[j].N2]=i;
										PART[i].N2++;
										PART[j].N2++;
									}
									if(dis<=CON.get_re3()*le)//表面張力re3
									{       
										PART[i].NEI3[PART[i].N3]=j;
										PART[j].NEI3[PART[j].N3]=i;
										PART[i].N3++;
										PART[j].N3++;
									}
									if(dis<=CON.get_re4()*le)
									{       
									    double r=CON.get_re4()*le;
										double w=kernel(r,dis);
										PND4[i]+=w;
										PND4[j]+=w;
									}
								}
							}
							else
							{
								double X=PART[j].r[A_X]-PART[i].r[A_X];
								double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
								double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
								double dis=sqrt(X*X+Y*Y+Z*Z);
							
								if(dis<=CON.get_re()*le)//勾配・発散
								{       
									double r=CON.get_re()*le;
									double w=kernel(r,dis);
									PART[i].PND+=w;
									PART[i].NEI[PART[i].N]=j;
									PART[i].N++;
								}
								if(dis<=CON.get_re2()*le)
								{       
									double r=CON.get_re2()*le;
									double w=kernel(r,dis);
									PART[i].PND2+=w;
									PART[i].NEI2[PART[i].N2]=j;
									PART[i].N2++;
								}
								if(dis<=CON.get_re3()*le)//表面張力re3
								{       
									PART[i].NEI3[PART[i].N3]=j;
									PART[i].N3++;
								}
								if(dis<=CON.get_re4()*le)
								{       
								    double r=CON.get_re4()*le;
									double w=kernel(r,dis);
									PND4[i]+=w;
								}
							}
						}
					}
				}
			}
		}
		flag1[i]=1;//検査終了

		if(PART[i].N3>450) cout<<"error in freeon over. Change more than "<<PART[i].N3<<" coods:"<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		////////////////////

		if(PART[i].type==ELASTIC)
//		if(PART[i].type==FLUID)
		{
			if(PND4[i]<n0_4*CON.get_beta())//β以下なら
			{
				PART[i].surface=ON;//表面粒子とする
				PART[i].P=0;//境界条件？
			}
			else PART[i].surface=OFF;
		}
		else if(PART[i].type==WALL)
//		else if(PART[i].type==INWALL)
		{
			if(PND4[i]<n0_4*CON.get_beta())//β以下なら
			{
				PART[i].surface=ON;//壁表面粒子とする
				PART[i].P=0;
			}
			else  PART[i].surface=OFF;
		}
		else if(PART[i].type==SOLID)
//		else if(PART[i].type==INWALL)
		{
			if(PND4[i]<n0_4*CON.get_beta())//β以下なら
			{
				PART[i].surface=ON;//壁表面粒子とする
				PART[i].P=0;
			}
			else  PART[i].surface=OFF;
		}
	}
	for(int i=out;i<particle_number;i++)
	{
		PART[i].N=0;
		PART[i].N2=0;
		PART[i].N3=0;
	}

	////最低粒子間距離をもとめる
	double min0=CON.get_distancebp();//最低粒子間距離
	int type1,type2,surface1,surface2;
	double X1,Y1,Z1;
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			if(dis<min0)
			{
				type1=PART[i].type;
				surface1=PART[i].surface;
				type2=PART[j].type;
				surface2=PART[j].surface;
				min0=dis;
				X1=PART[i].r[A_X];
				Y1=PART[i].r[A_Y];
				Z1=PART[i].r[A_Z];
			}
		}
	}/////最低粒子間距離がもとまった
	cout<<"mindis="<<min0<<" between "<<type1<<" "<<surface1<<" & "<<type2<<" "<<surface2<<endl;
	/////////*/
	*mindis=min0;

	delete [] PND4;
	delete [] flag1;
}

//表面判定関数ver.2
void surface_judge2(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number)
{
	//従来の表面判定に加えて、ここでさらにふるいにかける。
	//粒子iの法線ベクトルと、粒子j方向ﾍﾞｸﾄﾙとの角度がある値を超えていたら表面ではないと判断
	double *direct[DIMENSION];
    for(int D=0;D<DIMENSION;D++) direct[D]=new double [fluid_number];
		
	
    //////法線ﾍﾞｸﾄﾙ計算
    for(int i=0;i<fluid_number;i++)
    {
        if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
        else  for(int D=0;D<3;D++) direct[D][i]=0;
	}

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)//この表面粒子が本当に表面粒子足りうるか判定する
		{
			int flag=ON;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
					double dis=sqrt(X*X+Y*Y+Z*Z);
					double nx=X/dis;//i粒子からj粒子へと向かう単位ﾍﾞｸﾄﾙ
					double ny=Y/dis;
					double nz=Z/dis;
					double inp=nx*direct[A_X][i]+ny*direct[A_Y][i]+nz*direct[A_Z][i];//法線ﾍﾞｸﾄﾙと単位ﾍﾞｸﾄﾙとの内積
					if(inp<-0.5) flag=OFF;//内積が-0.5のものがひとつでもあれば内部流体粒子。内積が-0.5ということは、その角度が120度が許容範囲ということ
				}
			}
			PART[i].surface=flag;
		}
	}

	for(int D=0;D<DIMENSION;D++) delete [] direct[D];

}

//freeon関数 ver.3 粒子数密度のみ再計算。粒子―粒子関係は変化しないと仮定している
void freeon3(mpsconfig &CON,vector<mpselastic> &PART,int particle_number,int out)
{
	double d=2;
	if(CON.get_dimension()==3) d=3;
	double R=CON.get_re()*CON.get_distancebp();
	double R2=CON.get_re2()*CON.get_distancebp();
	for(int i=0;i<out;i++)
	{
		double W=0;
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			double X=PART[j].r[A_X]-PART[i].r[A_X];
			double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
			double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			double dis=sqrt(X*X+Y*Y+Z*Z);
			double w=kernel(R,dis);
			//double w=kernel2(R,dis,d);
			W+=w;
		}
		PART[i].PND=W;
	}
	if(CON.get_re()!=CON.get_re2())
	{
		for(int i=0;i<out;i++)
		{
			double W=0;
			for(int k=0;k<PART[i].N2;k++)
			{
				int j=PART[i].NEI2[k];
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
				double dis=sqrt(X*X+Y*Y+Z*Z);
				double w=kernel(R2,dis);
				//double w=kernel2(R2,dis,d);
				W+=w;
			}
			PART[i].PND2=W;
		}
	}
	else for(int i=0;i<out;i++) PART[i].PND2=PART[i].PND;
}