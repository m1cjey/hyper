#include "stdafx.h"		//主要なヘッダーファイルはまとめてこのなか。
//only main function
#include "Model.h"
#include "Check.h"
//#include "Rigidbody.h"
#include "Micro_AVS.h"
void file_initialization();
void Make_STL();//STL用足し合わせ関数
//////////////////////
int _tmain(int argc, _TCHAR* argv[])
{

//メモリリーク検出
_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
/*_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG);
_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);*/
_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_DELAY_FREE_MEM_DF | _CRTDBG_CHECK_CRT_DF);



	mpsconfig CON;				//解析条件格納クラス
	//CON.Set_initial_config();	//解析条件入力


	clock_t t1=clock();	//経過時間を秒で表現するには、CLOCKS_PER_SECで割る
	
    int particle_number=0;	//全粒子数
	int fluid_number=0;		//流体粒子数（弾性体の場合は弾性体粒子数）
	int out;				//fluid_number<=i<outがINWALL群で、out<=iがOUTWALL群になる。 なってないなら自動であとで並び替える
	int count;				//カウント用変数
	int order_sw=OFF;		//粒子を並び替える必要があるかどうかのフラグ
	int t=0;				//ステップ数
	int dimension=CON.get_dimension();
	int model_set=CON.get_model_set_way();
    double dt=CON.get_dt();
    int STEP=CON.get_step();
    double g[DIMENSION]={0,0,0};
    if(dimension==2) g[A_Y]=CON.get_g();
    if(dimension==3) g[A_Z]=CON.get_g();

    //MPSにおける定数計算
	double n0=initial_pnd(CON.get_re(), dimension, model_set);		//初期粒子密度n0
    double n0_4=initial_pnd(CON.get_re4(), dimension, model_set);	//freeon用初期粒子密度
    double TIME=0; //解析時間
	double Umax=0; //最大速度
	double mindis; //CFLの最低粒子間距離 minimum distance


	//超弾性
	int hyper_number=0;
	int magnetic_number=0;
	double Dt=CON.get_dt()*CON.get_interval();

	cout<<"初期粒子数密度 n0= "<<setprecision(10)<<n0<<endl;

	
	//FEM用のclass作成
	vector<point3D> NODE_FEM3D;
	vector<element3D> ELEM_FEM3D;
	int node_FEM3D=0;
	int nelm_FEM3D=0;	//全節点数,要素数

	//初期粒子配置書き込み　restart=ONの場合は粒子数読み込み
	initial_model_input(&CON, &particle_number, &TIME);
//	Model model;
//	particle_number=model.Model_set();

	

	//INDEX関係
    int *INDEX=new int[CON.get_number_of_mesh()];	//各格子に含まれる粒子数を格納(格子の数だけ配列が必要)
    cout<<"X_mesh="<<CON.get_X_mesh()<<" Y_mesh="<<CON.get_Y_mesh()<<" Z_mesh="<<CON.get_Z_mesh()<<endl;
    cout<<"number_of_mesh="<<CON.get_number_of_mesh()<<endl;
	Check check;
	check.Out_put_config();
	//vectorで確保するとアルゴリズムが利用できる
	//初期化・・・push_back()ではコンストラクタは動かない！！（コピーコンストラクタが動く！！）angは不定！！
	//これ以降はparticle_number=PART.size()でよい
	mpselastic PART0;
	vector<mpselastic> PART1;	//粒子配列をtyoe毎に並べ替える。一時保管
	vector<mpselastic> PART;

	cout<<"ベクトル作成完了"<<endl;

	cout<<"1"<<endl;
//	PART.reserve(20000);
	for(int i=0;i<particle_number;i++){
		PART1.push_back(PART0);
		PART.push_back(PART0);
	}



	////////////STL/////////
//	Make_STL();
	///////////////////////

	//粒子データをファイルから読み取り
	//PARTの初期値はファイルから読み込む・・・コンストラクタは動かない！！
	//ここで粒子番号の整理をしている MRE→シリコーン→その他の順番
	input_particle_data(&CON, PART1,PART, 1);	//最後の引数に1を渡したらinitialデータを読み取り、それ以外なら前ステップデータを読み取る
	Micro_AVS avs;
	for(int i=0;i<PART.size();i++){
		avs.make_list(PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z],0,0,0);
	}
	avs.Output_mgf_MicroAVS("check_particle",1);



	//各粒子数をカウント or 並び替え．※ほとんど意味がない．実行しなくてもOK(2012/02/21)
	calc_numbers_of_particles_and_change_the_order(&CON, PART, &fluid_number,&hyper_number,&magnetic_number, &out, &order_sw);


	check_initial_position(&CON, PART);


	//FEM3D
	double **F=new double*[DIMENSION];//F[D][i]となっているがこれでOK??
	double **old_F=new double*[DIMENSION];//F[D][i]となっているがこれでOK??

	for(int D=0;D<DIMENSION;D++)
	{
		F[D]=new double [(unsigned)PART.size()]; //各粒子に働く電磁力//particle_numberはinitial_model_input()で求まる．
		for(int i=0;i<PART.size();i++)	F[D][i]=0.0; //初期化
		old_F[D]=new double [(unsigned)PART.size()]; 
		for(int i=0;i<PART.size();i++)	old_F[D][i]=0.0; //初期化
	}

	//陽解析の前にreloadINDEX	
	//粒子ID更新・粒子数密度更新
	reload_INDEX(CON,PART, INDEX);//格子内の粒子数更新
	
	count=0;
	int **MESH0=new int *[CON.get_number_of_mesh()];
	for(int i=0;i<CON.get_number_of_mesh();i++)
	{       
		count+=INDEX[i];
		MESH0[i]=new int [INDEX[i]];
	}
	reload_INDEX2(&CON, PART, MESH0);

	//表面判定（弾性計算の場合最初だけやればOK！）

	freeon(CON, PART, n0_4, INDEX, MESH0, &mindis, t);
	

	for(int i=0;i<CON.get_number_of_mesh();i++) delete [] MESH0[i];
	delete [] MESH0;

	//初期化も同時に行うのでこの位置(ePND[i]=PART[i].PNDなのでこの位置でなければ)
	elastic ELAST(PART);
	for(int i=0;i<PART.size();i++) PART[i].initialize_particles(ELAST, t);
	//ELAST.set_ground_position(CON.get_ground_position());

	//HYPERクラスに粒子情報格納
	hyperelastic HYPER0;
	hyperelastic2 HYPER10;
	vector<hyperelastic> HYPER;
	vector<hyperelastic2> HYPER1;

		for(int i=0;i<hyper_number*hyper_number;i++)	HYPER1.push_back(HYPER10);
		for(int i=0;i<hyper_number;i++)	HYPER.push_back(HYPER0);
	//初期粒子配置で解析領域外に粒子がないかチェック


	//プリプロセス終了


	//file initialization
	file_initialization();
	/////////////////////
	ofstream pt("PART_model.dat");
	for(int i=0;i<PART.size();i++)	pt<<PART[i].r[A_Y]<<PART[i].r[A_Z]<<endl;
	pt.close();

	double L=0;
	double W=100;
	double P=0;
	int m=0;
	int wait=0;	//待ちステップ
	double limP=0;
	bool ff=0;

	
	//表面だけ表示//表面初期化
	/*
	if(bool cat=ON)
	{
		for(int i=0 ;i<PART.size();i++)
		{
			PART[i].surface=OFF;	
			if(PART[i].PND<=18)
			{
				PART[i].surface=ON;
			}
		}
	}*/

	//表面粒子の存在確認
/*	count_suf=0;
	for(int i=0 ;i<PART.size();i++)
	{
		if(PART[i].surface==ON) count_suf++;
	}
	cout<<count_suf<<endl;*/


	//計算スタート
	for(t=1;t<=STEP;t++)
	{	
		cout<<"陽解析 start:step="<<t<<", 粒子数="<<particle_number<<", dt="<<dt<<endl;

		//陽解析の前にreloadINDEX（粒子は移動するので毎ステップ実行する必要がある）
		reload_INDEX(CON, PART, INDEX);	//格子内の粒子数更新
		cout<<"reload_INDEX_OK"<<endl;

		//MESHはmpsconfigのメンバ関数にするのが良い。new/delete危険なのでshared_ptrかvectorを使うべき
		int **MESH = new int *[CON.get_number_of_mesh()];	//get_number_of_mesh() {return (int)((maxX-minX)/(distancebp*dx)*(maxY-minY)/(distancebp*dx)*(maxZ-minZ)/(distancebp*dx)+0.001);}//格子数：X_mesh*Y_mesh*Z_mesh
		count=0;
		for(int i=0;i<CON.get_number_of_mesh();i++)
		{       
			count+=INDEX[i];
			MESH[i]=new int [INDEX[i]];
		}
		if(count>PART.size()) cout<<"INDEX error 粒子数増加?"<<endl;
		if(count<PART.size()) cout<<"INDEX error 粒子数減少?"<<endl;
		//reload_INDEX2(&CON, PART, MESH);

		//陽解析
		unsigned timeA=GetTickCount();


		//影響半径内の粒子をカウント＋mindisの計算
		//弾性計算用の関数を使う
		cout<<"freeon_start"<<endl;
		if(CON.get_flag_ELAST()==ON)	freeon(ELAST, PART, n0_4, INDEX, MESH, &mindis, t); //表面判定

		cout<<"陽解析前の粒子依存関係計算終了　－－time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

		//CFL条件・拡散数条件によるdtの決定
		//シンプレティックスキームを使う場合は時間解像度可変は適用してはならない
		//if(OFF==ELAST.get_symp_flag()) courant_elastic(&CON, PART, fluid_number, t, &dt, mindis, Umax, g);
		//courant_elastic(&CON, PART, fluid_number, t, &dt, mindis, Umax, g);
		//弾性体計算が落ち着いてからdtを変更する

		//model1だと1311、model7だと731、model11だと1781でICCG法で誤差の計算値が発散し、解析が破綻する



		if(CON.get_FEM_flag()==true && t>wait)
		{
			///////HYPERの場合差分をFは足す & 初期化
			if(CON.get_flag_HYPER()==ON)
			{
				if(t==1)
				{
					for(int i=0;i<hyper_number;i++)
					{
						for(int D=0;D<DIMENSION;D++)
						{
							HYPER[i].p[D]=0;
							old_F[D][i]=0;
						}
					}
				}
				else{for(int i=0;i<hyper_number;i++)for(int D=0;D<DIMENSION;D++)old_F[D][i]=F[D][i];}
			}

			//ELASTからCONの関数呼び出せるので二系統要らない・・・
			if(ELAST.get_FEM_switch()==true && CON.get_mesh_input()!=2)
			{
				CON.change_step_size(); //こんな関数は混乱を招くだけなので要らない
				ELAST.change_step_size();	//必要？15/2/3
				if(t==1 || (t-1)%CON.get_EM_interval()==0)
				{
					//自作デローニ分割
					FEM3D(CON, PART, F, &node_FEM3D, &nelm_FEM3D, NODE_FEM3D, ELEM_FEM3D, t, TIME, fluid_number);	//HYPERELASTに対応できていない
				}
			}
			else if(CON.get_mesh_input()==2)
			{
				if(t==1 || (t-1)%CON.get_EM_interval()==0)
				{	
					//TetGenによるメッシュ分割
					TetGenInterface(CON, PART, F, fluid_number, dt, t, particle_number, n0, TIME);
				}	
			}

			////////超弾性体計算の場合は、ローレンツ力の変化分を運動量に追加する15/05/23
			if(CON.get_flag_HYPER()==ON)
			{
				double det_F[DIMENSION];
				for(int i=0;i<hyper_number;i++)
				{
					for(int D=0;D<DIMENSION;D++)
					{
						det_F[D]=0;
						det_F[D]=F[D][i]-old_F[D][i];
						HYPER[i].p[D]+=Dt*det_F[D];
					}
				}
			}
		}
		

		////粒子移動計算
		if(CON.get_flag_ELAST()==ON)	calc_elastic(PART, ELAST, t, F);	//if分の追加15/2/10

		cout<<"hyper_number="<<hyper_number<<endl;

		//超弾性計算
		if(CON.get_flag_HYPER()==ON)	calc_hyper(CON,PART,HYPER,HYPER1,hyper_number,t);//if分の追加15/2/10

		cout<<"hyper_calculation is ended."<<endl;


		cout<<"陽解析終了"<<endl;
//		cout<<"陽解析終了 umax="<<sqrt(Umax)<<"  limit U="<<0.2*mindis/dt<<endl;

		//粒子が動いたのでINDEX更新
		for(int i=0;i<CON.get_number_of_mesh();i++) delete [] MESH[i];
		delete [] MESH;

		double umax2=0.0;//最大速度
		for(int i=0;i<PART.size();i++)
		{ 
			double speed=0.0;//粒子速度
			for(int D=0;D<DIMENSION;D++) speed+=PART[i].u[D]*PART[i].u[D];
			if(speed>umax2) umax2=speed;
				
		}
		cout<<"umax="<<sqrt(umax2)<<endl;

		if(t>=wait){
		TIME+=dt;///時間更新
		CON.set_current_time(TIME);
		}
		cout<<"解析物理時間="<<TIME<<"/ "<<(CON.get_step()*CON.get_dt())<<endl;

		

		//ポスト処理：
		post_processing(CON, PART, ELAST, particle_number, particle_number, dt, Umax, t, TIME,F); //各物理量出力＆クーラン数によるdt改変&microAVS出力
		post_processing3(CON, PART, particle_number, particle_number, t, TIME); //restart用ファイル生成
//		order_sw=check_position(&CON, PART, particle_number, &particle_number); //領域外の粒子を検査 //領域外粒子を検知すれば、order_sw=ONになる
		clock_t t3=clock();
		cout<<"CPU time="<<(t3-t1)/CLOCKS_PER_SEC<<"[sec]"<<endl;
		ofstream t_log("time_log.dat", ios::app);

		t_log<<"step="<<t<<", time="<<(t3-t1)/CLOCKS_PER_SEC<<"[sec]"<<endl;
		t_log.close();

		check.Courant_condition(PART);

		cout<<endl;
		
	}
	//ループ終了

	for(int D=0;D<DIMENSION;D++){
		delete[] F[D];
		delete[]old_F[D];
//		delete [] laplacian[D];
	}

	delete [] INDEX;
	delete [] F;
	delete[]old_F;
	clock_t t2=clock();
	cout<<"CPU last time="<<(t2-t1)/CLOCKS_PER_SEC<<"[sec]"<<endl;
//	MessageBeep(MB_ICONEXCLAMATION);//作業の終了を知らせるBEEP音 マングリングエラーが起こるのでコメントアウト

	return 0;
}

//重み関数　定義がおかしい！！
double kernel(double r,double dis)
{
	return r/dis-1;
//	return r/dis-1;
	//return r*r/(dis*dis)-1;
	//return r*r*r/(dis*dis*dis)-1;
	//return (1-dis/r)*(1-dis/r);
	//return 1;
}

//重み関数２
double kernel2(double r,double dis,double d)
{
	return r/dis-1;
    //return r*r*r/(dis*dis*dis)-1;
	//return r*r*r*r/(dis*dis*dis*dis)-1;
	//return pow(r,d)/pow(dis,d);
}

//重み関数3　壁用
double kernel3(double r,double dis)
{
	double ndis=0;
	ndis=dis-(r*0.1);
	return (dis<r) ? 5*(r/ndis-1): 0;
}

//重み関数４　
double kernel4(double r,double dis)
{
	if(dis<r)	return	r/dis-1;
	else
	{
		return	0;
	}

}


//初期粒子密度の計算（90行）
double initial_pnd(double r,int dimension,int calc_type)
{
	int size = (int)(r+1);//計算領域
	double dis;//距離
	double pnd=0;
	int count=0;
	if(dimension==2)
	{
		if(calc_type==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis=sqrt((double)(i*i+j*j));
					if(dis!=0 && dis<=r )
					{
						pnd+=kernel(r,dis);
						count++;
					}			
				}
			}
		}
		if(calc_type==1)				//初期配置として細密六法をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
					dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis<=r )
					{
						pnd+=kernel(r,dis);
						count++;
					}			
				}
			}
		}
	}
	else if(dimension==3)
	{
		if(calc_type==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis=sqrt((double)(i*i+j*j+k*k));
						if(dis!=0 && dis<=r )
						{
							pnd+=kernel(r,dis);
							count++;
						}
					}			
				}
			}
		}
		if(calc_type==1)				//初期配置として細密六法をとった場合
		{
			for(int i=-2*size;i<=2*size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					for(int k=-2*size;k<=2*size;k++)
					{
						double ii=(double)i;
						double jj=j*sqrt(3.0)/2;
						double kk=k*sqrt(2.0/3);
						if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//kが奇数ならiiとjjをずらす
						dis=sqrt(ii*ii+jj*jj+kk*kk);
						if(dis!=0 && dis<=r )
						{
							pnd+=kernel(r,dis);
							count++;
						}
					}			
				}
			}
		}
	}
	cout<<"n0のcount="<<count<<endl;
	return pnd;  
}

//ラプラシアン用変数λ計算関数（108行）
double calclambda(mpsconfig &CON)
{
	
	int dimension=CON.get_dimension();	//解析次元
	int Ini_place=CON.get_model_set_way();	//初期粒子配置方法　0=正方 1=細密
	double R=CON.get_re2();			//ラプラシアン用影響半径
	int size = (int)(R+1);//計算領域
	int count=0;
	double dis;//距離
	double w;
	double pnd=0;
	double lam=0;
	if(dimension==2)
	{
		if(Ini_place==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					dis=sqrt((double)(i*i+j*j));
					if(dis!=0 && dis<=R)
					{
						double length=dis*CON.get_distancebp();
						w=kernel(R,dis);
						pnd+=w;
						lam+=length*length*w;
						count++;
					}
				}				
			}
		}
		else if(Ini_place==1)				//初期配置として細密六法をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					double jj=j*sqrt(3.0)/2;
					double ii=(double)i;
					if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
					dis=sqrt(ii*ii+jj*jj);
					if(dis!=0 && dis<=R )
					{
						double length=dis*CON.get_distancebp();
						w=kernel(R,dis);
						pnd+=w;
						lam+=length*length*w;
						count++;
					}			
				}
			}
		}
	}
	else if(dimension==3)
	{
		if(Ini_place==0)				//初期配置として正方格子をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-size;j<=size;j++)
				{
					for(int k=-size;k<=size;k++)
					{
						dis=sqrt((double)(i*i+j*j+k*k));
						if(dis!=0 && dis<=R)
						{
							double length=dis*CON.get_distancebp();
							w=kernel(R,dis);
							pnd+=w;
							lam+=length*length*w;
							count++;
						}
					}
				}				
			}
		}
		else if(Ini_place==1)				//初期配置として細密六法をとった場合
		{
			for(int i=-size;i<=size;i++)
			{
				for(int j=-2*size;j<=2*size;j++)
				{
					for(int k=-2*size;k<=2*size;k++)
					{
						double ii=(double)i;
						double jj=j*sqrt(3.0)/2;
						double kk=k*sqrt(2.0/3);
						if(j%2!=0) ii+=0.5;//jが奇数ならiiを0.5格子だけずらす
						if(k%2!=0) {ii+=0.5; jj+=sqrt(3.0)/6;}//kが奇数ならiiとjjをずらす
						dis=sqrt(ii*ii+jj*jj+kk*kk);
						if(dis!=0 && dis<=R )
						{
							double length=dis*CON.get_distancebp();
							w=kernel(R,dis);
							pnd+=w;
							lam+=length*length*w;
							count++;
						}
					}			
				}
			}
		}
	}
	lam/=pnd;
	cout<<"λのcount="<<count<<endl;
	return lam;  
}

//粒子データ読み取り関数（52行）
void input_particle_data(mpsconfig *CON, vector<mpselastic> &PART1, vector<mpselastic> &PART, int t)
{
	int p=0;
	if(t==1)//最初はinitial_input.datから読み込み
	{
		ifstream fin("initial_input.dat");
		if(!fin) cout<<"\"initial_input.dat\" cannot be opened."<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
		for(int i=0;i<PART1.size();i++)
		{       
			fin>>PART1[i].ID;
			for(int D=0;D<DIMENSION;D++) fin>>PART1[i].r[D];
			for(int D=0;D<DIMENSION;D++) fin>>PART1[i].u[D];
			for(int D=0;D<DIMENSION;D++) fin>>PART1[i].PAcc[D];
			fin>>PART1[i].P;
			fin>>PART1[i].h;
			fin>>PART1[i].val;
			fin>>PART1[i].type;
			fin>>PART1[i].materialID;	//これはなんの意味がある？？
			fin>>PART1[i].surface;
			fin>>PART1[i].toFEM;
			PART1[i].dir_Pem=0;			//初期化
			PART1[i].dir_Pst=0;
			for(int D=0;D<3;D++) PART1[i].ang[D]=0.0;
			PART1[i].ang[3]=1.0;
			for(int D=0;D<3;D++) PART1[i].ang_u[D]=0.0;

		}
		fin.close();	
	}
	///////////////////////*/
	cout<<"粒子並べ替え開始"<<endl;
	if(t!=1) //2STEP以降はmps_input.datから読み込み
	{
		ifstream fin("restart_input.dat");
		if(!fin) cout<<"cannot open restart_input.dat"<<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
		for(int i=0;i<PART.size();i++)
		{
			fin>>PART[i].ID;
			for(int D=0;D<DIMENSION;D++) fin>>PART[i].r[D];
			for(int D=0;D<DIMENSION;D++) fin>>PART[i].u[D];
			for(int D=0;D<DIMENSION;D++) fin>>PART1[i].PAcc[D];
			fin>>PART[i].P;
			fin>>PART[i].h;
			fin>>PART[i].val;
			fin>>PART[i].type;
			fin>>PART[i].materialID;
			fin>>PART[i].surface;
			fin>>PART[i].toFEM;
			PART[i].dir_Pem=0;
			PART[i].dir_Pst=0;
	//		for(int D=0;D<DIMENSION;D++) PART[i].u[D]=0; //位置だけ欲しいので速度は０ 

		}
		fin.close();
	}
	bool num=false;
	////粒子並べ替え/////
	if(num==true){
	for(int i=0;i<PART1.size();i++)
	{
	 if(PART1[i].type==MAGELAST) 
	 {
		 swap(PART[p],PART1[i]);
		 p++;
	 }
	}
	for(int i=0;i<PART1.size();i++)
	{
	 if(PART1[i].type==MAGELAST2) 
	 {
		 swap(PART[p],PART1[i]);
		 p++;
	 }
	}
	for(int i=0;i<PART1.size();i++)
	{
	 if(PART1[i].type==ELASTIC) 
	 {
		 swap(PART[p],PART1[i]);
		 p++;
	 }
	}
	for(int i=0;i<PART1.size();i++)
	{
	 if(PART1[i].type==TERMINAL1 || PART1[i].type==TERMINAL2 || PART1[i].type==WALL) 
	 {
		 swap(PART[p],PART1[i]);
		 p++;
		 
	 }
	}
	for(int i=0;i<PART1.size();i++)
	{
		if(PART1[i].type==HYPERELAST) 
		 {
			 swap(PART[p],PART1[i]);
			 p++;
		 }
	}

	cout<<"並べ替え終了"<<endl;
	//////////////////////////
	}
	else if(num==false){
		for(int i=0;i<PART1.size();i++)
		{
			PART[i]=PART1[i];
		}
	}
}

//粒子数カウント関数 & 並び替え（28行）
//基本的にパーツの順番はset_initial_placementで計算しているので並び替える必要はない
void calc_numbers_of_particles_and_change_the_order(mpsconfig *CON,vector <mpselastic> &PART, int *fluid_number,int *hyper_number,int *magnetic_number,int *out,int *order_sw)
{
	//各粒子数をカウント
	int elastic_num=0;						//エラストマー粒子
	int magelast_num=0;
	int wall_num=0;
	int solid_num=0;
	int out_num=0;
	int magelast_num2=0;
	int hyper_num=0;


	for(int i=0;i<PART.size();i++) //mapでOK
	{
		if(PART[i].type==ELASTIC) elastic_num++;
		else if(PART[i].type==MAGELAST) magelast_num++;
		else if(PART[i].type==WALL) wall_num++;
		else if(PART[i].type==TERMINAL1) wall_num++;
		else if(PART[i].type==TERMINAL2) wall_num++;
		else if(PART[i].type==MAGELAST2) magelast_num2++;
		else if(PART[i].type==HYPERELAST) hyper_num++;
//		else if(PART[i].type==FLUID) fluid_num++;
//		else if(PART[i].type==INWALL) inwall_num++;
	}
	
	cout<<"elastic: "<<elastic_num<<", magelast: "<<magelast_num+magelast_num2<<", solid: "<<wall_num<<", hyperelast:"<<hyper_num<<endl;
	out_num=elastic_num+magelast_num+magelast_num2+solid_num+wall_num+hyper_num;	//fluid_number<=i<outがINWALL群で、out<=iがOUTWALL群になる。
	*out=out_num;
	*fluid_number=elastic_num+magelast_num+magelast_num2+hyper_num;
	*hyper_number=hyper_num;
	*magnetic_number=hyper_num+magelast_num+magelast_num2;
	if(out_num!=PART.size()){
		cout<<"\n！粒子数合計エラー！"<<endl;
		exit(EXIT_FAILURE);
	}

	//並び替え
/*	mpselastic PART_temp;			//並び替え用粒子クラス
	for(int i=0;i<fluid_num;i++)
	{
		if(PART[i].type!=ELASTIC)//if(PART[i].type!=FLUID)
		{
			cout<<"並び替え必要あり"<<endl;
		}
	}
*/	*order_sw=OFF;
}

//INDEX更新関数（17行）
//INDEX: 各格子のなかに含まれる粒子数。
//INDEXを数える。格子番号は２次元では左下で0。X方向につれ＋１で、右上で最大(X_mesh*Y_mesh)。３次元ではＺ方向にも増えていく
void reload_INDEX(mpsconfig &CON, vector<mpselastic> &PART, int *INDEX)
{       	
	
	int X, Y, Z;	//X, Y, Z方向に何個目の格子か 
	int lattice_number;		//粒子iを含む格子の番号
	double width=CON.get_distancebp()*CON.get_dx();		//格子幅
//	cout<<"error_check"<<endl;
	for(int i=0;i<CON.get_number_of_mesh();i++) INDEX[i]=0; //
	for(int i=0;i<PART.size();i++)
	{
		//領域外粒子
		if(!(PART[i].r[A_X]>CON.get_minX() && PART[i].r[A_X]<CON.get_maxX())) cout<<"X="<<PART[i].r[A_X]<<", i="<<i<<endl;
		else if(!(PART[i].r[A_Y]>CON.get_minY() && PART[i].r[A_Y]<CON.get_maxY())) cout<<"Y="<<PART[i].r[A_Y]<<", i="<<i<<endl;
//		else if(!(PART[i].r[A_Z]>CON->get_minZ() && PART[i].r[A_Z]<CON->get_maxZ())) cout<<"Z="<<PART[i].r[A_Z]<<", i="<<i<<endl; 

		//////////////////

		X=(int)((PART[i].r[A_X]-CON.get_minX())/width);
		Y=(int)((PART[i].r[A_Y]-CON.get_minY())/width);
		Z=(int)((PART[i].r[A_Z]-CON.get_minZ())/width);

		lattice_number=(Z*CON.get_X_mesh()*CON.get_Y_mesh())+(Y*CON.get_X_mesh())+X;

        PART[i].index=lattice_number; //粒子iが何番目の格子に含まれているかを計算
		INDEX[lattice_number]++; //格子に含まれる粒子数を計算
	}
//	cout<<"error_check_end"<<endl;
}

//INDEX更新関数その２ MESHに粒子番号を格納する（15行）
//MESH: int **MESH = new int *[CON.get_number_of_mesh()];	
//get_number_of_mesh() {return (int)((maxX-minX)/(distancebp*dx)*(maxY-minY)/(distancebp*dx)*(maxZ-minZ)/(distancebp*dx)+0.001);} //格子数：X_mesh*Y_mesh*Z_mesh
void reload_INDEX2(mpsconfig *CON, vector<mpselastic> &PART, int **MESH)
{
	int number_of_mesh = CON->get_number_of_mesh();
	int *count = new int [number_of_mesh];
	for(int i=0;i<number_of_mesh;i++) count[i]=0;//初期化

	for(int i=0;i<PART.size();i++)
	{
		int lattice_number=PART[i].index; //粒子iが何番目の格子に含まれているかを取得
        
		MESH[lattice_number][count[lattice_number]]=i;
		count[lattice_number]++;
	}
	delete [] count;
}

//法線ベクトル作成関数（30行）
void direct_f(mpsconfig &CON,vector<mpselastic> &PART,int i,double *direct[DIMENSION])
{
	
	double R=CON.get_re3()*CON.get_distancebp();//法線ベクトル計算に利用する影響半径

	double px=PART[i].r[A_X]+CON.get_distancebp();//x+le
	double mx=PART[i].r[A_X]-CON.get_distancebp();//x-le
	double py=PART[i].r[A_Y]+CON.get_distancebp();//y+le
	double my=PART[i].r[A_Y]-CON.get_distancebp();//y-le
	double pz=PART[i].r[A_Z]+CON.get_distancebp();//z+le
	double mz=PART[i].r[A_Z]-CON.get_distancebp();//z-le
	
	double pnd_px=pnd_for_direct(CON,PART,px,PART[i].r[A_Y],PART[i].r[A_Z],R,i);
	double pnd_mx=pnd_for_direct(CON,PART,mx,PART[i].r[A_Y],PART[i].r[A_Z],R,i);
	double pnd_py=pnd_for_direct(CON,PART,PART[i].r[A_X],py,PART[i].r[A_Z],R,i);
	double pnd_my=pnd_for_direct(CON,PART,PART[i].r[A_X],my,PART[i].r[A_Z],R,i);
	double pnd_pz=pnd_for_direct(CON,PART,PART[i].r[A_X],PART[i].r[A_Y],pz,R,i);
	double pnd_mz=pnd_for_direct(CON,PART,PART[i].r[A_X],PART[i].r[A_Y],mz,R,i);

	direct[A_X][i]=(pnd_px-pnd_mx)/(2*CON.get_distancebp());
	direct[A_Y][i]=(pnd_py-pnd_my)/(2*CON.get_distancebp());
	direct[A_Z][i]=(pnd_pz-pnd_mz)/(2*CON.get_distancebp());
	
	double a=sqrt(direct[A_X][i]*direct[A_X][i]+direct[A_Y][i]*direct[A_Y][i]+direct[A_Z][i]*direct[A_Z][i]);
	if(a!=0)
	{ 
		direct[A_X][i]/=a;
		direct[A_Y][i]/=a;
		direct[A_Z][i]/=a;
	}
}

//direct用粒子数密度測定（25行）
double pnd_for_direct(mpsconfig &CON,vector<mpselastic> &PART,double x,double y,double z,double R,int i)
{
	//粒子iの位置とはずれるので、MESHを使用したほうがよい
	//教科書では個数をかぞえているが、ここでは重み関数を用いる。
	//R=CON->get_re3()*CON.get_distancebp();
	double spnd=0;

	for(int k=0;k<PART[i].N3;k++)
	{       
		int j=PART[i].NEI3[k];
		double X=PART[j].r[A_X]-x;
		double Y=PART[j].r[A_Y]-y;
		double Z=PART[j].r[A_Z]-z;
		double dis=sqrt(X*X+Y*Y+Z*Z);
		//if(dis<R) spnd++;   //教科書どうり
		if(dis<R)
		{
			double w=(1-dis/R)*(1-dis/R);
			spnd+=w;
		}
	}
	return spnd;
}

//クーラン数
void courant_elastic(vector<mpselastic> &PART, int fluid_number, int t, double *dt, double mindis, double Umax, double *g)
{
	mpsconfig CON;
	double CFL=CON.get_courant();	//クーラン条件数
	double CFL2=0;
	double le=mindis;				//最短粒子間距離
	double E=CON.get_E_m();		//MAGELASTのヤング率
	double v=CON.get_v_m();
	double lambda=v*E/((1.0+v)*(1.0-2.0*v));
	double mu=E/(2.0*(1.0+v));
//	double density=CON->Get_density();
	double factor=20.0;
	double uu;
	double cfl2;
	ofstream aaa("CFL.dat", ios::app);
	for(int i=0;i<PART.size();i++){
		if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2)CFL=sqrt((lambda+2.0*mu)/CON.Get_MRE_density());
		else if(PART[i].type==ELASTIC)CFL=sqrt((lambda+2.0*mu)/CON.Get_Silicone_density());
		uu=pow(pow(PART[i].u[A_X],2)+pow(PART[i].u[A_Y],2)+pow(PART[i].u[A_Z],2),0.5);
		cfl2=(*dt*uu)/0.001;
		if(cfl2>CFL2) CFL2=cfl2; 
	}
	aaa<<"le="<<le<<", dt= "<<*dt<<", CFL= "<<CFL<<", CFL2="<<CFL2<<", new dt="<<le/CFL<<endl;
	aaa.close();
/*	///////////クーラン数//////////////////
	if(CFL>0)
	{      
		double newdt=*dt;	//新しいdt
		CFL=sqrt((lambda+2.0*mu)/density);

	//	if(Umax!=0) newdt=CFL*le/Umax;
		if(Umax!=0) newdt=(le/CFL)/factor;
		if(newdt>CON->get_dt()) *dt=CON->get_dt();
		else *dt=newdt;
		
		if(*dt!=CON->get_dt()) cout<<"CFL条件によるdt更新 dt="<<*dt<<endl;
	}

	///拡散数の正確な定義を調べて書きなおせ
	if(CON->get_vis()!=0 && CON->get_vis_calc_type()==POSITIVE)
	{
		if(*dt>0.25*le*le/CON->get_vis()) cout<<"拡散数違反"<<endl;
	}
	/////////////////////////////*/
}

//CG法（90行）
void CG_method(mpsconfig *CON,double *r,double *P,double *AP,double *val,int *ind,int *ptr,int pn,double *X,int *countN,double EP)
{
	cout<<"CG法スタート------";
	//unsigned int timeCG=GetTickCount();
	int count=0;
	double rr=0;
	double E=1;//誤差
	double alp,beta;
	for(int n=0;n<pn;n++) rr+=r[n]*r[n];
	
	if(CON->get_omp_P()==OFF)//通常版
	{
		while(E>EP)// EP=CON->get_CGep();//収束判定(convergence test)
		{
			count++;
			//////////////alpを求める
			for(int n=0;n<pn;n++)
			{      
				AP[n]=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
			}
			double PAP=0;
			for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
			alp=rr/PAP;
		//	cout<<"alp="<<alp<<" rr="<<rr<<" PAP="<<PAP<<endl;
			//////////////////////
		
			//////////////// 解更新　X(k+1)=X(k)+alp*P
			for(int n=0;n<pn;n++) X[n]+=alp*P[n];
			//////////////////////////////
			
			//////////////// r=r-alp*AP
			for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
			/////////////////////////////
			
			///////////////////////beta
			beta=1.0/rr;
			rr=0;
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			beta=beta*rr;
			///////////////////////

			//////////////////誤差
			E=sqrt(rr);
			//cout<<"E="<<E<<endl;
			////////////////////////
			
			///////////////////// P=r+beta*P
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];
		}
	}
	else if(CON->get_omp_P()==ON)//openMPを使用する場合
	{
		while(E>EP)
		{
			count++;
			//////////////alpを求める
			double PAP=0;
			#pragma omp parallel for reduction(+:PAP)
			for(int n=0;n<pn;n++)
			{
				AP[n]=0;
				for(int j=ptr[n];j<ptr[n+1];j++) AP[n]+=val[j]*P[ind[j]];
				PAP+=P[n]*AP[n];
			}
			alp=rr/PAP;
			//////////////////////

			//////////////
			E=0;//誤差
			beta=1.0/rr;
			rr=0;
			#pragma omp parallel for reduction(+:rr)
			for(int n=0;n<pn;n++) 
			{
				X[n]+=alp*P[n];// 解更新　X(k+1)=X(k)+alp*P
				r[n]-=alp*AP[n];// r=r-alp*AP
				rr+=r[n]*r[n];
			}
			E=sqrt(rr);
			//cout<<"E="<<E<<endl;
		
			beta=beta*rr;///beta
			
			for(int n=0;n<pn;n++) P[n]=r[n]+beta*P[n];/// P=r+beta*P
		}
	}
	*countN=count;//反復回数を渡す
}

//ICCG法(230行)
void iccg(mpsconfig *CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X,double *r,double *P,double EP,int *count2)
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

	num2=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
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
	
	for(int n=0;n<pn;n++) P[n]=LDLt_r[n];

	cout<<"ICCG法:未知数="<<pn<<" ---";
	unsigned int time=GetTickCount();
	int count=0;
	double ep=EP;//収束判定
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
		//cout<<"alp="<<alp<<endl;
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
		//cout<<"E="<<E<<endl;
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
	}
	//cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"/";
		
	delete [] AP;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

	delete [] val2;
	delete [] ind2;
	delete [] ptr2;

	for(int i=0;i<pn;i++)
	{
		delete [] VAL[i];
		delete [] IND[i];
	}
	delete [] VAL;
	delete [] IND;
	delete [] NUM;
	*count2=count;//反復回数を格納して返す
}


//ガウスの消去法（24行）
//解は最終的にBのなかへ
void gauss(double *matrix,double *B,int N)
{
	for(int k=0;k<N;k++)
	{
		double akk=matrix[k*N+k];
		
		for(int i=0;i<N;i++)
		{
			if(i!=k)
			{
				double A=matrix[i*N+k]/akk;
				//for(int j=0;j<N;j++)
				for(int j=k;j<N;j++)
				{					
					matrix[i*N+j]-=A*matrix[k*N+j];
				}
				B[i]-=A*B[k];
			}
		}
	}
	for(int k=0;k<N;k++) B[k]/=matrix[k*N+k];

}




//仮の速度および位置決定（46行）
void renewal_u_and_r_in_positive(vector<mpselastic> &PART,int fluid_number,int t,double dt,double *Umax,double **potential,double **laplacian,double *g,double **previous_Un,double **F)
{
	mpsconfig CON;
	double U=0;						//最大速度
	double vis;//=CON->Get_vis();
	double mass=CON.get_particle_mass();	//粒子の質量
	int d=CON.get_dimension();
	int sw=CON.get_temporary_r();	//ＯＮなら仮の位置を計算する

	double *old_U[DIMENSION];
	for(int D=0;D<DIMENSION;D++) old_U[D]=new double [fluid_number];//変更前の速度を記憶しておく

	if(t==1) for(int i=0;i<fluid_number;i++)for(int D=0;D<DIMENSION;D++) previous_Un[D][i]=0;//t=1のときは初期化      
			

	//potential[D][i]を場合によってはゼロに初期化する
	if(CON.get_dir_for_P()==1 || CON.get_dir_for_P()==3) //表面粒子の表面張力は圧力値として計算されているので,ここでは考慮しないよう初期化する
	{
		for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) potential[D][i]=0;
	}
	//////////////////////*/

	/////////////速度更新
	for(int i=0;i<fluid_number;i++)
	{        
		double speed=0;//粒子速度
		for(int D=0;D<d;D++)
		{   
			if(PART[i].type==ELASTIC) vis=CON.Get_Silicone_vis();
			else if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2) vis=CON.Get_MRE_vis();
			old_U[D][i]=PART[i].u[D];
			
			PART[i].u[D]+=dt*(vis *laplacian[D][i]+potential[D][i]+g[D]+F[D][i]/mass);
			
			//PART[i].u[D]=previous_Un[D][i]+dt*(vis*laplacian[D][i]+potential[D][i]+g[D]);//蛙とび法
			speed+=PART[i].u[D]*PART[i].u[D];
		}
		if(speed>U) U=speed;
	}
	*Umax=U;

	//位置更新
	if(sw==ON) for(int i=0;i<fluid_number;i++) for(int D=0;D<d;D++) PART[i].r[D]+=dt*0.5*(PART[i].u[D]+old_U[D][i]);//台形則
	
}

//速度発散計算関数（27行）
double divergence(mpsconfig *CON,vector<mpselastic> &PART,int i,double n0)
{
    double W=0;										//粒子数密度
    double R=CON->get_distancebp()*CON->get_re();	//影響半径
    double div=0;									//発散の値

	for(int k=0;k<PART[i].N;k++)
    {    
        int j=PART[i].NEI[k]; 
        double X=PART[j].r[A_X]-PART[i].r[A_X];
		double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
		double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
		double dis=sqrt(X*X+Y*Y+Z*Z);
			       
		double w=kernel(R,dis);
		
		div+=(PART[j].u[A_X]-PART[i].u[A_X])*X*w/(dis*dis);
		div+=(PART[j].u[A_Y]-PART[i].u[A_Y])*Y*w/(dis*dis);
		div+=(PART[j].u[A_Z]-PART[i].u[A_Z])*Z*w/(dis*dis);
		W+=w;
    }
    if(W!=0)
	{
		div*=CON->get_dimension()/W;
	}
    return div;
}


///速度発散計算関数(WLSM法)（198行）
///WLSM=Weighed Least Square Method：重み付き最小二乗法
double divergence2(mpsconfig *CON,vector<mpselastic> &PART,int i)
{
    double div=0;//発散の値
	double le=CON->get_distancebp();
	
	double r=CON->get_re();
	double R=r*le;
	int d=CON->get_dimension();
	int N=0;					//係数行列の元
	int order=1;				//近似曲面のｵｰﾀﾞｰ。 1=線形 2=二次
    
	//係数行列の大きさの決定
	if(d==2)
	{
		if(order==1) N=2;
		else if(order==2) N=5;
	}
	else if(d==3)
	{
		if(order==1) N=4;
		else if(order==2) N=10;
	}
	////////////////////////////////

	double *matrix=new double [N*N];	//N×Nの係数行列
	double *B1=new double [N];			//Nの解行列
	double *B2=new double [N];			//Nの解行列
	double *B3=new double [N];			//Nの解行列

	for(int n=0;n<N*N;n++) matrix[n]=0;	//初期化
	for(int n=0;n<N;n++) {B1[n]=0;B2[n]=0;B3[n]=0;}

	if(d==2 && order==1)				//二次元
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			
			double X=(PART[j].r[A_X]-PART[i].r[A_X])/le;// leで割るのは打ち切り誤差防止
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y])/le;
			double U=(PART[j].u[A_X]-PART[i].u[A_X]);
			double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
			double dis=sqrt(X*X+Y*Y);
					
			double w=1;
			if(dis>1) w=1/(dis*dis*dis*dis);
					
			matrix[0]+=X*X*w;			//ΣXjwj
			matrix[1]+=X*Y*w;		//ΣXjYjwj
			matrix[3]+=Y*Y*w;			//ΣYjwj
				
			B1[0]+=U*X*w;//ΣujXjwj
			B1[1]+=U*Y*w;//ΣujYjwj
			B2[0]+=V*X*w;//ΣvjXjwj
			B2[1]+=V*Y*w;//ΣvjYjwj
		}
			
		matrix[2]=matrix[1];		//ΣXjYjwj
		for(int n=0;n<N;n++)
		{
			B1[n]/=le;//打ち切り誤差防止
			B2[n]/=le;//打ち切り誤差防止
		}

		double determinant=matrix[0]*matrix[3]-matrix[1]*matrix[2];//行列式
			
		double dudx=(B1[0]*matrix[3]-matrix[1]*B1[1])/determinant;
		double dvdy=(B2[1]*matrix[0]-matrix[2]*B2[0])/determinant;

		div=dudx+dvdy;		
	}
	if(d==2 && order==2)//二次元2次式
	{
		for(int k=0;k<PART[i].N;k++)
		{
			int j=PART[i].NEI[k];
			
			double X=(PART[j].r[A_X]-PART[i].r[A_X]);// leで割るのは打ち切り誤差防止
			double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
			double U=(PART[j].u[A_X]-PART[i].u[A_X]);
			double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
			double dis=sqrt(X*X+Y*Y);
					
			double w=1;
			//if(dis>1) w=r*r*r*r/(dis*dis*dis*dis);
			if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
					
			matrix[0]+=X*X*w;			//ΣXjwj
			matrix[1]+=X*Y*w;		//ΣXjYjwj
			matrix[2]+=X*X*X*w;			//ΣXj^3wj
			matrix[3]+=X*X*Y*w;			//ΣXj^2Yjwj
			matrix[4]+=X*Y*Y*w;			//ΣXjYj^2wj

			matrix[6]+=Y*Y*w;			//ΣYj^2wj
			matrix[9]+=Y*Y*Y*w;			//ΣYj^3wj

			matrix[12]+=X*X*X*X*w;			//ΣXj^4wj
			matrix[13]+=X*X*X*Y*w;			//ΣXj^3Yjwj
			matrix[14]+=X*X*Y*Y*w;			//ΣXj^2Yj^2wj
	
			matrix[19]+=X*Y*Y*Y*w;			//ΣXjYj^3wj

			matrix[24]+=Y*Y*Y*Y*w;			//ΣYj^4wj

				
			B1[0]+=U*X*w;//ΣujXjwj
			B1[1]+=U*Y*w;//ΣujYjwj
			B1[2]+=U*X*X*w;//ΣujXj^2wj
			B1[3]+=U*X*Y*w;//ΣujXjYjwj
			B1[4]+=U*Y*Y*w;//ΣujYj^2wj

			B2[0]+=V*X*w;//ΣvjXjwj
			B2[1]+=V*Y*w;//ΣvjYjwj
			B2[2]+=V*X*X*w;//ΣvjXj^2wj
			B2[3]+=V*X*Y*w;//ΣvjXjYjwj
			B2[4]+=V*Y*Y*w;//ΣvjYj^2wj
		}
			
		matrix[5]=matrix[1];		//ΣXjYjwj
		matrix[7]=matrix[3];		//ΣXj^2Yjwj
		matrix[8]=matrix[4];		//ΣXjYj^2wj
		matrix[10]=matrix[2];		//ΣXj^3Yjwj
		matrix[11]=matrix[3];		//ΣXj^2Yjwj
		matrix[15]=matrix[3];		//ΣXj^2Yjwj
		matrix[16]=matrix[4];		//ΣXjYj^2wj
		matrix[17]=matrix[13];		//ΣXj^3Yjwj
		matrix[18]=matrix[14];		//ΣXj^2Yj^2wj
		matrix[20]=matrix[4];		//ΣXjYj^2wj
		matrix[21]=matrix[9];		//ΣYj^3wj
		matrix[22]=matrix[14];		//ΣXj^2Yj^2wj
		matrix[23]=matrix[19];		//ΣXjYj^3wj

		///丸め誤差防止
		for(int n=0;n<N;n++)
		{
			for(int m=0;m<N;m++) matrix[n*N+m]*=1e7;
			B1[n]*=1e7;
			B2[n]*=1e7;
		}//*/

		double dudx=0;
		double dvdy=0;
		return_X_for5N(matrix,N,B1,B2,&dudx,&dvdy);//5次連立方程式の、解１，２を返す関数
		
		div=(dudx+dvdy);
			
	}
	else if(d==3 && order==1)//3次元1次式
	{
		if(PART[i].N>5)
		{
			//div=calc_WLSM_divu_D3_order1(CON,PART,matrix,B1,B2,B3,i,3);//最後の引数は未知数
			div=calc_WLSM_divu_D3_order1_2(CON,PART,matrix,B1,B2,B3,i,4);//最後の引数は未知数
		}
		else 
		{
			div=0;
			//cout<<"警告 近隣粒子数が5以下("<<PART[i].N<<")です"<<endl;
		}
	}
	else if(d==3 && order==2)//3次元2次式
	{
		//P=Pi+aΔx+bΔy+cΔz+dΔx2+eΔy2+fΔz2+gΔxΔy+hΔyΔz+iΔzΔxとおくと、
		///係数行列は
		///   ΣΔx2      ΣΔxΔy    ΣΔxΔz    ΣΔx3       ΣΔxΔy2    ΣΔxΔz2    ΣΔx2Δy     ΣΔxΔyΔz  ΣΔx2Δz     a = ΣΔxΔP  
		///   ΣΔxΔy    ΣΔy2      ΣΔyΔz    ΣΔx2Δy    ΣΔy3       ΣΔyΔz2    ΣΔxΔy2     ΣΔy2Δz    ΣΔxΔyΔz   b = ΣΔyΔP
		///   ΣΔxΔz    ΣΔyΔz    ΣΔz2      ΣΔx2Δz    ΣΔy2Δz    ΣΔz3       ΣΔxΔyΔz   ΣΔyΔz2    ΣΔxΔz2     c = ΣΔzΔP
		///   ΣΔx3      ΣΔx2Δy   ΣΔx2Δz   ΣΔx4       ΣΔx2Δy2   ΣΔx2Δz2   ΣΔx3Δy     ΣΔx2ΔyΔz ΣΔx3Δz     d = ΣΔx2ΔP
		///   ΣΔxΔy2   ΣΔy3      ΣΔy2Δz   ΣΔx2Δy2   ΣΔy4       ΣΔy2Δz2   ΣΔxΔy3     ΣΔy3Δz    ΣΔxΔy2Δz  e = ΣΔy2ΔP
		///   ΣΔxΔz2   ΣΔyΔz2   ΣΔz3      ΣΔx2Δz2   ΣΔy2Δz2   ΣΔz4       ΣΔxΔyΔz2  ΣΔyΔz3    ΣΔxΔz3     f = ΣΔz2ΔP
		///   ΣΔx2Δy   ΣΔxΔy2   ΣΔxΔyΔz ΣΔx3Δy    ΣΔxΔy3    ΣΔxΔyΔz2 ΣΔx2Δy2    ΣΔxΔy2Δz ΣΔx2ΔyΔz  g = ΣΔxΔyΔP
		///   ΣΔxΔyΔz ΣΔy2Δz   ΣΔyΔz2   ΣΔx2ΔyΔz ΣΔy3Δz    ΣΔyΔz3    ΣΔxΔy2Δz  ΣΔy2Δz2   ΣΔxΔyΔz2  h = ΣΔyΔzΔP
		///   ΣΔx2Δz   ΣΔxΔyΔz ΣΔxΔz2   ΣΔx3Δz    ΣΔxΔy2Δz  ΣΔxΔz3   ΣΔx2Δy Δz ΣΔxΔyΔz2 ΣΔx2Δz2    g = ΣΔxΔzΔP
		
		if(PART[i].N>8)
		{
			//div=calc_WLSM_divu_D3_order2(CON,PART,matrix,B1,B2,B3,i,9);
			div=calc_WLSM_divu_D3_order2_2(CON,PART,matrix,B1,B2,B3,i,10);
		}
		else if(PART[i].N>5)
		{
			//div=calc_WLSM_divu_D3_order1(CON,PART,matrix,B1,B2,B3,i,3);//この場合、最後の引数はN=3を渡すことに注意
			div=calc_WLSM_divu_D3_order1_2(CON,PART,matrix,B1,B2,B3,i,4);//最後の引数は未知数
		}
		else 
		{
			div=0;
			//cout<<"警告 近隣粒子数が9以下("<<PART[i].N<<")です"<<endl;
		}
	}

	delete [] matrix;
	delete [] B1;
	delete [] B2;
	delete [] B3;

    return div;
}

//5次の連立方程式の解1,2を返す関数（24行）
void return_X_for5N(double *matrix,int N,double *B1,double *B2,double *dudx,double *dudy)
{
	double a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55;
	double b1,b2,b3,b4,b5;
	double c1,c2,c3,c4,c5;

	a11=matrix[0];a12=matrix[1];a13=matrix[2];a14=matrix[3];a15=matrix[4];
	a21=matrix[5];a22=matrix[6];a23=matrix[7];a24=matrix[8];a25=matrix[9];
	a31=matrix[10];a32=matrix[11];a33=matrix[12];a34=matrix[13];a35=matrix[14];
	a41=matrix[15];a42=matrix[16];a43=matrix[17];a44=matrix[18];a45=matrix[19];
	a51=matrix[20];a52=matrix[21];a53=matrix[22];a54=matrix[23];a55=matrix[24];

	b1=B1[0];b2=B1[1];b3=B1[2];b4=B1[3];b5=B1[4];
	c1=B2[0];c2=B2[1];c3=B2[2];c4=B2[3];c5=B2[4];
	
	double determinant=(a11*a22*a33*a44*a55-a11*a22*a33*a45*a54-a11*a22*a34*a43*a55+a11*a22*a34*a45*a53+a11*a22*a35*a43*a54-a11*a22*a35*a44*a53-a11*a23*a32*a44*a55+a11*a23*a32*a45*a54+a11*a23*a34*a42*a55-a11*a23*a34*a45*a52-a11*a23*a35*a42*a54+a11*a23*a35*a44*a52+a11*a24*a32*a43*a55-a11*a24*a32*a45*a53-a11*a24*a33*a42*a55+a11*a24*a33*a45*a52+a11*a24*a35*a42*a53-a11*a24*a35*a43*a52-a11*a25*a32*a43*a54+a11*a25*a32*a44*a53+a11*a25*a33*a42*a54-a11*a25*a33*a44*a52-a11*a25*a34*a42*a53+a11*a25*a34*a43*a52-a12*a21*a33*a44*a55+a12*a21*a33*a45*a54+a12*a21*a34*a43*a55-a12*a21*a34*a45*a53-a12*a21*a35*a43*a54+a12*a21*a35*a44*a53+a12*a23*a31*a44*a55-a12*a23*a31*a45*a54-a12*a23*a34*a41*a55+a12*a23*a34*a45*a51+a12*a23*a35*a41*a54-a12*a23*a35*a44*a51-a12*a24*a31*a43*a55+a12*a24*a31*a45*a53+a12*a24*a33*a41*a55-a12*a24*a33*a45*a51-a12*a24*a35*a41*a53+a12*a24*a35*a43*a51+a12*a25*a31*a43*a54-a12*a25*a31*a44*a53-a12*a25*a33*a41*a54+a12*a25*a33*a44*a51+a12*a25*a34*a41*a53-a12*a25*a34*a43*a51+a13*a21*a32*a44*a55-a13*a21*a32*a45*a54-a13*a21*a34*a42*a55+a13*a21*a34*a45*a52+a13*a21*a35*a42*a54-a13*a21*a35*a44*a52-a13*a22*a31*a44*a55+a13*a22*a31*a45*a54+a13*a22*a34*a41*a55-a13*a22*a34*a45*a51-a13*a22*a35*a41*a54+a13*a22*a35*a44*a51+a13*a24*a31*a42*a55-a13*a24*a31*a45*a52-a13*a24*a32*a41*a55+a13*a24*a32*a45*a51+a13*a24*a35*a41*a52-a13*a24*a35*a42*a51-a13*a25*a31*a42*a54+a13*a25*a31*a44*a52+a13*a25*a32*a41*a54-a13*a25*a32*a44*a51-a13*a25*a34*a41*a52+a13*a25*a34*a42*a51-a14*a21*a32*a43*a55+a14*a21*a32*a45*a53+a14*a21*a33*a42*a55-a14*a21*a33*a45*a52-a14*a21*a35*a42*a53+a14*a21*a35*a43*a52+a14*a22*a31*a43*a55-a14*a22*a31*a45*a53-a14*a22*a33*a41*a55+a14*a22*a33*a45*a51+a14*a22*a35*a41*a53-a14*a22*a35*a43*a51-a14*a23*a31*a42*a55+a14*a23*a31*a45*a52+a14*a23*a32*a41*a55-a14*a23*a32*a45*a51-a14*a23*a35*a41*a52+a14*a23*a35*a42*a51+a14*a25*a31*a42*a53-a14*a25*a31*a43*a52-a14*a25*a32*a41*a53+a14*a25*a32*a43*a51+a14*a25*a33*a41*a52-a14*a25*a33*a42*a51+a15*a21*a32*a43*a54-a15*a21*a32*a44*a53-a15*a21*a33*a42*a54+a15*a21*a33*a44*a52+a15*a21*a34*a42*a53-a15*a21*a34*a43*a52-a15*a22*a31*a43*a54+a15*a22*a31*a44*a53+a15*a22*a33*a41*a54-a15*a22*a33*a44*a51-a15*a22*a34*a41*a53+a15*a22*a34*a43*a51+a15*a23*a31*a42*a54-a15*a23*a31*a44*a52-a15*a23*a32*a41*a54+a15*a23*a32*a44*a51+a15*a23*a34*a41*a52-a15*a23*a34*a42*a51-a15*a24*a31*a42*a53+a15*a24*a31*a43*a52+a15*a24*a32*a41*a53-a15*a24*a32*a43*a51-a15*a24*a33*a41*a52+a15*a24*a33*a42*a51);
	
	*dudx=(b1*a22*a33*a44*a55-b1*a22*a33*a45*a54-b1*a22*a34*a43*a55+b1*a22*a34*a45*a53+b1*a22*a35*a43*a54-b1*a22*a35*a44*a53-b1*a23*a32*a44*a55+b1*a23*a32*a45*a54+b1*a23*a34*a42*a55-b1*a23*a34*a45*a52-b1*a23*a35*a42*a54+b1*a23*a35*a44*a52+b1*a24*a32*a43*a55-b1*a24*a32*a45*a53-b1*a24*a33*a42*a55+b1*a24*a33*a45*a52+b1*a24*a35*a42*a53-b1*a24*a35*a43*a52-b1*a25*a32*a43*a54+b1*a25*a32*a44*a53+b1*a25*a33*a42*a54-b1*a25*a33*a44*a52-b1*a25*a34*a42*a53+b1*a25*a34*a43*a52-a12*b2*a33*a44*a55+a12*b2*a33*a45*a54+a12*b2*a34*a43*a55-a12*b2*a34*a45*a53-a12*b2*a35*a43*a54+a12*b2*a35*a44*a53+a12*a23*b3*a44*a55-a12*a23*b3*a45*a54-a12*a23*a34*b4*a55+a12*a23*a34*a45*b5+a12*a23*a35*b4*a54-a12*a23*a35*a44*b5-a12*a24*b3*a43*a55+a12*a24*b3*a45*a53+a12*a24*a33*b4*a55-a12*a24*a33*a45*b5-a12*a24*a35*b4*a53+a12*a24*a35*a43*b5+a12*a25*b3*a43*a54-a12*a25*b3*a44*a53-a12*a25*a33*b4*a54+a12*a25*a33*a44*b5+a12*a25*a34*b4*a53-a12*a25*a34*a43*b5+a13*b2*a32*a44*a55-a13*b2*a32*a45*a54-a13*b2*a34*a42*a55+a13*b2*a34*a45*a52+a13*b2*a35*a42*a54-a13*b2*a35*a44*a52-a13*a22*b3*a44*a55+a13*a22*b3*a45*a54+a13*a22*a34*b4*a55-a13*a22*a34*a45*b5-a13*a22*a35*b4*a54+a13*a22*a35*a44*b5+a13*a24*b3*a42*a55-a13*a24*b3*a45*a52-a13*a24*a32*b4*a55+a13*a24*a32*a45*b5+a13*a24*a35*b4*a52-a13*a24*a35*a42*b5-a13*a25*b3*a42*a54+a13*a25*b3*a44*a52+a13*a25*a32*b4*a54-a13*a25*a32*a44*b5-a13*a25*a34*b4*a52+a13*a25*a34*a42*b5-a14*b2*a32*a43*a55+a14*b2*a32*a45*a53+a14*b2*a33*a42*a55-a14*b2*a33*a45*a52-a14*b2*a35*a42*a53+a14*b2*a35*a43*a52+a14*a22*b3*a43*a55-a14*a22*b3*a45*a53-a14*a22*a33*b4*a55+a14*a22*a33*a45*b5+a14*a22*a35*b4*a53-a14*a22*a35*a43*b5-a14*a23*b3*a42*a55+a14*a23*b3*a45*a52+a14*a23*a32*b4*a55-a14*a23*a32*a45*b5-a14*a23*a35*b4*a52+a14*a23*a35*a42*b5+a14*a25*b3*a42*a53-a14*a25*b3*a43*a52-a14*a25*a32*b4*a53+a14*a25*a32*a43*b5+a14*a25*a33*b4*a52-a14*a25*a33*a42*b5+a15*b2*a32*a43*a54-a15*b2*a32*a44*a53-a15*b2*a33*a42*a54+a15*b2*a33*a44*a52+a15*b2*a34*a42*a53-a15*b2*a34*a43*a52-a15*a22*b3*a43*a54+a15*a22*b3*a44*a53+a15*a22*a33*b4*a54-a15*a22*a33*a44*b5-a15*a22*a34*b4*a53+a15*a22*a34*a43*b5+a15*a23*b3*a42*a54-a15*a23*b3*a44*a52-a15*a23*a32*b4*a54+a15*a23*a32*a44*b5+a15*a23*a34*b4*a52-a15*a23*a34*a42*b5-a15*a24*b3*a42*a53+a15*a24*b3*a43*a52+a15*a24*a32*b4*a53-a15*a24*a32*a43*b5-a15*a24*a33*b4*a52+a15*a24*a33*a42*b5)/determinant;
	*dudy=(a11*c2*a33*a44*a55-a11*c2*a33*a45*a54-a11*c2*a34*a43*a55+a11*c2*a34*a45*a53+a11*c2*a35*a43*a54-a11*c2*a35*a44*a53-a11*a23*c3*a44*a55+a11*a23*c3*a45*a54+a11*a23*a34*c4*a55-a11*a23*a34*a45*c5-a11*a23*a35*c4*a54+a11*a23*a35*a44*c5+a11*a24*c3*a43*a55-a11*a24*c3*a45*a53-a11*a24*a33*c4*a55+a11*a24*a33*a45*c5+a11*a24*a35*c4*a53-a11*a24*a35*a43*c5-a11*a25*c3*a43*a54+a11*a25*c3*a44*a53+a11*a25*a33*c4*a54-a11*a25*a33*a44*c5-a11*a25*a34*c4*a53+a11*a25*a34*a43*c5-c1*a21*a33*a44*a55+c1*a21*a33*a45*a54+c1*a21*a34*a43*a55-c1*a21*a34*a45*a53-c1*a21*a35*a43*a54+c1*a21*a35*a44*a53+c1*a23*a31*a44*a55-c1*a23*a31*a45*a54-c1*a23*a34*a41*a55+c1*a23*a34*a45*a51+c1*a23*a35*a41*a54-c1*a23*a35*a44*a51-c1*a24*a31*a43*a55+c1*a24*a31*a45*a53+c1*a24*a33*a41*a55-c1*a24*a33*a45*a51-c1*a24*a35*a41*a53+c1*a24*a35*a43*a51+c1*a25*a31*a43*a54-c1*a25*a31*a44*a53-c1*a25*a33*a41*a54+c1*a25*a33*a44*a51+c1*a25*a34*a41*a53-c1*a25*a34*a43*a51+a13*a21*c3*a44*a55-a13*a21*c3*a45*a54-a13*a21*a34*c4*a55+a13*a21*a34*a45*c5+a13*a21*a35*c4*a54-a13*a21*a35*a44*c5-a13*c2*a31*a44*a55+a13*c2*a31*a45*a54+a13*c2*a34*a41*a55-a13*c2*a34*a45*a51-a13*c2*a35*a41*a54+a13*c2*a35*a44*a51+a13*a24*a31*c4*a55-a13*a24*a31*a45*c5-a13*a24*c3*a41*a55+a13*a24*c3*a45*a51+a13*a24*a35*a41*c5-a13*a24*a35*c4*a51-a13*a25*a31*c4*a54+a13*a25*a31*a44*c5+a13*a25*c3*a41*a54-a13*a25*c3*a44*a51-a13*a25*a34*a41*c5+a13*a25*a34*c4*a51-a14*a21*c3*a43*a55+a14*a21*c3*a45*a53+a14*a21*a33*c4*a55-a14*a21*a33*a45*c5-a14*a21*a35*c4*a53+a14*a21*a35*a43*c5+a14*c2*a31*a43*a55-a14*c2*a31*a45*a53-a14*c2*a33*a41*a55+a14*c2*a33*a45*a51+a14*c2*a35*a41*a53-a14*c2*a35*a43*a51-a14*a23*a31*c4*a55+a14*a23*a31*a45*c5+a14*a23*c3*a41*a55-a14*a23*c3*a45*a51-a14*a23*a35*a41*c5+a14*a23*a35*c4*a51+a14*a25*a31*c4*a53-a14*a25*a31*a43*c5-a14*a25*c3*a41*a53+a14*a25*c3*a43*a51+a14*a25*a33*a41*c5-a14*a25*a33*c4*a51+a15*a21*c3*a43*a54-a15*a21*c3*a44*a53-a15*a21*a33*c4*a54+a15*a21*a33*a44*c5+a15*a21*a34*c4*a53-a15*a21*a34*a43*c5-a15*c2*a31*a43*a54+a15*c2*a31*a44*a53+a15*c2*a33*a41*a54-a15*c2*a33*a44*a51-a15*c2*a34*a41*a53+a15*c2*a34*a43*a51+a15*a23*a31*c4*a54-a15*a23*a31*a44*c5-a15*a23*c3*a41*a54+a15*a23*c3*a44*a51+a15*a23*a34*a41*c5-a15*a23*a34*c4*a51-a15*a24*a31*c4*a53+a15*a24*a31*a43*c5+a15*a24*c3*a41*a53-a15*a24*c3*a43*a51-a15*a24*a33*a41*c5+a15*a24*a33*c4*a51)/determinant;
	
}

//divergence2における、3次元1次近似を行う関数（169行）
double calc_WLSM_divu_D3_order1(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	///係数行列は
	///   ΣΔx2    ΣΔxΔy  ΣΔxΔz  a = ΣΔxΔf  
	///  ΣΔxΔy    ΣΔy2   ΣΔyΔz  b = ΣΔyΔf 
	///  ΣΔxΔz   ΣΔyΔz   ΣΔz2   c = ΣΔzΔf 

	double le=CON->get_distancebp();
	double matrix_val[9];			//matrixの値は一度gauss()で使用すると値が変わるので、変更前のmatrixを保存する
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
			
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);// leで割るのは打ち切り誤差防止
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
					
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
		weight[k]=w;
					
		matrix[0]+=X*X*w;			//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
					
		matrix[4]+=Y*Y*w;			//ΣYj^2wj
		matrix[5]+=Y*Z*w;		//ΣYjZjwj

		matrix[8]+=Z*Z*w;			//ΣZj^2wj
			
		B1[0]+=U*X*w;//ΣfjXjwj
		B1[1]+=U*Y*w;//ΣfjYjwj
		B1[2]+=U*Z*w;//ΣfjZjwj

		B2[0]+=V*X*w;//ΣfjXjwj
		B2[1]+=V*Y*w;//ΣfjYjwj
		B2[2]+=V*Z*w;//ΣfjZjwj

		B3[0]+=W*X*w;//ΣfjXjwj
		B3[1]+=W*Y*w;//ΣfjYjwj
		B3[2]+=W*Z*w;//ΣfjZjwj
	}
			
	matrix[3]=matrix[1];		//ΣXjYjwj
	matrix[6]=matrix[2];		//ΣXjZjwj
	matrix[7]=matrix[5];		//ΣYjZjwj

	for(int L=0;L<9;L++) matrix_val[L]=matrix[L];//行列の値を保存

	/*double dudx=0;//こっちのほうが若干早い。けど誤差評価したいならガウス
	double dvdy=0;
	double dwdz=0;
	double determinant=(matrix[0]*matrix[4]*matrix[8]-matrix[0]*matrix[5]*matrix[7]-matrix[1]*matrix[3]*matrix[8]+matrix[1]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*matrix[7]-matrix[2]*matrix[4]*matrix[6]);//行列式
			
	dudx=(B1[0]*matrix[4]*matrix[8]-B1[0]*matrix[5]*matrix[7]-matrix[1]*B1[1]*matrix[8]+matrix[1]*matrix[5]*B1[2]+matrix[2]*B1[1]*matrix[7]-matrix[2]*matrix[4]*B1[2])/determinant;
	dvdy=(matrix[0]*B2[1]*matrix[8]-matrix[0]*matrix[5]*B2[2]-B2[0]*matrix[3]*matrix[8]+B2[0]*matrix[5]*matrix[6]+matrix[2]*matrix[3]*B2[2]-matrix[2]*B2[1]*matrix[6])/determinant;
	dwdz=(matrix[0]*matrix[4]*B3[2]-matrix[0]*B3[1]*matrix[7]-matrix[1]*matrix[3]*B3[2]+matrix[1]*B3[1]*matrix[6]+B3[0]*matrix[3]*matrix[7]-B3[0]*matrix[4]*matrix[6])/determinant;
	*/	

	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Yflag=OFF;
	int Zflag=OFF;
	for(int k=0;k<3;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B2[k]!=0) Yflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B3[k]!=0) Zflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<9;L++) matrix[L]=matrix_val[L];//matrixを変更前に戻す
	}
	else B1[0]=0; //flagがOFFならどのみちdudxはゼロ。
	
	if(Yflag==ON)
	{
		gauss(matrix,B2,N);
		for(int L=0;L<9;L++) matrix[L]=matrix_val[L];//matrixを変更前に戻す
	}
	else B2[1]=0;	//flagがOFFならどのみちdvdyはゼロ。
	
	if(Zflag==ON) gauss(matrix,B3,N);
	else B3[2]=0;	//flagがOFFならどのみちdwdzはゼロ。

	//計算終了

	//誤差を調査
	double Q[3]={0,0,0};//各方向の誤差
	double W=0;//重みの総和
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	/////*/

	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];

	double div=(dudx+dvdy+dwdz);

	delete [] weight;

	return div;
}

//divergence2における、3次元1次近似を行う関数ver.2（256行）
double calc_WLSM_divu_D3_order1_2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	///係数行列は
	///   ΣΔx2    ΣΔxΔy  ΣΔxΔz ΣΔx a = ΣΔxfj  
	///  ΣΔxΔy    ΣΔy2   ΣΔyΔz ΣΔy b = ΣΔyfj 
	///  ΣΔxΔz   ΣΔyΔz  ΣΔz2  ΣΔz  c = ΣΔzfj 
	///  ΣΔx      ΣΔy     ΣΔz     Σ1  d = Σfj

	double le=CON->get_distancebp();
	double matrix_val[16];			//matrixの値は一度gauss()で使用すると値が変わるので、変更前のmatrixを保存する
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
			
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);// leで割るのは打ち切り誤差防止
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		//double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		//double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		//double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double U=(PART[j].u[A_X]);
		double V=(PART[j].u[A_Y]);
		double W=(PART[j].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
					
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);
		weight[k]=w;
					
		matrix[0]+=X*X*w;			//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
		matrix[3]+=X*w;
					
		matrix[5]+=Y*Y*w;			//ΣYj^2wj
		matrix[6]+=Y*Z*w;		//ΣYjZjwj
		matrix[7]+=Y*w;

		matrix[10]+=Z*Z*w;			//ΣZj^2wj
		matrix[11]+=Z*w;

		matrix[15]+=w;
			
		B1[0]+=U*X*w;//ΣfjXjwj
		B1[1]+=U*Y*w;//ΣfjYjwj
		B1[2]+=U*Z*w;//ΣfjZjwj
		B1[3]+=U*w;//Σfjwj

		B2[0]+=V*X*w;//ΣfjXjwj
		B2[1]+=V*Y*w;//ΣfjYjwj
		B2[2]+=V*Z*w;//ΣfjZjwj
		B2[3]+=V*w;//Σfjwj

		B3[0]+=W*X*w;//ΣfjXjwj
		B3[1]+=W*Y*w;//ΣfjYjwj
		B3[2]+=W*Z*w;//ΣfjZjwj
		B3[3]+=W*w;//Σfjwj
	}
			
	matrix[4]=matrix[1];	
	matrix[8]=matrix[2];		
	matrix[9]=matrix[6];
	matrix[12]=matrix[3];
	matrix[13]=matrix[7];
	matrix[14]=matrix[11];

	matrix[15]+=1;//自分自身
	B1[3]+=PART[i].u[A_X];
	B2[3]+=PART[i].u[A_Y];
	B3[3]+=PART[i].u[A_Z];

	for(int L=0;L<16;L++) matrix_val[L]=matrix[L];//行列の値を保存

	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Yflag=OFF;
	int Zflag=OFF;
	for(int k=0;k<4;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B2[k]!=0) Yflag=ON;//B2[k]1がすべてゼロならXflagはOFFのまま。
		if(B3[k]!=0) Zflag=ON;//B3[k]1がすべてゼロならXflagはOFFのまま。
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<16;L++) matrix[L]=matrix_val[L];//matrixを変更前に戻す
	}
	else for(int k=0;k<4;k++) B1[k]=0; //flagがOFFならどのみちdudxはゼロ。
	
	if(Yflag==ON)
	{
		gauss(matrix,B2,N);
		for(int L=0;L<16;L++) matrix[L]=matrix_val[L];//matrixを変更前に戻す
	}
	else for(int k=0;k<4;k++) B2[k]=0;	//flagがOFFならどのみちdvdyはゼロ。
	
	if(Zflag==ON) gauss(matrix,B3,N);
	else for(int k=0;k<4;k++) B3[k]=0;	//flagがOFFならどのみちdwdzはゼロ。

	//計算終了

	//誤差を調査
	double Q[3]={0,0,0};//各方向の誤差
	double err[3];
	double W=1;//重みの総和
	err[A_X]=B1[3]-PART[i].u[A_X];//自身の誤差
	err[A_Y]=B2[3]-PART[i].u[A_Y];
	err[A_Z]=B3[3]-PART[i].u[A_Z];
	Q[A_X]+=err[A_X]*err[A_X];
	Q[A_Y]+=err[A_Y]*err[A_Y];
	Q[A_Z]+=err[A_Z]*err[A_Z];
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]);
		double V=(PART[j].u[A_Y]);
		double W=(PART[j].u[A_Z]);
		double w=weight[k];
		W+=w;
		
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]-U;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]-V;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]-W;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	/////*/

	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];

	double div=(dudx+dvdy+dwdz);

	delete [] weight;

	return div;
}

//divergence2における、3次元2次近似を行う関数（227行）
double calc_WLSM_divu_D3_order2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	double le=CON->get_distancebp();
	double matrix_val[81];			//matrixの値は一度gauss()で使用すると値が変わるので、変更前のmatrixを保存する
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
						
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
		weight[k]=w;

		matrix[0]+=X*X*w;		//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
		matrix[3]+=X*X*X*w;		//ΣXj^3wj
		matrix[4]+=X*Y*Y*w;		//ΣXjYj^2wj
		matrix[5]+=X*Z*Z*w;		//ΣXjZj^2wj
		matrix[6]+=X*X*Y*w;		//ΣXj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//ΣXjYjZjwj
		matrix[8]+=X*X*Z*w;		//ΣXj^2Zjwj
	
		matrix[10]+=Y*Y*w;		
		matrix[11]+=Y*Z*w;		
		matrix[12]+=X*X*Y*w;		
		matrix[13]+=Y*Y*Y*w;		
		matrix[14]+=Y*Z*Z*w;
		matrix[15]+=X*Y*Y*w;
		matrix[16]+=Y*Y*Z*w;
						
		matrix[20]+=Z*Z*w;			
		matrix[23]+=Z*Z*Z*w;		
					
		matrix[30]+=X*X*X*X*w;
		matrix[31]+=X*X*Y*Y*w;
		matrix[32]+=X*X*Z*Z*w;	
		matrix[33]+=X*X*X*Y*w;	
		matrix[34]+=X*X*Y*Z*w;	
		matrix[35]+=X*X*X*Z*w;	
					
		matrix[40]+=Y*Y*Y*Y*w;
		matrix[41]+=Y*Y*Z*Z*w;
		matrix[42]+=X*Y*Y*Y*w;
		matrix[43]+=Y*Y*Y*Z*w;
		matrix[44]+=X*Y*Y*Z*w;

		matrix[50]+=Z*Z*Z*Z*w;	//6行目
		matrix[51]+=X*Y*Z*Z*w;
		matrix[52]+=Y*Z*Z*Z*w;
		matrix[53]+=X*Z*Z*Z*w;

		//7～9行目はすべて既存の要素から転用が可能


		B1[0]+=U*X*w;		//a
		B1[1]+=U*Y*w;		//b
		B1[2]+=U*Z*w;		//c
		B1[3]+=U*X*X*w;		//d
		B1[4]+=U*Y*Y*w;		//e
		B1[5]+=U*Z*Z*w;		//f
		B1[6]+=U*X*Y*w;		//g
		B1[7]+=U*Y*Z*w;		//h
		B1[8]+=U*X*Z*w;		//i

		B2[0]+=V*X*w;		//a
		B2[1]+=V*Y*w;		//b
		B2[2]+=V*Z*w;		//c
		B2[3]+=V*X*X*w;		//d
		B2[4]+=V*Y*Y*w;		//e
		B2[5]+=V*Z*Z*w;		//f
		B2[6]+=V*X*Y*w;		//g
		B2[7]+=V*Y*Z*w;		//h
		B2[8]+=V*X*Z*w;		//i

		B3[0]+=W*X*w;		//a
		B3[1]+=W*Y*w;		//b
		B3[2]+=W*Z*w;		//c
		B3[3]+=W*X*X*w;		//d
		B3[4]+=W*Y*Y*w;		//e
		B3[5]+=W*Z*Z*w;		//f
		B3[6]+=W*X*Y*w;		//g
		B3[7]+=W*Y*Z*w;		//h
		B3[8]+=W*X*Z*w;		//i
		
	}
	matrix[9]=matrix[1];		//ΣXjYjwj
	matrix[17]=matrix[7];

	matrix[18]=matrix[2];
	matrix[19]=matrix[11];
	matrix[21]=matrix[8];
	matrix[22]=matrix[16];
	matrix[24]=matrix[7];
	matrix[25]=matrix[14];
	matrix[26]=matrix[5];

	matrix[27]=matrix[3];
	matrix[28]=matrix[12];
	matrix[29]=matrix[21];

	matrix[36]=matrix[4];
	matrix[37]=matrix[13];
	matrix[38]=matrix[22];
	matrix[39]=matrix[31];

	matrix[45]=matrix[5];
	matrix[46]=matrix[14];
	matrix[47]=matrix[23];
	matrix[48]=matrix[32];
	matrix[49]=matrix[41];

	matrix[54]=matrix[6];
	matrix[55]=matrix[15];
	matrix[56]=matrix[24];
	matrix[57]=matrix[33];
	matrix[58]=matrix[42];
	matrix[59]=matrix[51];
	matrix[60]=matrix[31];
	matrix[61]=matrix[44];
	matrix[62]=matrix[34];

	matrix[63]=matrix[7];
	matrix[64]=matrix[16];
	matrix[65]=matrix[25];
	matrix[66]=matrix[34];
	matrix[67]=matrix[43];
	matrix[68]=matrix[52];
	matrix[69]=matrix[61];
	matrix[70]=matrix[41];
	matrix[71]=matrix[51];

	matrix[72]=matrix[8];
	matrix[73]=matrix[17];
	matrix[74]=matrix[26];
	matrix[75]=matrix[35];
	matrix[76]=matrix[44];
	matrix[77]=matrix[53];
	matrix[78]=matrix[62];
	matrix[79]=matrix[71];
	matrix[80]=matrix[32];

	for(int L=0;L<81;L++) matrix_val[L]=matrix[L];//行列の値を保存
			
	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;
	int Yflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Zflag=OFF;
	for(int k=0;k<9;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B2[k]!=0) Yflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B3[k]!=0) Zflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_val[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) B1[k]=0;} //flagがOFFならどのみちdudxはゼロ。

	if(Yflag==ON) 
	{
		gauss(matrix,B2,N);
		for(int L=0;L<81;L++) matrix[L]=matrix_val[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) B2[k]=0;}	//flagがOFFならどのみちdvdyはゼロ。

	if(Zflag==ON) gauss(matrix,B3,N);
	else {for(int k=0;k<N;k++) B3[k]=0;}	//flagがOFFならどのみちdwdzはゼロ。

	////計算終了

	//誤差を調査
	double Q[3]={0,0,0};//各方向の標準偏差
	double W=0;		//重みの総和
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]*X*X+B1[4]*Y*Y+B1[5]*Z*Z+B1[6]*X*Y+B1[7]*Y*Z+B1[8]*X*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]*X*X+B2[4]*Y*Y+B2[5]*Z*Z+B2[6]*X*Y+B2[7]*Y*Z+B2[8]*X*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]*X*X+B3[4]*Y*Y+B3[5]*Z*Z+B3[6]*X*Y+B3[7]*Y*Z+B3[8]*X*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	//////


	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];
			
	double div=dudx+dvdy+dwdz;

	delete [] weight;

	return div;
}

//divergence2における、3次元2次近似を行う関数ver.2 未知数がひとつ多い（233行）
double calc_WLSM_divu_D3_order2_2(mpsconfig *CON,vector<mpselastic> &PART,double *matrix, double *B1,double *B2,double *B3,int i,int N)
{
	//P=aΔx+bΔy+cΔz+dΔx2+eΔy2+fΔz2+gΔxΔy+hΔyΔz+iΔzΔx+Pとおくと、
	///係数行列は
	///   ΣΔx2      ΣΔxΔy    ΣΔxΔz    ΣΔx3       ΣΔxΔy2    ΣΔxΔz2    ΣΔx2Δy     ΣΔxΔyΔz  ΣΔx2Δz    ΣΔx    a = ΣΔxf  
	///   ΣΔxΔy    ΣΔy2      ΣΔyΔz    ΣΔx2Δy    ΣΔy3       ΣΔyΔz2    ΣΔxΔy2     ΣΔy2Δz    ΣΔxΔyΔz  ΣΔy    b = ΣΔyf
	///   ΣΔxΔz    ΣΔyΔz    ΣΔz2      ΣΔx2Δz    ΣΔy2Δz    ΣΔz3       ΣΔxΔyΔz   ΣΔyΔz2    ΣΔxΔz2    ΣΔz    c = ΣΔzf
	///   ΣΔx3      ΣΔx2Δy   ΣΔx2Δz   ΣΔx4       ΣΔx2Δy2   ΣΔx2Δz2   ΣΔx3Δy     ΣΔx2ΔyΔz ΣΔx3Δz    ΣΔx2   d = ΣΔx2f
	///   ΣΔxΔy2   ΣΔy3      ΣΔy2Δz   ΣΔx2Δy2   ΣΔy4       ΣΔy2Δz2   ΣΔxΔy3     ΣΔy3Δz    ΣΔxΔy2Δz ΣΔy2   e = ΣΔy2f
	///   ΣΔxΔz2   ΣΔyΔz2   ΣΔz3      ΣΔx2Δz2   ΣΔy2Δz2   ΣΔz4       ΣΔxΔyΔz2  ΣΔyΔz3    ΣΔxΔz3	 ΣΔz2   f = ΣΔz2f
	///   ΣΔx2Δy   ΣΔxΔy2   ΣΔxΔyΔz ΣΔx3Δy    ΣΔxΔy3    ΣΔxΔyΔz2 ΣΔx2Δy2    ΣΔxΔy2Δz ΣΔx2ΔyΔz ΣΔxΔy g = ΣΔxΔyf
	///   ΣΔxΔyΔz ΣΔy2Δz   ΣΔyΔz2   ΣΔx2ΔyΔz ΣΔy3Δz    ΣΔyΔz3    ΣΔxΔy2Δz  ΣΔy2Δz2   ΣΔxΔyΔz2 ΣΔyΔz h = ΣΔyΔzf
	///   ΣΔx2Δz   ΣΔxΔyΔz ΣΔxΔz2   ΣΔx3Δz    ΣΔxΔy2Δz  ΣΔxΔz3   ΣΔx2Δy Δz ΣΔxΔyΔz2 ΣΔx2Δz2   ΣΔzΔx i = ΣΔxΔzf
	///   ΣΔx       ΣΔy ΣΔz ΣΔx2      ΣΔy2       ΣΔz        ΣΔxΔy     ΣΔyΔz      ΣΔxΔz     ΣΔzΔx     Σ1      P = Σfj

	double le=CON->get_distancebp();
	double matrix_val[100];			//matrixの値は一度gauss()で使用すると値が変わるので、変更前のmatrixを保存する
	double *weight=new double [PART[i].N];	//周辺粒子の重みを格納する。あとで誤差評価のとき再計算が不要になる。

	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double U=(PART[j].u[A_X]-PART[i].u[A_X]);
		double V=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double W=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double dis=sqrt(X*X+Y*Y+Z*Z);
						
		double w=1;
		if(dis>le) w=le*le*le*le/(dis*dis*dis*dis);//この式だとdis=Rのとき重みがゼロに近くなる
		weight[k]=w;

		matrix[0]+=X*X*w;		//ΣXj^2wj
		matrix[1]+=X*Y*w;		//ΣXjYjwj
		matrix[2]+=X*Z*w;		//ΣXjZjwj
		matrix[3]+=X*X*X*w;		//ΣXj^3wj
		matrix[4]+=X*Y*Y*w;		//ΣXjYj^2wj
		matrix[5]+=X*Z*Z*w;		//ΣXjZj^2wj
		matrix[6]+=X*X*Y*w;		//ΣXj^2Yjwj
		matrix[7]+=X*Y*Z*w;		//ΣXjYjZjwj
		matrix[8]+=X*X*Z*w;		//ΣXj^2Zjwj
		matrix[9]+=X*w;		//ΣXj^2Zjwj
	
		matrix[11]+=Y*Y*w;		
		matrix[12]+=Y*Z*w;		
		matrix[13]+=X*X*Y*w;		
		matrix[14]+=Y*Y*Y*w;		
		matrix[15]+=Y*Z*Z*w;
		matrix[16]+=X*Y*Y*w;
		matrix[17]+=Y*Y*Z*w;
		matrix[19]+=Y*w;
					
		matrix[22]+=Z*Z*w;			
		matrix[23]+=X*X*Z*w;
		matrix[24]+=Y*Y*Z*w;
		matrix[25]+=Y*Y*Y*w;
		matrix[29]+=Z*w;
					
		matrix[33]+=X*X*X*X*w;
		matrix[34]+=X*X*Y*Y*w;
		matrix[35]+=X*X*Z*Z*w;	
		matrix[36]+=X*X*X*Y*w;	
		matrix[37]+=X*X*Y*Z*w;	
		matrix[38]+=X*X*X*Z*w;	
					
		matrix[44]+=Y*Y*Y*Y*w;
		matrix[45]+=Y*Y*Z*Z*w;
		matrix[46]+=X*Y*Y*Y*w;
		matrix[47]+=Y*Y*Y*Z*w;
		matrix[48]+=X*Y*Y*Z*w;

		matrix[55]+=Z*Z*Z*Z*w;	//6行目
		matrix[56]+=X*Y*Z*Z*w;
		matrix[57]+=Y*Z*Z*Z*w;
		matrix[58]+=X*Z*Z*Z*w;

		matrix[99]+=w;
		//7～9行目はすべて既存の要素から転用が可能

		B1[0]+=U*X*w;		//a
		B1[1]+=U*Y*w;		//b
		B1[2]+=U*Z*w;		//c
		B1[3]+=U*X*X*w;		//d
		B1[4]+=U*Y*Y*w;		//e
		B1[5]+=U*Z*Z*w;		//f
		B1[6]+=U*X*Y*w;		//g
		B1[7]+=U*Y*Z*w;		//h
		B1[8]+=U*X*Z*w;		//i
		B1[9]+=U*w;

		B2[0]+=V*X*w;		//a
		B2[1]+=V*Y*w;		//b
		B2[2]+=V*Z*w;		//c
		B2[3]+=V*X*X*w;		//d
		B2[4]+=V*Y*Y*w;		//e
		B2[5]+=V*Z*Z*w;		//f
		B2[6]+=V*X*Y*w;		//g
		B2[7]+=V*Y*Z*w;		//h
		B2[8]+=V*X*Z*w;		//i
		B2[9]+=V*w;	

		B3[0]+=W*X*w;		//a
		B3[1]+=W*Y*w;		//b
		B3[2]+=W*Z*w;		//c
		B3[3]+=W*X*X*w;		//d
		B3[4]+=W*Y*Y*w;		//e
		B3[5]+=W*Z*Z*w;		//f
		B3[6]+=W*X*Y*w;		//g
		B3[7]+=W*Y*Z*w;		//h
		B3[8]+=W*X*Z*w;		//i
		B3[9]+=W*w;
	}
	matrix[10]=matrix[1];
	matrix[18]=matrix[7];

	matrix[20]=matrix[2];
	matrix[21]=matrix[12];
	matrix[24]=matrix[16];
	matrix[26]=matrix[7];
	matrix[27]=matrix[15];
	matrix[28]=matrix[5];

	for(int k=0;k<=2;k++) matrix[30+k]=matrix[3+10*k];//30～32要素
	matrix[39]=matrix[0];

	for(int k=0;k<=3;k++) matrix[40+k]=matrix[4+10*k];//40～43要素
	matrix[49]=matrix[11];

	for(int k=0;k<=4;k++) matrix[50+k]=matrix[5+10*k];//50～54要素
	matrix[59]=matrix[22];

	for(int k=0;k<=5;k++) matrix[60+k]=matrix[6+10*k];//60～65要素
	matrix[66]=matrix[34];
	matrix[67]=matrix[48];
	matrix[68]=matrix[37];
	matrix[69]=matrix[1];

	for(int k=0;k<=6;k++) matrix[70+k]=matrix[7+10*k];//70～76要素
	matrix[77]=matrix[54];
	matrix[78]=matrix[56];
	matrix[79]=matrix[12];
	
	for(int k=0;k<=7;k++) matrix[80+k]=matrix[8+10*k];//80～87要素
	matrix[88]=matrix[35];
	matrix[89]=matrix[20];

	for(int k=0;k<=8;k++) matrix[90+k]=matrix[9+10*k];//90～98要素

	matrix[99]+=1;//自身
	B1[9]+=PART[i].u[A_X];
	B2[9]+=PART[i].u[A_Y];
	B3[9]+=PART[i].u[A_Z];

	for(int L=0;L<100;L++) matrix_val[L]=matrix[L];//行列の値を保存
			
	//行列をガウスの消去法で解く　解はBに格納される
	int Xflag=OFF;
	int Yflag=OFF;//ガウスの消去法をするかしないか。解行列がゼロならガウスの消去法したらだめ。
	int Zflag=OFF;
	for(int k=0;k<10;k++)
	{
		if(B1[k]!=0) Xflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B2[k]!=0) Yflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
		if(B3[k]!=0) Zflag=ON;//B1[k]1がすべてゼロならXflagはOFFのまま。
	}

	if(Xflag==ON)
	{
		gauss(matrix,B1,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_val[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) B1[k]=0;} //flagがOFFならどのみちdudxはゼロ。

	if(Yflag==ON) 
	{
		gauss(matrix,B2,N);
		for(int L=0;L<100;L++) matrix[L]=matrix_val[L];//行列の値戻す
	}
	else {for(int k=0;k<N;k++) B2[k]=0;}	//flagがOFFならどのみちdvdyはゼロ。

	if(Zflag==ON) gauss(matrix,B3,N);
	else {for(int k=0;k<N;k++) B3[k]=0;}	//flagがOFFならどのみちdwdzはゼロ。

	////計算終了

	//誤差を調査
	double Q[3]={0,0,0};//各方向の標準偏差
	double W=0;		//重みの総和
	for(int k=0;k<PART[i].N;k++)
	{
		int j=PART[i].NEI[k];
		double X=(PART[j].r[A_X]-PART[i].r[A_X]);
		double Y=(PART[j].r[A_Y]-PART[i].r[A_Y]);
		double Z=(PART[j].r[A_Z]-PART[i].r[A_Z]);
		double dU=(PART[j].u[A_X]-PART[i].u[A_X]);
		double dV=(PART[j].u[A_Y]-PART[i].u[A_Y]);
		double dW=(PART[j].u[A_Z]-PART[i].u[A_Z]);
		double w=weight[k];
		W+=w;
		double err[3];
		err[A_X]=B1[0]*X+B1[1]*Y+B1[2]*Z+B1[3]*X*X+B1[4]*Y*Y+B1[5]*Z*Z+B1[6]*X*Y+B1[7]*Y*Z+B1[8]*X*Z-dU;
		err[A_Y]=B2[0]*X+B2[1]*Y+B2[2]*Z+B2[3]*X*X+B2[4]*Y*Y+B2[5]*Z*Z+B2[6]*X*Y+B2[7]*Y*Z+B2[8]*X*Z-dV;
		err[A_Z]=B3[0]*X+B3[1]*Y+B3[2]*Z+B3[3]*X*X+B3[4]*Y*Y+B3[5]*Z*Z+B3[6]*X*Y+B3[7]*Y*Z+B3[8]*X*Z-dW;
		Q[A_X]+=err[A_X]*err[A_X]*w;
		Q[A_Y]+=err[A_Y]*err[A_Y]*w;
		Q[A_Z]+=err[A_Z]*err[A_Z]*w;
	}
	for(int D=0;D<3;D++)
	{
		if(W>0) Q[D]/=W;
		Q[D]=sqrt(Q[D]);
		double speed=sqrt(PART[i].u[D]*PART[i].u[D]);
		if(speed>1e-6) Q[D]/=speed;
		else Q[D]=0;
	}
	//if(Q[A_X]>10 ||Q[A_Y]>10||Q[A_Z]>10) cout<<"Q("<<i<<")="<<Q[A_X]<<" "<<Q[A_Y]<<" "<<Q[A_Z]<<endl;
	//////


	double dudx=B1[0];
	double dvdy=B2[1];
	double dwdz=B3[2];
			
	double div=dudx+dvdy+dwdz;

	delete [] weight;

	return div;
}

//粒子が解析領域の外にでていないかチェック（70行）
int check_position(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,int *particle_number)
{
	int sw=OFF;	//本関数で返す値。OFFなら粒子数に変化なし。ONなら変化したという印
	int num=0;	//消滅する粒子数
	double le=CON->get_distancebp();
	int *flag=new int [fluid_number];	//ONなら領域内 OFFなら領域外

	double Xmax=CON->get_maxX(); double Xmin=CON->get_minX();
	double Ymax=CON->get_maxY(); double Ymin=CON->get_minY();
	double Zmax=CON->get_maxZ(); double Zmin=CON->get_minZ();

	double dx=CON->get_dx()*le;	//格子幅

	Xmax-=dx; Ymax-=dx; Zmax-=2*dx;		//保険をかけて1格子分内側に境界をとる。これより外側なら粒子を消す
	Xmin+=dx; Ymin+=dx; Zmin+=2*dx;

	vector<mpselastic>::iterator p,p0;//反復子
	p0=PART.begin();
	
	for(int i=0;i<fluid_number;i++) 
	{
		flag[i]=ON;
		if(PART[i].r[A_X]<Xmin || PART[i].r[A_X]>Xmax) flag[i]=OFF;
		else if(PART[i].r[A_Y]<Ymin || PART[i].r[A_Y]>Ymax) flag[i]=OFF;
		else if(PART[i].r[A_Z]<Zmin || PART[i].r[A_Z]>Zmax) flag[i]=OFF;
	}//flag[i]が求まった
	
	for(int i=0;i<fluid_number;i++) if(PART[i].N<=4) flag[i]=OFF;//周辺粒子数が少ない粒子も削除?

	for(int i=0;i<fluid_number;i++) if(flag[i]==OFF) num++;


	if(num>0)//領域外粒子を検知したなら
	{
		sw=ON;
		int *erase_id=new int[num];
		int count=0;
		for(int i=0;i<fluid_number;i++)
		{
			if(flag[i]==OFF)
			{
				erase_id[count]=i;//消すべきidを記憶
				count++;
			}
			
		}
	
		for(int i=0;i<num;i++)
		{
			p=PART.begin();
			p+=erase_id[i];
			
			PART.erase(p);
			for(int j=i+1;j<num;j++)
			{
				erase_id[j]=erase_id[j]-1;//1つ値をさげる
			}
		}

		delete [] erase_id;

		//idがずれたから戻す
		for(int i=0;i<(int) PART.size();i++) if(PART[i].ID!=i) PART[i].ID=i; 
	}
	
	if(sw==ON)
	{
		cout<<"領域外粒子を探知 "<<num<<"個の粒子を消去--現在の粒子数="<<PART.size()<<endl;
		*particle_number=*particle_number-num;
	}
	
	delete [] flag;

	return sw;
}

//粒子が近づきすぎるのを防ぐ関数（63行）
void modify_position(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double dt)
{
	double le=CON->get_distancebp();

	for(int i=0;i<fluid_number;i++)
	{
		//if(PART[i].surface==ON)
		{
			double mindis=le;
			//mindis=100;			//距離の短いものを遠ざけるだけでなく、遠いものをleにしたいときはこっち
			int J=i;			//最近接粒子
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				//if(PART[j].type==FLUID)// && j>i)
				{
					//if(PART[j].surface==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
						if(dis<mindis)
						{
							mindis=dis;
							J=j;
						}
					}
				}
			}
			if(J!=i && (PART[J].type==FLUID || PART[J].type==ELASTIC || PART[J].type==MAGELAST || PART[J].type==MAGELAST2))//leより近接している流体粒子があったなら
			{
				double L=le-mindis;//開くべき距離
				double dL[DIMENSION];
				for(int D=0;D<DIMENSION;D++) dL[D]=PART[J].r[D]-PART[i].r[D];

				for(int D=0;D<DIMENSION;D++)
				{
					double dU=0.5*L/dt;	//変化すべき速度
				//	PART[J].u[D]+=dL[D]/mindis*dU;
					PART[J].r[D]+=dL[D]/mindis*dU*dt;

				//	PART[i].u[D]-=dL[D]/mindis*dU;
					PART[i].r[D]-=dL[D]/mindis*dU*dt;
				}				
			}
			else if(J!=i && (PART[J].type!=FLUID || PART[J].type!=ELASTIC || PART[J].type!=MAGELAST || PART[J].type!=MAGELAST2))//leより近接している壁粒子があったなら
			{
				double L=le-mindis;//開くべき距離
				double dL[DIMENSION];
				for(int D=0;D<DIMENSION;D++) dL[D]=PART[J].r[D]-PART[i].r[D];

				for(int D=0;D<DIMENSION;D++)
				{
					double dU=L/dt;	//変化すべき速度
				//	PART[i].u[D]-=dL[D]/mindis*dU;
					PART[i].r[D]-=dL[D]/mindis*dU*dt;
				}				
			}
		}
	}
}


//粒子体積計算関数（18行）
double get_volume(mpsconfig *CON)
{
	double V=0;//体積
	double le=CON->get_distancebp();
	if(CON->get_model_set_way()==0)	//正方格子のとき
	{
		if(CON->get_dimension()==2){V=le*le;}
		else	{V=le*le*le;}
	}
	else if(CON->get_model_set_way()==1)	//細密格子のとき
	{
		if(CON->get_dimension()==2){V=sqrt(3.0)/2*le*le;}
		else V=le*le*le/sqrt(2.0);
	}	
	else cout<<"モデルの積み方が不定です 体積を計算できません"<<endl;
	return V;
}

//粒子数密度のずれ測定関数（18行）
void output_particle_density(mpsconfig *CON,vector<mpselastic> &PART,int fluid_number,double n0,int particle_number,int t)
{
	///初期粒子数密度分布を出力する

	ofstream fp("initial_n0.dat");
	double le=CON->get_distancebp();

	if(CON->get_dimension()==2) {for(int i=0;i<particle_number;i++) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].PND<<endl;}
	else if(CON->get_dimension()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].r[A_Y]>-0.5*le && PART[i].r[A_Y]<0.5*le) fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<PART[i].PND<<endl;
		}
	}
	fp.close();
}


//初期配置と解析領域の関係調査関数（16行）
void check_initial_position(mpsconfig *CON,vector<mpselastic> &PART)
{
	double region[3][2];
	region[A_X][0]=CON->get_minX(); region[A_X][1]=CON->get_maxX(); 
	region[A_Y][0]=CON->get_minY(); region[A_Y][1]=CON->get_maxY(); 
	region[A_Z][0]=CON->get_minZ(); region[A_Z][1]=CON->get_maxZ(); 

	for(int i=0;i<PART.size();i++) 
	{
		int flag=OFF;
		for(int D=0;D<3;D++) if(PART[i].r[D]>=region[D][1] || PART[i].r[D]<=region[D][0]) flag=ON;
		
		if(flag==ON)	cout<<"error:解析領域外に初期粒子(id="<<i<<")が存在します x,y,z="<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
	}
}


//25行
void FEM3D(mpsconfig &CON, vector<mpselastic> &PART, double **F, int *static_node, int *static_nelm, vector<point3D> &static_NODE, vector<element3D> &static_ELEM, int t, double TIME, int fluid_number)
{
	size_t particle_number=PART.size();
	double dt=CON.get_dt();
	if(CON.get_EM_interval()==1 || t==1 || (t-1)%CON.get_EM_interval()==0){
		FEM3D_calculation(CON, static_node, static_nelm, static_NODE, static_ELEM, t, TIME, PART, fluid_number, particle_number, dt, F);
	}else{ //何ステップかに一度だけ電磁場計算を行う場合
	
		if(CON.get_dir_for_P()!=2 && CON.get_dir_for_P()!=3)//電磁力が圧力ディリクレ値のときはここでは計算しない
		{
			double F_val;
			ifstream fin5("FEM_interval.dat");
			if(!fin5) cout<<"cannot open FEM_interval.dat"<<endl;
			fin5.unsetf(ifstream::dec);
			fin5.setf(ifstream::skipws);
			////電磁力読み取り
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].type!=WALL)
				{
					for(int D=0;D<3;D++)
					{
						fin5>>F_val;
						F[D][i]+=F_val;//F[D][i]は宣言した直後にゼロセットしてある
					}
				}
			}
			fin5.close();
		}
	}
}


void initial_model_input(mpsconfig *CON, int *particle_number, double *TIME)
{
	if(CON->get_restart()==OFF)
	{
		//モデル作成と全粒子数の計算。
		set_initial_placement_using_MD(CON, particle_number);//set_initial_placement_using_MD(&CON,&particle_number);//分子動力学によるモデルセット。計算時間はかかるが、表面をきれいに表現
	}
    else if(CON->get_restart()==ON)
    {
		ifstream fin("number.dat");
		if(!fin) cout<<"number.dat cannot be opened" <<endl;
		fin.unsetf(ifstream::dec);
		fin.setf(ifstream::skipws);
    
		fin>>*particle_number;
		TIME=0;//fin>>*TIME;
		fin.close();
		cout<<"前回における解析を引き継ぎ"<<endl;
    }
	
    cout<<"粒子数は"<<*particle_number<<endl;
}

//TetGenによるメッシュ生成
void TetGenInterface(mpsconfig &CON, vector<mpselastic> &PART, double **F, int fluid_number, double dt, int t, int particle_number, double n0, double TIME)
{
	
	vector<int>	flag_hyper;

	//HYPERELAST適用のための工夫
	for(int i=0;i<particle_number;i++)
	{
		if(PART[i].type==HYPERELAST)
		{
			PART[i].type=MAGELAST;
			flag_hyper.push_back(ON);
		}
		else
		{
			flag_hyper.push_back(OFF);
		}
	}

	/////////////////////////////////////////////////
	if(CON.get_EM_interval()==1 || t==1 || (t-1)%CON.get_EM_interval()==0){

			usingTetGen(CON,PART,F,fluid_number,particle_number,t,TIME); //delaun3D関係の変数は全て不要

			if(CON.get_Ferror()==ON){ //エラーなら
				double F_val;
			ifstream fin5("FEM_interval.dat");
			if(!fin5) cout<<"cannot open FEM_interval.dat"<<endl;
			fin5.unsetf(ifstream::dec);
			fin5.setf(ifstream::skipws);
			////電磁力読み取り
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].type==MAGELAST)
				{//MAGELAST
						for(int D=0;D<3;D++)
						{
							fin5>>F_val;
							F[D][i]+=F_val;//F[D][i]は宣言した直後にゼロセットしてある
						}
					}
				}
			fin5.close();			
			}
	}
	else
	{ //何ステップかに一度だけ電磁場計算を行う場合
	
		if(CON.get_dir_for_P()!=2 && CON.get_dir_for_P()!=3)//電磁力が圧力ディリクレ値のときはここでは計算しない
		{
			double F_val;
			ifstream fin5("FEM_interval.dat");
			if(!fin5) cout<<"cannot open FEM_interval.dat"<<endl;
			fin5.unsetf(ifstream::dec);
			fin5.setf(ifstream::skipws);
			////電磁力読み取り
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].type==MAGELAST)
				{//MAGELAST
					for(int D=0;D<3;D++)
					{
						fin5>>F_val;
						F[D][i]+=F_val;//F[D][i]は宣言した直後にゼロセットしてある
					}
				}
			}
			fin5.close();
		}
	}
	//////////////////////////////////////////////////*/
/*	//電磁力プロット
	plot_F(CON, PART, fluid_number,t);
	if(CON.get_F_interval()>0)
	{
		if(t==1 || t%CON.get_F_interval()==0) plot_F_log(CON,PART,fluid_number,t);
	}

	//電磁力AVSファイル出力
	if(CON.get_avs_eforce_interval()>0)
	{
		if(t==1 || t%CON.get_avs_eforce_interval()==0) plot_avs_eforce(CON,PART,fluid_number,t);
	}//*/


	//HYPERELAST適用のための工夫
	for(int i=0;i<particle_number;i++)		if(flag_hyper[i]==ON) 	PART[i].type=HYPERELAST;

}

void file_initialization()
{
	//////////////////////////records//////////////////////////////////////
	//hamiltonian file
	ofstream init1("./Elastic/hamiltonian.dat", ios::trunc);
	ofstream init2("./Elastic/kinetic_energy.dat", ios::trunc);
	ofstream init3("./Elastic/elastic_energy.dat", ios::trunc);
	ofstream init4("./Elastic/potential_energy.dat", ios::trunc);
	
	ofstream init5("aveave_P_history.txt", ios::trunc);
	ofstream init6("node.dat", ios::trunc);
	ofstream init7("time_log.dat", ios::trunc);
	ofstream init8("A-t.dat", ios::trunc);
	ofstream init9("longZ.dat", ios::trunc);
	system("rmdir /s /q Mesh");
	system("rmdir /s /q Lorentz");
	system("rmdir /s /q Speed");
	system("rmdir /s /q FluxContour");
	system("rmdir /s /q FluxAVS");
	
	system("rmdir /s /q Residual");
	system("rmdir /s /q Pressure");
	system("rmdir /s /q Current");
	/////////////////////////make file///////////////////////////////////////
	system("mkdir Mesh");
	system("mkdir Lorentz");
	system("mkdir Speed");
	system("mkdir FluxContour");
	system("mkdir FluxAVS");

	system("mkdir Residual");
	system("mkdir Pressure");
	system("mkdir Current");

	//超弾性
	ofstream init10("sum_P.csv", ios::trunc);
	ofstream init11("P_X.csv", ios::trunc);
	ofstream init12("P_Y.csv", ios::trunc);
	ofstream init13("P_Z.csv", ios::trunc);
	ofstream init14("P_corner.csv", ios::trunc);
	ofstream init15("d_P_X.csv", ios::trunc);
	ofstream init16("d_P_Y.csv", ios::trunc);
	ofstream init17("d_P_Z.csv", ios::trunc);
	ofstream init18("h_P_X.csv", ios::trunc);
	ofstream init19("h_P_Y.csv", ios::trunc);
	ofstream init20("h_P_Z.csv", ios::trunc);
	ofstream init21("lambda.csv", ios::trunc);
	ofstream init22("pnd.csv", ios::trunc);
	ofstream init23("E.csv", ios::trunc);
	ofstream init24("E_T.csv", ios::trunc);
	ofstream init25("E_g.csv", ios ::trunc);
	ofstream init26("E_W.csv", ios::trunc);
	ofstream init27("E_lam.csv", ios::trunc);
	system("mkdir Position");
	system("mkdir Hy_stress");
	system("mkdir Fi");
	system("mkdir Newton_raphson");
	system("mkdir Wall");

	//close file
	init1.close();
	init2.close();
	init3.close();
	init4.close();
	init5.close();
	init6.close();
	init7.close();
	init8.close();
	init9.close();
	init10.close();
	init11.close();
	init12.close();
	init13.close();
	init14.close();
	init15.close();
	init16.close();
	init17.close();
	init18.close();
	init19.close();
	init20.close();
	init21.close();
	init22.close();
	init23.close();
	init24.close();
	init25.close();
	init26.close();
	init27.close();
}

void Make_STL(){
			double normal[368][3];
			double cooda[368][3];
			double coodb[368][3];
			double coodc[368][3];
			double normal_no[368][3];
			double cood1[368][3];
			double cood2[368][3];
			double cood3[368][3];
		/////////////////coordinate file///////////////
		ifstream fin1("cood.txt");
		if(!fin1) cout<<"cood.txt cannot be opened."<<endl;
		fin1.unsetf(ifstream::dec);
		fin1.setf(ifstream::skipws);

		for(int i=0;i<368;i++)
		{       
			for(int D=0;D<3;D++) fin1>>normal[i][D];
			for(int D=0;D<3;D++) fin1>>cooda[i][D];
			for(int D=0;D<3;D++) fin1>>coodb[i][D];
			for(int D=0;D<3;D++) fin1>>coodc[i][D];
			
		}
		fin1.close();	
		//////////////////////////////////////////////////
		/////////////////non coordinate file/////////////
		ifstream fin2("noncood.txt");
		if(!fin2) cout<<"noncood.txt cannot be opened."<<endl;
		fin2.unsetf(ifstream::dec);
		fin2.setf(ifstream::skipws);
		for(int i=0;i<368;i++)
		{       
			for(int D=0;D<3;D++) fin2>>normal_no[i][D];
			for(int D=0;D<3;D++) fin2>>cood1[i][D];
			for(int D=0;D<3;D++) fin2>>cood2[i][D];
			for(int D=0;D<3;D++) fin2>>cood3[i][D];
			
		}
		fin2.close();	
		//////////////////////////////////////////////////

		///////////////足し合わせ/////////////////////////
		for(int i=0;i<368;i++){
			for(int D=0;D<3;D++) cood1[i][D]+=cooda[i][D];
			for(int D=0;D<3;D++) cood2[i][D]+=coodb[i][D];
			for(int D=0;D<3;D++) cood3[i][D]+=coodc[i][D];
		}
		/////////////////////////////////////////////////

		////////////////STLファイル出力//////////////////
		ofstream fout1("COIL.STL");
		if(!fout1) cout<<"COIL.STL cannot be opened."<<endl;
		fout1.precision(6);
	
		fout1<<"COIL"<<endl;
		for(int i=0;i<368;i++)
		{       
			fout1<<"   facet normal";
			for(int D=0;D<3;D++) fout1<<" "<<normal[i][D];
			fout1<<endl<<"      outer loop"<<endl;
			fout1<<"         vertex";
			for(int D=0;D<3;D++) fout1<<" "<<cooda[i][D];
			fout1<<endl<<"         vertex";
			for(int D=0;D<3;D++) fout1<<" "<<coodb[i][D];
			fout1<<endl<<"         vertex";
			for(int D=0;D<3;D++) fout1<<" "<<coodc[i][D];
			fout1<<endl<<"      endloop"<<endl;
			fout1<<"   endfacet"<<endl;
			
		}
		fout1.close();	
		/////////////////////////////////////////////////
}