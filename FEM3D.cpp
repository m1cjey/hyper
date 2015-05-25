#include "stdafx.h"	
#include "PostProcess.h"

#define FULL 3
#define REMESH 4
#define FULL_IMPORT 5

using namespace std;

/*********用語解説
静的要素：リメッシュ時に再配置しない要素（初期状態のものを使い続ける）
動的要素：リメッシュ時に再配置する要素
*****************/
//要素の深さ決定関数
void set_depth(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *depth,int KTE);
void non_linear_Avector3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *V,int *jnb,int *branch_num,double **current,double *RP,double II,int *depth,double **Be,int t);
//電流ベクトル
void carrent_vector(mpsconfig &CON, vector<point3D> &NODE, vector<element3D> &ELEM, int nelm, double **B,int t);
void FEM3D_calculation(mpsconfig &CON, int *static_node, int *static_nelm, vector<point3D> &static_NODE, vector<element3D> &static_ELEM, int t,double TIME, vector<mpselastic> &PART, int fluid_number, int particle_number, double dt, double **F)
{
	int delaun_flag;//デローニ分割を行うか、行わないか

	if(CON.get_mesh_input()==0)//MPSTOFEMによる節点、要素生成
	{
		if(CON.get_remesh_sw()==OFF) delaun_flag=FULL;	//remesh領域を想定せず、常にすべてをデローニ分割
		else if(CON.get_remesh_sw()==ON)
		{
			if(t==1) delaun_flag=FULL;	//全モデルをデローニ分割
			else delaun_flag=REMESH;	//remesh領域のみデローニ分割
		}
	}

	if(delaun_flag==FULL) cout<<"FULL デローニ分割実行"<<endl;
	else if(delaun_flag==REMESH) cout<<"remesh領域のみデローニ分割実行"<<endl;
	else if(delaun_flag==FULL_IMPORT) cout<<"Magnet生成ファイルより要素情報等読み込み"<<endl;

	vector<point3D> NODE;
	vector<element3D> ELEM;
	int node=0;					//全節点数
	int nelm=0;					//全要素数
	int KTJ;					//最大節点数
	int KTE;					//最大要素数　3次元式もとめよ
	const double err=1.0e-14;			//誤差判定のしきい値・・・1e-14マシンイプシロンより厳しいのでは？
	
	if(delaun_flag==FULL)
	{
		//粒子配置より節点配置を入手
		MPS_TO_FEM3D(CON, &node, NODE, PART, fluid_number, particle_number);

		KTJ=node;
		if(CON.get_fine()!=OFF) KTJ+=CON.get_add_points();
		KTE=12*KTJ;
	}else if(delaun_flag==REMESH){
		node=(int)static_NODE.size()-1;			//静的節点数
		nelm=(int)static_ELEM.size()-1;
		KTJ=node+particle_number;				//このあと動的節点(流体)を格納しないといけないから、KTJを増加
		if(CON.get_fine()!=OFF) KTJ+=CON.get_add_points();
		KTE=12*KTJ;
	}

	if(delaun_flag==FULL)//すべてをデローニ分割
	{
		point3D NODE0;
		element3D ELEM0;
		for(int i=KTJ+8;i>node;i--) NODE.push_back(NODE0);//10^5個とかintの最大値超えるのでは？？？	//これは空気節点情報の追加？15/2/4
		for(int i=0;i<KTE;i++) ELEM.push_back(ELEM0);//配列を確保
		int FINE_sw=CON.get_fine();//自動再分割スイッチ

		//デローニ分割
		delaun3D_main(CON,NODE,ELEM, KTJ, KTE,&node,&nelm, FINE_sw);
		cout<<"要素数="<<nelm<<" 節点数＝"<<node<<endl;

		//メッシュ生成のポスト処理
		double *val=new double[KTJ+1];
		for(int i=1;i<=node;i++) val[i]=NODE[i].material;
//		if(t==1 || t%(CON.get_EM_interval()*CON.get_mesh_output_interval())==0) 
		{
			data_avs(node,nelm,NODE,ELEM,KTJ,val,CON);
			data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t);//断面図
			data_avs3(node,nelm,NODE,ELEM,CON);//材質	
		}

		if(CON.get_remesh_sw()==ON)//remesh領域を想定するなら、静的節点・要素情報を記憶しておく
		{
			//NODE,ELEMのうち、動かない要素、節点だけをstatic_NODE,staticELEMに格納する
			static_ELEM.clear();
			static_NODE.clear();
			memorize_static_NODE_and_ELEM(CON,NODE,ELEM,static_NODE,static_ELEM, node, nelm);

			/*/チェック
			int snode=(int) static_NODE.size()-1;
			int snelm=(int) static_ELEM.size()-1;
			for(int i=1;i<=snode;i++) val[i]=static_NODE[i].remesh;
			data_avs2(CON,snode,snelm,static_NODE,static_ELEM,KTJ,t);//断面図*/
		}
		delete [] val;
	}

//メッシュが切れてから力の計算まで

	//節点-要素関係
    int *jnb=new int[node+1];//各節点に隣接する要素数格納
    set_jnb3D(NODE,ELEM,node,nelm,jnb);

	int **nei=new int* [node+1];//各節点の周辺要素番号格納
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	for(int i=1;i<=node;i++)
	{
		if(jnb[i]==0)
		{
			//cout<<"jnb=0 i="<<i<<" material="<<NODE[i].material<<" particle="<<NODE[i].particleID<<endl;
			NODE[i].boundary_condition=1;//境界条件をディリクレ型にすることで、ICCGに参加させない
			//if(NODE[i].material==FLUID) NODE[i].particleID=-1;//流体節点が消失する場合、必要な値が計算されないため、粒子にフィードバックできない(してはいけない)
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC) if(NODE[i].particleID>=0) cout<<"suf="<<PART[NODE[i].particleID].surface<<endl; 
		}
	}

	//FEM
	if(CON.get_EM_calc_type()==1 || CON.get_EM_calc_type()==4) potential_calculation(CON,NODE,ELEM, node, nelm,jnb, TIME,PART, fluid_number,nei,F);
	if(CON.get_EM_calc_type()==2) calc_static_magnetic_field(CON, node, nelm,NODE,ELEM,jnb, dt,TIME, t,nei, KTE,PART,fluid_number,F,KTJ);
	if(CON.get_EM_calc_type()==3) calc_transitional_EM_field(CON, node, nelm,NODE,ELEM,jnb, dt, TIME,t,nei, KTE,PART,fluid_number,F,KTJ);
	if(CON.get_EM_calc_type()==5) calc_variable_magnetic_field(CON, node, nelm,NODE,ELEM,jnb, dt,TIME, t,nei, KTE,PART,fluid_number,F);

	//FEM_intervalによる計算時間短縮の場合
    if(CON.get_EM_interval()>1)
    {
		//ofstream bb("FEM_interval.dat");///電磁力出力　次ステップはこれを読み取る
        FILE *b=fopen("FEM_interval.dat","w");///電磁力出力　次ステップはこれを読み取る
		for(int i=0;i<fluid_number;i++)
		{
			//bb<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
			fprintf(b,"%1.15lf %1.15lf %1.15lf\n",F[A_X][i],F[A_Y][i],F[A_Z][i]);
		}
		fclose(b);
		//bb.close();
    }

	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;
}

//渦電流を考慮にいれた節点要素のベクトルポテンシャル計算関数ver.2 クーロンゲージを適用して左辺をラプラシアンとする
void Avector3D_node_eddy2(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double **A,int *jnb,int *depth,int **nei2,int *branch_num,double **old_A,double II,double dt,int mps_num,double *V,int t)
{	cout<<"渦電流を考慮にいれた節点要素によるﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始"<<endl;

	double u0=4*PI*0.0000001;			//空気の透磁率
    double v0=1/u0;						//磁気抵抗率
	double j0x,j0y,j0z;					//電流密度
	double sigma=CON.get_ele_conduc();	//電気伝導率
	
	double *current[3];					//各要素の電流密度[A/m3]格納
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];

	//////////////////////ｺｲﾙがなければ電流計算を行う必要がない。そこでｺｲﾙ節点数をかぞえる
	int coil_node_num=0;	//ｺｲﾙ節点数
	for(int i=1;i<=node;i++) if(NODE[i].material==COIL) coil_node_num++;
	////////////////////////*/

	int *save_bound=new int [node+1];//本関数は何度も呼び出されるので、電流に関する境界条件を完全に消去できない。そこで保存する
	
	if(coil_node_num>0)///ｺｲﾙ節点があるなら電流密度計算
	{
		for(int i=1;i<=node;i++) save_bound[i]=NODE[i].boundary_condition;//境界条件を保存

		if(II!=0)//電流値が非ゼロなら電流計算
		{
			calc_current(CON,NODE,ELEM,SIDE,node,nelm,side_num,jnb,branch_num,current,depth,II);
			cout<<"電流計算完了"<<endl;
		}
		else	//ゼロなら初期化
		{
			for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) for(int D=0;D<3;D++) current[D][i]=0;
		}
		///電流境界条件の初期化（もう電流はもとまっているからいらない)
		for(int i=1;i<=side_num;i++) if(SIDE[i].boundary_condition>=10) SIDE[i].boundary_condition=0;
		for(int i=1;i<=node;i++) if(NODE[i].boundary_condition>=10) NODE[i].boundary_condition=0;
	}
	carrent_vector(CON, NODE, ELEM, nelm, current,t);

    int NN=0;//ディリクレ型境界節点数
    int *dn=new int [node+1]; //各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
	double *PHAT[3];
	for(int D=0;D<3;D++) PHAT[D]=new double [CON.get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
    for(int i=1;i<=node;i++)
    {
        if(jnb[i]==0)
		//if(NODE[i].boundary_condition==3)
		{    
			//cout<<"Jnb["<<i<<"]=0"<<endl;
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==2)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else if(NODE[i].boundary_condition==1)
		{    
	        dn[i]=NN;
	        PHAT[A_X][NN]=0;
			PHAT[A_Y][NN]=0;
			PHAT[A_Z][NN]=0;
			A[A_X][i]=0;
			A[A_Y][i]=0;
			A[A_Z][i]=0;
	        NN++;
		}
		else
		{
			dn[i]=node+1;
		}
    }/////////////*/
    cout<<"ﾃﾞｨﾘｸﾚ数＝"<<3*NN<<endl;//実際に計算しなくていいﾃﾞｨﾘｸﾚ型節点数は3*NNである
	    
	///Ａφ法には導体(流体)節点数の数だけ電位に関する未知数が存在する。それを数える

	int conducter_num=0;//導体(流体)節点数
	for(int i=1;i<=node;i++) if( NODE[i].material==WATER && NODE[i].boundary_condition==0 &&jnb[i]!=0) conducter_num++;
	//////////////
    int pn=3*(node-NN)+conducter_num;///未知数
    int *ppn=new int [pn];		///行列でn番目の節点は節点番号ppn[n]
    int *npp=new int [node+1];	///各節点が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int i=1;i<=node;i++)//行列はA1x,A1y,A1z,A1φ,A2x,A2y・・・の順に格納
    {
        if(NODE[i].boundary_condition==0 && jnb[i]!=0)
		//if(NODE[i].boundary_condition==0)
		{
			ppn[num]=i;
			ppn[num+1]=-1;
			ppn[num+2]=-2;
			npp[i]=num;
			num+=3;
			if( NODE[i].material==WATER)
			{
				//cout<<"節点番号"<<i<<endl;
				ppn[num]=-3;
				num++;
			}
		}
		else npp[i]=pn+1;
    }
    
    ////行列の最大幅計算  もっとも幅の大きいのはφの要素だろうと仮定している。確実ではないことに注意
    int mat_w=0;
	for(int i=1;i<=node;i++) if(branch_num[i]+1>mat_w && NODE[i].material==WATER) mat_w=branch_num[i]+1;
	mat_w=4*mat_w;//X,Y,Z,φ成分とあるので×4 //左辺をﾗﾌﾟﾗｼｱﾝにしても、φに関しては４成分のままに注意
    //////
//	mat_w=200;//空芯の場合渦電流が流れない
    
    ////配列確保
    double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    ////
    
    /////////全体行列を作成する
    unsigned int time=GetTickCount();
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int J,J2,J3,flag,flag2,flag3;
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
	
		delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
		
		if(ELEM[je].material==COIL)
		{
			j0x=current[A_X][je];
			j0y=current[A_Y][je];
			j0z=current[A_Z][je];
		}
		else
		{
			j0x=0;
			j0y=0;
			j0z=0;
		}

		///比透磁率
		double v=v0;
		if(ELEM[je].material==WATER) v/=CON.get_RP();
//		if(ELEM[je].material==IRON) v/=5000;
		double rp=u0;
		if(ELEM[je].material==WATER) rp*=CON.get_RP();
//		if(ELEM[je].material==IRON) rp*=5000;
		
		for(int n=1;n<=4;n++)
		{
			if(NODE[N[n]].boundary_condition==0)
			{
				/////X方向

				int I=npp[N[n]]+1;
				//各Ｂはここでまとめて計算する
				B[I-1]+=delta/4*j0x;//x,y,zは１つずつずれる
				B[I]+=delta/4*j0y;//
				B[I+1]+=delta/4*j0z;
				
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
							    flag=1;
							}
						}
						if(flag==0)//相当するＪが存在しなかったら作る
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aixの項
						    ROW[I][H]=J;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_X][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}

				/////Y方向
				I=npp[N[n]]+1+1;/////Y方向
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						flag=0;
						flag2=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J2==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aiyの項
							    flag2=1;
							}
						}
						if(flag2==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=J2;	
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Y][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
				/////*/

				/////Z方向
				I=npp[N[n]]+2+1;
				for(int m=1;m<=4;m++)
				{	
					if(NODE[N[m]].boundary_condition==0)
					{
						J=npp[N[m]]+1;	//Aixの項
						J2=J+1;			//Aiyの項
						J3=J2+1;		//Aizの項
						flag3=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J3==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;//Aizの項
							    flag3=1;
							}
						}
						if(flag3==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
							G[I][H]+=(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
						    ROW[I][H]=J3;
						}
					}
					else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
					{
					    int NN=dn[N[m]];
					    B[I-1]-=PHAT[A_Z][NN]*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m])*delta6/rp;
					}
				}
			}
		}
	}
	
	//////渦電流項計算
	int J4,flag4;
	
	for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==WATER)//MAGELASTには渦電流が流れない
		{
			
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
			double Xs=0;//重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
		
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
			//ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
	
			delta6=1/delta6/6;//計算に必要なのは1/(36V)なのでここで反転しておく(1/36V^2)*V=1/36V
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1.0;
					c[i]*=-1.0;
					d[i]*=-1.0;
					e[i]*=-1.0;
				}
			}
			/////////
	
			double Sx=0;//渦電流項のうち、１ｽﾃｯﾌﾟ前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙを考慮することで得られる、定ﾍﾞｸﾄﾙ
			double Sy=0;
			double Sz=0;

			double co=delta/20.0/dt*sigma;//σV/(20dt)　何度も計算することになるので係数化する
			
			double co2=delta6*sigma;//   σ/(36V)
			double co3=sigma*dt*delta6;

			for(int n=1;n<=4;n++)
			{
				if(NODE[N[n]].boundary_condition==0)
				{
					/////X方向

					int I=npp[N[n]]+1;
					Sx=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J4=J+3;			//φの項
							flag=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)//Aixの項
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==J) G[I][h]+=co*2;
									else G[I][h]+=co;
									flag=1;
								}
							
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
									flag4=1;
								}///
							}
							if(flag==0)//相当するＪが存在しなかったら作る(ここではそんなはずないけど)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
						    
								if(I==J) G[I][H]+=co*2;
								else G[I][H]+=co;
								ROW[I][H]=J;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*c[m];
								ROW[I][H]=J4;
							}

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							if(I==J) Sx+=2*Ax;
							else Sx+=Ax;
							
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							//int NN=dn[N[m]];
							//B[I-1]-=(d[n]*d[m]+e[n]*e[m])*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-(d[n]*c[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*c[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}
					B[I-1]+=co*Sx;//支配方程式のX成分から得られるBの値

					/////Y方向

					I=npp[N[n]]+1+1;/////Y方向
					Sy=0;

					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J4=J+3;			//φの項
							flag2=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==J2) G[I][h]+=co*2;//Aiyの項
									else G[I][h]+=co;
									flag2=1;
								}

								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
									flag4=1;
								}
							}
							if(flag2==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								if(I==J2) G[I][H]+=co*2;
								else G[I][H]+=co;
								ROW[I][H]=J2;
							}
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*d[m];
								ROW[I][H]=J4;
							}//*/

							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							if(I==J2) Sy+=2*Ay;
							else Sy+=Ay;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*d[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+e[n]*e[m])*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=-(e[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					B[I-1]+=co*Sy;//支配方程式のX成分から得られるBの値

					/////Z方向
					Sz=0;
					I=npp[N[n]]+2+1;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(I==J3) G[I][h]+=co*2;
									else G[I][h]+=co;
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
									flag4=1;
								}
							}
							
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    if(I==J3) G[I][H]+=co*2;
								else G[I][H]+=co;
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co2*(b[n]+Xs*c[n]+Ys*d[n]+Zs*e[n])*e[m];
								ROW[I][H]=J4;
							}//*/
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							if(I==J3) Sz+=2*Az;
							else Sz+=Az;
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}////////*/
					B[I-1]+=co*Sz;//支配方程式のX成分から得られるBの値

					/////φ
					
					I=npp[N[n]]+3+1;
					Sx=0;
					Sy=0;
					Sz=0;
					for(int m=1;m<=4;m++)
					{	
						if(NODE[N[m]].boundary_condition==0)
						{
							J=npp[N[m]]+1;	//Aixの項
							J2=J+1;			//Aiyの項
							J3=J2+1;		//Aizの項
							J4=J+3;			//φの項
							flag=0;
							flag2=0;
							flag3=0;
							flag4=0;
							for(int h=1;h<=NUM[I];h++)
							{
								//Aixの項
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									///co2*c[n]*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])というふうにc[n]を前に持ってくると、打ち切り誤差の関係で非対称と認識されるので、上の式と同様に最後にもってきた
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
									flag=1;
								}
								//Aiyの項
								if(J2==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
									flag2=1;
								}
								//Aizの項
								if(J3==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
									flag3=1;
								}
								//φの項
								if(J4==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
									flag4=1;
								}
							}
							if(flag==0)
							{   
								//cout<<I<<" "<<J<<" co2="<<co2<<"  c[m]="<<c[m]<<" "<<co2*c[m]<<" B  n,m="<<n<<" "<<m<<endl;
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    //G[I][H]+=co2*c[m];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*c[n];
							    ROW[I][H]=J;
							}
							if(flag2==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
								G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*d[n];
							    ROW[I][H]=J2;
							}
							if(flag3==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    G[I][H]+=co2*(b[m]+Xs*c[m]+Ys*d[m]+Zs*e[m])*e[n];
							    ROW[I][H]=J3;
							}
							
							if(flag4==0)//相当するＪが存在しなかったら作る
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								G[I][H]+=co3*(c[n]*c[m]+d[n]*d[m]+e[n]*e[m]);
								ROW[I][H]=J4;
							}
							///Bの計算 //仮にN[m]が固定境界上だった場合もSxの計算に組み込まないといけないが、ディリクレ値は０と仮定して計算を除外
							double Ax=old_A[A_X][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							double Ay=old_A[A_Y][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							double Az=old_A[A_Z][N[m]];//節点N[m]の１ｽﾃｯﾌﾟ 前のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙX成分
							Sx+=Ax*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*c[n];
							Sy+=Ay*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*d[n];
							Sz+=Az*(b[m]+c[m]*Xs+d[m]*Ys+e[m]*Zs)*e[n];
						}
						else //N[m]がﾃﾞｨﾘｸﾚ型境界節点なら
						{
							int NN=dn[N[m]];
							//B[I-1]-=-c[n]*e[m]*delta6*PHAT[A_X][NN]/rp;
							//B[I-1]-=-d[n]*e[m]*delta6*PHAT[A_Y][NN]/rp;
							//B[I-1]-=(c[n]*c[m]+d[n]*d[m])*delta6*PHAT[A_Z][NN]/rp;
						}
					}//////*/
					B[I-1]+=co2*(Sx+Sy+Sz);
				}
			}
		}
	}//////*/
	
    for(int D=0;D<3;D++) delete [] current[D];
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
		    for(int m=1;m<j;m++)
		    {
		        if(ROW[i][j]<ROW[i][m])
				{
				    tempG=G[i][m];
				    tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }//////////*/

	/*for(int i=1;i<=node;i++)
	{
		if(NODE[i].boundary_condition==0)
		{
			//if(branch_num[i]==0) cout<<i<<endl;
			int I=npp[i]+1;
			//if((branch_num[i]+1)*4!=NUM[I]) cout<<branch_num[i]<<" "<<NUM[I]<<" i="<<i<<" I="<<I<<endl;
		}
	}
	//for(int k=1;k<=NUM[4];k++) if(ppn[ROW[4][k]-1]>0) cout<<ppn[ROW[4][k]-1]<<" ";
	//cout<<endl;
	*/

	///行列の実際の最大幅を求める
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	cout<<"最大幅："<<maxN<<"/"<<mat_w;
	
	ofstream fout("matrix.dat");
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			fout<<G[i][j]<<" ";
			//cout<<ROW[i][j]<<" ";
		}
		fout<<endl;
		//cout<<endl;
	}
	fout.close();/////////*/

	///対称性チェック
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					if(G[i][j]!=G[J][k]) 
					{
						cout<<"非対称"<<G[i][j]<<" "<<G[J][k]<<endl;
						G[i][j]=G[J][k];
					}
				}
			}
		}
	}///////////*/
	
    
    double *val = new double [number];
    int *ind = new int [number];//非ゼロ要素の列番号格納配列
    int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
    
    /////////////////////val,ind ,ptrに値を格納
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
    cout<<"  行列作成終了  time="<<GetTickCount()-time<<endl;
    
	//CG法実行
	double *XX=new double [pn];//行列の答え格納
//    if(CON.get_FEMCG()==0) CG3D_Avector(val,ind,ptr,pn,ppn,B,A);//CG法実行
//	else 
//	{
		if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
		else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
		else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
		else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
		for(int n=0;n<pn;n++)
		{
			int i=ppn[n];
			if(i>0) A[A_X][i]=XX[n];
			else if(i==-1) A[A_Y][ppn[n-1]]=XX[n];
			else if(i==-2) A[A_Z][ppn[n-2]]=XX[n];
			else if(i==-3) V[ppn[n-3]]=XX[n];//電位φ
		}	
//	}
	delete [] XX;



	/*//old_Aの更新 //渦電流計算のところでdA/dtを計算したいから、ここではold_Aはそのままにしておく
	for(int i=1;i<=mps_num;i++)
	{
		for(int D=0;D<3;D++) old_A[D][i]=A[D][i];
	}///*/

	///old_Aのﾌｧｲﾙ出力(流体粒子のみ出力。最初のmps_num番目までが流体粒子)
	ofstream g("old_A.dat");
	for(int i=1;i<=mps_num;i++) g<<old_A[A_X][i]<<" "<<old_A[A_Y][i]<<" "<<old_A[A_Z][i]<<endl;
	g.close();
	
	if(coil_node_num>0)//境界条件をもとに戻す(次のｽﾃｯﾌﾟで再び電流を計算するかもしれないから)
	{	
		for(int i=1;i<=node;i++) NODE[i].boundary_condition=save_bound[i];
	}

    delete [] dn;
    for(int D=0;D<3;D++) delete [] PHAT[D];
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] save_bound;
	
}

//電位、磁位などのポテンシャル計算関数
void potential_calculation(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,double TIME,vector<mpselastic> &PART,int fluid_number,int **nei,double **F)
{
	double *V=new double [node+1];	//potential
    
    double *Ee[3];					//要素内電界 or 要素内磁界
    for(int D=0;D<3;D++) Ee[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//各要素の透磁率または誘電率格納（現在では誘電率は使っていない）

	double fluid_rp=CON.get_r_perm();		//誘電率
	if(CON.get_EM_calc_type()==4) fluid_rp=CON.get_RP();	//透磁率
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2) RP[i]=fluid_rp;
		else if(ELEM[i].material==ELASTIC) RP[i]=1;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL)
		{
			if(CON.get_EM_calc_type()==1) RP[i]=1000;//本当は導体だから∞
			else if(CON.get_EM_calc_type()==4) RP[i]=1;//透磁率
		}
		else cout<<"error:材質に対し、誘電率あるいは透磁率が不定"<<endl;
	}

	//電位解決
	VOLT3D(CON,NODE,ELEM, node, nelm,V,jnb, TIME,PART, fluid_number,nei,RP);

	//電界計算
	ELECTRO3D(CON,NODE,ELEM, node, nelm,V,Ee);

	//静電力計算関数 表面に働く応力はP=0.5epE^2 これは誘電体を導体近似したもの
	//conductive_approximation(CON,NODE,ELEM, node, nelm,Ee, jnb, nei, RP,PART,F,fluid_number);

	

	delete [] V;
    for(int D=0;D<3;D++) delete [] Ee[D];
	delete [] RP;
}

///電位計算関数
void VOLT3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,int *jnb,double TIME,vector<mpselastic> &PART,int fluid_number,int **nei,double *RP)
{
    double V1=0;
    double V2=CON.get_V();			//電位
	double u0=0.000001256;			//真空の透磁率
	double ep0=8.854e-12;			//真空の誘電率
	double le=CON.get_distancebp();//粒子間距離
	unsigned timeA=GetTickCount();	//計算開始時刻

	////////電位決定
	if(CON.get_EM_calc_type()==1)//電場計算なら
	{
		double tau=CON.get_dt()*CON.get_V_step();//時定数
		if(CON.get_V_con()==2)
		{
			V2=CON.get_V()*(1-exp(-TIME/tau));//電位を指数関数的に増加させる
		}
		else if(CON.get_V_con()==1)
		{
			V2=(CON.get_V()-CON.get_initial_V())*TIME/(CON.get_dt()*CON.get_V_step())+CON.get_initial_V();//電位を直線的に増加させる
		}
		if(V2==0) V2=1;//0だとエラーになるから1にする
		else if(V2>CON.get_V()) V2=CON.get_V();
		cout<<"電位計算開始 V="<<V2<<" ";
	}
	else if(CON.get_EM_calc_type()==4)//磁位計算
	{
		double R=1;//比率
		V1=0;
		if(CON.get_uniform_B_sw()==OFF) V2=CON.get_magnet_B()/u0*CON.get_magnet_H();
		if(CON.get_uniform_B_sw()==ON)  V2=CON.get_uniform_B()/u0*(CON.get_ZU()-CON.get_ZD());//一様磁場
		if(CON.get_V_con()==1) 
		{
			V2=(V2)*TIME/(CON.get_dt()*CON.get_V_step());
			R=TIME/(CON.get_dt()*CON.get_V_step());
			if(V2>CON.get_magnet_B()/u0*CON.get_magnet_H()) {V2=CON.get_magnet_B()/u0*CON.get_magnet_H();R=1;}
		}
		if(V2==0) V2=1;//0だとエラーになるから1にする
		
		cout<<"磁位計算開始 V2="<<V2<<"("<<R*CON.get_magnet_B()<<"T)"<<endl;
	}
	////////////////////////////

	///磁位計算の場合、解析領域の境界条件を自由境界条件にする。（ここで、MPSTOFEMの段階でそのようにしたらFINE3Dがうまくいかない）
	if(CON.get_EM_calc_type()==4)//磁位計算
	{
		if(CON.get_uniform_B_sw()==OFF)//通常はこっち
		{
			for(int i=1;i<=nelm;i++)
			{
				if(ELEM[i].material==AIR)
				{
					for(int k=1;k<=4;k++)
					{
						int kelm=ELEM[i].elm[k];
						if(kelm==0)
						{
							int j=k%4+1;//ielmとjelmの接する三角形を構成する節点番号 
							int m=4-(k-1)/2*2;
							int n=3-(k/2%2)*2;
							NODE[ELEM[i].node[j]].boundary_condition=0;//境界条件の抹消
							NODE[ELEM[i].node[m]].boundary_condition=0;
							NODE[ELEM[i].node[n]].boundary_condition=0;
						}
					}
				}
			}
		}
		if(CON.get_uniform_B_sw()==ON)//一様磁場の場合
		{
			double ZU=CON.get_ZU();		//解析領域上端
			double ZD=CON.get_ZD();		//解析領域下端
			double err=1e-14;
			for(int i=1;i<=node;i++)
			{
				if(NODE[i].r[A_Z]>ZU-err) NODE[i].boundary_condition=2;//上端
				else if(NODE[i].r[A_Z]<ZD+err) NODE[i].boundary_condition=1;//下端
				else NODE[i].boundary_condition=0;
			}
		}
	}////////*/


    int NN=0;					//ディリクレ型境界節点数
    int *dn=new int [node+1];	//各節点がディリクレ型境界の何番目か。ディリクレ型境界上でないならnode+1を格納
    double *PHAT=new double [CON.get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
    for(int i=1;i<=node;i++)
    {
		if(NODE[i].boundary_condition==1)//上で境界条件を書き換えてるから、この文はelse ifにすること
		{    
	        dn[i]=NN;
	        PHAT[NN]=V1;
	        V[i]=V1;
	        NN++;
		}
		else if(NODE[i].boundary_condition==2)
		{   
	        dn[i]=NN;
	        PHAT[NN]=V2;
	        V[i]=V2;
	        NN++;
		}
		else dn[i]=node+1;
    }/////////////*/
    cout<<"ディリクレ数＝"<<NN<<" ";
	    
    int pn=node-NN;///未知数
    int *ppn=new int [pn];		//行列でn番目の未知数は節点番号ppn[n]に相当
    int *npp=new int [node+1];	//i番目の節点は行列のnpp[i]番目に相当。ディリクレ型の場合は行列に入っていないのでpn+1を格納
    int num=0; 
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].boundary_condition==0)	//未知数なら
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
    
    ////行列の最大幅計算
    int mat_w=0;
	mat_w=calc_matrix_width(CON,NODE,ELEM,node,nelm,jnb,nei);//メモリをきっちりもとめる。だけど若干の計算コストになる。
    //////
    
    ////配列確保
    double **G=new double *[pn+1];						///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1];							///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1];							///各行の、非ゼロ要素数
    
    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    ////
	
	/*//解行列Ｂに節点の電荷を代入
	if(CON.get_charge()==1 && CON.get_FEM_calc_type()==1)
	{
		for(int k=1;k<=WATER_N;k++)//初めのWATER_N個の節点はTRANS[k]の粒子に相当
		{
			if(NODE[k].boundary_condition==0)//未知数ならBにその節点の領域があるから代入
			{
				int i=TRANS[k];//節点kに相当する粒子番号
				for(int n=1;n<=jnb[k];n++)
				{
					int jelm=nei[k][n];//節点kの隣接する要素
					int N[5];
					for(int j=1;j<=4;j++) N[j]=ELEM[jelm].node[j];
		
					///デローニ分割の際に求めた体積は、スーパーボックス内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
					ELEM[jelm].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
					B[k]+=PART[i].angle/ep0/4*ELEM[jelm].volume;
					//if(PART[i].angle!=0)cout<<"EE"<<endl;
				}
				
			}
		}
	}///////////*/
	
	double *charge=new double [nelm+1];
	for(int n=1;n<=nelm;n++) charge[n]=0;
	
    /////////全体行列を作成する
    
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
    for(int je=1;je<=nelm;je++)
    {   
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
		
		///デローニ分割の際に求めた体積は、スーパーボックス内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
		
		double delta6=ELEM[je].volume;//体積の6倍
		
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]))*delta6;//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]))*delta6;//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]))*delta6;//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
	
		/////比誘電率を定義
		double ep=RP[je];
		/*if(ELEM[je].material==FLUID)
		{
			if(CON.get_EM_calc_type()==1) ep=CON.get_r_perm();
			else if(CON.get_EM_calc_type()==4) ep=CON.get_RP();//磁位計算だから比透磁率
		}*/
	
		////要素マトリクス作成開始
		for(int i=1;i<=4;i++)
		{
			if(NODE[N[i]].boundary_condition==0)///未知なら
			{   
				int I=npp[N[i]]+1;///節点N[i]は行列のI番目
				for(int j=1;j<=4;j++)
				{					
					int J=npp[N[j]]+1;///節点N[j]は行列のJ番目
					if(NODE[N[j]].boundary_condition==0)///未知なら
					{   
						int flag=0;
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
								//if(I==85839 && h==1) cout<<h<<" "<<G[I][h]<<" "<<delta<<" "<<(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])<<endl;
								flag=1;
							}
						}
						if(flag==0)
						{   
							NUM[I]=NUM[I]+1;
							int H=NUM[I];
			    
							G[I][H]+=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta;
							ROW[I][H]=J;
						}
					}
					else //N[j]がディリクレ型境界節点なら
					{
						int n=dn[N[j]];
						B[I-1]-=ep*(c[i]*c[j]+d[i]*d[j]+e[i]*e[j])*delta*PHAT[n];
					}
				}
				//if(CON.get_charge()==ON) B[I-1]+=charge[je]/4*delta;
			}
		}   
    }
    ///////////////////////*/
    
	///実際の最大幅を計算
	double realwidth=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>realwidth) realwidth=NUM[i];
    
    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];

    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//非ゼロ要素の列番号格納配列
    int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
    
    /////////////////////val,ind ,ptrに値を格納
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
	cout<<"行列作成終了 幅:"<<realwidth<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    
	///////////////////////行列計算開始
	
	double *XX=new double [pn];//行列の答え格納
	if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
	else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		V[i]=XX[n];//電位φ
	}
	delete [] XX;
	////////////////////////////
    
    ///////電位をファイル出力
	ofstream fp("V.dat");
    for(int i=1;i<=node;i++)
    {
        if(NODE[i].r[A_Y]<2*le && NODE[i].r[A_Y]>-2*le) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<V[i]<<endl;
    }
	fp.close(); 
    
    ////////////////////////

    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;

	delete [] charge;
}

//行列の幅計算関数
int calc_matrix_width(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei)
{
	//係数行列のなかで、i行目に含まれる非ぜろ要素数は、計算点iからのびる辺の数(すなわち近隣計算点数)+1(自身)である。これの最大値を求める
	///配列確保
		
	int *temp_check=new int[node+1];	//一時的なﾁｪｯｸ配列
		
	///初期化
	for(int i=1;i<=node;i++) temp_check[i]=0;
	/////////////////////

	int maxwidth=0;

	//if(CON.get_FEM_calc_type()==1 || CON.get_FEM_calc_type()==4 || CON.get_iteration_count()==0)
	{
		for(int i=1;i<=node;i++)
		{
			if(NODE[i].boundary_condition==0)
			{
				int width=0;
				vector<int> NEI2;//各節点の隣接する節点番号格納(neiは節点-要素であるが、NEI2は節点-節点)
				for(int k=1;k<=jnb[i];k++)
				{
					int jelm=nei[i][k];//節点iに隣接する要素番号
					for(int j=1;j<=4;j++)
					{
						int p=ELEM[jelm].node[j];
						if(p!=i && temp_check[p]==0 && NODE[p].boundary_condition==0)
						{	
							width++;
							NEI2.push_back(p);
							temp_check[p]=1;//もう見た
						}
					}
				}
				for(size_t k=0;k<NEI2.size();k++) temp_check[NEI2[k]]=0;//初期化
				if(width>maxwidth) maxwidth=width;
			}
		}
	}

	delete [] temp_check;

	return maxwidth+1;//widthは自分から伸びる辺の数に等しい。行列の幅は自分自身もいれてwidth+1
}

//行列並び替え関数
void arrange_matrix(int pn,int *NUM,int **ROW,double **G)
{
	///G[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
			for(int m=1;m<j;m++)
			{
				if(ROW[i][j]<ROW[i][m])
				{
					tempG=G[i][m];
					tempR=ROW[i][m];
					G[i][m]=G[i][j];


					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////
}

//ICCG法
void ICCG3D2(mpsconfig &CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	int count;						//数え上げ変数
	double accel=CON.get_CGaccl();	//加速ファクタ
	
	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
	
	double *matrix=new double[num2];//係数行列を保存(対角成分を含む下三角行列) 非ゼロ要素だけを1次配列として保存
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
				matrix[num2]=val[m];
				ind2[num2]=ind[m];
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];			//列方向にみた、各列の要素数
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];						//ゼロ要素の値 VAL[i][k]はi列のk番目の非ゼロ要素
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];							//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg法
	double alp,beta;
	double rLDLt_r;
	double E=1;//誤差
    double *r=new double[pn];
	double *AP = new double [pn];
	double *P = new double [pn];
	double *y=new double [pn];
	double *LDLt_r= new double [pn];
	double *L=new double[num2];//不完全コレスキー分解後の下三角行列L格納
	double *D1 = new double [pn];//D行列
	
	/////不完全コレスキー分解
	Incomplete_Cholesky_Decomposition(CON,L,D1,matrix,ptr2,ind2,pn,num2);//LとD1に値が書き込まれる
	
	delete [] matrix;

	///列を基準にした配列に値を代入
	for(int k=0;k<pn;k++) NUM[k]=0;
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			VAL[J][NUM[J]]=L[m];
			IND[J][NUM[J]]=k;
			NUM[J]=NUM[J]+1;
		}
	}////////*/

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0;
		r[n]=B[n];
	}

	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
	for(int i=0;i<pn;i++)
	{
		if(i==0) y[0]=r[0]/L[0]; //式（3.77） 
		else
		{
		    double sum=0;
			for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=L[m]*y[ind2[m]];//式（3.78）
		    int m=ptr2[i+1]-1;
			y[i]=(r[i]-sum)/L[m];
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

    cout<<"ICCG法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	count=0;
	double ep=CON.get_EMCGep();//収束判定
	while(E>ep)
	{
		count++;
		if(count==pn) cout<<"count=pn"<<endl;
		//////////////alpを求める
		rLDLt_r=0;

		for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];
		
		for(int n=0;n<pn;n++)
		{
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
		}
		double PAP=0;
		for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];
		alp=rLDLt_r/PAP;
		//cout<<"alp="<<alp<<" PAP="<<PAP<<endl;
		//////////////////////
	
		//////////////// X(k+1)=X(k)+alp*P
		for(int n=0;n<pn;n++) X[n]+=alp*P[n];
		
		//////////////// r=r-alp*AP
		for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
		
		//////////////////誤差
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];
		
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		////////////////////////
		
		///////////////////////beta
		beta=1.0/rLDLt_r;
		rLDLt_r=0;
		
        /////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
		for(int i=0;i<pn;i++)
		{
			if(i==0) y[0]=r[0]/L[0]; //式（3.77） 新
			else
			{
			    double sum=0;
			    for(int m=ptr2[i];m<ptr2[i+1]-1;m++)//対角成分は除くからptr[i+1]-1
			    {
					sum+=L[m]*y[ind2[m]];//式（3.78）
			    }
			    int m=ptr2[i+1]-1;
				y[i]=(r[i]-sum)/L[m];
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
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	
    delete [] r;
	delete [] AP;
	delete [] P;

	delete [] y;
	delete [] LDLt_r;
	delete [] D1;

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

	delete [] L;
}



//不完全コレスキー分解
void Incomplete_Cholesky_Decomposition(mpsconfig &CON,double *L,double *D1,double *matrix,int *ptr2,int *ind2,int pn,int num2)
{
	//matrix[i]:係数行列のうち、対角を含む下三角が格納されている
	//ptr2[k]:matrix[i]に対するアクセス配列
	//L[i],D[i]:不完全コレスキー分解格納配列
	//pn:未知数
	//num2:matrix[i]の要素数

	cout<<"不完全コレスキー分解開始---";
	unsigned int timeC=GetTickCount();
	double accel=CON.get_CGaccl();		//加速係数
	double accel2=0;					//自動計算された加速係数
	double maxvalue=0;					//matrix中の最大値
	int Line=1;							//最大値を持つ行番号
	int Column=1;						//最大値を持つ列番号
	
	//1回目:まずは普通に不完全コレスキー分解を行う。そこでD1[i]<0なら加速係数を修正し、2回目を行う
	int onemoreflag=OFF;
	L[0]=matrix[0]*accel;
	D1[0]=1/L[0];			//L[0]は係数行列の値と変わらないので、ゼロでもなければ負でもないと予想される
	if(D1[0]<0) cout<<"error in 不完全コレスキー分解 D1[0]<0"<<endl;
	maxvalue=D1[0]*D1[0];
	for(int k=1;k<pn;k++)	//k=0はもう済んだ。k>=1からloop
	{	
		for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0からk-1までの要素
		{
			int i=ind2[m];//列番号
		
			double sum=0;
			for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i行目のなかから列の一致するものを探している。少し手間か？
			L[m]=matrix[m]-sum;//i=0のときはsum=0となり自動的にL[m]=matrix[m]となる
			if(L[m]*L[m]>maxvalue)
			{
				maxvalue=L[m]*L[m];
				Line=k; Column=i;
			}
		}
		int m=ptr2[k+1]-1;				//このmは必ず対角成分に該当する
				
		double akk=matrix[m];			//係数行列の対角成分.値を保存.
		double sum=0;
		for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
		L[m]=matrix[m]*accel-sum;		//加速係数をかける
		if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkkが非常に小さいとき、このようにする。ここで、どのみちLkk<0なら収束しづらいから、Lkkが負の極小値でも正の値にしてしまってよいのでは？
		D1[k]=1/L[m];
		
		if(L[m]<0)						//L[m]<0のときD1[m]<0となり、収束が極端に遅くなる.これを防ぐために加速係数を大きくする
		{
			accel2=sum/akk*1.1;
			if(accel2>accel) accel=accel2;
			onemoreflag=ON;				//このスイッチがONなら、不完全コレスキー分解をもう一度行う。
		}
		if(L[m]*L[m]>maxvalue)
		{
			maxvalue=L[m]*L[m];
			Line=k; Column=k;
		}
	}
	if(Line!=Column) cout<<"value="<<maxvalue<<" Line="<<Line<<" Column="<<Column<<endl;
	if(onemoreflag==ON)
	{
		L[0]=matrix[0]*accel;
		D1[0]=1/L[0];			//L[0]は係数行列の値と変わらないので、ゼロでもなければ負でもないと予想される
		for(int k=1;k<pn;k++)	//k=0はもう済んだ。k>=1からloop
		{	
			for(int m=ptr2[k];m<ptr2[k+1]-1;m++)///0からk-1までの要素
			{
				int i=ind2[m];//列番号
		
				double sum=0;		//accelの変更によりLとD1の値が変わったので、sumも計算しなおし。
				for(int j=ptr2[k];j<m;j++) for(int J=ptr2[i];J<ptr2[i+1];J++) if(ind2[J]==ind2[j]) sum+=L[j]*L[J]*D1[ind2[j]];//i行目のなかから列の一致するものを探している。少し手間か？
				L[m]=matrix[m]-sum;//i=0のときはsum=0となり自動的にL[m]=matrix[m]となる
			}
			int m=ptr2[k+1]-1;				//このmは必ず対角成分に該当する
				
			double akk=matrix[m];			//係数行列の対角成分.値を保存.
			double sum=0;
			for(int j=ptr2[k];j<m;j++) sum+=L[j]*L[j]*D1[ind2[j]];
			L[m]=matrix[m]*accel-sum;		//加速係数をかける
			if(L[m]*L[m]<1e-8)L[m]=0.0001;	//Lkkが非常に小さいとき、このようにする。ここで、どのみちLkk<0なら収束しづらいから、Lkkが負の極小値でも正の値にしてしまってよいのでは？
			D1[k]=1/L[m];
		}
	}	
	cout<<"加速係数="<<accel<<" time="<<(GetTickCount()-timeC)*0.001<<"[sec]"<<endl;
	///不完全コレスキー分解完了/////////*/
}

//行列の対称性検査関数
void check_matrix_symmetry(int pn,int *NUM,int **ROW,double **G)
{
	ofstream sym("symmetry_out.csv");
	for(int i=1;i<=pn;i++)
	{
		int flag_sym=0;
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			int flag=0;
			for(int k=1;k<=NUM[J];k++)
			{
				if(ROW[J][k]==i)
				{
					flag=1;
					if(G[i][j]!=G[J][k])
					{
						flag_sym=ON;
						cout<<"対称性エラー "<<i<<" "<<G[i][j]<<" "<<G[J][k]<<endl;							
					}
				}
			}
			if(flag==0)cout<<"DDD"<<endl;
		}
		if(flag_sym==ON)	sym<<i<<endl;
	}
	sym.close();
}

//電界計算関数
void ELECTRO3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,double **Ee)
{
    //　Ee:要素の電界 or 要素の磁界B
    
	if(CON.get_EM_calc_type()==1) cout<<"電界計算-------";
	if(CON.get_EM_calc_type()==4) cout<<"磁界計算(磁位)------";
	unsigned timeA=GetTickCount();		//計算開始時刻
	double ep0=8.854e-12;				//真空の誘電率
	double u0=12.56e-7;					//真空の透磁率
    
    int N[4+1];							//要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];

	//double maxE=0;
	//double Xc[3];
    ///要素の電界・磁界をもとめる
    for(int je=1;je<=nelm;je++)
    {
        for(int j=1;j<=4;j++)
		{
			N[j]=ELEM[je].node[j];
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
		}
	
		double delta6=ELEM[je].volume;//体積の6倍  体積は電位求めるときに計算してある
	
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
		
		///係数c,d,e計算
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]))*delta6;//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]))*delta6;//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]))*delta6;//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}   
		}
		/////////
    
		Ee[A_X][je]=c[1]*V[N[1]]+c[2]*V[N[2]]+c[3]*V[N[3]]+c[4]*V[N[4]];
        Ee[A_Y][je]=d[1]*V[N[1]]+d[2]*V[N[2]]+d[3]*V[N[3]]+d[4]*V[N[4]];
		Ee[A_Z][je]=e[1]*V[N[1]]+e[2]*V[N[2]]+e[3]*V[N[3]]+e[4]*V[N[4]];
	
		for(int D=0;D<3;D++) Ee[D][je]*=-1;//電界は電位勾配のマイナス

		double EE=sqrt(Ee[A_X][je]*Ee[A_X][je]+Ee[A_Y][je]*Ee[A_Y][je]+Ee[A_Z][je]*Ee[A_Z][je]);
		/*if(EE>maxE) 
		{
			maxE=EE;
			Xc[A_X]=(X[1]+X[2]+X[3]+X[4])*0.25;
			Xc[A_Y]=(Y[1]+Y[2]+Y[3]+Y[4])*0.25;
			Xc[A_Z]=(Z[1]+Z[2]+Z[3]+Z[4])*0.25;
		}*/
    }///電界計算終了

	//cout<<"MAX="<<maxE<<" x,y,z="<<Xc[A_X]<<" "<<Xc[A_Y]<<" "<<Xc[A_Z]<<endl;
    
	/*if(CON.get_FEM_calc_type()==4)//磁位による磁場解析の場合
	{
		///E[D][je]に格納されているのはHの値なので、これに透磁率をかける
		
		for(int je=1;je<=nelm;je++)
		{
			if(ELEM[je].material==WATER) for(int D=0;D<3;D++) Ee[D][je]*=u0*CON.get_RP();
			else for(int D=0;D<3;D++) Ee[D][je]*=u0;
		}
	}///////////////////////*/

    /////ファイル出力
	if(CON.get_EM_calc_type()==1 && CON.get_E_times()>0)
	{
		ofstream fp("E.dat");
		double times=CON.get_E_times();
		double le=CON.get_distancebp();
		double Xmin=CON.get_XL()+le; double Xmax=CON.get_XR()-le;//解析領域
		double Zmin=CON.get_ZD()+le; double Zmax=CON.get_ZU()-le;
		double dx=3*le;
		int Nx=(int)((Xmax-Xmin)/dx);//各方向の分割数
		int Nz=(int)((Zmax-Zmin)/dx);
		int search=nelm;//locate関数で最初に探索する要素番号

		if(CON.get_plot_E_type()==1)
		{
			for(int n=0;n<Nx;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=dx*n+Xmin;//出力する点の座標
					double yp=0;
					double zp=dx*m+Zmin;
					int loc=locate3D(NODE,ELEM,search,xp,yp,zp);

					fp<<xp<<" "<<zp<<" "<<-Ee[A_X][loc]*times<<" "<<-Ee[A_Z][loc]*times<<endl;///見ずらいからわざと＋－反転して出力
					search=loc;
				}
			}
		}
		else if(CON.get_plot_E_type()==2)//スカラー
		{
			for(int n=0;n<Nx;n++)
			{
				for(int m=0;m<Nz;m++)
				{
					double xp=dx*n+Xmin;//出力する点の座標
					double yp=0;
					double zp=dx*m+Zmin;
					int loc=locate3D(NODE,ELEM,search,xp,yp,zp);

					fp<<xp<<" "<<zp<<" "<<sqrt(Ee[A_X][loc]*Ee[A_X][loc]+Ee[A_Z][loc]*Ee[A_Z][loc])<<endl;
					
					//fp<<xp<<" "<<zp<<" "<<sqrt(Ee[A_X][loc]*Ee[A_X][loc]+Ee[A_Y][loc]*Ee[A_Y][loc]+Ee[A_Z][loc]*Ee[A_Z][loc])<<endl;
					search=loc;
				}
			}
		}
		fp.close();
	}
	if(CON.get_EM_calc_type()==4)
	{
		ofstream fp("H.dat");
		double times=CON.get_B_times()*u0;//u0: 透磁率
		double Le=CON.get_distancebp();
		int plot_type=CON.get_plot_B_type();	//出力タイプ　1=ベクトル　2=スカラー
		for(int i=1;i<=nelm;i++)
		{
			if(ELEM[i].material!=MAGNET)
			{
				double X=0;
				double Y=0;
				double Z=0;
				for(int j=1;j<=4;j++)
				{
					X+=NODE[ELEM[i].node[j]].r[A_X]*0.25;
					Y+=NODE[ELEM[i].node[j]].r[A_Y]*0.25;
					Z+=NODE[ELEM[i].node[j]].r[A_Z]*0.25;
				}
				
				if(Y>-2*Le && Y<2*Le)
				{
					if(plot_type==1) fp<<X<<" "<<Z<<" "<<Ee[A_X][i]*times<<" "<<Ee[A_Z][i]*times<<endl;
					if(plot_type==2) fp<<X<<" "<<Z<<" "<<sqrt(Ee[A_X][i]*Ee[A_X][i]+Ee[A_Y][i]*Ee[A_Y][i]+Ee[A_Z][i]*Ee[A_Z][i])<<endl;
				}
			}
		}
		fp.close();
	}
    ////////*/
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

//静磁場計算関数
void calc_static_magnetic_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F,int KTJ)
{
	Post post;
	sides *SIDE=new sides[KTE];			//辺クラス
	int *branch_num=new int[node+1];	//各節点が隣接する節点数(各節点から延びる辺の数)
	int max=1000;						//節点に隣接する最大節点数・・・vectorにするべき
	int **nei2=new int* [node+1];		//各節点の隣接する節点番号格納(neiは節点-要素であるが、nei2は節点-節点)
	for(int i=1;i<=node;i++) nei2[i]=new int [max];
	int side_num=0;						//全辺数格納
		
//辺要素生成 (節点要素を使用する場合でも、電流密度を求めるときに辺要素がほしい)
//sideかedgeに用語を統一すること
	side_num=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,SIDE,branch_num,nei2,KTE);

	double *current[3];//各要素の電流密度[A/m^2]
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];

	double *Be[3];//要素内電界 or 要素内磁界
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//各要素の透磁率

//透磁率決定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2) RP[i]=CON.get_RP();
		else if(ELEM[i].material==ELASTIC) RP[i]=1;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==MAGNET) RP[i]=1;
		else if(ELEM[i].material==IRON) RP[i]=2500;
		else cout<<"error:材質に対し、透磁率が不定 材質="<<ELEM[i].material<<endl;
	}

//電流密度決定
	if(CON.get_J_input_way()==2)//電流密度を他のソフトから読み込む
	{
		import_J0_density(CON, node, nelm,NODE,ELEM,current);
		check_J0(CON, node, nelm,NODE,ELEM,current);
	}
	else if(CON.get_J_input_way()==1)
	{
		//自分で
	}

	else if(CON.get_J_input_way()==0) for(int i=0;i<=nelm;i++) for(int D=0;D<3;D++) current[D][i]=0.0;//強制電流が存在しないから初期化

	if(CON.get_FEM_elm_type()==1)//辺要素
	{
		double *A=new double [side_num+1];//ベクトルポテンシャル
		
		//ベクトルポテンシャルを計算する
		if(CON.get_EM_calc_type()==2) Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP);	

		//磁束密度を計算する
		Bflux3D_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,Be,RP);

		//磁束密度をプロット
		if(post.get_plot_B()==true){
			plot_magnetic_flux_density(CON, PART, NODE, ELEM, nelm, Be, t);
			C_Fluix_data_avs(node,nelm,NODE,ELEM,KTJ,Be,CON,t);
		}
		
		//力の計算
		if(CON.get_m_force()==1) NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number);
		if(CON.get_m_force()==2) kelvin_force3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);
		if(CON.get_m_force()==4) direct_divT3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

		//スムージングと力のプロット
		if(post.get_plot_NF()==true){
			smoothingF3D(CON,PART,fluid_number,F,t);
		}

		delete [] A;
	}

	delete [] RP;
	for(int D=0;D<3;D++)
	{
		delete [] current[D];
		delete [] Be[D];
	}
	for(int i=1;i<=node;i++) delete [] nei2[i];
	delete [] nei2;
	delete [] branch_num;
	delete [] SIDE;
}

//静磁場計算関数
void calc_variable_magnetic_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F)
{
	Post post;
	sides *SIDE=new sides[KTE];			//辺クラス
	int *branch_num=new int[node+1];	//各節点が隣接する節点数(各節点から延びる辺の数)
	int max=1000;						//節点に隣接する最大節点数・・・vectorにするべき
	int **nei2=new int* [node+1];		//各節点の隣接する節点番号格納(neiは節点-要素であるが、nei2は節点-節点)
	for(int i=1;i<=node;i++) nei2[i]=new int [max];
	int side_num=0;						//全辺数格納
		
//辺要素生成 (節点要素を使用する場合でも、電流密度を求めるときに辺要素がほしい)
//sideかedgeに用語を統一すること
	side_num=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,SIDE,branch_num,nei2,KTE);

	double *current[3];//各要素の電流密度[A/m^2]
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];

	double *Be[3];//要素内電界 or 要素内磁界
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//各要素の透磁率

//透磁率決定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2) RP[i]=CON.get_RP();
		else if(ELEM[i].material==ELASTIC) RP[i]=1;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==MAGNET) RP[i]=1;
		else if(ELEM[i].material==IRON) RP[i]=2500;
		else cout<<"error:材質に対し、透磁率が不定 要追加"<<endl;
	}

//電流密度決定
	if(CON.get_J_input_way()==2)//電流密度を他のソフトから読み込む
	{
		import_J0_density(CON, node, nelm,NODE,ELEM,current);
		check_J0(CON, node, nelm,NODE,ELEM,current);
	}
	else if(CON.get_J_input_way()==1)
	{
		//自分で
	}

	else if(CON.get_J_input_way()==0) for(int i=0;i<=nelm;i++) for(int D=0;D<3;D++) current[D][i]=0.0;//強制電流が存在しないから初期化

	if(CON.get_FEM_elm_type()==1)//辺要素
	{
		double *A=new double [side_num+1];//ベクトルポテンシャル
		
		//ベクトルポテンシャルを計算する
		if(CON.get_EM_calc_type()==2) Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP);

		//磁束密度を計算する
		Bflux3D_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,Be,RP);

		//磁束密度をプロット
		if(post.get_plot_B()==true){
			plot_magnetic_flux_density(CON, PART, NODE, ELEM, nelm, Be, t);
		}

		//力の計算
		if(CON.get_m_force()==1) NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number);
		if(CON.get_m_force()==2) kelvin_force3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);
		if(CON.get_m_force()==4) direct_divT3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

		//スムージングと力のプロット
		if(post.get_plot_NF()==true){
			smoothingF3D(CON,PART,fluid_number,F,t);
		}

		delete [] A;
	}

	delete [] RP;
	for(int D=0;D<3;D++)
	{
		delete [] current[D];
		delete [] Be[D];
	}
	for(int i=1;i<=node;i++) delete [] nei2[i];
	delete [] nei2;
	delete [] branch_num;
	delete [] SIDE;
}

//電磁力スムージング関数・・・影響半径内にある周辺粒子の電磁力の平均値を与えている
void smoothingF3D(mpsconfig &CON, vector<mpselastic> &PART,int fluid_number,double *F[3],int t)
{

	double le=CON.get_distancebp();
    double ep=8.854e-12;///真空の誘電率。
    double *newF[3];
    for(int D=0;D<3;D++) newF[D]=new double [PART.size()];

	if(CON.get_FEM_smn()>0)
	{
		for(int n=0;n<CON.get_FEM_smn();n++)//これを１以上にする必要ある？？？
		{
		    for(size_t i=0;i<PART.size();i++) 
		    {  
				if(PART[i].type!=WALL && PART[i].type==MAGELAST)
				{
					for(int D=0;D<3;D++) newF[D][i]=F[D][i];
					int num=1; //自分自身をカウントするから1
					for(int k=0;k<PART[i].N;k++)
					{       
						int j=PART[i].NEI[k];
						if(PART[j].type==FLUID || PART[j].type==MAGELAST)
						{
							num++;
							for(int D=0;D<3;D++) newF[D][i]+=F[D][j];
						}
					}
					for(int D=0;D<3;D++) newF[D][i]/=num;
				}			
			} 
			for(size_t i=0;i<PART.size();i++)
			{
				if(PART[i].type!=WALL && PART[i].type==MAGELAST)	for(int D=0;D<3;D++) F[D][i]=newF[D][i];
			}
		}
	}
	else if(CON.get_FEM_smn()<0)//表面のみでスムージング
	{
		int N=-1*CON.get_FEM_smn();
		for(int n=0;n<N;n++)
		{
			for(size_t i=0;i<PART.size();i++) 
			{  
				if(PART[i].type!=WALL)
				{
					for(int D=0;D<3;D++) newF[D][i]=F[D][i];
					if(PART[i].surface==ON)
					{
						int num=1; //自分自身をカウントするから1
						for(int k=0;k<PART[i].N;k++)
						{       
							int j=PART[i].NEI[k];
							if(PART[j].surface==ON && (PART[j].type==FLUID || PART[j].type==MAGELAST))
							{
								num++;
								for(int D=0;D<3;D++) newF[D][i]+=F[D][j];
							}
						}
						for(int D=0;D<3;D++) newF[D][i]/=num;
					}
				}
			} 
			for(size_t i=0;i<PART.size();i++){
				if(PART[i].type!=WALL){
				for(int D=0;D<3;D++) F[D][i]=newF[D][i];
				}
			}
		}
	}

    for(int D=0;D<3;D++) delete [] newF[D];
    /////////////////////////////////////*/
    
    //電磁力をプロット, GNUplot用

	ofstream fp("F.dat");
	double xmax=-100;//出力粒子の最大横座標
	double ymax=-100;//出力粒子の最大縦座標
	double mass=CON.get_particle_mass();
	double times=CON.get_times()/mass*le*le*CON.get_FEMtimes();
    for(size_t i=0;i<PART.size();i++)
    {
		if(PART[i].type!=WALL){
		if(PART[i].r[A_Y]>-le*0.5 && PART[i].r[A_Y]<+le*0.5)
		{
			fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<F[A_X][i]*times<<"\t"<<F[A_Z][i]*times<<endl;
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		}
		}
    }
	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	if(CON.get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力
	fp.close();

	int BOnum=0;
	for(size_t i=0;i<PART.size();i++) if(PART[i].surface==ON && PART[i].type!=WALL) BOnum++;

	//構造格子データファイル
	int timestep=CON.get_current_step();
	stringstream sstr;
	sstr<<"./Lorentz/Lorentz"<<timestep<<".fld";
	string NodalForce=sstr.str();

	ofstream fout2(NodalForce);
	if(fout2.fail()){
		system("mkdir Lorentz");
		ofstream fout2(NodalForce);
		if(fout2.fail()){
			cout<<"./Lorentzフォルダを開けませんでした"<<endl;
			exit(1);
		}
	}

	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./Lorentz"<<timestep<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	//データファイル
	stringstream ss;
	ss<<"./Lorentz/Lorentz"<<timestep;
	string filename=ss.str();

	ofstream fout(filename);
	if(fout.fail()){
		cout<<"./Lorentzフォルダを開けませんでした"<<endl;
		exit(1);
	}

	fout<<"e-x e-y e-z x y z"<<endl;
	for(size_t i=0;i<PART.size();i++)
    {
		if(PART[i].type!=WALL){
//		fout<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		fout<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		}
	}
	fout.close();
}

//辺要素生成関数
int make_edge_element(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,sides *SIDE,int *branch_num,int **nei2,int KTE)
{
	unsigned timeA=GetTickCount();		//計算開始時刻

	//辺要素生成
	int *check=new int [node+1];		//各節点がﾁｪｯｸしおわったかどうか
	int *temp_check=new int[node+1];	//一時的なﾁｪｯｸ配列
	int max=1000;						//節点に隣接する最大節点数
	int *ele_side_num=new int[nelm+1];	//各要素の第何辺までがもとまっているか
	int side_num=0;						//全辺数格納
		
	///初期化
	for(int i=1;i<=node;i++) 
	{
		check[i]=0;
		temp_check[i]=0;
		branch_num[i]=0;//初期化
	}
	for(int i=1;i<=nelm;i++) ele_side_num[i]=0;
	/////////////////////

	//辺番号と辺-節点情報生成
	for(int i=1;i<=node;i++)
	{
		for(int k=1;k<=jnb[i];k++)
		{
			int jelm=nei[i][k];//節点iに隣接する要素番号
			for(int j=1;j<=4;j++)
			{
				int p=ELEM[jelm].node[j];//pは何？？
				if(p!=i && temp_check[p]==0)
				{
					branch_num[i]=branch_num[i]+1;//節点iの隣接する節点数をプラス
					temp_check[p]=1;//もう見た
					nei2[i][branch_num[i]]=p;

					if(check[p]==0)
					{	
						side_num++;//辺要素数更新
						SIDE[side_num].node[1]=i;//このアルゴリズム下においては、常にi<pである.なぜならiより小さな番号の節点はすでにchek=1だから
						SIDE[side_num].node[2]=p;
						///節点ベースの境界条件を辺ベースに拡張
						if(NODE[i].boundary_condition==NODE[p].boundary_condition){
							SIDE[side_num].boundary_condition=NODE[i].boundary_condition;//両端が同じ境界条件なら、その辺もその境界条件に従う。仮に両方未知数でも問題ない
						}else{
							SIDE[side_num].boundary_condition=0;//それ以外は未知数
						}
					}
				}
			}
		}
		check[i]=1;
		for(int k=1;k<=branch_num[i];k++) temp_check[nei2[i][k]]=0;//初期化
//		cout<<"for in i= "<<i<<"/"<<node<<"node, NODE["<<i<<"].attr="<<NODE[i].material<<endl; //なぜかmaterialにFLUID=0が入っている・・・
	}

	//要素-辺情報生成
	cout<<"ELEM[].sides格納"<<endl;	
	for(int i=1;i<=side_num;i++)
	{
		int ia=SIDE[i].node[1];//辺iを構成する節点番号(ia<ib)
		int ib=SIDE[i].node[2];
		for(int k=1;k<=jnb[ia];k++)
		{
			int jelm=nei[ia][k];//節点iaの隣接する要素番号
			for(int j=1;j<=4;j++)
			{
				if(ELEM[jelm].node[j]==ib)
				{
					ele_side_num[jelm]=ele_side_num[jelm]+1;
					int a=ele_side_num[jelm];
					ELEM[jelm].sides[a]=i;
				}
			}
		}
	}////////////

	for(int i=1;i<=side_num;i++) SIDE[i].boundary_condition=0;		//初期化・・・この位置はダメでは？

	delete [] check;
	delete [] temp_check;
	delete [] ele_side_num;

	std::cout<<"辺数="<<side_num<<" 最大辺数="<<KTE<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	return side_num;
}

//動磁場計算関数
void calc_transitional_EM_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F,int KTJ)
{
	sides *SIDE=new sides[KTE];			//辺クラス
	int *branch_num=new int[node+1];	//各節点が隣接する節点数(各節点から延びる辺の数)
	int max=1000;						//節点に隣接する最大節点数
	int **nei2=new int* [node+1];		//各節点の隣接する節点番号格納(neiは節点-要素であるが、nei2は節点-節点)
		for(int i=1;i<=node;i++) nei2[i]=new int [max];
		int side_num=0;						//全辺数格納
	int *depth=new int [KTE+1];//各要素の深さ格納

	Post post;
	double *V=new double[node+1];
	for(int i=1;i<=node;i++) V[i]=0;				//初期化　Vは渦電流の電位格納に使う
	double *current[3];								//各要素の電流密度[A/m2]格納
	for(int D=0;D<3;D++) current[D]=new double [nelm+1];
	double *Be[3];					//要素内電界 or 要素内磁界
    for(int D=0;D<3;D++) Be[D]=new double [nelm+1];
	double *RP=new double [nelm+1];	//各要素の透磁率

	side_num=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,SIDE,branch_num,nei2,KTE);
	//要素深さ決定
	set_depth(CON,NODE,ELEM,node,nelm,depth,KTE);

	//透磁率決定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==MAGELAST || ELEM[i].material==MAGELAST2) RP[i]=CON.get_RP();	
		else if(ELEM[i].material==ELASTIC) RP[i]=1;
		else if(ELEM[i].material==AIR) RP[i]=1;
		else if(ELEM[i].material==COIL) RP[i]=1;
		else if(ELEM[i].material==IRON) RP[i]=2500;
		else cout<<"error:材質に対し、透磁率が不定	材質="<<ELEM[i].material<<endl;
	}

	if(CON.get_J_input_way()==2)					//電流密度を他のソフトから読み込む
	{
		import_J0_density(CON, node, nelm,NODE,ELEM,current);
		check_J0(CON, node, nelm,NODE,ELEM,current);
	}
	else if(CON.get_J_input_way()==1)
	{
		//自分で
	}
	//////////////////ｺｲﾙ電流値を決定する
				int nm=450;			//巻きすう
				double I0=CON.get_I0();			//電流の振幅
				double f=CON.get_Hz();				//交流の周波数[Hz]
			//	double II=I0*cos(2*3.1415*f*TIME)+900;		//電流波形
				double II=I0*(-cos(2*3.1415*f*TIME))+450;		//0から
			//	double II=I0*sin(2*3.1415*f*TIME);		//SIN波電流波形
			//	double II=I0+450;//直流
/*				/////////////////////////////
				double II;
				if(TIME<(1/100.0)/2.0){
					II=I0*(-cos(2*3.1415*100*TIME))+450;
				}
				else II=I0+450;
				if(II==0)II=0.00000001;
			//	cout<<"Current="<<II<<endl;
				/////////////////////////////*/
		///////////////////////////////////////////
				ofstream Hz("A-t.dat", ios::app);
				Hz<<TIME<<" "<<II<<endl;
				Hz.close();
	int check=1;
	if(CON.get_FEM_calc_type()==3)//動磁場
	{
		if(check==0){//接点要素
			double *A[3];									//節点におけるベクトルポテンシャル
			for(int D=0;D<3;D++) A[D]=new double [node+1];
			double *old_A[3];		
			for(int D=0;D<3;D++) old_A[D]=new double [node+1]; //1step前のベクトルポテンシャル
				
			/////old_Aに値を格納
			if(t==1 && CON.get_restart()==OFF)				//最初のステップはゼロで初期化
			{
				for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old_A[D][i]=0;	//初期化
			}
			else//それ以外はファイルから読み込み
			{
				ifstream old("old_A.dat");
				if(!old) cout<<"cannot open old_A.dat"<<endl;
				old.unsetf(ifstream::dec);
				old.setf(ifstream::skipws);

				for(int i=1;i<=node;i++) for(int D=0;D<3;D++) old>>old_A[D][i];
					
				old.close();
			}
		
			Avector3D_node_eddy2(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,depth,nei2,branch_num,old_A,II,dt,fluid_number,V,t);
			//磁束密度解決
			Bflux3D_node(CON,NODE,ELEM,node,nelm,A,Be);
			/////////////////////////////////////////////////////
			//磁束密度をプロット
			if(post.get_plot_B()==true){
				plot_magnetic_flux_density(CON, PART, NODE, ELEM, nelm, Be, t);
				C_Fluix_data_avs(node,nelm,NODE,ELEM,KTJ,Be,CON,t);
			
			}

			//力の計算
			if(CON.get_m_force()==1) NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number);
			if(CON.get_m_force()==2) kelvin_force3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);
			if(CON.get_m_force()==4) direct_divT3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

			//スムージングと力のプロット
			if(post.get_plot_NF()==true){
				smoothingF3D(CON,PART,fluid_number,F,t);
			}
			/////////////////////////////////////////////////////

			for(int D=0;D<3;D++)
			{
				delete [] A[D];
				delete [] old_A[D];
			}
		}
		else if(check==1){//辺要素

			double *A=new double [side_num+1];//ベクトルポテンシャル
		
			//ベクトルポテンシャルを計算する
			Avector3D2(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP,II,depth,t);//辺要素
//			non_linear_Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP, II,depth,Be, t);
			//磁束密度を計算する
			Bflux3D_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,Be,RP);

			//磁束密度をプロット
			if(post.get_plot_B()==true){
				plot_magnetic_flux_density(CON, PART, NODE, ELEM, nelm, Be, t);
				C_Fluix_data_avs(node,nelm,NODE,ELEM,KTJ,Be,CON,t);
			}
		
			//力の計算
			if(CON.get_m_force()==1) NODE_F3D(CON,NODE,ELEM,node,nelm,Be,jnb,nei,RP,PART,F,fluid_number);
			if(CON.get_m_force()==2) kelvin_force3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);
			if(CON.get_m_force()==4) direct_divT3D(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

			//スムージングと力のプロット
			if(post.get_plot_NF()==true){
				smoothingF3D(CON,PART,fluid_number,F,t);
			}

			delete [] A;
		
		}
	}

	delete [] V;
	delete [] RP;
	delete [] depth;
	delete [] SIDE;
	delete [] branch_num;
	for(int i=1;i<=node;i++)delete [] nei2[i];
	delete [] nei2;
	for(int D=0;D<3;D++)
	{
		delete [] current[D];
		delete [] Be[D];
	}
}

//電流密度読み込み関数
void import_J0_density(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	cout<<"current_density.txtより強制電流分布を読み込み--";
	int id;

	ifstream fp("current_density.dat");
	if(!fp) cout<<"cannot open current_density.dat"<<endl;
	fp.unsetf(ifstream::dec);
	fp.setf(ifstream::skipws);

	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==COIL)						//current_density.datには、コイル要素のみ出力されている仕様にしておくこと
		{												//その場合、コイル要素は静的要素でなければならない。動的なら他のソフトから読み込めない
			fp>>id;
			for(int D=0;D<3;D++) fp>>current[D][i];		//
		}
		else for(int D=0;D<3;D++) current[D][i]=0;
	}
	fp.close();
	cout<<"ok"<<endl;
}

//電流密度出力関数
void check_J0(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current)
{
	int coil_num=0;							//コイル要素数
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) coil_num++;

	ofstream fout2("J0.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << coil_num <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./J0 filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./J0 filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./J0 filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./J0 filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./J0 filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./J0 filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("J0");

	fout<<"e-x e-y e-z x y z"<<endl;
	double times=1e-12;
	for(int i=1;i<=nelm;i++)
    {
		if(ELEM[i].material==COIL)
		{
			double r[3]={0,0,0};
			for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
			fout<<current[A_X][i]*times<<" "<<current[A_Y][i]*times<<" "<<current[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
			
		}
	}
	fout.close();
}

//節点要素磁束密度計算関数
void Bflux3D_node(mpsconfig &CON, vector<point3D> &NODE, vector<element3D> &ELEM, int node,int nelm, double **A, double **B)
{
	cout<<"磁束密度計算----------";

	int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];

	double times=CON.get_B_times();
	double le=CON.get_distancebp();
	int plot_type=CON.get_plot_B_type();//1:ベクトル　2:スカラー
	//ofstream fp2("test.dat");
    for(int je=1;je<=nelm;je++)
    {   
		for(int D=0;D<3;D++) B[D][je]=0;
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
    
		double delta6=ELEM[je].volume;//体積の6倍

		delta6=1/delta6;

		double Xs=0;
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];Ys+=Y[j];Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		//係数作成
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////

		for(int i=1;i<=4;i++)
		{
			B[A_X][je]+=d[i]*A[A_Z][N[i]]-e[i]*A[A_Y][N[i]];
			B[A_Y][je]+=e[i]*A[A_X][N[i]]-c[i]*A[A_Z][N[i]];
			B[A_Z][je]+=c[i]*A[A_Y][N[i]]-d[i]*A[A_X][N[i]];
		}

		for(int D=0;D<3;D++) B[D][je]*=delta6;

		//if(Ys>-5*le && Ys<5*le) fp2<<Xs<<" "<<Zs<<" "<<B[A_X][je]*100<<" "<<B[A_Z][je]*100<<endl;

	}
	cout<<"ok"<<endl;
//	fp2.close();

	

	//磁束密度出力
	ofstream fp("Bflux.dat");
	
	double Xmin=CON.get_XL()+le; double Xmax=CON.get_XR()-le;//解析領域
	double Zmin=CON.get_ZD()+le; double Zmax=CON.get_ZU()-le;
	double Rmax=CON.get_RU()-le;
	
	double dx=5*le;
	int Nx=(int)((Xmax-Xmin)/dx);//各方向の分割数
	int Nr=(int)(2*Rmax/dx);
	int Nz=(int)((Zmax-Zmin)/dx);
	int search=nelm;//locate関数で最初に探索する要素番号
	
	if(CON.get_region_shape()==1) //円筒領域のときは変数を書き換えて処理
	{
		Nx=Nr;
		Xmin=-Rmax;
	}


	if(plot_type==1)		//ベクトル表示
	{
		for(int n=0;n<Nx;n++)
		{
			for(int m=0;m<Nz;m++)
			{
				double xp=dx*n+Xmin;//出力する点の座標
				double yp=0;
				double zp=dx*m+Zmin;
				int loc=locate3D(NODE,ELEM,search,xp,yp,zp);
				fp<<xp<<" "<<zp<<" "<<B[A_X][loc]*times<<" "<<B[A_Z][loc]*times<<endl;
				search=loc;
				//if(loc==0) cout<<"EE"<<endl;
			}
		}
	}
	else if(plot_type==2)	//スカラー表示
	{
		for(int n=0;n<Nx;n++)
		{
			for(int m=0;m<Nz;m++)
			{
				double xp=dx*n+Xmin;//出力する点の座標
				double yp=0;
				double zp=dx*m+Zmin;
				int loc=locate3D(NODE,ELEM,search,xp,yp,zp);
				double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc]);
				fp<<xp<<" "<<zp<<" "<<BB<<endl;
				search=loc;
			}
		}
	}
	fp.close();//*/

	//流体付近だけピックアップしたeforce
	int count=0;
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		if(r[A_Y]<=0 &&  r[A_Y]>-0.03 && -0.03<r[A_X] && r[A_X]<0.03 &&  0.12<r[A_Z] && r[A_Z]<0.18) 
		count++;
	}


	ofstream fout2("eforce.fld");
	fout2 << "# AVS field file" << endl;
	fout2 << "ndim=1" << endl;
	//fout2 << "dim1=" << fluid_number <<endl;
	fout2 << "dim1=" << count/*nelm*/ <<endl;
	fout2 << "nspace=3" << endl;
	fout2 << "veclen=3" << endl;
	fout2 << "data=float" << endl;
	fout2 << "field=irregular" << endl;
	fout2 << "label=e-x e-y e-z" << endl << endl;
	fout2 << "variable 1 file=./eforce filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout2 << "variable 2 file=./eforce filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout2 << "variable 3 file=./eforce filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout2 << "coord    1 file=./eforce filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout2 << "coord    2 file=./eforce filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout2 << "coord    3 file=./eforce filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout2.close();

	ofstream fout("eforce");
	fout<<"e-x e-y e-z x y z"<<endl;
	for(int i=1;i<=nelm;i++)
    {
		double r[3]={0,0,0};
		for(int D=0;D<3;D++) for(int j=1;j<=4;j++) r[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		if(r[A_Y]<=0 &&  r[A_Y]>-0.03 && -0.03<r[A_X] && r[A_X]<0.03 &&  0.12<r[A_Z] && r[A_Z]<0.18)
		{
			fout<<B[A_X][i]*times<<" "<<B[A_Y][i]*times<<" "<<B[A_Z][i]*times<<" "<<r[A_X]<<" "<<r[A_Y]<<" "<<r[A_Z]<<endl;
		}
	}
	fout.close();

}

//表面力スムージング関数
void smoothingFs3D(mpsconfig &CON, vector<mpselastic> &PART,int fluid_number,double *Fs)
{
	double le=CON.get_distancebp();
    double *newFs=new double [fluid_number];

	int N=abs(CON.get_FEM_smn());//反復回数
	for(int n=0;n<N;n++)
	{
		for(int i=0;i<fluid_number;i++) 
		{  
			newFs[i]=Fs[i];
			if(PART[i].surface==ON)
			{
				int num=1; //自分自身をカウントするから1
				for(int k=0;k<PART[i].N;k++)
				{       
					int j=PART[i].NEI[k];
					if(PART[j].surface==ON && (PART[j].type==FLUID || PART[j].type==ELASTIC))
					{
						num++;
						newFs[i]+=Fs[j];
					}
				}
				newFs[i]/=num;
			}
		} 
		for(int i=0;i<fluid_number;i++)  Fs[i]=newFs[i];
	}

    delete [] newFs;
}

//ベクトルポテンシャル計算関数(辺要素用)
void Avector3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP)
{
	cout<<"ベクトルポテンシャル計算開始---"<<endl;

	double u0=4*PI*0.0000001;	//真空の透磁率
    double v0=1/u0;				//磁気抵抗率
	double j0x,j0y,j0z;			//電流密度[A/m^3]
	unsigned timeA=GetTickCount();

	//磁石の着磁方向を決定
	double MA=CON.get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON.get_magnet_B()*magnet_direction[A_X];
	double My=CON.get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON.get_magnet_B()*magnet_direction[A_Z];

	//周囲壁に固定境界条件を設定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					//第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
					int p=ELEM[i].node[j];//境界三角形に属さない節点
					for(int k=1;k<=6;k++)
					{
						int iside=ELEM[i].sides[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//節点pを含まない辺は境界辺
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
			
		}
	}
	////////////////*/

    int NN=0;//ディリクレ型境界辺数
    int *dn=new int [side_num+1]; //各辺がディリクレ型境界の何番目か。ディリクレ型境界上でないならside_num+1を格納
    double *PHAT=new double [CON.get_max_DN()];//ディリクレ型境値
    
    //ディリクレ型境界条件入力・・・ここでなぜかNNが更新されていない
	//NNも計算する
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"ディリクレ数="<<NN<<endl;
	//*/
    
    int pn=side_num-NN;				//未知数
    int *ppn=new int [pn];			//行列のn番目は辺ppn[n]
    int *npp=new int [side_num+1];	//各辺が行列の何番目にあたるか。ディリクレ型の場合はpn+1を格納
    int num=0;

    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//未知数
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
   
	cout<<"未知数: "<<pn<<endl;

    //行列の最大幅計算
    int mat_w=0;
	int *nume=new int [side_num+1];				//SIDE[i]を含む要素数
	int *wid= new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].sides[j];
			nume[edge]=nume[edge]+1;
		}
	}						
	//nume[i]がもとまった

	for(int i=1;i<=side_num;i++)
	{	//考え方:ある辺周りの要素群を考える。まず、その内の一つを空間においた段階で、辺の数は6つ。以後、要素を加えるごとに3つ辺が増える。よって次式となる
		int width=6+3*(nume[i]-1);//場合によってはこれより小さい値になるかもしれない。けど多くメモリを確保するぶんには問題ない
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}
	delete [] nume;
	///////////////////////////////////////////

//	cout<<"配列確保の前まで成功, pn: "<<pn<<endl;
	//for(int i=0;i<=pn;i++) cout<<"wid["<<i<<"]="<<wid[i]<<endl;

    //配列確保
    double **G=new double*[pn+1];//全体行列
	for(int i=1;i<=pn;i++){
		G[i]=new double[wid[i]+1];
	}
//	cout<<"配列確保成功　G"<<endl;

    int **ROW=new int *[pn+1]; //各行の、非ゼロ要素の列番号記憶
	for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
//	cout<<"配列確保成功　ROW"<<endl;

    int *NUM=new int [pn+1]; //各行の、非ゼロ要素数
    
//	cout<<"配列確保成功　NUM"<<endl;

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        //for(int j=1;j<=mat_w;j++)
		for(int j=1;j<=wid[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    //

   
//全体行列を作成する
    
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	//cout<<"行列作成開始 ";
    for(int je=1;je<=nelm;je++)
    {   
		//辺－節点テーブル作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].sides[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//要素の重心座標
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		////////////////////////////

		//デローニ分割の際に求めた体積は、スーパーボックス内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
		
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;//i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////	

		double rp=RP[je];
	//	if(ELEM[je].material==FLUID) rp=CON.get_RP();//比透磁率
		//要素マトリクス作成開始
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].sides[i];//要素jeの辺番号
			if(SIDE[iside].boundary_condition==0)//ディリクレでない。つまり未知なら
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				int I1=SIDE[iside].node[1];//isideを構成する2点
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].sides[j];
					
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)//ディリクレでない。つまり未知なら
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jsideを構成する2点
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがディリクレ型境界節点なら
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}

				///B[I-1]を計算する（解行列）
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)	B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
			}
		}   	
    }
    ///////////////////////*/
	

    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    //行列の実際の最大幅を求める
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	
	//このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//非ゼロ要素の列番号格納配列
    int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
    
    /////////////////////val,ind ,ptrに値を格納
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;

	delete [] wid;
    
	cout<<"行列作成 幅："<<maxN<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    

//行列計算開始
	double *XX=new double [pn];//行列の答え格納
	if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
	else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		A[i]=XX[n];
		//cout<<A[i]<<endl;
	}
	delete [] XX;
	///////////////////////////*/
    
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}
//ベクトルポテンシャル計算関数(辺要素用)
void Avector3D2(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP,double II,int *depth,int t)
{
	cout<<"ベクトルポテンシャル計算開始---"<<endl;

	double u0=4*PI*0.0000001;	//真空の透磁率
    double v0=1/u0;				//磁気抵抗率
	double j0x,j0y,j0z;			//電流密度[A/m^3]
	unsigned timeA=GetTickCount();

	//磁石の着磁方向を決定
	double MA=CON.get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON.get_magnet_B()*magnet_direction[A_X];
	double My=CON.get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON.get_magnet_B()*magnet_direction[A_Z];
	
	//////////////////////////////////////////////////////*/
	//////////////////////ｺｲﾙがなければ電流計算を行う必要がない。そこでｺｲﾙ節点数をかぞえる
	int coil_node_num=0;	//ｺｲﾙ節点数
	for(int i=1;i<=node;i++) if(NODE[i].material==COIL) coil_node_num++;
	////////////////////////*/

	int *save_bound=new int [node+1];//本関数は何度も呼び出されるので、電流に関する境界条件を完全に消去できない。そこで保存する
	
	if(coil_node_num>0)///ｺｲﾙ節点があるなら電流密度計算
	{
		for(int i=1;i<=node;i++) save_bound[i]=NODE[i].boundary_condition;//境界条件を保存

		if(II!=0)//電流値が非ゼロなら電流計算
		{
			calc_current(CON,NODE,ELEM,SIDE,node,nelm,side_num,jnb,branch_num,current,depth,II);
			cout<<"電流計算完了"<<endl;
		}
		else	//ゼロなら初期化
		{
			for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) for(int D=0;D<3;D++) current[D][i]=0;
		}
		///電流境界条件の初期化（もう電流はもとまっているからいらない)
		for(int i=1;i<=side_num;i++) if(SIDE[i].boundary_condition>=10) SIDE[i].boundary_condition=0;
		for(int i=1;i<=node;i++) if(NODE[i].boundary_condition>=10) NODE[i].boundary_condition=0;
	}
	carrent_vector(CON, NODE, ELEM, nelm, current,t);
	//周囲壁に固定境界条件を設定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					//第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
					int p=ELEM[i].node[j];//境界三角形に属さない節点
					for(int k=1;k<=6;k++)
					{
						int iside=ELEM[i].sides[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//節点pを含まない辺は境界辺
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
			
		}
	}
	////////////////*/

    int NN=0;//ディリクレ型境界辺数
    int *dn=new int [side_num+1]; //各辺がディリクレ型境界の何番目か。ディリクレ型境界上でないならside_num+1を格納
    double *PHAT=new double [CON.get_max_DN()];//ディリクレ型境値
    
    //ディリクレ型境界条件入力・・・ここでなぜかNNが更新されていない
	//NNも計算する
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"ディリクレ数="<<NN<<endl;
	//*/
    
    int pn=side_num-NN;				//未知数
    int *ppn=new int [pn];			//行列のn番目は辺ppn[n]
    int *npp=new int [side_num+1];	//各辺が行列の何番目にあたるか。ディリクレ型の場合はpn+1を格納
    int num=0;

    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//未知数
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
   
	cout<<"未知数: "<<pn<<endl;

    //行列の最大幅計算
    int mat_w=0;
	int *nume=new int [side_num+1];				//SIDE[i]を含む要素数
	int *wid= new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].sides[j];
			nume[edge]=nume[edge]+1;
		}
	}						
	//nume[i]がもとまった

	for(int i=1;i<=side_num;i++)
	{	//考え方:ある辺周りの要素群を考える。まず、その内の一つを空間においた段階で、辺の数は4つ。以後、要素を加えるごとに3つ辺が増える。よって次式となる
		int width=6+3*(nume[i]-1);//場合によってはこれより小さい値になるかもしれない。けど多くメモリを確保するぶんには問題ない
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}
	delete [] nume;
	///////////////////////////////////////////

//	cout<<"配列確保の前まで成功, pn: "<<pn<<endl;
	//for(int i=0;i<=pn;i++) cout<<"wid["<<i<<"]="<<wid[i]<<endl;

    //配列確保
    double **G=new double*[pn+1];//全体行列
	for(int i=1;i<=pn;i++){
		G[i]=new double[wid[i]+1];
	}
//	cout<<"配列確保成功　G"<<endl;

    int **ROW=new int *[pn+1]; //各行の、非ゼロ要素の列番号記憶
	for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
//	cout<<"配列確保成功　ROW"<<endl;

    int *NUM=new int [pn+1]; //各行の、非ゼロ要素数
    
//	cout<<"配列確保成功　NUM"<<endl;

    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        //for(int j=1;j<=mat_w;j++)
		for(int j=1;j<=wid[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    //

   
//全体行列を作成する
    
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	//cout<<"行列作成開始 ";
    for(int je=1;je<=nelm;je++)
    {   
		//辺－節点テーブル作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].sides[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//要素の重心座標
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		////////////////////////////

		//デローニ分割の際に求めた体積は、スーパーボックス内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
		
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;//i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////	

		double rp=RP[je];
	//	if(ELEM[je].material==FLUID) rp=CON.get_RP();//比透磁率
		//要素マトリクス作成開始
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].sides[i];//要素jeの辺番号
			if(SIDE[iside].boundary_condition==0)//ディリクレでない。つまり未知なら
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				int I1=SIDE[iside].node[1];//isideを構成する2点
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].sides[j];
					
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)//ディリクレでない。つまり未知なら
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jsideを構成する2点
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがディリクレ型境界節点なら
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}

				///B[I-1]を計算する（解行列）
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)
				{
					B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				}
			}
		}   	
    }
    ///////////////////////*/
	

    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    //行列の実際の最大幅を求める
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	
	//このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//非ゼロ要素の列番号格納配列
    int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
    
    /////////////////////val,ind ,ptrに値を格納
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;

	delete [] wid;
    
	cout<<"行列作成 幅："<<maxN<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    

//行列計算開始
	double *XX=new double [pn];//行列の答え格納
	if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
	else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		A[i]=XX[n];
		//cout<<A[i]<<endl;
	}
	delete [] XX;
	///////////////////////////*/
    
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}
void Avector3D_OK(mpsconfig &CON,vector <point3D> &NODE,vector <element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,int *depth)
{
	cout<<"ベクトルポテンシャル計算開始---"<<endl;

	double u0=4*PI*0.0000001;	//真空の透磁率
    double v0=1/u0;				//磁気抵抗率
	double j0x,j0y,j0z;			//電流密度[A/m^3]
	unsigned timeA=GetTickCount();

	//磁石の着磁方向を決定
	double MA=CON.get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON.get_magnet_B()*magnet_direction[A_X];
	double My=CON.get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON.get_magnet_B()*magnet_direction[A_Z];
	

	//////////////////////////////////////////////////////*/

	///周囲壁に固定境界条件を設定
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==AIR)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm==0)
				{
					///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
					int p=ELEM[i].node[j];//境界三角形に属さない節点
					for(int k=1;k<=6;k++)
					{
						int iside=ELEM[i].sides[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//節点pを含まない辺は境界辺
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
			
		}
	}
	////////////////*/

    int NN=0;//ディリクレ型境界辺数
    int *dn=new int [side_num+1]; //各辺がディリクレ型境界の何番目か。ディリクレ型境界上でないならside_num+1を格納
    double *PHAT=new double [CON.get_max_DN()];//ディリクレ型境値
    
    ///ディリクレ型境界条件入力
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"ディリクレ数="<<NN;
	/////////////*/
    
	    
    int pn=side_num-NN;				///未知数
    int *ppn=new int [pn];			//行列のn番目は辺ppn[n]
    int *npp=new int [side_num+1];	///各辺が行列の何番目にあたるか。ディリクレ型の場合はpn+1を格納
    int num=0; 
    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//未知数
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
   
    ////行列の最大幅計算
    int mat_w=0;
	int *nume=new int [side_num+1];				//SIDE[i]を含む要素数
	int *wid= new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].sides[j];
			nume[edge]=nume[edge]+1;
		}
	}											///nume[i]がもとまった
	for(int i=1;i<=side_num;i++)
	{	//考え方:ある辺周りの要素群を考える。まず、その内の一つを空間においた段階で、辺の数は4つ。以後、要素を加えるごとに3つ辺が増える。よって次式となる
		int width=6+3*(nume[i]-1);//場合によってはこれより小さい値になるかもしれない。けど多くメモリを確保するぶんには問題ない
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}
	delete [] nume;
	///////////////////////////////////////////

	
    ////配列確保
    double **G=new double *[pn+1];///全体行列
	for(int i=1;i<=pn;i++) G[i]=new double [wid[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
	for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        //for(int j=1;j<=mat_w;j++)
		for(int j=1;j<=wid[i];j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    ////
    
    /////////全体行列を作成する
    
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	//cout<<"行列作成開始 ";
    for(int je=1;je<=nelm;je++)
    {   
		//if(ELEM[je].material!=7) cout<<ELEM[je].material<<endl;
		//辺－節点テーブル作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].sides[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//要素の重心座標
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j];
			Ys+=Y[j];
			Zs+=Z[j];
		}
		Xs/=4;Ys/=4;Zs/=4;
		////////////////////////////

		///デローニ分割の際に求めた体積は、スーパーボックス内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
		double delta6=ELEM[je].volume;//体積の6倍
		
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
				b[i]*=-1;
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////
	
		
		double rp=CON.get_RP();
	//	if(ELEM[je].material==FLUID) rp=CON.get_RP();//比透磁率
		////要素マトリクス作成開始
		for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].sides[i];//要素jeの辺番号
			if(SIDE[iside].boundary_condition==0)///ディリクレでない。つまり未知なら
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				int I1=SIDE[iside].node[1];//isideを構成する2点
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].sides[j];
					
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)///ディリクレでない。つまり未知なら
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jsideを構成する2点
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						    G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがディリクレ型境界節点なら
					{
					    int n=dn[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////*/
				}
				double* current[3]={0,0,0};
				///B[I-1]を計算する
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)
				{
					//cout<<Mz<<endl;
					B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				}
			}
		}   	
    }
    ///////////////////////*/
	

    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
    
    ///行列の実際の最大幅を求める
	int maxN=0;
	for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
	

	///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
	arrange_matrix(pn,NUM,ROW,G);

	//対称性チェック
	check_matrix_symmetry(pn,NUM,ROW,G);

    double *val = new double [number];
    int *ind = new int [number];//非ゼロ要素の列番号格納配列
    int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
    
    /////////////////////val,ind ,ptrに値を格納
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;

	delete [] wid;
    
	cout<<"行列作成 幅："<<maxN<<"/"<<mat_w<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
    

	///////////////////////行列計算開始
	double *XX=new double [pn];//行列の答え格納
	if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
	else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	for(int n=0;n<pn;n++)
	{
		int i=ppn[n];
		A[i]=XX[n];
		//cout<<A[i]<<endl;
	}
	delete [] XX;
	///////////////////////////*/
    
    
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}

//境界条件適用関数(３Ｄ辺要素用)
void set_boundary_condition3D_edge(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,int *dn,int *NN,double *PHAT,double *A)
{
	int N=0;	//数え上げ変数

	//ディリクレ型境界条件入力
	if(CON.get_uniform_B_sw()==OFF)
	{
		for(int i=1;i<=side_num;i++)
		{
		    if(SIDE[i].boundary_condition==1)
			{
		        dn[i]=N;
		        PHAT[N]=0;
		        A[i]=0;
		        N++;
			}
			else if(SIDE[i].boundary_condition==2)
			{   
				/*int ia=SIDE[i].node[1];
				int ib=SIDE[i].node[2];
				A[i]=(NODE[ib].r[A_Z]-NODE[ia].r[A_Z])*0;
		        dn[i]=NN;
		        PHAT[NN]=(NODE[ib].r[A_Z]-NODE[ia].r[A_Z])*0;
		        //A[i]=0;
		        NN++;*/
				dn[i]=N;
		        PHAT[N]=0.0;
		        A[i]=0.0;
		        N++;
			}
			else
			{
				dn[i]=side_num+1;
			}
		}
	}
	if(CON.get_uniform_B_sw()==ON)	//解析領域全体に一様磁場を与える場合
	{
		double B=CON.get_uniform_B();//一様磁場の大きさ[Ｔ］
		double R[3];					//辺の長さ格納(X,Y,Z方向)
		double r[3];					//辺の中点座標格納
		double err=1e-14;
		
		for(int i=1;i<=side_num;i++)
		{
		    if(SIDE[i].boundary_condition==1)
			{   
				int ia=SIDE[i].node[1];
				int ib=SIDE[i].node[2];
				for(int D=0;D<3;D++)
				{
					R[D]=NODE[ib].r[D]-NODE[ia].r[D];
					r[D]=(NODE[ib].r[D]+NODE[ia].r[D])*0.5;//中点
				}

			//	if(r[A_Z]<CON.get_ZU()-err && r[A_Z]>CON.get_ZD()+err)
				{

					double L=sqrt(R[A_X]*R[A_X]+R[A_Y]*R[A_Y]+R[A_Z]*R[A_Z]);//辺の長さ
				
					A[i]=0.5*(-B*r[A_Y]*R[A_X]+B*r[A_X]*R[A_Y]);	//なぜこうなるのかはストークスの定理を利用すればわかる。
					dn[i]=N;
					PHAT[N]=0.5*(-B*r[A_Y]*R[A_X]+B*r[A_X]*R[A_Y]);
					N++;
				}
				//else 
				//{
				////	SIDE[i].boundary_condition=0;//自然境界条件
				//	dn[i]=side_num+1;
				//}
			}
			else dn[i]=side_num+1;
		}
	}

	*NN=N;
}

//辺要素磁束密度計算関数
void Bflux3D_side(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,double **B,double *RP)
{
	cout<<"磁束密度計算開始----";
	unsigned timeA=GetTickCount();
	int plot_type=CON.get_plot_B_type();	//ファイル出力形式　1=ベクトル 2=スカラー
	double times=CON.get_B_times();		//ファイル出力時の倍率
	double u0=4*PI*1e-7;					//真空の透磁率

	double *Xg=new double [nelm+1];			//要素の重心座標
	double *Yg=new double [nelm+1];
	double *Zg=new double [nelm+1];


	//#pragma omp parallel for
	for(int je=1;je<=nelm;je++)
    {   
		//cout<<je<<" "<<omp_get_thread_num()<<endl;//各iの計算を担当しているスレッド番号出力
		int N[4+1];								//要素の各節点番号格納
		double X[4+1];
		double Y[4+1];
		double Z[4+1];
		double c[4+1];
		double d[4+1];
		double e[4+1];
		int table[6+1][2+1];					//要素を構成する辺とその要素内節点番号関係

		//辺－節点テーブル作成
		for(int i=1;i<=6;i++)
		{
			int iside=ELEM[je].sides[i];
			int ia=SIDE[iside].node[1];
			int ib=SIDE[iside].node[2];
			for(int j=1;j<=4;j++)
			{
				if(ELEM[je].node[j]==ia) table[i][1]=j;
				else if(ELEM[je].node[j]==ib) table[i][2]=j;
			}
		}////////////////
	
		for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
		double Xs=0;//要素の重心座標
		double Ys=0;
		double Zs=0;
		for(int j=1;j<=4;j++)
		{
			X[j]=NODE[N[j]].r[A_X];
			Y[j]=NODE[N[j]].r[A_Y];
			Z[j]=NODE[N[j]].r[A_Z];
			Xs+=X[j]*0.25;
			Ys+=Y[j]*0.25;
			Zs+=Z[j]*0.25;
		}
		Xg[je]=Xs; Yg[je]=Ys; Zg[je]=Zs;	//重心代入
		////////////////////////////

		double delta6=ELEM[je].volume;//体積の6倍(正しい体積の値はすでにベクトルポテンシャルを求める際に計算している)
	
		delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
		double delta=ELEM[je].volume/6;//本当の体積
	
		for(int i=1;i<=4;i++)
		{
			int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
			int m=j%4+1;
			int n=m%4+1;
	    
			c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
			d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
			e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
			if(i%2!=0)//iが奇数なら
			{
			    c[i]*=-1;
				d[i]*=-1;
				e[i]*=-1;
			}
		}
		/////////

		for(int D=0;D<3;D++) B[D][je]=0;//初期化

		for(int i=1;i<=6;i++)
		{
			int s=ELEM[je].sides[i];
			int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
			int k2=table[i][2];

			B[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*A[s]; //c=x, d=y, e=zに対応, B=rotAを計算している
			B[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*A[s];
			B[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*A[s];
			//cout<<A[s]<<endl;
		}

		double BB=0;
		for(int D=0;D<3;D++)
		{
			B[D][je]*=delta6*delta6*2; //式(3.9)
		//	BB+=B[D][je]*B[D][je];
			//cout<<B[D][je]<<endl;
		}
		//cout<<sqrt(BB)<<"[T]"<<endl;
	}

	//磁束密度出力
	string filename[4];
	filename[0]="vectorBflux";
	filename[1]="vectorField";
	filename[2]="scalarBflux";
	filename[3]="scalarField";

	for(int suffix=0;suffix<4;suffix++)
	{
		stringstream ss;

		if(suffix%2==0){
			ss<<"./Bflux/"<<filename[suffix]<<CON.get_current_step()<<".dat";
		}else{
			ss<<"./Field/"<<filename[suffix]<<CON.get_current_step()<<".dat";
		}
		filename[suffix]=ss.str();
	}

	ofstream vectorflux(filename[0]);
	ofstream vectorfield(filename[1]);
	ofstream scalarflux(filename[2]);
	ofstream scalarfield(filename[3]);

	double le=CON.get_distancebp();
	double Xmin=CON.get_XL()+le; double Xmax=CON.get_XR()-le; //解析領域(境界ちょうどだと、その点を含む要素が打ち切り誤差で見つからないときがあるのでleだけ保険かける)
	double Zmin=CON.get_ZD()+le; double Zmax=CON.get_ZU()-le;

	if(CON.get_region_shape()==1)//円筒形の解析領域なら
	{
		Xmin=-1*CON.get_RU()+le; 
		Xmax=CON.get_RU()-le;
	}
	
	double dx=3*le;
	int Nx=(int)((Xmax-Xmin)/dx);								//各方向の分割数
	int Nz=(int)((Zmax-Zmin)/dx);
	int search=nelm;												//locate関数で最初に探索する要素番号
	
//	if(plot_type==1)	//ベクトル表示
	{
		for(int n=0;n<Nx;n++)
		{
			for(int m=0;m<Nz;m++)
			{
				double xp=dx*n+Xmin;							//出力する点の座標
				double yp=0;
				double zp=dx*m+Zmin;
				int loc=locate3D(NODE,ELEM,search,xp,yp,zp);
				vectorflux<<xp<<"\t"<<zp<<"\t"<<B[A_X][loc]*times<<"\t"<<B[A_Z][loc]*times<<endl;
				vectorfield<<xp<<"\t"<<zp<<"\t"<<B[A_X][loc]*times/(u0*RP[loc])<<"\t"<<B[A_X][loc]*times/(u0*RP[loc])<<endl;
				search=loc;
			}
		}
	}
//	else if(plot_type==2)										//スカラー表示
	{
		for(int n=0;n<Nx;n++)
		{
			for(int m=0;m<Nz;m++)
			{
				double xp=dx*n+Xmin;							//出力する点の座標
				double yp=0;
				double zp=dx*m+Zmin;
				int loc=locate3D(NODE,ELEM,search,xp,yp,zp);
				double BB=sqrt(B[A_X][loc]*B[A_X][loc]+B[A_Y][loc]*B[A_Y][loc]+B[A_Z][loc]*B[A_Z][loc]);
				scalarflux<<xp<<"\t"<<zp<<"\t"<<BB<<endl;
				scalarfield<<xp<<"\t"<<zp<<"\t"<<BB/(u0*RP[loc])<<endl;
				search=loc;
			}
		}
	}

	vectorflux.close();
	scalarflux.close();
	vectorfield.close();
	scalarfield.close();

	delete [] Xg;
	delete [] Yg;
	delete [] Zg;
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
		
}

//マクスウェルテンソル体積力計算関数
void direct_divT3D(mpsconfig &CON, vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int *jnb,int **nei)
{
    cout<<"div(T)の直接評価による電磁力計算--";
	unsigned int timeA=GetTickCount();
    double le=CON.get_distancebp();
	double R=CON.get_re()*le;
    double V=CON.get_particle_volume();//粒子の体積
    double u0=1.257e-6;				//真空の透磁率
	int N=4;
	int order=1;
	double Fz=0;
    for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) F[D][i]=0;//初期化
	double *matrix=new double [N*N];	//N×Nの係数行列
	double *matrix_val=new double [N*N];	//matrixの保存用
	double *Bxx=new double [N];	//Nの解行列
	double *Bxy=new double [N];	//Nの解行列
	double *Bxz=new double [N];	//Nの解行列
	double *Byx=new double [N];	//Nの解行列
	double *Byy=new double [N];	//Nの解行列
	double *Byz=new double [N];	//Nの解行列
	double *Bzx=new double [N];	//Nの解行列
	double *Bzy=new double [N];	//Nの解行列
	double *Bzz=new double [N];	//Nの解行列

	

	if(order==1)//3次元1次式
	{
		///係数行列は
		///   ΣΔx2    ΣΔxΔy  ΣΔxΔz ΣΔx  a = ΣΔxΔf  
		///  ΣΔxΔy   ΣΔy2    ΣΔyΔz ΣΔy  b = ΣΔyΔf 
		///  ΣΔxΔz   ΣΔyΔz  ΣΔz2   ΣΔz  c = ΣΔzΔf 
		///  ΣΔx      ΣΔy     ΣΔz    Σ1    d = ΣΔf

		double Ri[3];//節点iの座標格納
		double Rg[3];	//要素の重心格納
		//for(int i=1;i<=NN;i++)//流体節点
		for(int i=1;i<=node;i++)//流体節点
		{
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC || NODE[i].material==MAGELAST)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;		//初期化
				for(int n=0;n<N;n++) 
				{
					Bxx[n]=0;  Bxy[n]=0;  Bxz[n]=0;
					Byx[n]=0;  Byy[n]=0;  Byz[n]=0;
					Bzx[n]=0;  Bzy[n]=0;  Bzz[n]=0;
				}
				for(int D=0;D<3;D++) Ri[D]=NODE[i].r[D];
				//if(jnb[i]<3) cout<<i<<" "<<jnb[i]<<endl;
				if(jnb[i]>3)
				{
					for(int j=1;j<=jnb[i];j++)//節点iはjnb[n]個の要素に隣接している
					{
						int jelm=nei[i][j];
						double u=RP[jelm]*u0;
						///マクスウェルの応力テンソル
						double Txx=(Be[A_X][jelm]*Be[A_X][jelm]-Be[A_Y][jelm]*Be[A_Y][jelm]-Be[A_Z][jelm]*Be[A_Z][jelm])/(2*u);
						double Txy=2*Be[A_X][jelm]*Be[A_Y][jelm]/(2*u);
						double Txz=2*Be[A_X][jelm]*Be[A_Z][jelm]/(2*u);
						double Tyx=Txy;
						double Tyy=(-Be[A_X][jelm]*Be[A_X][jelm]+Be[A_Y][jelm]*Be[A_Y][jelm]-Be[A_Z][jelm]*Be[A_Z][jelm])/(2*u);
						double Tyz=2*Be[A_Y][jelm]*Be[A_Z][jelm]/(2*u);
						double Tzx=Txz;
						double Tzy=Tyz;
						double Tzz=(-Be[A_X][jelm]*Be[A_X][jelm]-Be[A_Y][jelm]*Be[A_Y][jelm]+Be[A_Z][jelm]*Be[A_Z][jelm])/(2*u);
						//////////
						for(int D=0;D<3;D++) Rg[D]=0;
						for(int k=1;k<=4;k++)
						{
							int node=ELEM[jelm].node[k];//jelmを構成する節点番号
							for(int D=0;D<3;D++) Rg[D]+=NODE[node].r[D]*0.25;
						}//Rg[D]がもとまった。
	
						double X=Rg[A_X]-Ri[A_X];
						double Y=Rg[A_Y]-Ri[A_Y];
						double Z=Rg[A_Z]-Ri[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						double w=1;
						if(dis>le) w=le*le/(dis*dis);
		
						matrix[0]+=X*X*w;		
						matrix[1]+=X*Y*w;	
						matrix[2]+=X*Z*w;	
						matrix[3]+=X*w;
	
						matrix[5]+=Y*Y*w;		
						matrix[6]+=Y*Z*w;
						matrix[7]+=Y*w;
	
						matrix[10]+=Z*Z*w;
						matrix[11]+=Z*w;
	
						matrix[15]+=w;
	
						Bxx[0]+=X*w*Txx;  Bxy[0]+=X*w*Txy;  Bxz[0]+=X*w*Txz;
						Byx[0]+=X*w*Tyx;  Byy[0]+=X*w*Tyy;  Byz[0]+=X*w*Tyz;
						Bzx[0]+=X*w*Tzx;  Bzy[0]+=X*w*Tzy;  Bzz[0]+=X*w*Tzz;
	
						Bxx[1]+=Y*w*Txx;  Bxy[1]+=Y*w*Txy;  Bxz[1]+=Y*w*Txz;
						Byx[1]+=Y*w*Tyx;  Byy[1]+=Y*w*Tyy;  Byz[1]+=Y*w*Tyz;
						Bzx[1]+=Y*w*Tzx;  Bzy[1]+=Y*w*Tzy;  Bzz[1]+=Y*w*Tzz;
	
						Bxx[2]+=Z*w*Txx;  Bxy[2]+=Z*w*Txy;  Bxz[2]+=Z*w*Txz;
						Byx[2]+=Z*w*Tyx;  Byy[2]+=Z*w*Tyy;  Byz[2]+=Z*w*Tyz;
						Bzx[2]+=Z*w*Tzx;  Bzy[2]+=Z*w*Tzy;  Bzz[2]+=Z*w*Tzz;
	
						Bxx[3]+=w*Txx;  Bxy[3]+=w*Txy;  Bxz[3]+=w*Txz;
						Byx[3]+=w*Tyx;  Byy[3]+=w*Tyy;  Byz[3]+=w*Tyz;
						Bzx[3]+=w*Tzx;  Bzy[3]+=w*Tzy;  Bzz[3]+=w*Tzz;
					
					}	
					matrix[4]=matrix[1];
					matrix[8]=matrix[2];
					matrix[9]=matrix[6];
					matrix[12]=matrix[3];
					matrix[13]=matrix[7];
					matrix[14]=matrix[11];

					for(int n=0;n<N*N;n++) matrix_val[n]=matrix[n];//値を保存
	
					//ガウスの消去法で解く
					gauss(matrix,Bxx,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
					gauss(matrix,Bxy,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
					gauss(matrix,Bxz,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
	
					gauss(matrix,Byx,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
					gauss(matrix,Byy,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
					gauss(matrix,Byz,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
	
					gauss(matrix,Bzx,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
					gauss(matrix,Bzy,N);
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
					gauss(matrix,Bzz,N);
	
					double fx=Bxx[0]+Bxy[1]+Bxz[2];//体積力fx:(divT)のX成分
					double fy=Byx[0]+Byy[1]+Byz[2];//体積力fy:(divT)のX成分
					double fz=Bzx[0]+Bzy[1]+Bzz[2];//体積力fz:(divT)のX成分
	
					int p=NODE[i].particleID;//対応する粒子番号
					if(p>=0)
					{
						F[A_X][p]=fx*V;
						F[A_Y][p]=fy*V;//単位は[N]
						F[A_Z][p]=fz*V;
		
						
					}
					Fz+=fz*V;
				}
			}

		}
	}
	int *check=new int[fluid_number];//check[i]=ONならその粒子は対応する節点が存在するということ
	for(int i=0;i<fluid_number;i++) check[i]=OFF;
	for(int i=1;i<=node;i++)
	{
		if(NODE[i].material==FLUID || NODE[i].material==ELASTIC || NODE[i].material==MAGELAST)//流体節点
		{
			int p=NODE[i].particleID;//粒子番号
			if(p>=0) check[p]=ON;
		}
	}
	for(int i=0;i<fluid_number;i++)
	{
		if(check[i]==OFF)//対応する節点がないので、力がゼロのままな粒子
		{
			double W=0;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					if(check[j]==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);//粒子間の距離
						double w=R/dis-1;//重み関数
						W+=w;
						for(int D=0;D<3;D++) F[D][i]+=F[D][j]*w;
					}
				}
			}
			if(W!=0) for(int D=0;D<3;D++) F[D][i]/=W;
		}
	}

	delete [] matrix;
	delete [] matrix_val;
	delete [] Bxx;
	delete [] Bxy;
	delete [] Bxz;
	delete [] Byx;
	delete [] Byy;
	delete [] Byz;
	delete [] Bzx;
	delete [] Bzy;
	delete [] Bzz;

	delete [] check;

	cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
} 

//ケルビン力計算関数
void kelvin_force3D(mpsconfig &CON, vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int*jnb,int **nei)
{
	cout<<"kelvin力による電磁力計算--";
	//kelvin力は∫u0(M・∇)H dv
	unsigned int timeA=GetTickCount();
    double le=CON.get_distancebp();
	double R=CON.get_re()*le;
    double V=CON.get_particle_volume();//粒子の体積
    double u0=1.257e-6;				//真空の透磁率
	double kai=CON.get_RP()-1;
	double Ms=14700;					//飽和磁化[A/m]
	double kai0=1.172;					//磁気感受率
	double gamma=3*kai0/Ms;

	int N=4;						//未知数
	int order=1;					//近似精度
	double Fz=0;
	double Fx=0;
    for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) F[D][i]=0;//初期化

	double *matrix=new double [N*N];	//N×Nの係数行列
	double *matrix_val=new double [N*N];	//matrixの保存用
	double *Bx=new double [N];	//Nの解行列
	double *By=new double [N];	//Nの解行列		//たとえばByはHyの微分を求めるのに使う
	double *Bz=new double [N];	//Nの解行列
	double *Bh=new double [N];	//Nの解行列		//Hの大きさの勾配を求めるのに使う

	double H[3];

	double *direct[DIMENSION];
    for(int D=0;D<DIMENSION;D++) direct[D]=new double [fluid_number];
	double *fs=new double [fluid_number];//表面節点の表面力格納
	int *num=new int [fluid_number];//表面節点の表面力格納
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)  direct_f(CON,PART,i,direct);
		else  for(int D=0;D<DIMENSION;D++) direct[D][i]=0;
		fs[i]=0;
		num[i]=0;
	}

	
	//Hの勾配をWLSMで求める
	if(order==1)//3次元1次式
	{
		///係数行列は
		///   ΣΔx2    ΣΔxΔy  ΣΔxΔz ΣΔx  a = ΣΔxΔf  
		///  ΣΔxΔy   ΣΔy2    ΣΔyΔz ΣΔy  b = ΣΔyΔf 
		///  ΣΔxΔz   ΣΔyΔz  ΣΔz2   ΣΔz  c = ΣΔzΔf 
		///  ΣΔx      ΣΔy     ΣΔz    Σ1    d = ΣΔf

		double Ri[3];//節点iの座標格納
		double Rg[3];	//要素の重心格納
		for(int i=1;i<=node;i++)
		{
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;		//初期化
				for(int n=0;n<N;n++)  {Bx[n]=0; By[n]=0; Bz[n]=0; Bh[n]=0;}
				
				for(int D=0;D<3;D++) Ri[D]=NODE[i].r[D];
				if(jnb[i]>5)
				{
					double Hi[3]={0,0,0};//節点iにおけるH格納
					double HHi=0;		////節点iにおけるHの大きさ格納
					for(int j=1;j<=jnb[i];j++)//節点iはjnb[n]個の要素に隣接している
					{
						int jelm=nei[i][j];
						double u=RP[jelm]*u0;
						double HHj=0;			//jelmにおける磁場の大きさ
						for(int D=0;D<3;D++) 
						{
							H[D]=Be[D][jelm]/u;//jelmにおけるトータルな磁場　本当は外部磁場じゃないといけないらしいけど？
							Hi[D]+=H[D];
							HHj+=H[D]*H[D];
						}
						HHj=sqrt(HHj);
						HHi+=HHj;
						
						for(int D=0;D<3;D++) Rg[D]=0;
						for(int k=1;k<=4;k++)
						{
							int node=ELEM[jelm].node[k];//jelmを構成する節点番号
							for(int D=0;D<3;D++) Rg[D]+=NODE[node].r[D]*0.25;
						}//jelmの重心Rg[D]がもとまった。
	
						double X=Rg[A_X]-Ri[A_X];
						double Y=Rg[A_Y]-Ri[A_Y];
						double Z=Rg[A_Z]-Ri[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);
					
						double w=1;
						if(dis>le) w=le*le/(dis*dis);
		
						matrix[0]+=X*X*w;		
						matrix[1]+=X*Y*w;	
						matrix[2]+=X*Z*w;	
						matrix[3]+=X*w;
	
						matrix[5]+=Y*Y*w;		
						matrix[6]+=Y*Z*w;
						matrix[7]+=Y*w;
	
						matrix[10]+=Z*Z*w;
						matrix[11]+=Z*w;
	
						matrix[15]+=w;
	
						Bx[0]+=X*w*H[A_X];
						By[0]+=X*w*H[A_Y];
						Bz[0]+=X*w*H[A_Z]; 
						Bh[0]+=X*w*HHj;

						Bx[1]+=Y*w*H[A_X];
						By[1]+=Y*w*H[A_Y]; 
						Bz[1]+=Y*w*H[A_Z];
						Bh[1]+=Y*w*HHj;
	
						Bx[2]+=Z*w*H[A_X]; 
						By[2]+=Z*w*H[A_Y]; 
						Bz[2]+=Z*w*H[A_Z]; 
						Bh[2]+=Z*w*HHj;

						Bx[3]+=w*H[A_X];
						By[3]+=w*H[A_Y]; 
						Bz[3]+=w*H[A_Z]; 
						Bh[3]+=w*HHj;	
					}
					matrix[4]=matrix[1];
					matrix[8]=matrix[2];
					matrix[9]=matrix[6];
					matrix[12]=matrix[3];
					matrix[13]=matrix[7];
					matrix[14]=matrix[11];

					for(int D=0;D<3;D++) Hi[D]/=jnb[i];//節点iにおけるH
					HHi/=jnb[i];
				
					matrix[15]+=1;//自身
					Bx[3]+=Hi[A_X];
					By[3]+=Hi[A_Y];
					Bz[3]+=Hi[A_Z];
					Bh[3]+=HHi;

					for(int n=0;n<N*N;n++) matrix_val[n]=matrix[n];//値を保存

					//ガウスの消去法で解く
					gauss(matrix,Bx,N);
				
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
				
					gauss(matrix,By,N);
				
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
				
					gauss(matrix,Bz,N);

					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
				
					gauss(matrix,Bh,N);
				
					double Hix=Bx[3];//節点iにおけるHx
					double Hiy=By[3];//節点iにおけるHy
					double Hiz=Bz[3];//節点iにおけるHz
					HHi=Bh[3];		//WLSMにより求まったHを格納しなおす

					double Mx,My,Mz,MM;
					if(CON.get_NLMH()==OFF)
					{
						MM=kai*HHi;
						Mx=kai*Hix;
						My=kai*Hiy;
						Mz=kai*Hiz;
					}
					else if(CON.get_NLMH()==ON)
					{
						double HH2=sqrt(Hix*Hix+Hiy*Hiy+Hiz*Hiz);//Hの大きさ
						MM=Ms*(cosh(gamma*HH2)/sinh(gamma*HH2)-1/(gamma*HH2));//ランジュバンの式によるMの大きさ
						Mx=MM*Hix/HH2;
						My=MM*Hiy/HH2;
						Mz=MM*Hiz/HH2;
					}

				//	double fx=u0*(Mx*Bx[0]+My*Bx[1]+Mz*Bx[2]);//体積力fx
				//	double fy=u0*(Mx*By[0]+My*By[1]+Mz*By[2]);//体積力fy
				//	double fz=u0*(Mx*Bz[0]+My*Bz[1]+Mz*Bz[2]);//体積力fz
	
					double fx=u0*MM*Bh[0];//体積力fx
					double fy=u0*MM*Bh[1];//体積力fy
					double fz=u0*MM*Bh[2];//体積力fz
					
					int p=NODE[i].particleID;//対応する粒子番号
				
					if(p>=0)//対応する粒子があるならF[D][p]に代入
					{
						if(PART[p].surface==OFF)//内部流体粒子のみ力を計算
						{
							F[A_X][p]=fx*V;
							F[A_Y][p]=fy*V;//単位は[N]
							F[A_Z][p]=fz*V;
						}
					}
					Fz+=fz*V;
				}
			}
		}
	}///////////*/

	//Pm=u0∫MdHの勾配をWLSMで求める
	if(order==1)//3次元1次式
	{
		///係数行列は
		///   ΣΔx2    ΣΔxΔy  ΣΔxΔz ΣΔx  a = ΣΔxΔf  
		///  ΣΔxΔy   ΣΔy2    ΣΔyΔz ΣΔy  b = ΣΔyΔf 
		///  ΣΔxΔz   ΣΔyΔz  ΣΔz2   ΣΔz  c = ΣΔzΔf 
		///  ΣΔx      ΣΔy     ΣΔz    Σ1    d = ΣΔf

		double Ri[3];//節点iの座標格納
		double Rg[3];	//要素の重心格納
		for(int i=1;i<=node;i++)
		{
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC)
			{
				for(int n=0;n<N*N;n++) matrix[n]=0;		//初期化
				for(int n=0;n<N;n++)  {Bx[n]=0; By[n]=0; Bz[n]=0;}
				
				for(int D=0;D<3;D++) Ri[D]=NODE[i].r[D];
				if(jnb[i]>5)
				{
					double Hi[3]={0,0,0};//節点iにおけるH格納
					int num_nei=0;			//計算に寄与する要素数カウント(空気は除く)
					for(int j=1;j<=jnb[i];j++)//節点iはjnb[n]個の要素に隣接している
					{
						int jelm=nei[i][j];
						//if(ELEM[jelm].material==WATER)
						{
							num_nei++;
							double u=RP[jelm]*u0;
							for(int D=0;D<3;D++) 
							{
								H[D]=Be[D][jelm]/u;//jelmにおけるトータルな磁場　本当は外部磁場じゃないといけないらしいけど？
								Hi[D]+=H[D];
							}
							double HHj=sqrt(H[A_X]*H[A_X]+H[A_Y]*H[A_Y]+H[A_Z]*H[A_Z]);//jelmのHの大きさ
							double integ_of_M=0;//u0∫MdH
							//if(ELEM[jelm].material==WATER)
							{
								if(CON.get_NLMH()==ON)	integ_of_M=u0*Ms/gamma*(log(sinh(gamma*HHj))-log(gamma*HHj));//ランジュバンの式
								else integ_of_M=u0*kai*HHj*HHj*0.5;
							}
							
							for(int D=0;D<3;D++) Rg[D]=0;
							for(int k=1;k<=4;k++)
							{
								int node=ELEM[jelm].node[k];//jelmを構成する節点番号
								for(int D=0;D<3;D++) Rg[D]+=NODE[node].r[D]*0.25;
							}//jelmの重心Rg[D]がもとまった。
		
							double X=Rg[A_X]-Ri[A_X];
							double Y=Rg[A_Y]-Ri[A_Y];
							double Z=Rg[A_Z]-Ri[A_Z];
							double dis=sqrt(X*X+Y*Y+Z*Z);
						
							double w=1;
							if(dis>le) w=le*le/(dis*dis);
			
							matrix[0]+=X*X*w;		
							matrix[1]+=X*Y*w;	
							matrix[2]+=X*Z*w;	
							matrix[3]+=X*w;
		
							matrix[5]+=Y*Y*w;		
							matrix[6]+=Y*Z*w;
							matrix[7]+=Y*w;
	
							matrix[10]+=Z*Z*w;
							matrix[11]+=Z*w;
		
							matrix[15]+=w;
		
							Bx[0]+=X*w*integ_of_M;
		
							Bx[1]+=Y*w*integ_of_M;
		
							Bx[2]+=Z*w*integ_of_M; 
		
							Bx[3]+=w*integ_of_M;
						}
					}
					matrix[4]=matrix[1];
					matrix[8]=matrix[2];
					matrix[9]=matrix[6];
					matrix[12]=matrix[3];
					matrix[13]=matrix[7];
					matrix[14]=matrix[11];
		
					if(num_nei>0) for(int D=0;D<3;D++) Hi[D]/=num_nei;//節点iにおけるH
					double HHi=sqrt(Hi[A_X]*Hi[A_X]+Hi[A_Y]*Hi[A_Y]+Hi[A_Z]*Hi[A_Z]);//iのHの大きさ
					double integ_of_Mi;//iのu0∫MdH
					if(CON.get_NLMH()==ON) integ_of_Mi=u0*Ms/gamma*(log(sinh(gamma*HHi))-log(gamma*HHi));//ランジュバンの式
					else integ_of_Mi=0.5*u0*HHi*HHi*kai;
	
					matrix[15]+=1;//自身
					Bx[3]+=integ_of_Mi;
	
					for(int n=0;n<N*N;n++) matrix_val[n]=matrix[n];//値を保存
	
					//ガウスの消去法で解く
					gauss(matrix,Bx,N);
					
					for(int n=0;n<N*N;n++) matrix[n]=matrix_val[n];//値を戻す
					
					integ_of_Mi=Bx[3];//節点iにおけるu0∫MdH
	
					double fx=Bx[0];//体積力fx
					double fy=Bx[1];//体積力fy
					double fz=Bx[2];//体積力fz
				
					int p=NODE[i].particleID;
					
					if(p>=0)//対応する粒子があるならF[D][p]に代入
					{
						if(PART[p].surface==OFF)//内部流体粒子のみ力を計算
						{
							F[A_X][p]-=fx*V;
							F[A_Y][p]-=fy*V;//単位は[N]
							F[A_Z][p]-=fz*V;
						}
					}
					Fz+=fz*V;
				}
			}
		}
	}///////////////////*/

	//表面積分実行
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==FLUID || ELEM[i].material==ELASTIC)
		{
			double u=RP[i]*u0;
			for(int j=1;j<=4;j++)
			{
				int jelm=ELEM[i].elm[j];
				if(ELEM[jelm].material==AIR)
				{
					int ia=ELEM[i].node[j%4+1];//頂点ｊの向いに位置する三角形の頂点 ia→ib→icはｊからみて半時計まわり
					int ib=ELEM[i].node[4-(j-1)/2*2];
					int ic=ELEM[i].node[3-(j/2%2)*2];

					double iaic[3];//ia→icのベクトル成分格納
					double iaib[3];//ia→ibのベクトル成分格納
					for(int D=0;D<3;D++)
					{
						iaic[D]=NODE[ic].r[D]-NODE[ia].r[D];
						iaib[D]=NODE[ib].r[D]-NODE[ia].r[D];
					}
					///0.5(iaic×iaib)は三角形ia,ib,icの面積の大きさをもち、方向は外向きのベクトルとなる
					double S[3];//上記のベクトル成分格納
					S[A_X]=0.5*(iaic[A_Y]*iaib[A_Z]-iaic[A_Z]*iaib[A_Y]);
					S[A_Y]=0.5*(iaic[A_Z]*iaib[A_X]-iaic[A_X]*iaib[A_Z]);
					S[A_Z]=0.5*(iaic[A_X]*iaib[A_Y]-iaic[A_Y]*iaib[A_X]);
					double SS=sqrt(S[A_X]*S[A_X]+S[A_Y]*S[A_Y]+S[A_Z]*S[A_Z]);
					////面積Sがもとまった

					double n[3]={S[A_X]/SS,S[A_Y]/SS,S[A_Z]/SS};//法線ベクトル
					double H[3]={Be[A_X][i]/u,Be[A_Y][i]/u,Be[A_Z][i]/u};//H
					double ave_B[3]={(Be[A_X][i]+Be[A_X][jelm])*0.5,(Be[A_Y][i]+Be[A_Y][jelm])*0.5,(Be[A_Z][i]+Be[A_Z][jelm])*0.5};//界面のB
					double H2[3]={ave_B[A_X]/u0,ave_B[A_Y]/u0,ave_B[A_Z]/u0};//空気側のH
					double H1[3]={ave_B[A_X]/u,ave_B[A_Y]/u,ave_B[A_Z]/u};//流体側のH
					double Hn2=H2[A_X]*n[A_X]+H2[A_Y]*n[A_Y]+H2[A_Z]*n[A_Z];//H2の法線成分
					double Hn1=H1[A_X]*n[A_X]+H1[A_Y]*n[A_Y]+H1[A_Z]*n[A_Z];
					double Mn=Hn2-Hn1;//ferrohydrodynamicsの式5.21b参照
					double BB=sqrt(ave_B[A_X]*ave_B[A_X]+ave_B[A_Y]*ave_B[A_Y]+ave_B[A_Z]*ave_B[A_Z]);//界面のBの大きさ
					double Bn=ave_B[A_X]*n[A_X]+ave_B[A_Y]*n[A_Y]+ave_B[A_Z]*n[A_Z];
					
					double M[3];//磁化
					for(int D=0;D<3;D++) M[D]=0.5*(H1[D]+H2[D])*kai;
					double MM=sqrt(M[A_X]*M[A_X]+M[A_Y]*M[A_Y]+M[A_Z]*M[A_Z]);//磁化の大きさ
					double HH=MM/kai;
					/*if(CON.get_NLMH()==OFF)
					{
						MM=kai*HH;
						for(int D=0;D<3;D++) M[D]=H[D]*kai;
					}
					else
					{
						MM=Ms*(cosh(gamma*HH)/sinh(gamma*HH)-1/(gamma*HH));//ランジュバンの式によるMの大きさ
						for(int D=0;D<3;D++) M[D]=MM*H[D]/HH;
					}*/
					//double Mn=0;
					//for(int D=0;D<3;D++) Mn+=M[D]*n[D];//法線方向の磁化
					//cout<<Mn<<endl;*/
					double Fs[3]={0,0,0};
					for(int D=0;D<3;D++) Fs[D]=0.5*u0*Mn*Mn*SS*n[D];//表面に働く力 //単位は[N]
					//for(int D=0;D<3;D++) Fs[D]+=0.5*u0*MM*MM*SS*n[D];//???
					//for(int D=0;D<3;D++) Fs[D]+=0.5*u0*kai*HH*HH*SS*n[D];//???

					double integ_of_Mi;//u0∫MdH
					if(CON.get_NLMH()==ON) integ_of_Mi=u0*Ms/gamma*(log(sinh(gamma*HH))-log(gamma*HH));//ランジュバンの式
					else integ_of_Mi=0.5*u0*HH*HH*kai;
					for(int D=0;D<3;D++) Fs[D]+=integ_of_Mi*SS*n[D];//表面に働く力 //単位は[N]

					Fz+=Fs[A_Z];

					//得られたFsを頂点に割り当てる
					int node[3]={ia,ib,ic};//頂点の番号を記憶
					for(int k=0;k<3;k++)
					{
						int node_ID=node[k];//節点番号
						if(NODE[node_ID].particleID>=0)//流体節点なら
						{
							int p=NODE[node_ID].particleID;
							for(int D=0;D<3;D++) F[D][p]+=Fs[D]/3;
							for(int D=0;D<3;D++) fs[p]+=Fs[D]/SS*n[D];
							num[p]++;
						}
					}
				}
			}
		}
	}///*/

	//磁気双曲子間の粒子間力を計算
	//coroid_pole(CON,PART,NODE,ELEM, node, nelm,Be, fluid_number,F,RP,jnb,nei);

	int *check=new int[fluid_number];//check[i]=ONならその粒子は対応する節点が存在するということ
	for(int i=0;i<fluid_number;i++) check[i]=OFF;
	for(int I=1;I<=node;I++)
	{
		if(NODE[I].material==FLUID || NODE[I].material==ELASTIC)
		{	
			int i=NODE[I].particleID;//粒子番号
			if(i>=0) check[i]=ON;
		}
	}
	
	for(int i=0;i<fluid_number;i++)
	{
		if(check[i]==OFF)//対応する節点がないので、力がゼロのままな粒子
		{
			double W=0;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(j<fluid_number)
				{
					if(check[j]==ON)
					{
						double X=PART[j].r[A_X]-PART[i].r[A_X];
						double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
						double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
						double dis=sqrt(X*X+Y*Y+Z*Z);//粒子間の距離
						double w=R/dis-1;//重み関数
						W+=w;
						for(int D=0;D<3;D++) F[D][i]+=F[D][j]*w;
					}
				}
			}
			if(W!=0) for(int D=0;D<3;D++) F[D][i]/=W;
		}
	}

	if(CON.get_dir_for_P()==2 ||CON.get_dir_for_P()==3 )
    {
		for(int i=0;i<fluid_number;i++)
		{
			double val=0;//表面力
			if(PART[i].surface==ON)
			{
				if(num[i]>0) val=fs[i]/num[i];
				else	//デローニ分割でjnb=0となった粒子はここではnum[i]=0となってしまいディリクレ値が決まらない。そういうときは周辺粒子より値を補間する
				{
					double count=0;
					for(int k=0;k<PART[i].N;k++)
					{
						int j=PART[i].NEI[k];
						if(j<fluid_number && PART[j].surface==ON)
						{
							if(num[j]>0)
							{
								val+=fs[j];
								count+=1;
							}
						}
					}
					if(count>0) val/=count;		//ここでcount=0となるようならその粒子は完全に孤立粒子なので、ディリクレ値はゼロを格納させておけばよい
				}
				for(int D=0;D<3;D++) F[D][i]=0;//圧力デﾞｨﾘｸﾚとして電磁力を使用するのでここでは初期化
			}
			PART[i].dir_Pem=-val;
        }
	}
	ofstream fg("Fss.dat");
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			fg<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].dir_Pem*-1<<endl;
		}
	}

	delete [] matrix;
	delete [] matrix_val;
	delete [] Bx;
	delete [] By;
	delete [] Bz;
	delete [] Bh;

	delete [] check;
	for(int D=0;D<DIMENSION;D++) delete [] direct[D];
	delete [] fs;
	delete [] num;

	cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
	cout<<"Fx="<<Fx<<endl;
}

void coroid_pole(mpsconfig &CON, vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int*jnb,int **nei)
{
	double le=CON.get_distancebp();
    double V=CON.get_particle_volume();//粒子の体積
	double u0=1.257e-6;	
	double u=CON.get_RP()*u0;
	double kai=CON.get_RP()-1;

	double *M[DIMENSION];
	for(int D=0;D<DIMENSION;D++) M[D]=new double [fluid_number];//各粒子位置でのM
	double *F_di[DIMENSION];
	for(int D=0;D<DIMENSION;D++) F_di[D]=new double [fluid_number];//この関数で計算されるF
	double *B_di[DIMENSION];
	for(int D=0;D<DIMENSION;D++) B_di[D]=new double [fluid_number];//この関数で計算されるB
	int *check=new int[fluid_number];//check=ONの粒子は対応する節点が存在する。ゆえに対応するMが手に入る。check=OFFなら周囲からMの値を補間しなければならない

	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<DIMENSION;D++)
		{
			M[D][i]=0;//初期化
			F_di[D][i]=0;//初期化
			B_di[D][i]=0;//初期化
		}
		check[i]=OFF;
	}

	for(int i=1;i<=node;i++)
	{
		if(NODE[i].particleID>=0)//対応する粒子があるなら
		{
			int ip=NODE[i].particleID;
			int num=0;
			check[ip]=ON;				//対応する節点が存在する
			for(int j=1;j<=jnb[i];j++)
			{
				int jelm=nei[i][j];//隣接する要素番号
				if(ELEM[jelm].material==FLUID || ELEM[jelm].material==ELASTIC)
				{
					double H[3];
					for(int D=0;D<3;D++) H[D]=Be[D][jelm]/u;
					for(int D=0;D<3;D++) M[D][ip]+=H[D]*kai;
					num++;
				}
			}
			if(num>0) for(int D=0;D<3;D++) M[D][ip]/=num;
		}
	}

	for(int i=0;i<fluid_number;i++)
	{
		if(check[i]==OFF)
		{
			int num=0;
			for(int k=0;k<PART[i].N;k++)
			{
				int j=PART[i].NEI[k];
				if(PART[j].type==FLUID || PART[j].type==ELASTIC) if(check[j]==ON) for(int D=0;D<3;D++) M[D][i]+=M[D][j];
			}
			if(num>0) for(int D=0;D<3;D++) M[D][i]/=num;
			check[i]=ON;
		}
	}//全粒子のMが求まった

	for(int i=0;i<fluid_number;i++) for(int D=0;D<3;D++) M[D][i]*=V*0.32;//磁気モーメントmに置き換え

	//力を計算
	double co=u0/(4*PI);//よく使う係数
	double f[3]={0,0,0};
	double L=le/10;
	for(int i=0;i<fluid_number;i++)
	{
		for(int k=0;k<PART[i].N3;k++)
		{
			int j=PART[i].NEI3[k];
			if(PART[j].type==FLUID || PART[j].type==ELASTIC)
			{	
				/*double dX[3];
				for(int D=0;D<3;D++) dX[D]=PART[i].r[D]-PART[j].r[D];//jからiへ向かうベクトル
				double r=sqrt(dX[A_X]*dX[A_X]+dX[A_Y]*dX[A_Y]+dX[A_Z]*dX[A_Z]);//iとjの距離
			
				double MMi=sqrt(M[A_X][i]*M[A_X][i]+M[A_Y][i]*M[A_Y][i]+M[A_Z][i]*M[A_Z][i]);//Miの大きさ
				double MMj=sqrt(M[A_X][j]*M[A_X][j]+M[A_Y][j]*M[A_Y][j]+M[A_Z][j]*M[A_Z][j]);//Mjの大きさ
				double cos1=(dX[A_X]*M[A_X][i]+dX[A_Y]*M[A_Y][i]+dX[A_Z]*M[A_Z][i])/(r*MMi);		//rとMiのなす角度
				double cos2=(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])/(r*MMj);		//rとMjのなす角度
				double cos3=(M[A_X][i]*M[A_X][j]+M[A_Y][i]*M[A_Y][j]+M[A_Z][i]*M[A_Z][j])/(MMj*MMi);		//MjとMiのなす角度
				for(int D=0;D<3;D++) F_di[D][i]+=-co*3*MMi*MMj/pow(r,4.0)*(3*cos1*cos2-cos3)*dX[D]/r;
				//for(int D=0;D<3;D++) F_di[D][i]+=co*3*MMi*MMj/pow(r,4.0)*(cos1+cos2-5*cos1*cos2+cos3)*dX[D]/r;
				for(int D=0;D<3;D++) f[D]+=-co*3*MMi*MMj/pow(r,4.0)*(3*cos1*cos2-cos3)*dX[D]/r;*/

				double Xnp[6];
				double Ynp[6];
				double Znp[6];
				double Bx[6],By[6],Bz[6];
				Xnp[0]=PART[i].r[A_X]-L; Xnp[1]=PART[i].r[A_X]+L; Xnp[2]=PART[i].r[A_X]; Xnp[3]=PART[i].r[A_X]; Xnp[4]=PART[i].r[A_X]; Xnp[5]=PART[i].r[A_X];
				Ynp[0]=PART[i].r[A_Y]; Ynp[1]=PART[i].r[A_Y]; Ynp[2]=PART[i].r[A_Y]-L; Ynp[3]=PART[i].r[A_Y]+L; Ynp[4]=PART[i].r[A_Y]; Ynp[5]=PART[i].r[A_Y]; 
				Znp[0]=PART[i].r[A_Z]; Znp[1]=PART[i].r[A_Z]; Znp[2]=PART[i].r[A_Z]; Znp[3]=PART[i].r[A_Z]; Znp[4]=PART[i].r[A_Z]-L; Znp[5]=PART[i].r[A_Z]+L; 

				for(int k=0;k<6;k++)
				{
					double dX[3];
					dX[A_X]=Xnp[k]-PART[j].r[A_X];//jからiへ向かうベクトル
					dX[A_Y]=Ynp[k]-PART[j].r[A_Y];//jからiへ向かうベクトル
					dX[A_Z]=Znp[k]-PART[j].r[A_Z];//jからiへ向かうベクトル
					double r=sqrt(dX[A_X]*dX[A_X]+dX[A_Y]*dX[A_Y]+dX[A_Z]*dX[A_Z]);//iとjの距離
				
					Bx[k]=co*3*(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])*dX[A_X]/pow(r,5.0)-co*M[A_X][j]/pow(r,3.0);
					By[k]=co*3*(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])*dX[A_Y]/pow(r,5.0)-co*M[A_Y][j]/pow(r,3.0);
					Bz[k]=co*3*(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])*dX[A_Z]/pow(r,5.0)-co*M[A_Z][j]/pow(r,3.0);
				}

				F_di[A_X][i]+=M[A_X][i]*(Bx[1]-Bx[0])/(2*L)+M[A_Y][i]*(Bx[3]-Bx[2])/(2*L)+M[A_Z][i]*(Bx[5]-Bx[4])/(2*L);
				F_di[A_Y][i]+=M[A_X][i]*(By[1]-By[0])/(2*L)+M[A_Y][i]*(By[3]-By[2])/(2*L)+M[A_Z][i]*(By[5]-By[4])/(2*L);
				F_di[A_Z][i]+=M[A_X][i]*(Bz[1]-Bz[0])/(2*L)+M[A_Y][i]*(Bz[3]-Bz[2])/(2*L)+M[A_Z][i]*(Bz[5]-Bz[4])/(2*L);

				//Bを求める
				double dX[3];
				for(int D=0;D<3;D++) dX[D]=PART[i].r[D]-PART[j].r[D];//jからiへ向かうベクトル
				double r=sqrt(dX[A_X]*dX[A_X]+dX[A_Y]*dX[A_Y]+dX[A_Z]*dX[A_Z]);//iとjの距離
				
				for(int D=0;D<3;D++) B_di[D][i]+=co*3*(dX[A_X]*M[A_X][j]+dX[A_Y]*M[A_Y][j]+dX[A_Z]*M[A_Z][j])*dX[D]/pow(r,5.0)-co*M[D][j]/pow(r,3.0);

			}
		}
	}
	//F[D][i]に加える
	for(int i=0;i<fluid_number;i++)
	{
		for(int D=0;D<3;D++)
		{
			//F[D][i]+= F_di[D][i];
			f[D]+=F_di[D][i];
		}
	}

	cout<<f[A_X]<<" "<<f[A_Y]<<" "<<f[A_Z]<<endl;

	ofstream fp("M.dat");
	ofstream fq("F_di.dat");
	ofstream fr("B_di.dat");
	double times=1000;
	double timesF=CON.get_times()/CON.get_density()/le*CON.get_FEMtimes();
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le)
		{
			fp<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<M[A_X][i]*times<<" "<<M[A_Z][i]*times<<endl;
			fq<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<F_di[A_X][i]*times<<" "<<F_di[A_Z][i]*times<<endl;
			fr<<PART[i].r[A_X]<<" "<<PART[i].r[A_Z]<<" "<<B_di[A_X][i]*CON.get_B_times()<<" "<<B_di[A_Z][i]*CON.get_B_times()<<endl;
		}
	}
	fp.close();
	fq.close();
	fr.close();

	for(int D=0;D<DIMENSION;D++) delete [] M[D];
	for(int D=0;D<DIMENSION;D++) delete [] F_di[D];
	for(int D=0;D<DIMENSION;D++) delete [] B_di[D];
	delete [] check;
}

//節点力法計算関数
void NODE_F3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpselastic> &PART,double **F,int fluid_number)
{
	cout<<"節点力法による電磁力計算--------";
	double ep0=8.854e-12;	//真空の誘電率。
	double u0=12.5e-7;		//真空の透磁率
	int N[4+1];				//要素の各節点番号格納
	double X[4+1];
	double Y[4+1];
	double Z[4+1];
	double Fz=0;			//Z方向の合力[N]
	double Fx=0;
	unsigned timeA=GetTickCount();//計算開始時刻

	double *Fn[3];
	for(int D=0;D<3;D++) Fn[D]=new double [node+1];//

	for(int i=0;i<=node;i++) for(int D=0;D<3;D++) Fn[D][i]=0;//初期化

	if(CON.get_EM_calc_type()==1)//静電力計算
	{
		for(int I=1;I<=node;I++)
		{
			if(NODE[I].material==FLUID || NODE[I].material==MAGELAST  || NODE[I].material==ELASTIC || NODE[I].material==IRON)
			{			
				for(int k=1;k<=jnb[I];k++)
				{
					int jelm=nei[I][k];//節点Iが隣接する要素番号
				
					double ep=ep0;
					if(ELEM[jelm].material==FLUID) ep*=CON.get_r_perm();

					//マクスウェルの応力テンソル
					double Txx=0.5*ep*(Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Txy=ep*Ee[A_X][jelm]*Ee[A_Y][jelm];
					double Txz=ep*Ee[A_X][jelm]*Ee[A_Z][jelm];
					double Tyx=Txy;
					double Tyy=0.5*ep*(-Ee[A_X][jelm]*Ee[A_X][jelm]+Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Tyz=ep*Ee[A_Y][jelm]*Ee[A_Z][jelm];
					double Tzx=Txz;
					double Tzy=Tyz;
					double Tzz=0.5*ep*(-Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]+Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					//////////
		   
					/////係数c,d,e計算
					for(int j=1;j<=4;j++)
					{
						N[j]=ELEM[jelm].node[j];
						X[j]=NODE[N[j]].r[A_X];
						Y[j]=NODE[N[j]].r[A_Y];
						Z[j]=NODE[N[j]].r[A_Z];
					}

					int i=0;///節点iは要素jelmの第j番目の節点
					for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==I) i=j;
					if(i==0) cout<<"ERROR IN Fn"<<endl;

					int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
					int m=j%4+1;
					int n=m%4+1;
					//delta6は相殺されるのでいらない
					double c=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
					double d=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
					double e=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
					if(i%2!=0)//iが奇数なら
					{
						c*=-1;
						d*=-1;
						e*=-1;
					}
					///////////////
		
					Fn[A_X][I]+=Txx*c+Txy*d+Txz*e;
					Fn[A_Y][I]+=Tyx*c+Tyy*d+Tyz*e;
					Fn[A_Z][I]+=Tzx*c+Tzy*d+Tzz*e;
				}
				for(int D=0;D<3;D++) Fn[D][I]*=-1.00000000000/6.00000000000;
				Fz+=Fn[A_Z][I];
				
			}
		}
	}
	else if(CON.get_EM_calc_type()==4)//磁位解析 EeにはＨが格納されていることに注意
	{
		for(int I=1;I<=node;I++)
		{
			if(NODE[I].material==FLUID || NODE[I].material==MAGELAST || NODE[I].material==ELASTIC || NODE[I].material==IRON)
			{
				for(int k=1;k<=jnb[I];k++)
				{
					int jelm=nei[I][k];//節点iが隣接する要素番号
					//if(ELEM[jelm].material==AIR)
					{
						double u=u0*RP[jelm];

						//マクスウェルの応力テンソル
						double Txx=0.5*u*(Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
						double Txy=u*Ee[A_X][jelm]*Ee[A_Y][jelm];
						double Txz=u*Ee[A_X][jelm]*Ee[A_Z][jelm];
						double Tyx=Txy;
						double Tyy=0.5*u*(-Ee[A_X][jelm]*Ee[A_X][jelm]+Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
						double Tyz=u*Ee[A_Y][jelm]*Ee[A_Z][jelm];
						double Tzx=Txz;
						double Tzy=Tyz;
						double Tzz=0.5*u*(-Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]+Ee[A_Z][jelm]*Ee[A_Z][jelm]);
						//*/
						
					//係数c,d,e計算
						for(int j=1;j<=4;j++)
						{
							N[j]=ELEM[jelm].node[j];
							X[j]=NODE[N[j]].r[A_X];
							Y[j]=NODE[N[j]].r[A_Y];
							Z[j]=NODE[N[j]].r[A_Z];
						}
						int i=0;///節点iは要素jelmの第J番目の節点
						for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==I) i=j;

						//if(i==0) cout<<"ERROR IN Fn"<<endl;
						int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
						int m=j%4+1;
						int n=m%4+1;
						//delta6は相殺されるのでいらない
						double c=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
						double d=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
						double e=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
						if( i & 1 )//iが奇数なら
						{
							c*=-1;
							d*=-1;
							e*=-1;
						}
						///////////////
			
						
						
						Fn[A_X][I]+=(Txx*c+Txy*d+Txz*e);
						Fn[A_Y][I]+=(Tyx*c+Tyy*d+Tyz*e);
						Fn[A_Z][I]+=(Tzx*c+Tzy*d+Tzz*e);
						//cout<<(Tzx*c+Tzy*d+Tzz*e)<<" "<<c<<" "<<d<<" "<<e<<endl;
					}
				}
				for(int D=0;D<3;D++) Fn[D][I]*=-1.00000000000/6.00000000000;
				//if(Fn[A_Z][I]<0) Fn[A_Z][I]=0;
				Fz+=Fn[A_Z][I];
			}
		}
	}
	else //磁場解析
	{
	for(int I=1;I<=node;I++)
	{
		if(NODE[I].material==FLUID || NODE[I].material==ELASTIC || NODE[I].material==MAGELAST || NODE[I].material==IRON)
			{
				for(int k=1;k<=jnb[I];k++)
				{
					int jelm=nei[I][k];//節点iが隣接する要素番号
						
					//if(ELEM[jelm].material==AIR){
					//マクスウェルの応力テンソル
					double Txx=(Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Txy=2*Ee[A_X][jelm]*Ee[A_Y][jelm];
					double Txz=2*Ee[A_X][jelm]*Ee[A_Z][jelm];
					double Tyx=Txy;
					double Tyy=(-Ee[A_X][jelm]*Ee[A_X][jelm]+Ee[A_Y][jelm]*Ee[A_Y][jelm]-Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					double Tyz=2*Ee[A_Y][jelm]*Ee[A_Z][jelm];
					double Tzx=Txz;
					double Tzy=Tyz;
					double Tzz=(-Ee[A_X][jelm]*Ee[A_X][jelm]-Ee[A_Y][jelm]*Ee[A_Y][jelm]+Ee[A_Z][jelm]*Ee[A_Z][jelm]);
					//////////
					
					//係数c,d,e計算
					for(int j=1;j<=4;j++)
					{
						N[j]=ELEM[jelm].node[j];
						X[j]=NODE[N[j]].r[A_X];
						Y[j]=NODE[N[j]].r[A_Y];
						Z[j]=NODE[N[j]].r[A_Z];
					}

					int i=0;//節点iは要素jelmの第j番目の節点
					for(int j=1;j<=4;j++) if(ELEM[jelm].node[j]==I) i=j;
					int j=i%4+1;//i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
					int m=j%4+1;
					int n=m%4+1;

				//delta6は相殺されるのでいらない
					double c=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
					double d=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
					double e=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
					if(i & 1)//iが奇数なら
					{
						c*=-1;
						d*=-1;
						e*=-1;
					}
					///////////////
		
					double u=RP[jelm];
						
					Fn[A_X][I]+=(Txx*c+Txy*d+Txz*e)/(2*u);//u0はこの下でかけてる
					Fn[A_Y][I]+=(Tyx*c+Tyy*d+Tyz*e)/(2*u);
					Fn[A_Z][I]+=(Tzx*c+Tzy*d+Tzz*e)/(2*u);
					//}
				}
				for(int D=0;D<3;D++) Fn[D][I]*=-1.00000000000/(6.00000000000*u0);
				//if(Fn[A_Z][I]<0) Fn[A_Z][I]=0;

				Fz+=Fn[A_Z][I];
				if(NODE[I].r[A_X]<0) Fx+= Fn[A_X][I];//これは何？？

			}//if(NODE[I].material==FLUID || NODE[I].material==ELASTIC)
		}
	}
	//cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	int *bound=new int[node+1];
	for(int n=1;n<=node;n++)
	{
		bound[n]=OFF;
		if(NODE[n].material==FLUID || NODE[n].material==ELASTIC || NODE[n].material==MAGELAST|| NODE[n].material==IRON)
		{
			for(int k=1;k<=jnb[n];k++)
			{
				int jelm=nei[n][k];//節点iが隣接する要素番号
				if(ELEM[jelm].material==AIR) bound[n]=ON;
			}
		}
	}
	
	

	

	//ファイル出力
	ofstream fp("Fn.dat");
	double le=CON.get_distancebp();
	double times=CON.get_times()/CON.get_density()/le*CON.get_FEMtimes();

	for(int i=1;i<=node;i++)//流体節点のみ出力
	{
		if(NODE[i].material==FLUID || NODE[i].material==ELASTIC || NODE[i].material==IRON) if(NODE[i].r[A_Y]>-le*0.5&& NODE[i].r[A_Y]<le*0.5) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Z]<<" "<<Fn[A_X][i]*times<<" "<<Fn[A_Z][i]*times<<endl;	
	}
	fp.close();//


	for(int D=0;D<3;D++) delete [] Fn[D];
	cout<<"OK Fz="<<Fz<<"[N] time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	if(2*Fz!=Fz+Fz){			//エラーがでたらFerrorをON
		CON.set_Ferror(ON);
		cout<<"Fzでエラー、前回の電磁力を使用します";
	}
	else {
		cout<<"not error";
		CON.set_Ferror(OFF);	//エラーがないなら
		///////////F更新////////////
		for(int n=1;n<=node;n++)
		{
			if(bound[n]==ON)
			{
				int i=NODE[n].particleID;
				if(i>=0) for(int D=0;D<3;D++) F[D][i]=Fn[D][n];
			}
		}
		//節点力法の力の単位は[N]なので、ここからさらに粒子≠節点の粒子に関しても力をもとめてやる必要はない
		cout<<"!"<<endl;
	}
	delete [] bound;
	cout<<"ok1"<<endl;
}

//磁束密度をプロット(加嶋作成)
void plot_magnetic_flux_density(mpsconfig &CON, vector<mpselastic> &PART, vector<point3D> &NODE, vector<element3D> &ELEM, int nelm, double **B, int t)
{

	//磁束密度スムージング関数
	double le=CON.get_distancebp();
//とりあえずスムージングかけずにやる

	/*  double *newB[3];
	for(int D=0;D<3;D++) newB[D]=new double [nelm+1];
	if(CON.get_FEM_smn()>0)
	{
		for(int n=0;n<CON.get_FEM_smn();n++)
		{
		    for(int i=1;i<=nelm;i++) 
		    {  
		        for(int D=0;D<3;D++) newB[D][i]=B[D][i];
				int num=1; //自分自身をカウントするから1
				for(int k=0;k<PART[i].N;k++)
				{       
					int j=PART[i].NEI[k];
					if(PART[j].type==FLUID || PART[j].type==ELASTIC)
					{
						num++;
						for(int D=0;D<3;D++) newB[D][i]+=B[D][j];
					}
				}
				for(int D=0;D<3;D++) newB[D][i]/=num;
		    } 
		    for(int i=1;i<nelm;i++) for(int D=0;D<3;D++) B[D][i]=newB[D][i];
		}
	}
	else if(CON.get_FEM_smn()<0)//表面のみでスムージング
	{
		int N=-1*CON.get_FEM_smn();
		for(int n=0;n<N;n++)
		{
			for(int i=1;i<=nelm;i++) 
			{  
			    for(int D=0;D<3;D++) newB[D][i]=B[D][i];
				if(PART[i].surface==ON)
				{
					int num=1; //自分自身をカウントするから1
					for(int k=0;k<PART[i].N;k++)
					{       
						int j=PART[i].NEI[k];
						if(PART[j].surface==ON && (PART[j].type==FLUID || PART[j].type==ELASTIC))
						{
							num++;
							for(int D=0;D<3;D++) newB[D][i]+=B[D][j];
						}
					}
					for(int D=0;D<3;D++) newB[D][i]/=num;
				}
			} 
			for(int i=1;i<=nelm;i++) for(int D=0;D<3;D++) B[D][i]=newB[D][i];
		}
	}
    for(int D=0;D<3;D++) delete [] newB[D];
*/
	ofstream fp("./FluxAVS/Flux.dat");
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	double times=1.0; //=CON.get_times()/CON.get_particle_mass();*le*le*CON.get_FEMtimes(); //?
	
	double **center_of_element=new double* [nelm+1];
	for(int i=1;i<=nelm;i++) center_of_element[i]=new double[3];

	//重心の計算・・・nelm全て出力
	for(int i=1;i<=nelm;i++)
	{
		for(int D=0;D<3;D++)
		{
			center_of_element[i][D]=0.0;
			for(int j=1;j<=4;j++) center_of_element[i][D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		}
	}

	for(int i=1;i<=nelm;i++)
    {
		if(center_of_element[i][A_Y]>-le*0.5 && center_of_element[i][A_Y]<+le*0.5)
		{
			fp<<center_of_element[i][A_X]<<"\t"<<center_of_element[i][A_Z]<<"\t"<<B[A_X][i]*times<<"\t"<<B[A_Z][i]*times<<endl;
			if(center_of_element[i][A_X]>xmax) xmax=center_of_element[i][A_X];
			if(center_of_element[i][A_Z]>ymax) ymax=center_of_element[i][A_Z];
		}
    }
	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	if(CON.get_B_times()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_B_times()*times<<" "<<0*times<<endl;//最後に凡例出力
	fp.close();////

	//////////////////////////////////////
//	if(t==1 || t%(CON.get_EM_interval()*2)==0)
//	{
	int timestep=CON.get_current_step();
	int half_nelm=0;
	for(int i=1;i<=nelm;i++){				//X-Z面に表示する断面までの要素数
		if(center_of_element[i][A_Y]<0.0) half_nelm++;
	}
	///////////////////////////////////////磁束密度ベクトルファイル///////////////////////////////////////////////
		stringstream ssfd;
		ssfd<<"./FluxAVS/FluxDensity"<<t<<".fld";
		string FluxDensity=ssfd.str();	

		ofstream ffd(FluxDensity);
		if(ffd.fail()){
			system("mkdir FluxAVS");
		ofstream ffd(FluxDensity);
		if(ffd.fail()){
			cout<<"ファイルを開けませんでした。FluxAVSフォルダがあるかどうか確認してください"<<endl;
			exit(1);
			}
		}

		ffd << "# AVS field file" << endl;
		ffd << "ndim=1" << endl;
		ffd << "dim1=" << half_nelm <<endl;
		ffd << "nspace=3" << endl;
		ffd << "veclen=3" << endl;
		ffd << "data=float" << endl;
		ffd << "field=irregular" << endl;
		ffd << "label=e-x e-y e-z" << endl << endl;
		ffd << "variable 1 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		ffd << "variable 2 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		ffd << "variable 3 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		ffd << "coord    1 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		ffd << "coord    2 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		ffd << "coord    3 file=./FluxDensity"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;

		ffd.close();

		stringstream ssfdd;
		ssfdd<<"./FluxAVS/FluxDensity"<<t;
		string FluxDensitydata=ssfdd.str();

		ofstream fout(FluxDensitydata);
		if(fout.fail()){
			cout<<"データファイルを開けませんでした。"<<endl;
			exit(1);
		}

		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=1;i<=nelm;i++)
		{
			if(center_of_element[i][A_Y]<0){
				fout<<B[A_X][i]<<" "<<B[A_Y][i]<<" "<<B[A_Z][i]<<" "<<center_of_element[i][A_X]<<" "<<center_of_element[i][A_Y]<<" "<<center_of_element[i][A_Z]<<endl;
			}
		}
		fout.close();
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////磁束密度ベクトル断面図///////////////////////////////////////
		//読み込み用データファイルの作成
		stringstream sscfd;
		sscfd<<"FluxAVS/C_FluxDensity"<<t;
		string CutFluxDensity=sscfd.str();

		ofstream flux_numbered(CutFluxDensity);
		if(flux_numbered.fail()){
			cout<<"ファイルを開けませんでした。FluxAVSフォルダがあるかどうか確認してください2"<<endl;
			exit(1);
		}

		flux_numbered<<"e-x e-y e-z x y z"<<endl;

		int counter=0;//データファイルのデータ点数
		for(int i=1;i<=nelm;i++)
		{	
			if(center_of_element[i][A_Y]>-le*0.5 && center_of_element[i][A_Y]<le*0.5){
				flux_numbered<<B[A_X][i]<<" "<<B[A_Y][i]<<" "<<B[A_Z][i]<<" "<<center_of_element[i][A_X]<<" "<<center_of_element[i][A_Y]<<" "<<center_of_element[i][A_Z]<<endl;
				counter++;
			}
		}
		flux_numbered.close();

		//表示用構造格子型データファイルの作成
		sscfd.clear();
		string filename=CutFluxDensity+".fld";
		ofstream flux(filename);

		flux << "# AVS field file" << endl;
		flux << "ndim=1" << endl;
		flux << "dim1=" << counter <<endl;
		flux << "nspace=3" << endl;
		flux << "veclen=3" << endl;
		flux << "data=float" << endl;
		flux << "field=irregular" << endl;
		flux << "label=e-x e-y e-z" << endl << endl;
		flux << "variable 1 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		flux << "variable 2 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		flux << "variable 3 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		flux << "coord    1 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		flux << "coord	  2 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		flux << "coord    3 file=./C_FluxDensity"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;

		flux.close();
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
/*	//構造格子データファイル
	///////////////////////////////////磁束密度コンター図///////////////////////////////////////////////////////
	stringstream sstr;
	sstr<<"./FluxContour/FluxContour"<<t<<".fld";
	string FluxContour=sstr.str();
	
	ofstream ffc(FluxContour);
	if(ffc.fail()){
		system("mkdir FluxContour");
		ofstream ffc(FluxContour);
		if(ffc.fail()){
		cout<<"./FluxContourフォルダを開けませんでした"<<endl;
		exit(1);
		}
	}
	
	ffc << "# AVS field file" << endl;
	ffc << "ndim=1" << endl;
	ffc << "dim1=" << half_nelm <<endl;
	ffc << "nspace=3" << endl;
	ffc << "veclen=1" << endl;
	ffc << "data=float" << endl;
	ffc << "field=irregular" << endl;
	ffc << "label=FluxContour" << endl << endl; //e-x e-y e-z
	ffc << "variable 1 file=./FluxContour"<<t<<" filetype=ascii skip=1 offset=0 stride=4" << endl;
	ffc << "coord    1 file=./FluxContour"<<t<<" filetype=ascii skip=1 offset=1 stride=4" << endl;
	ffc << "coord    2 file=./FluxContour"<<t<<" filetype=ascii skip=1 offset=2 stride=4" << endl;
	ffc << "coord    3 file=./FluxContour"<<t<<" filetype=ascii skip=1 offset=3 stride=4" << endl;
	ffc.close();
	//データファイル
	stringstream ssfc;
	ssfc<<"./FluxContour/FluxContour"<<t;
	string datafluxcontour=ssfc.str();
	ofstream fout3(datafluxcontour);
	if(fout3.fail()){
		cout<<"./FluxContourフォルダを開けませんでした"<<endl;
		exit(1);
	}
	fout3<<"e-x e-y e-z x y z"<<endl;
	for(int i=1;i<=nelm;i++)
    {
		if(center_of_element[i][A_X]<0){
//		fout<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		fout3<<sqrt(B[A_X][i]*B[A_X][i]+B[A_Y][i]*B[A_Y][i]+B[A_Z][i]*B[A_Z][i])<<" "<<center_of_element[i][A_X]<<" "<<center_of_element[i][A_Y]<<" "<<center_of_element[i][A_Z]<<endl;
		}
	}
	fout3.close();*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	}
	for(int i=1;i<=nelm;i++) delete [] center_of_element[i];
	delete [] center_of_element;
}

//TetGen用 三次元有限要素法
void usingTetGen(mpsconfig &CON,vector<mpselastic> &PART, double **F, int fluid_number,int particle_number,int t,double TIME)
{
	double dt=CON.get_dt();
//	double TIME=t*dt;
	//////////////////////////////////////////////////////////////////////////////////////////
	//TetGen用配列作成
	vector<tetgen_node> NODEall;
	vector<tetgen_facet> FACEall;
	vector<tetgen_element> ELEMall;

	//TRANS[i]には、節点番号iに対応する粒子番号を格納する。
	//FEM3D.cpp では、節点番号が1から始まるので、TRANS[0]には宣言後に-1を入れる。
	vector<int> TRANS;
	TRANS.push_back(-1);	//TRANS[0]には-1

    ///TetGenによるメッシュ生成////////////////////////////////////////////////////

	tetgen_function TETFUNC;
    TETFUNC.call_TetGen(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS);
	
	///////////////////////////////////////////////////////////////////////////////

	int N=static_cast<int>(TRANS.size())-1;	//FEM節点に含まれる粒子数(0番目を除くTRANS[]の長さ)
	int node=static_cast<int>(NODEall.size()); //節点数
	int nelm=static_cast<int>(ELEMall.size()); //要素数
	int KTJ=node; //最大節点数
	int KTE=nelm*6; //最大要素数
	
	int *depth=new int [KTE+1];			//各要素の深さ格納
//	vector<int> depth(KTE+1);

	//配列確保
	vector<point3D> NODE(NODEall.size()+1);
//	NODE.reserve(NODEall.size()+1);
	vector<element3D> ELEM(ELEMall.size()+1);
//	ELEM.reserve(ELEMall.size()+1);
	
	cout<<"NODE.size(): "<<NODE.size()<<", ELEM.size(): "<<ELEM.size()<<endl;
	//TetGenの節点・要素データを取得	節点番号を1つずらす
	//節点データ
//	int num_magnet_attr=0;
    for(int i=0;i<node;i++)
    {
		NODE[i+1].r[A_X]=NODEall[i].r[A_X];
		NODE[i+1].r[A_Y]=NODEall[i].r[A_Y];
		NODE[i+1].r[A_Z]=NODEall[i].r[A_Z];
		NODE[i+1].material=static_cast<int>(NODEall[i].attribute);
		NODE[i+1].boundary_condition=static_cast<int>(NODEall[i].boundary);
		NODE[i].particleID=-1; //-1で初期化
//		if(static_cast<int>(NODEall[i].attribute==MAGNET)) num_magnet_attr++;
    }

	////////////粒子と対応する接点はここで粒子番号を入れる///////////////
	for(int i=0;i<TRANS.size();i++){
		NODE[i].particleID=TRANS[i];
	}
	//////////////////////////////////////////////////////////////////////

//	cout<<"num_magnet_attr: "<<num_magnet_attr<<endl;
	
	//要素データ
    for(int i=0;i<nelm;i++)
	{
		//材質
		ELEM[i+1].material=ELEMall[i].attribute;

		//構成節点
		for(int n=0;n<4;n++) ELEM[i+1].node[n+1]=ELEMall[i].node[n]+1;
		
		//要素-要素関係
		for(int n=0;n<4;n++) 
		{
			if(ELEMall[i].nei_elem[n]==-1){
				ELEM[i+1].elm[n+1]=0;	//近隣要素なしの場合TetGenは-1を返してくるため0に修正する
			}else{
				ELEM[i+1].elm[n+1]=ELEMall[i].nei_elem[n]+1;
			}
		}
	}
	//メッシュ生成のポスト処理
		double *val=new double[KTJ+1];
		for(int i=1;i<=node;i++) val[i]=NODE[i].material;
//		if(t==1 || t%(CON.get_EM_interval()*CON.get_mesh_output_interval())==0) 
		{
//			data_avs(node,nelm,NODE,ELEM,KTJ,val,CON);
			data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t);//断面図
			data_avs3(node,nelm,NODE,ELEM,CON);//材質
		}

		delete [] val;
	//double min_volume=10;
/*	//体積および外接球パラメータ計算
    for(int i=1;i<=nelm;i++)
	{
		//4つの節点
		int ia=ELEM[i].node[1];
		int ib=ELEM[i].node[2];
		int ic=ELEM[i].node[3];
		int ip=ELEM[i].node[4];
		//要素体積計算
		ELEM[i].volume=volume3D(NODE,ia,ib,ic,ip);
		
		//外接球中心座標および半径計算
		sphere3D(NODE,ELEM,ia,ib,ic,ip,i);
	}
	//cout<<"最小体積="<<min_volume<<endl;*/
	////////////////////////////////////////////////////////////////
	//メッシュが切れてから力の計算まで

	//節点-要素関係
    int *jnb=new int[node+1];//各節点に隣接する要素数格納
    set_jnb3D(NODE,ELEM,node,nelm,jnb);

	int **nei=new int* [node+1];//各節点の周辺要素番号格納
    for(int i=1;i<=node;i++) nei[i]=new int [jnb[i]+1];
    set_nei3D(NODE,ELEM,node,nelm,jnb,nei);

	for(int i=1;i<=node;i++)
	{
		if(jnb[i]==0)
		{
			//cout<<"jnb=0 i="<<i<<" material="<<NODE[i].material<<" particle="<<NODE[i].particleID<<endl;
			NODE[i].boundary_condition=1;//境界条件をディリクレ型にすることで、ICCGに参加させない
			//if(NODE[i].material==FLUID) NODE[i].particleID=-1;//流体節点が消失する場合、必要な値が計算されないため、粒子にフィードバックできない(してはいけない)
			if(NODE[i].material==FLUID || NODE[i].material==ELASTIC) if(NODE[i].particleID>=0) cout<<"suf="<<PART[NODE[i].particleID].surface<<endl;
		}
	}

	//FEM
	if(CON.get_EM_calc_type()==1 || CON.get_EM_calc_type()==4) potential_calculation(CON,NODE,ELEM, node, nelm,jnb, TIME,PART, fluid_number,nei,F);
	if(CON.get_EM_calc_type()==2) calc_static_magnetic_field(CON, node, nelm,NODE,ELEM,jnb, dt,TIME, t,nei, KTE,PART,fluid_number,F,KTJ);
	if(CON.get_EM_calc_type()==3) calc_transitional_EM_field(CON, node, nelm,NODE,ELEM,jnb, dt, TIME,t,nei, KTE,PART,fluid_number,F,KTJ);
	if(CON.get_EM_calc_type()==5) calc_variable_magnetic_field(CON, node, nelm,NODE,ELEM,jnb, dt,TIME, t,nei, KTE,PART,fluid_number,F);

	//FEM_intervalによる計算時間短縮の場合
    if(CON.get_EM_interval()>1 && CON.get_Ferror()==OFF) //えらーがないなら電磁力をファイルに書き込み
    {
		//ofstream bb("FEM_interval.dat");///電磁力出力　次ステップはこれを読み取る
        FILE *b=fopen("FEM_interval.dat","w");///電磁力出力　次ステップはこれを読み取る
		for(int i=0;i<fluid_number;i++)
		{
			//bb<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<endl;
			fprintf(b,"%1.15lf %1.15lf %1.15lf\n",F[A_X][i],F[A_Y][i],F[A_Z][i]);
		}
		fclose(b);
		//bb.close();
    }

	delete [] jnb;
    for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;

/*	///////////////////////////////////////////////////////////////
    //節点-要素関係
	int *jnb=new int[node+1]; //各節点に隣接する要素数格納
	set_jnb3D(NODE, ELEM, node, nelm, jnb); //delaun3Dにset_jnb3Dのバグ対策が書いてある
	int **nei=new int*[node+1]; //各節点の周辺要素番号格納
	for(int i=1;i<=node;i++) nei[i]=new int[jnb[i]+1]; 
    set_nei3D(NODE, ELEM, node, nelm, jnb, nei);
	cout<<"要素数="<<nelm<<" 節点数="<<node<<endl;
//メッシュ生成完了
	//表面メッシュ出力
	if(t==1 || t%CON.get_avs_mesh3_interval()==0)
	{
		t_data_avs3(node,nelm,NODE,ELEM,CON,t);//材質
	}
    //ここまではOK, 2012-10-30, 14:14
	//モデル毎の設定(境界条件など)
	//静電霧化
	if(CON.get_model_number()==14)
	{
	    for(int i=0;i<node;i++)
		{
			//iが1ずれていることに注意
			if(NODEall[i].boundary==ELECTRODE1)	NODE[i+1].boundary_condition=1;	//円柱電極および土台
			else if(NODEall[i].boundary==ELECTRODE2) NODE[i+1].boundary_condition=2;	//平板電極
			else NODE[i+1].boundary_condition=0;	//その他の節点は未知数
		}
		//要素-要素関係による飛散判定
	//	if(CON.get_fly_judge()==2)	fly_judge_FEMelement(CON,PART,NODE,ELEM,node,nelm,TRANS);
	}
	//磁性エラストマーの境界条件指定・・・実際の境界条件は？
	//AIRは？
	if(CON.get_model_number()==5)
	{
	    for(int i=0;i<node;i++)
		{
			//iが1ずれていることに注意			
			if(NODEall[i].boundary==MAGNET){
				NODE[i+1].boundary_condition=1;//1;	//磁石
				NODE[i+1].particleID=-1;
			}else{
				NODE[i+1].boundary_condition=0;	//その他の節点は未知数・・・AIRと
			}
		}
	}
	//set_material()＠delaun3D.cpp　これがないと要素材質定義できない？？？
//有限要素法計算開始
	//配列確保
	//磁場計算の準備
	double *V=new double[node+1];
	//vector<double> V(node+1);	//電位
	//vector<vector<double> > Eb(3,vector<double>(nelm+1));//要素内電界 or 要素内磁界	
    double *RP=new double[nelm+1];//各要素の透磁率または誘電率格納（現在では誘電率は使っていない）
	
	//FEMが失敗したときでもメッシュがみれるから、最初だけここで出力
	if(t==1){
		//data_avs(node,nelm,NODE,ELEM,KTJ,V,CON);
		//data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t);
	}
	//磁場解析
	if(CON.get_FEM_calc_type()==2 || CON.get_FEM_calc_type()==3 || CON.get_FEM_calc_type()==5)
	{
		sides* SIDE=new sides[KTE]; //辺クラス[KTE+1]ではない
		int* branch_num=new int[node+1]; //各節点が隣接する節点数(各節点から延びる辺の数)		
		int max=1000; //節点に隣接する最大節点数		
		int** nei2=new int*[node+1];
		for(int i=1;i<=node;i++) nei2[i]=new int[max];//イコールに注意
		int side_num=0;	//全辺数格納
		//辺要素生成 (節点要素を使用する場合でも、電流密度を求めるときに辺要素がほしい)
		side_num=make_edge_element(CON,NODE,ELEM,node,nelm,jnb,nei,SIDE,branch_num,nei2,KTE);
		cout<<"辺要素生成終了"<<endl;
		//磁場計算の準備
		double *B[3]; //要素内磁界
		for(int D=0;D<3;D++) B[D]=new double [nelm+1];
		double** current=new double*[DIMENSION];
		for(int D=0;D<3;D++) current[D]=new double[nelm+1];
		if(CON.get_J_input_way()==0){
			for(int D=0;D<3;D++){
				for(int i=0;i<=nelm;i++){
					current[D][i]=0.0;//強制電流が存在しないから初期化
				}
			}
		}
		for(int n=1;n<=nelm;n++)
		{
			RP[n]=1;
			if(ELEM[n].material==ELASTIC) RP[n]=CON.get_RP();
		}
		if(CON.get_ele_type()==0){
			cout<<"節点要素による磁場計算は未定義！"<<endl;
			exit(EXIT_FAILURE);
		}
		if(CON.get_ele_type()==1)//辺要素
		{
			double *A=new double[side_num+1];//ベクトルポテンシャル
		
			//ベクトルポテンシャル解決
			if(CON.get_FEM_calc_type()==2) Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,current,RP);
			if(CON.get_FEM_calc_type()==2) Avector3D_OK(CON, NODE, ELEM, SIDE, node, nelm, side_num, A, jnb, branch_num, depth);
			//else if(CON.get_FEM_calc_type()==5) non_linear_Avector3D(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,jnb,branch_num,RP,Eb);//非線形
			
			//磁束密度解決
			Bflux3D_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,A,B,RP);
	
			//電磁力解決
			if(CON.get_m_force()==0) NODE_F3D(CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,F,fluid_number);//節点力法
//			else if(CON.get_m_force()==1) VolumeForce3D(&CON,PART,NODE,ELEM,nelm,B,fluid_number,RP,jnb,nei,N,TRANS,particle_number);
//			else if(CON.get_m_force()==2) kelvin_force3D(&CON,PART,NODE,ELEM,nelm,B,fluid_number,RP,N,TRANS,jnb,nei);
//			else if(CON.get_m_force()==3) integral_surface_F3D(&CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,N,TRANS,fluid_number,depth);
//			else if(CON.get_m_force()==4) direct_divT3D(&CON,PART,NODE,ELEM,nelm,B,fluid_number,RP,jnb,nei,N,TRANS);
//			else if(CON.get_m_force()==5) virtual_air_gap_F3D(&CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,N,TRANS,fluid_number,depth);
//			else if(CON.get_m_force()==6) NODE_F3D_with_integral(&CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,N,TRANS,fluid_number);//節点力法
//			else if(CON.get_m_force()==7) Magnetic_Charge_Method_F3D(&CON,NODE,ELEM,node,nelm,B,jnb,nei,RP,PART,N,TRANS,fluid_number);
			delete [] A;
		}
		delete [] SIDE;
		delete [] branch_num;
		for(int i=0;i<=node;i++) delete [] nei2[i];
		delete [] nei2;
		for(int D=0;D<3;D++) delete [] current[D];
		delete [] current;
		for(int D=0;D<3;D++) delete [] B[D];
	}
	//電磁力スムージング
	smoothingF3D(CON,PART,fluid_number,F,t);
	//断面メッシュ出力
	if(t==1 || t%CON.get_avs_mesh2_interval()==0)
	{
		if(CON.get_FEM_calc_type()==1) data_avs(node,nelm,NODE,ELEM,KTJ,V,CON);
		//data_avs2(CON,node,nelm,NODE,ELEM,KTJ,V,t);
		data_avs2(CON,node,nelm,NODE,ELEM,KTJ,t); //物理量Vなどを指定すれば値が確認できる
		//data_avs2_movie(CON,node,nelm,NODE,ELEM,KTJ,V,t);
	}
	//表面メッシュ出力
	if(t==1 || t%CON.get_avs_mesh3_interval()==0)
	{
		if(CON.get_FEM_calc_type()==1) data_avs(node,nelm,NODE,ELEM,KTJ,V,CON);
//		data_avs3(node,nelm,NODE,ELEM,CON,t);//材質
	}
 
	delete [] V;
	delete [] RP;
    delete [] jnb;
	for(int i=1;i<=node;i++) delete [] nei[i];
    delete [] nei;*/
	delete [] depth;
	cout<<"ok2"<<endl;
}

//電磁力プロット関数
void plot_F(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number, int t)
{
	double le=CON.get_distancebp();
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	double times=CON.get_times()/CON.get_density()/le*CON.get_FEMtimes();
	double times_Pa=CON.get_times_Pa()*CON.get_FEMtimes();

	double Fz_sum=0;						//電磁力のZ方向成分の合計
	double Fall_sum=0;						//電磁力ベクトルの長さの合計
	
	ofstream fp("F.dat");

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]>-le*0.5 && PART[i].r[A_Y]<le*0.5)
		{
			if(CON.get_plot_F_type()==0) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].eforce[A_X]*times<<"\t"<<PART[i].eforce[A_Z]*times<<endl;//[N]表示
//			if(CON.get_plot_F_type()==1) fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].dirP_em*PART[i].n[A_X]*times_Pa<<"\t"<<PART[i].dirP_em*PART[i].n[A_Z]*times_Pa<<endl;//[Pa]表示
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		}
//		if(PART[i].fly==TOUCH)	Fz_sum+=PART[i].eforce[A_Z];	//Z方向成分を足していく
//		if(PART[i].type==BOFLUID)	Fall_sum+=sqrt(PART[i].eforce[A_X]*PART[i].eforce[A_X]+PART[i].eforce[A_Y]*PART[i].eforce[A_Y]+PART[i].eforce[A_Z]*PART[i].eforce[A_Z]);
	}
	
	//凡例   出力位置を保険をかけて少し斜めに移動
	xmax+=4*le;
	ymax+=4*le;
	if(CON.get_plot_F_type()==0)	if(CON.get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力 [N]表示
	if(CON.get_plot_F_type()==1)	if(CON.get_legend_F_Pa()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F_Pa()*times_Pa<<" "<<0*times_Pa<<endl;//最後に凡例出力 [Pa]表示
	fp.close();/////


	//Z方向の合力をプロット
	if(t==1)//1ステップ目はファイルをリセット
	{
		ofstream freset("Fz.dat");
		freset.close();
	}
	ofstream fout("Fz.dat",ios::app);
	fout<<t*CON.get_dt()<<setprecision(SIGNIFY)<<scientific<<"\t"<<Fz_sum<<endl;
	fout.close();

	//Z方向の合力をプロット
	if(t==1)//1ステップ目はファイルをリセット
	{
		ofstream allreset("Fall.dat");
		allreset.close();
	}
	ofstream all("Fall.dat",ios::app);
	all<<t*CON.get_dt()<<setprecision(SIGNIFY)<<scientific<<"\t"<<Fall_sum<<endl;
	all.close();
	//確認
	cout<<"静電力の大きさの総和 Fall="<<Fall_sum<<endl;

}

//電磁力プロット関数(ログ)
void plot_F_log(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number,int t)
{
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	double le=CON.get_distancebp();
	double times=CON.get_times()/CON.get_density()/le*CON.get_FEMtimes();
	double times_Pa=CON.get_times_Pa();
	
	char filename[20];
	sprintf_s(filename,"F%d.dat", t);
	ofstream fp(filename);

	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].r[A_Y]>-le*0.5 && PART[i].r[A_Y]<le*0.5)
		{
			if(CON.get_plot_F_type()==0)	fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].eforce[A_X]*times<<"\t"<<PART[i].eforce[A_Z]*times<<endl;//[N]表示
			if(CON.get_plot_F_type()==1)	fp<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].dirP_em*PART[i].n[A_X]*times_Pa<<"\t"<<PART[i].dirP_em*PART[i].n[A_Z]*times_Pa<<endl;//[Pa]表示
			if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
			if(PART[i].r[A_Z]>ymax) ymax=PART[i].r[A_Z];
		}
    }
	
	//凡例   出力位置を保険をかけて少し斜めに移動
	xmax+=4*le;
	ymax+=4*le;
	if(CON.get_plot_F_type()==0)	if(CON.get_legend_F()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F()*times<<" "<<0*times<<endl;//最後に凡例出力 [N]表示
	if(CON.get_plot_F_type()==1)	if(CON.get_legend_F_Pa()>0) fp<<xmax<<" "<<ymax<<" "<<CON.get_legend_F_Pa()*times_Pa<<" "<<0*times_Pa<<endl;//最後に凡例出力 [Pa]表示
	fp.close();/////
}


//OpenMPによるICCG法
void parallel_ICCG3D2(mpsconfig &CON,double *val,int *ind,int *ptr,int pn,double *B,int number,double *X)
{
	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解

	double accel=CON.get_CGaccl();//加速ファクタ
	
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
    double *r=new double[pn];
	for(int n=0;n<pn;n++) X[n]=0;
	
	double *AP = new double [pn];
	double *P = new double [pn];
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

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0;
		r[n]=B[n];
	}

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

    cout<<"parallel_ICCG法スタート  -----未知数="<<pn<<"  ---"<<endl;
	ofstream Eee("PICCG.dat", ios::trunc);
	Eee.close();
	unsigned int time=GetTickCount();
	int count=0;
	double lowE=100;
	double ep=CON.get_FEMCGep();//収束判定
	while(E>ep)
	{
//		if(E<=lowE) lowE=E;
		count++;
		if(count==pn){
			cout<<"count>pn E="<<E<<"lowE="<<lowE<<endl;
//			ep=lowE;
		}
		//////////////alpを求める
		rLDLt_r=0;
		double PAP=0;
		//for(int n=0;n<pn;n++) rLDLt_r+=r[n]*LDLt_r[n];//sirial
		#pragma omp parallel for reduction(+:rLDLt_r) reduction(+:PAP)
		for(int n=0;n<pn;n++)
		{	//printf("%d\n",omp_get_thread_num());
			rLDLt_r+=r[n]*LDLt_r[n];
			AP[n]=0;
			for(int m=ptr[n];m<ptr[n+1];m++) AP[n]+=val[m]*P[ind[m]];
			PAP+=P[n]*AP[n];
		}
		
		//for(int n=0;n<pn;n++)  PAP+=P[n]*AP[n];//sirial
		alp=rLDLt_r/PAP;
		
//		cout<<"alp="<<alp<<endl;
		//////////////////////
	
		//////////////// X(k+1)=X(k)+alp*P
		E=0;
		#pragma omp parallel for reduction(+:E)
		for(int n=0;n<pn;n++) 
		{
			X[n]+=alp*P[n];	// X(k+1)=X(k)+alp*P
			r[n]-=alp*AP[n];// r=r-alp*AP
			E+=r[n]*r[n];	//誤差
		}
		E=sqrt(E);
		//////////////////////////////
		
		//////////////// r=r-alp*AP
		//for(int n=0;n<pn;n++) r[n]-=alp*AP[n];
		/////////////////////////////
		
		/*/////////////////誤差
		E=0;
		for(int n=0;n<pn;n++) E+=r[n]*r[n];		
		E=sqrt(E);
		//cout<<"E="<<E<<endl;
		printf("残差 E=%lf          \r",E);
		///////////////////////*/
		ofstream Ee("PICCG.dat", ios::app);
		Ee<<count<<" "<<E<<endl;
		Ee.close();
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
	cout<<"反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	
    delete [] r;
	delete [] AP;
	delete [] P;

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
}

///電流密度計算関数(辺要素用)
void calc_current(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,int *jnb,int *branch_num,double **current,int *depth,double II)
{
	//II:電流[A]

	cout<<"電流密度計算開始"<<endl;
	
	double p=1.68e-8;//銅の電気抵抗率[Ωm]

	int side_num2=0;//コイルを構成する辺数
	for(int i=1;i<=side_num;i++)
	{
		int ia=SIDE[i].node[1];
		int ib=SIDE[i].node[2];
		if(NODE[ia].material==COIL && NODE[ib].material==COIL) side_num2++;
	}///side_num2がもとまった

	int *side_id=new int [side_num2+1];//コイルを構成する辺番号格納
	side_num2=0;
	for(int i=1;i<=side_num;i++)
	{
		int ia=SIDE[i].node[1];
		int ib=SIDE[i].node[2];
		if(NODE[ia].material==COIL && NODE[ib].material==COIL)
		{
			side_num2++;
			side_id[side_num2]=i;//辺番号iを格納
		}
	}///コイルを構成する辺番号を管理
	//////////////////////////ｺｲﾙ辺の電流に関する境界条件を設定
	for(int i=1;i<=nelm;i++)
	{
		if(depth[i]==1)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm!=0)
				{
					if(ELEM[kelm].material==COIL) 
					{
						///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
						int p=ELEM[i].node[j];//境界三角に属さない節点
						for(int k=1;k<=6;k++)
						{
							
							int iside=ELEM[i].sides[k];
							//cout<<"iside="<<iside<<endl;
							int ia=SIDE[iside].node[1];
							int ib=SIDE[iside].node[2];
							//cout<<"ia="<<ia<<"ib="<<ib<<endl;
							if(ib<ia)
							{
								int temp=ia;
								ia=ib;
								ib=temp;
							}///これで必ずia<ibとなった
							
							if(ia!=p && ib!=p)//辺isideはコイル境界上ということ
							{
								if(NODE[ia].boundary_condition==11 || NODE[ib].boundary_condition==11)
								{
									SIDE[iside].boundary_condition=11;//11をひとつでも含んだ辺は自然境界条件
									//cout<<(NODE[ia].r[A_X]+NODE[ib].r[A_X])/2<<" "<<(NODE[ia].r[A_Y]+NODE[ib].r[A_Y])/2<<" "<<(NODE[ia].r[A_Z]+NODE[ib].r[A_Z])/2<<endl;
								}
								else
								{
									SIDE[iside].boundary_condition=10;//T=0となる辺
								}
								//T=0の面および自然境界条件面に固定境界辺を設置する
								if((NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==22) ||(NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==23)||(NODE[ia].boundary_condition==23 && NODE[ib].boundary_condition==21)) SIDE[iside].boundary_condition=21;//辺は21→22の方向
								else if((NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==21) ||(NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==23)||(NODE[ia].boundary_condition==23 && NODE[ib].boundary_condition==22)) SIDE[iside].boundary_condition=22;//辺は22→21の方向
								
							}
						}
					}
				}
			}
		}
		if(depth[i]==5)//物体表面に隣接する空気要素
		{
			for(int j=1;j<=4;j++)
			{
				int kelm=ELEM[i].elm[j];
				if(kelm!=0)
				{
					if(ELEM[kelm].material==COIL) 
					{
						///第j面に隣接する三角形と向かい合っている頂点は、第j番節点である
						int p=ELEM[i].node[j];//境界三角に属さない節点
						for(int k=1;k<=6;k++)
						{
							
							int iside=ELEM[i].sides[k];
							//cout<<"iside="<<iside<<endl;
							int ia=SIDE[iside].node[1];
							int ib=SIDE[iside].node[2];
							//cout<<"ia="<<ia<<"ib="<<ib<<endl;
							if(ib<ia)
							{
								int temp=ia;
								ia=ib;
								ib=temp;
							}///これで必ずia<ibとなった
							
							if(ia!=p && ib!=p)//辺isideはコイル境界上ということ
							{
								if(NODE[ia].boundary_condition==11 || NODE[ib].boundary_condition==11)
								{
									SIDE[iside].boundary_condition=11;//11をひとつでも含んだ辺は自然境界条件
									//cout<<(NODE[ia].r[A_X]+NODE[ib].r[A_X])/2<<" "<<(NODE[ia].r[A_Y]+NODE[ib].r[A_Y])/2<<" "<<(NODE[ia].r[A_Z]+NODE[ib].r[A_Z])/2<<endl;
								}
								else
								{
									SIDE[iside].boundary_condition=10;//T=0となる辺
								}
								//T=0の面および自然境界条件面に固定境界辺を設置する
								if((NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==22) ||(NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==23)||(NODE[ia].boundary_condition==23 && NODE[ib].boundary_condition==21)) SIDE[iside].boundary_condition=23;//辺は21→22の方向
								else if((NODE[ia].boundary_condition==22 && NODE[ib].boundary_condition==21) ||(NODE[ia].boundary_condition==21 && NODE[ib].boundary_condition==23)||(NODE[ia].boundary_condition==23 && NODE[ib].boundary_condition==22)) SIDE[iside].boundary_condition=24;//辺は22→21の方向
								
							}
						}
					}
				}
			}
		}


	}//コイル端面に自然境界、側面にT=0の固定境界を敷いた。
	////////////////*/

	///境界条件出力　うまくいかないときにみる
	//data_avs_J_boundary(node,nelm,NODE,ELEM);

	for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
		//自然境界条件を未知数あつかいする
        if(SIDE[i].boundary_condition==11) SIDE[i].boundary_condition=0;
	}
	
	double *T=new double [side_num+1];//電流ベクトルﾎﾟﾃﾝｼｬﾙ
    int NN=0;//ディリクレ型境界辺数
    int *dn=new int [side_num+1]; //各辺がディリクレ型境界の何番目か。ディリクレ型境界上でないならside_num2+1を格納
    double *PHAT=new double [CON.get_max_DN()];//ディリクレ型境値
    double J1=II/0.0000896;//89.6
	for(int k=1;k<=side_num;k++)T[k]=0;
    ///ディリクレ型境界条件入力
    for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
        if(SIDE[i].boundary_condition==10)
		{
			dn[i]=NN;//i番目の辺はNN番目のディリクレ境界辺
	        PHAT[NN]=0;
	        T[i]=0;
	        NN++;
		}
		else if(SIDE[i].boundary_condition==21)
		{    
			int ia=SIDE[i].node[1];
			int ib=SIDE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//辺の長さ
	        dn[i]=NN;
	        PHAT[NN]=II;
	        T[i]=II;
	        NN++;
		}
		else if(SIDE[i].boundary_condition==22)
		{   
			int ia=SIDE[i].node[1];
			int ib=SIDE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//辺の長さ
	        dn[i]=NN;
	        PHAT[NN]=-II;
	        T[i]=-II;
	        NN++;
		}
		else if(SIDE[i].boundary_condition==23)
		{    
			int ia=SIDE[i].node[1];
			int ib=SIDE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//辺の長さ
	        dn[i]=NN;
	        PHAT[NN]=II;
	        T[i]=II;
	        NN++;
		}
		else if(SIDE[i].boundary_condition==24)
		{   
			int ia=SIDE[i].node[1];
			int ib=SIDE[i].node[2];
			double X=NODE[ia].r[A_X]-NODE[ib].r[A_X];
			double Y=NODE[ia].r[A_Y]-NODE[ib].r[A_Y];
			double Z=NODE[ia].r[A_Z]-NODE[ib].r[A_Z];
			double L=sqrt(X*X+Y*Y+Z*Z);//辺の長さ
	        dn[i]=NN;
	        PHAT[NN]=-II;
	        T[i]=-II;
	        NN++;
		}
		else dn[i]=side_num2+1;
    }
	cout<<"ﾃﾞｨﾘｸﾚ数＝"<<NN<<endl;
	/////////////*/

	
    //////int pn=side_num-NN;				///未知数
	int pn=side_num2-NN;				///未知数
    int *ppn=new int [pn];			//行列のn番目は辺ppn[n]
    int *npp=new int [side_num+1];	///各辺が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int k=1;k<=side_num2;k++)
    {
		int i=side_id[k];
        if(SIDE[i].boundary_condition==0)//未知数 
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
		}
		else npp[i]=pn+1;
    }
    
    ////行列の最大幅計算
    int mat_w=0;
	for(int k=1;k<=side_num2;k++)
	{
		int i=side_id[k];
		int ia=SIDE[i].node[1];
		int ib=SIDE[i].node[2];
		int width=branch_num[ia]+branch_num[ib];//行列の幅
		if(width>mat_w) mat_w=width;
	}
	//mat_w*=5;
	////////////
	
    ////配列確保
    double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [mat_w+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [mat_w+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    
    for(int i=1;i<=pn;i++)//初期化
    {
        NUM[i]=0;
        for(int j=1;j<=mat_w;j++)
		{
			G[i][j]=0;
			ROW[i][j]=0;
		}
    }
    double *B=new double [pn];//解行列
    for(int i=0;i<pn;i++) B[i]=0;//初期化
    ////
    
    /////////全体行列を作成する
    
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係
	
    for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==COIL)
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].sides[i];
				int ia=SIDE[iside].node[1];
				int ib=SIDE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
		
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
			
			double Xs=0;//要素の重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j];
				Ys+=Y[j];
				Zs+=Z[j];
			}
			Xs/=4;Ys/=4;Zs/=4;
			////////////////////////////
	
			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。
	        ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
		
			double delta6=ELEM[je].volume;//体積の6倍
		
			delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
		
			double delta=ELEM[je].volume/6;//本当の体積
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
		    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
		
			////要素ﾏﾄﾘｸｽ作成開始
			for(int i=1;i<=6;i++)
			{	
				int iside=ELEM[je].sides[i];//要素jeの辺番号
				if(SIDE[iside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
				{   
					int I=npp[iside]+1;///辺isideは行列のI番目
					int I1=SIDE[iside].node[1];//isideを構成する2点
					int I2=SIDE[iside].node[2];
					int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int k2=table[i][2];
				    for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
						
						int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[j][2];
							
						if(SIDE[jside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
						{   
							int J=npp[jside]+1;///辺jsideは行列のJ番目
							int flag=0;
							//if(J<=I){
							int J1=SIDE[jside].node[1];//jsideを構成する2点
							int J2=SIDE[jside].node[2];
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6;
								    flag=1;
								}
							}
							if(flag==0)
							{   
							    NUM[I]=NUM[I]+1;
							    int H=NUM[I];
							    
								G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6;
							    ROW[I][H]=J;
							}
							//}
						}
						////
						else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら
						{
						    int n=dn[jside];
						    B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*p*delta6*delta6*delta6*PHAT[n];
						}//////////*/
					}
				}
			}
		}
    }
    ///////////////////////*/
	

    int number=0;//行列の非ゼロ要素数
    for(int i=1;i<=pn;i++) number+=NUM[i];
   
    ///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
    for(int i=1;i<=pn;i++)
    {
        double tempG;
		int tempR;
		for(int j=2;j<=NUM[i];j++)
		{
		    for(int m=1;m<j;m++)
		    {
		        if(ROW[i][j]<ROW[i][m])
				{
				    tempG=G[i][m];
				    tempR=ROW[i][m];
					G[i][m]=G[i][j];
					ROW[i][m]=ROW[i][j];
					G[i][j]=tempG;
					ROW[i][j]=tempR;
				}
			}
		}
    }///////////

	///対称性チェック
	for(int i=1;i<=pn;i++)
	{
		for(int j=1;j<=NUM[i];j++)
		{
			int J=ROW[i][j];
			for(int k=1;k<=NUM[J];k++) if(ROW[J][k]==i) if(G[i][j]!=G[J][k])
			{
				cout<<"matrix isn't symmetric   "<<G[i][j]<<" "<<G[J][k]<<endl;
			}
		}
	}///////////*/

    double *val = new double [number];
    int *ind = new int [number];//非ゼロ要素の列番号格納配列
    int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
    
    /////////////////////val,ind ,ptrに値を格納
    int index=0;
    for(int n=0;n<pn;n++)
    {
        ptr[n]=index;
		for(int m=1;m<=NUM[n+1];m++)
		{
			val[index]=G[n+1][m];
			ind[index]=ROW[n+1][m]-1;
			index++;
		}
    }	    
    ptr[pn]=number;
    ////////////////////*/
    
    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
    
    cout<<"行列作成終了  "<<endl;
    
    double *XX=new double [pn];//行列の答え格納
	//CG3D(val,ind,ptr,pn,ppn,B,T);//CG法実行
	//ICCG3D(val,ind,ptr,pn,ppn,B,T,number);//ICCG法実行
	ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
	//MRTR(CON,B,pn,XX,val,ind,ptr);
    ///////////
	for(int n=0;n<pn;n++)
			{
				int i=ppn[n];
				T[i]=XX[n];
			}
			delete [] XX;
	denryu_side(CON,NODE,ELEM,SIDE,node,nelm,side_num,T,current,J1);
    
	delete [] side_id;
	delete [] T;
    ///////////////////////*/
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
    
    delete [] val;
    delete [] ind;
    delete [] ptr;
}
////辺要素電流密度計算関数
void denryu_side(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *T,double **current,double J1)
{
	int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];//要素を構成する辺とその要素内節点番号関係

	ofstream fp("j.dat");
	for(int je=1;je<=nelm;je++)
    {   
		if(ELEM[je].material==COIL)
		{
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].sides[i];
				int ia=SIDE[iside].node[1];
				int ib=SIDE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double Xs=0;//要素の重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j];
				Ys+=Y[j];
				Zs+=Z[j];
			}
			Xs/=4;
			Ys/=4;
			Zs/=4;
			////////////////////////////
	
			double delta6=ELEM[je].volume;//体積の6倍(正しい体積の値はすでにベクトルポテンシャルを求める際に計算している)
		
			delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
		
			double delta=ELEM[je].volume/6;//本当の体積
		
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
		    
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
				    c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			for(int D=0;D<3;D++) current[D][je]=0;//初期化
	
			for(int i=1;i<=6;i++)
			{
				int s=ELEM[je].sides[i];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];
	
				current[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*T[s];
				current[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*T[s];
				current[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*T[s];
				
			}
	
			for(int D=0;D<3;D++) current[D][je]*=delta6*delta6*2;
			///////////////////////////////////////////////////////////
			double Along=sqrt(pow(current[A_X][je],2)+pow(current[A_Y][je],2)+pow(current[A_Z][je],2));
			current[A_X][je]=current[A_X][je]*J1/Along;
			current[A_Y][je]=current[A_Y][je]*J1/Along;
			current[A_Z][je]=current[A_Z][je]*J1/Along;

			//////////////////////////////////////////////////////////*/
			//if(Zs>0 && Zs<0.0001) fp<<Xs<<" "<<Ys<<" "<<current[A_X][je]*1e-12<<" "<<current[A_Y][je]*1e-12<<endl;
			if(Zs>0.0005 && Zs<0.001)
			fp<<Xs<<" "<<Ys<<" "<<current[A_X][je]*1e-12<<" "<<current[A_Y][je]*1e-12<<endl;		
		}
		
	}	
	fp.close();
	
}

//要素の深さ決定関数
void set_depth(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *depth,int KTE)
{
	//for(int i=1;i<=KTE;i++) depth[i]=0;
	for(int i=1;i<=nelm;i++) depth[i]=0;

	for(int i=1;i<=nelm;i++)//物質節点をひとつでも含む空気要素は境界要素(深さ1)と定義する
	{
		if(ELEM[i].material==AIR)//要素iが空気
		{
			for(int j=1;j<=4;j++)
			{
				int ia=ELEM[i].node[j];
				if(NODE[ia].material!=AIR)
				{
					depth[i]=1;//物質要素に接する空気要素のdepthを１と定義
				}
			}
		}
		else if(ELEM[i].material==IRON)	//コイル用
		{
			{
			for(int j=1;j<=4;j++)
			{
				int ia=ELEM[i].node[j];
				if(NODE[ia].material==COIL)
				{
					depth[i]=5;//物質要素に接する空気要素のdepthを１と定義
				}
			}
		}
		}
	}///深さ１の要素がもとまった(ﾏｸｽｳｪﾙ積分面は深さ1の要素が深さ2と接する面とすればよい)
	
	int count=10;
	int A=1;//現在注目する深さ
	while(count>0)
	{
		count=0;//新しく深さがもとまった要素数
		for(int i=1;i<=nelm;i++)
		{
			if(ELEM[i].material==AIR)
			{
				if(depth[i]==A)
				{
					for(int j=1;j<=4;j++)
					{
						int jelm=ELEM[i].elm[j];//iの隣接する要素
						if(jelm!=0)
						{
							if(ELEM[jelm].material==AIR && depth[jelm]==0)
							{
								depth[jelm]=A+1;
								count++;
							}
						}
					}
				}
			}
		}
		A++;
	}////２以降の深さに関しては、ここでの定義とのちに自動分割したときの定義では少し異なるので注意。
	///check
	for(int i=1;i<=nelm;i++) if(ELEM[i].material==AIR) if(depth[i]==0) cout<<"要素深さerror"<<endl;
}

void carrent_vector(mpsconfig &CON, vector<point3D> &NODE, vector<element3D> &ELEM, int nelm, double **B,int t)
{
	double **center_of_element=new double* [nelm+1];
	for(int i=1;i<=nelm;i++) center_of_element[i]=new double[3];
	//重心の計算・・・nelm全て出力
	for(int i=1;i<=nelm;i++)
	{
		for(int D=0;D<3;D++)
		{
			center_of_element[i][D]=0.0;
			for(int j=1;j<=4;j++) center_of_element[i][D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		}
	}

	int timestep=CON.get_current_step();
	int half_nelm=0;
	for(int i=1;i<=nelm;i++){				//Y-Z面に表示する断面までの要素数
		if(center_of_element[i][A_X]<0.0) half_nelm++;
	}
//	if(t==0 || t%10==0){
		stringstream sscv;
		sscv<<"./Current/Cvector"<<t<<".fld";
		string Cvector=sscv.str();	
		ofstream fcv(Cvector);
		if(fcv.fail()){
			system("mkdir Current");
		ofstream fcv(Cvector);
		if(fcv.fail()){
			cout<<"ファイルを開けませんでした。Currentフォルダがあるかどうか確認してください"<<endl;
			exit(1);
			}
		}

		fcv << "# AVS field file" << endl;
		fcv << "ndim=1" << endl;
		fcv << "dim1=" << half_nelm <<endl;
		fcv << "nspace=3" << endl;
		fcv << "veclen=3" << endl;
		fcv << "data=float" << endl;
		fcv << "field=irregular" << endl;
		fcv << "label=e-x e-y e-z" << endl << endl;
		fcv << "variable 1 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
		fcv << "variable 2 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
		fcv << "variable 3 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
		fcv << "coord    1 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
		fcv << "coord    2 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
		fcv << "coord    3 file=./Cvector"<<t<<" filetype=ascii skip=1 offset=5 stride=6" << endl;

		fcv.close();

		stringstream sscvd;
		sscvd<<"./Current/Cvector"<<t;
		string Cvectordata=sscvd.str();

		ofstream fout(Cvectordata);
		if(fout.fail()){
			cout<<"データファイルを開けませんでした。"<<endl;
			exit(1);
		}

		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=1;i<=nelm;i++)
		{
			if(center_of_element[i][A_X]<0){
				fout<<B[A_X][i]/1000000<<" "<<B[A_Y][i]/1000000<<" "<<B[A_Z][i]/1000000<<" "<<center_of_element[i][A_X]<<" "<<center_of_element[i][A_Y]<<" "<<center_of_element[i][A_Z]<<endl;
			}
		}
		fout.close();
//	}//*/
		for(int i=1;i<=nelm;i++) delete [] center_of_element[i];
		delete [] center_of_element;
}
///非線形ﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算関数(辺要素用)(村尾)
void non_linear_Avector3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *V,int *jnb,int *branch_num,double **current,double *RP,double II,int *depth,double **Be,int t)
{
	//////////////////////////////////電流密度ベクトル計算///////////////////////////////////////////////
	//////////////////////ｺｲﾙがなければ電流計算を行う必要がない。そこでｺｲﾙ節点数をかぞえる
	int coil_node_num=0;	//ｺｲﾙ節点数
	for(int i=1;i<=node;i++) if(NODE[i].material==COIL) coil_node_num++;
	////////////////////////*/

	int *save_bound=new int [node+1];//本関数は何度も呼び出されるので、電流に関する境界条件を完全に消去できない。そこで保存する
	
	if(coil_node_num>0)///ｺｲﾙ節点があるなら電流密度計算
	{
		for(int i=1;i<=node;i++) save_bound[i]=NODE[i].boundary_condition;//境界条件を保存

		if(II!=0)//電流値が非ゼロなら電流計算
		{
			calc_current(CON,NODE,ELEM,SIDE,node,nelm,side_num,jnb,branch_num,current,depth,II);
			cout<<"電流計算完了"<<endl;
		}
		else	//ゼロなら初期化
		{
			for(int i=1;i<=nelm;i++) if(ELEM[i].material==COIL) for(int D=0;D<3;D++) current[D][i]=0;
		}
		///電流境界条件の初期化（もう電流はもとまっているからいらない)
		for(int i=1;i<=side_num;i++) if(SIDE[i].boundary_condition>=10) SIDE[i].boundary_condition=0;
		for(int i=1;i<=node;i++) if(NODE[i].boundary_condition>=10) NODE[i].boundary_condition=0;
	}
	carrent_vector(CON, NODE, ELEM, nelm, current,t);
	///////////////////////////////////////////////////////////////////////////////////////////////////


	//ここでは、V[i]は各辺のﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙの微小変化δAを表すことに注意。
	cout<<"非線形ﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ計算開始"<<endl;
	double u0=4*PI*0.0000001;	//真空の透磁率
    double v0=1/u0;///磁気抵抗率
	double j0x, j0y, j0z;//電流密度[A/m^3]
	double MA=CON.get_magnet_angle();
	double magnet_direction[3]={-sin(MA*2*PI/360),0,cos(MA*2*PI/360)};
	double Mx=CON.get_magnet_B()*magnet_direction[A_X];
	double My=CON.get_magnet_B()*magnet_direction[A_Y];
	double Mz=CON.get_magnet_B()*magnet_direction[A_Z];
	int graphtype=2;

	//周囲壁に固定境界条件を設定
	for(int i=1;i<=nelm; i++){
		if(ELEM[i].material==AIR){//物体表面に隣接する空気要素
			for(int j=1; j<=4; j++){
				int kelm=ELEM[i].elm[j];
				if(kelm==0){
					//第j面に隣接する三角形と向かい合っている頂点は，第j番接点である
					int p=ELEM[i].node[j];//境界三角形に属さない接点
					for(int k=1;k<=6;k++){
						int iside=ELEM[i].sides[k];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						if(ia!=p && ib!=p) SIDE[iside].boundary_condition=1;//接点pを含まない辺は境界辺
						else SIDE[iside].boundary_condition=0;
					}
				}
			}
		}
	}

	//////////////////////////////////////アキマ補完によるB-H曲線作成/////////////////////////////////
	////////////////////非線形情報入力
	double H[14];//Nn+4
	double M[14];
	double Bflux[14];
	int Nn=10;		//データ数
	if(graphtype==0)//サイト
	{
		H[2]=0.0;
		H[3]=17391;
		H[4]=43478;
		H[5]=78260;
		H[6]=173913;
		H[7]=391304;
		M[2]=0.0;
		M[3]=0.01;
		M[4]=0.015;
		M[5]=0.02;
		M[6]=0.025;
		M[7]=0.03;
		for(int n=2;n<Nn+2;n++) Bflux[n]=u0*H[n]+M[n];
	}
	else if(graphtype==1)//ランジュバン
	{
		double Ms=14700;			//飽和磁化[A/m]
		double kai0=1.172;			//磁気感受率
		double gam=3*kai0/Ms;
		H[2]=0.0;
		H[3]=2000;
		H[4]=4000;
		H[5]=6000;
		H[6]=8000;
		H[7]=12000;
		M[2]=0.0;
		for(int k=3;k<=7;k++) M[k]=Ms*(1/tanh(gam*H[k])-1/(gam*H[k]));
		for(int n=2;n<Nn+2;n++) Bflux[n]=u0*H[n]+u0*M[n];
	}
	else if(graphtype==2)
	{
		H[2]=29.083;
		H[3]=26119.17;
		H[4]=88160.77;
		H[5]=212216.5;
		H[6]=428597.7;
		H[7]=800722.5;
		H[8]=1020144;
		H[9]=1212705;
		H[10]=1483964;
		H[11]=1556321;
		Bflux[2]=0.0;
		Bflux[3]=0.102173;
		Bflux[4]=0.389118;
		Bflux[5]=0.859961;
		Bflux[6]=1.40583;
		Bflux[7]=1.930318;
		Bflux[8]=2.268766;
		Bflux[9]=2.523716;
		Bflux[10]=2.881576;
		Bflux[11]=2.976717;
	}

	///両端に2点追加
	Bflux[0]=Bflux[2]-(Bflux[4]-Bflux[2]);
	Bflux[1]=Bflux[2]-(Bflux[3]-Bflux[2]);
	H[0]=H[2]*(Bflux[0]-Bflux[3])*(Bflux[0]-Bflux[4])/((Bflux[2]-Bflux[3])*(Bflux[2]-Bflux[4]))+H[3]*(Bflux[0]-Bflux[2])*(Bflux[0]-Bflux[4])/((Bflux[3]-Bflux[2])*(Bflux[3]-Bflux[4]))+H[4]*(Bflux[0]-Bflux[2])*(Bflux[0]-Bflux[3])/((Bflux[4]-Bflux[2])*(Bflux[4]-Bflux[3]));
	H[1]=H[2]*(Bflux[1]-Bflux[3])*(Bflux[1]-Bflux[4])/((Bflux[2]-Bflux[3])*(Bflux[2]-Bflux[4]))+H[3]*(Bflux[1]-Bflux[2])*(Bflux[1]-Bflux[4])/((Bflux[3]-Bflux[2])*(Bflux[3]-Bflux[4]))+H[4]*(Bflux[1]-Bflux[2])*(Bflux[1]-Bflux[3])/((Bflux[4]-Bflux[2])*(Bflux[4]-Bflux[3]));

	Bflux[Nn+2]=Bflux[Nn+1]+(Bflux[Nn+1]-Bflux[Nn]);
	Bflux[Nn+3]=Bflux[Nn+1]+(Bflux[Nn+1]-Bflux[Nn-1]);
	H[Nn+2]=H[Nn-1]*(Bflux[Nn+2]-Bflux[Nn])*(Bflux[Nn+2]-Bflux[Nn+1])/((Bflux[Nn-1]-Bflux[Nn])*(Bflux[Nn-1]-Bflux[Nn+1]))+H[Nn]*(Bflux[Nn+2]-Bflux[Nn-1])*(Bflux[Nn+2]-Bflux[Nn+1])/((Bflux[Nn]-Bflux[Nn-1])*(Bflux[Nn]-Bflux[Nn+1]))+H[Nn+1]*(Bflux[Nn+2]-Bflux[Nn-1])*(Bflux[Nn+2]-Bflux[Nn])/((Bflux[Nn+1]-Bflux[Nn-1])*(Bflux[Nn+1]-Bflux[Nn]));
	H[Nn+3]=H[Nn-1]*(Bflux[Nn+3]-Bflux[Nn])*(Bflux[Nn+3]-Bflux[Nn+1])/((Bflux[Nn-1]-Bflux[Nn])*(Bflux[Nn-1]-Bflux[Nn+1]))+H[Nn]*(Bflux[Nn+3]-Bflux[Nn-1])*(Bflux[Nn+3]-Bflux[Nn+1])/((Bflux[Nn]-Bflux[Nn-1])*(Bflux[Nn]-Bflux[Nn+1]))+H[Nn+1]*(Bflux[Nn+3]-Bflux[Nn-1])*(Bflux[Nn+3]-Bflux[Nn])/((Bflux[Nn+1]-Bflux[Nn-1])*(Bflux[Nn+1]-Bflux[Nn]));
	////////////

	/*/////////////////対称性を考えて補正
	H[0]=-H[4];
	H[1]=-H[3];
	Bflux[0]=-Bflux[4];
	Bflux[1]=-Bflux[3];
	//////////////////////*/

	/////ユーザーの入力した点(B-H曲線をプロット
	ofstream fp3("B-H.dat");
	for(int i=0;i<Nn+4;i++) fp3<<H[i]<<" "<<Bflux[i]<<endl;
	fp3.close();
	//////////////////*/

	///各区間のH-B曲線の傾きを計算
	double *mm=new double [Nn+3];//各区間の傾き。
	for(int i=0;i<Nn+3;i++) mm[i]=(H[i+1]-H[i])/(Bflux[i+1]-Bflux[i]);

	///ユーザーの入力点におけるH-B曲線の傾き,つまり磁気抵抗率計算
	double *vm=new double[Nn+4];//傾き
	for(int i=2;i<Nn+2;i++)
	{
		double a1=sqrt((mm[i+1]-mm[i])*(mm[i+1]-mm[i]));
		double a2=sqrt((mm[i-1]-mm[i-2])*(mm[i-1]-mm[i-2]));
		double A0=a1+a2;

		if(A0<0.0001) vm[i]=0.5*(mm[i-1]+mm[i]);
		else vm[i]=(a1*mm[i-1]+a2*mm[i])/A0;
		
	}///////////

	ofstream fout("H-B_akima.dat");//アキマ補間によるH-B曲線ﾌﾟﾛｯﾄ
	for(double i=Bflux[2];i<=Bflux[Nn+1];i+=0.001)
	{
		for(int n=2;n<=Nn+1;n++)
		{
			if(i>=Bflux[n] && i<Bflux[n+1])
			{
				double h=Bflux[n+1]-Bflux[n];
				double aj0=H[n];
				double aj1=vm[n];
				double aj2=(3*mm[n]-2*vm[n]-vm[n+1])/h;
				double aj3=(vm[n]+vm[n+1]-2*mm[n])/(h*h);
				double val=aj0+aj1*(i-Bflux[n])+aj2*(i-Bflux[n])*(i-Bflux[n])+aj3*(i-Bflux[n])*(i-Bflux[n])*(i-Bflux[n]);
				double val2=aj1+2*aj2*(i-Bflux[n])+3*aj3*(i-Bflux[n])*(i-Bflux[n]);
				fout<<i<<" "<<val<<endl;
			}
		}
	}
	fout.close();
	///////////////////////////*/

	//両端2点のvmを求める
	vm[0]=vm[2]*(Bflux[0]-Bflux[3])*(Bflux[0]-Bflux[4])/((Bflux[2]-Bflux[3])*(Bflux[2]-Bflux[4]))+vm[3]*(Bflux[0]-Bflux[2])*(Bflux[0]-Bflux[4])/((Bflux[3]-Bflux[2])*(Bflux[3]-Bflux[4]))+vm[4]*(Bflux[0]-Bflux[2])*(Bflux[0]-Bflux[3])/((Bflux[4]-Bflux[2])*(Bflux[4]-Bflux[3]));
	vm[1]=vm[2]*(Bflux[1]-Bflux[3])*(Bflux[1]-Bflux[4])/((Bflux[2]-Bflux[3])*(Bflux[2]-Bflux[4]))+vm[3]*(Bflux[1]-Bflux[2])*(Bflux[1]-Bflux[4])/((Bflux[3]-Bflux[2])*(Bflux[3]-Bflux[4]))+vm[4]*(Bflux[1]-Bflux[2])*(Bflux[1]-Bflux[3])/((Bflux[4]-Bflux[2])*(Bflux[4]-Bflux[3]));
	vm[Nn+2]=vm[Nn-1]*(Bflux[Nn+2]-Bflux[Nn])*(Bflux[Nn+2]-Bflux[Nn+1])/((Bflux[Nn-1]-Bflux[Nn])*(Bflux[Nn-1]-Bflux[Nn+1]))+vm[Nn]*(Bflux[Nn+2]-Bflux[Nn-1])*(Bflux[Nn+2]-Bflux[Nn+1])/((Bflux[Nn]-Bflux[Nn-1])*(Bflux[Nn]-Bflux[Nn+1]))+vm[Nn+1]*(Bflux[Nn+2]-Bflux[Nn-1])*(Bflux[Nn+2]-Bflux[Nn])/((Bflux[Nn+1]-Bflux[Nn-1])*(Bflux[Nn+1]-Bflux[Nn]));
	vm[Nn+3]=vm[Nn-1]*(Bflux[Nn+3]-Bflux[Nn])*(Bflux[Nn+3]-Bflux[Nn+1])/((Bflux[Nn-1]-Bflux[Nn])*(Bflux[Nn-1]-Bflux[Nn+1]))+vm[Nn]*(Bflux[Nn+3]-Bflux[Nn-1])*(Bflux[Nn+3]-Bflux[Nn+1])/((Bflux[Nn]-Bflux[Nn-1])*(Bflux[Nn]-Bflux[Nn+1]))+vm[Nn+1]*(Bflux[Nn+3]-Bflux[Nn-1])*(Bflux[Nn+3]-Bflux[Nn])/((Bflux[Nn+1]-Bflux[Nn-1])*(Bflux[Nn+1]-Bflux[Nn]));
	
	//ユーザーの入力点におけるv-B曲線ﾌﾟﾛｯﾄ
	ofstream fp("v-B.dat");
	for(int i=0;i<Nn+4;i++) fp<<Bflux[i]<<" "<<vm[i]<<endl;
	fp.close();
	/////////////*/

	///各区間のv-B曲線の傾きを計算
	double *mm2=new double [Nn+3];//各区間の傾き。
	for(int i=0;i<Nn+3;i++) mm2[i]=(vm[i+1]-vm[i])/(Bflux[i+1]-Bflux[i]);

	///ユーザーの入力点におけるv-B曲線の傾き計算
	double *t2=new double[Nn+4];//傾き
	for(int i=2;i<Nn+2;i++)
	{
		double a1=sqrt((mm2[i+1]-mm2[i])*(mm2[i+1]-mm2[i]));
		double a2=sqrt((mm2[i-1]-mm2[i-2])*(mm2[i-1]-mm2[i-2]));
		double A=a1+a2;

		if(A<0.0001) t2[i]=0.5*(mm2[i-1]+mm2[i]);
		else t2[i]=(a1*mm2[i-1]+a2*mm2[i])/A;
	}///////////

	/////アキマ補間によるv-B曲線ﾌﾟﾛｯﾄ
	ofstream fout2("v-B_akima.dat");
	for(double i=Bflux[2];i<=Bflux[Nn+1];i+=0.001)
	{
		for(int n=2;n<=Nn+1;n++)
		{
			if(i>=Bflux[n] && i<Bflux[n+1])
			{
				double h=Bflux[n+1]-Bflux[n];
				double aj0=vm[n];
				double aj1=t2[n];
				double aj2=(3*mm2[n]-2*t2[n]-t2[n+1])/h;
				double aj3=(t2[n]+t2[n+1]-2*mm2[n])/(h*h);
				double val=aj0+aj1*(i-Bflux[n])+aj2*(i-Bflux[n])*(i-Bflux[n])+aj3*(i-Bflux[n])*(i-Bflux[n])*(i-Bflux[n]);
				fout2<<i<<" "<<val<<endl;
			}
		}
	}
	fout2.close();
	/////////*/

	double *Bflux2=new double [Nn+4];//Ｂ^2格納
	for(int i=0;i<Nn+4;i++) Bflux2[i]=Bflux[i]*Bflux[i];

	///各区間のv-B^2曲線の傾きを計算
	double *mm3=new double [Nn+3];//各区間の傾き。
	for(int i=0;i<Nn+3;i++) mm3[i]=(vm[i+1]-vm[i])/(Bflux2[i+1]-Bflux2[i]);

	///ユーザーの入力点におけるv-B^2曲線の傾き∂v/∂B^2計算
	double *dvdB2=new double[Nn+4];//傾き
	for(int i=2;i<Nn+2;i++)
	{
		double a1=sqrt((mm3[i+1]-mm3[i])*(mm3[i+1]-mm3[i]));
		double a2=sqrt((mm3[i-1]-mm3[i-2])*(mm3[i-1]-mm3[i-2]));
		double A=a1+a2;

		if(A<0.0001) dvdB2[i]=0.5*(mm3[i-1]+mm3[i]);
		else dvdB2[i]=(a1*mm3[i-1]+a2*mm3[i])/A;
	}///////////

	/////アキマ補間によるv-B^2曲線ﾌﾟﾛｯﾄ
	ofstream fout3("v-B2_akima.dat");
	for(double i=Bflux2[2];i<=Bflux2[Nn+1];i+=0.0001)
	{
		for(int n=2;n<=Nn+1;n++)
		{
			if(i>=Bflux2[n] && i<Bflux2[n+1])
			{
				double h=Bflux2[n+1]-Bflux2[n];
				double aj0=vm[n];
				double aj1=dvdB2[n];
				double aj2=(3*mm3[n]-2*dvdB2[n]-dvdB2[n+1])/h;
				double aj3=(dvdB2[n]+dvdB2[n+1]-2*mm3[n])/(h*h);
				double val=aj0+aj1*(i-Bflux2[n])+aj2*(i-Bflux2[n])*(i-Bflux2[n])+aj3*(i-Bflux2[n])*(i-Bflux2[n])*(i-Bflux2[n]);
				fout3<<i<<" "<<val<<endl;
			}
		}
	}
	fout3.close();
	/////////*/

	//両端2点のdvdB2を求める
	dvdB2[0]=dvdB2[2]*(Bflux2[0]-Bflux2[3])*(Bflux2[0]-Bflux2[4])/((Bflux2[2]-Bflux2[3])*(Bflux2[2]-Bflux2[4]))+dvdB2[3]*(Bflux2[0]-Bflux2[2])*(Bflux2[0]-Bflux2[4])/((Bflux2[3]-Bflux2[2])*(Bflux2[3]-Bflux2[4]))+dvdB2[4]*(Bflux2[0]-Bflux2[2])*(Bflux2[0]-Bflux2[3])/((Bflux2[4]-Bflux2[2])*(Bflux2[4]-Bflux2[3]));
	dvdB2[1]=dvdB2[2]*(Bflux2[1]-Bflux2[3])*(Bflux2[1]-Bflux2[4])/((Bflux2[2]-Bflux2[3])*(Bflux2[2]-Bflux2[4]))+dvdB2[3]*(Bflux2[1]-Bflux2[2])*(Bflux2[1]-Bflux2[4])/((Bflux2[3]-Bflux2[2])*(Bflux2[3]-Bflux2[4]))+dvdB2[4]*(Bflux2[1]-Bflux2[2])*(Bflux2[1]-Bflux2[3])/((Bflux2[4]-Bflux2[2])*(Bflux2[4]-Bflux2[3]));
	dvdB2[Nn+2]=dvdB2[Nn-1]*(Bflux2[Nn+2]-Bflux2[Nn])*(Bflux2[Nn+2]-Bflux2[Nn+1])/((Bflux2[Nn-1]-Bflux2[Nn])*(Bflux2[Nn-1]-Bflux2[Nn+1]))+dvdB2[Nn]*(Bflux2[Nn+2]-Bflux2[Nn-1])*(Bflux2[Nn+2]-Bflux2[Nn+1])/((Bflux2[Nn]-Bflux2[Nn-1])*(Bflux2[Nn]-Bflux2[Nn+1]))+dvdB2[Nn+1]*(Bflux2[Nn+2]-Bflux2[Nn-1])*(Bflux2[Nn+2]-Bflux2[Nn])/((Bflux2[Nn+1]-Bflux2[Nn-1])*(Bflux2[Nn+1]-Bflux2[Nn]));
	dvdB2[Nn+3]=dvdB2[Nn-1]*(Bflux2[Nn+3]-Bflux2[Nn])*(Bflux2[Nn+3]-Bflux2[Nn+1])/((Bflux2[Nn-1]-Bflux2[Nn])*(Bflux2[Nn-1]-Bflux2[Nn+1]))+dvdB2[Nn]*(Bflux2[Nn+3]-Bflux2[Nn-1])*(Bflux2[Nn+3]-Bflux2[Nn+1])/((Bflux2[Nn]-Bflux2[Nn-1])*(Bflux2[Nn]-Bflux2[Nn+1]))+dvdB2[Nn+1]*(Bflux2[Nn+3]-Bflux2[Nn-1])*(Bflux2[Nn+3]-Bflux2[Nn])/((Bflux2[Nn+1]-Bflux2[Nn-1])*(Bflux2[Nn+1]-Bflux2[Nn]));
	
	///各区間の∂v/∂B^2-B^2曲線の傾きを計算
	double *mm4=new double [Nn+3];//各区間の傾き。
	for(int i=0;i<Nn+3;i++) mm4[i]=(dvdB2[i+1]-dvdB2[i])/(Bflux2[i+1]-Bflux2[i]);

	///ユーザーの入力点における∂v/∂B^2-B^2曲線の傾き計算
	double *t4=new double[Nn+4];//傾き
	for(int i=2;i<Nn+2;i++)
	{
		double a1=sqrt((mm4[i+1]-mm4[i])*(mm4[i+1]-mm4[i]));
		double a2=sqrt((mm4[i-1]-mm4[i-2])*(mm4[i-1]-mm4[i-2]));
		double A=a1+a2;

		if(A<0.0001) t4[i]=0.5*(mm4[i-1]+mm4[i]);
		else t4[i]=(a1*mm4[i-1]+a2*mm4[i])/A;
	}///////////

	/////アキマ補間による∂v/∂B^2-B^2曲線ﾌﾟﾛｯﾄ
	ofstream fout4("dvdB2-B2_akima.dat");
	for(double i=Bflux2[2];i<=Bflux2[Nn+1];i+=0.001)
	{
		for(int n=2;n<=Nn+1;n++)
		{
			if(i>=Bflux2[n] && i<Bflux2[n+1])
			{
				double h=Bflux2[n+1]-Bflux2[n];
				double aj0=dvdB2[n];
				double aj1=t4[n];
				double aj2=(3*mm4[n]-2*t4[n]-t4[n+1])/h;
				double aj3=(t4[n]+t4[n+1]-2*mm4[n])/(h*h);
				double val=aj0+aj1*(i-Bflux2[n])+aj2*(i-Bflux2[n])*(i-Bflux2[n])+aj3*(i-Bflux2[n])*(i-Bflux2[n])*(i-Bflux2[n]);
				fout4<<i<<" "<<val<<endl;
			}
		}
	}
	fout4.close();
	/////////*/
	double *A=new double [side_num+1];
    int NN=0;//ディリクレ型境界辺数
    int *dn=new int [side_num+1]; //各辺がディリクレ型境界の何番目か。ディリクレ型境界上でないならside_num+1を格納
    double *PHAT=new double [CON.get_max_DN()];//ディリクレ型境値
	//NNも計算する
	set_boundary_condition3D_edge(CON,NODE,ELEM,SIDE,node,nelm,side_num,dn,&NN,PHAT,A);

	cout<<"ﾃﾞｨﾘｸﾚ数＝"<<NN<<endl;
	/////////////*/
    
	    
    int pn=side_num-NN;				///未知数
    int *ppn=new int [pn];			//行列のn番目は辺ppn[n]
    int *npp=new int [side_num+1];	///各辺が行列の何番目にあたるか。ﾃﾞｨﾘｸﾚ型の場合はpn+1を格納
    int num=0; 
    for(int i=1;i<=side_num;i++)
    {
        if(SIDE[i].boundary_condition==0)//未知数
		{
			ppn[num]=i;
			npp[i]=num;
			num++;
			A[i]=0;//初期化
		}
		else npp[i]=pn+1;
    }
	cout<<"未知数: "<<pn<<endl;
    cout<<"行列幅決定　";
    ////行列の最大幅計算
    int mat_w=0;
	int *nume=new int[side_num+1];
	int *wid=new int [pn+1];

	for(int i=1;i<=side_num;i++) nume[i]=0;
	for(int i=1;i<=nelm;i++)
	{
		for(int j=1;j<=6;j++)
		{
			int edge=ELEM[i].sides[j];
			nume[edge]=nume[edge]+1;
		}
	}
	for(int i=1;i<=side_num;i++)
	{
		int width=6+3*(nume[i]-1);
		if(width>mat_w) mat_w=width;
		if(npp[i]<pn) wid[npp[i]+1]=width;
	}

	delete [] nume;	
	//mat_w*=5;
	////////////
	
    ////配列確保
    double **G=new double *[pn+1];///全体行列
    for(int i=1;i<=pn;i++) G[i]=new double [wid[i]+1];
    int **ROW=new int *[pn+1]; ///各行の、非ゼロ要素の列番号記憶
    for(int i=1;i<=pn;i++) ROW[i]=new int [wid[i]+1];
    int *NUM=new int [pn+1]; ///各行の、非ゼロ要素数
    double *B=new double [pn];//解行列
    ////

	double *vB2=new double [nelm+1];//各要素のdv/dB^2の値格納
	for(int je=1;je<=nelm;je++) 
	{
		vB2[je]=0;	//初期化
		RP[je]=1;	//初期化
		if(ELEM[je].material==IRON) RP[je]=2500;
		else if(ELEM[je].material==MAGELAST) RP[je]=CON.get_RP();
	}
    
	
    /////////全体行列を作成する
    
    int N[4+1]; //要素の各節点番号格納
    double X[4+1];
    double Y[4+1];
    double Z[4+1];
	double b[4+1];
    double c[4+1];
    double d[4+1];
    double e[4+1];
	int table[6+1][2+1];		//要素を構成する辺とその要素内節点番号関係
	cout<<"行列作成開始    ";
	int ENDFLAG=OFF;
	int loop_count=0;			//反復回数。ある一定以上になったら計算を切り上げる

	
	double err2;				//非線形ループ収束判定誤差
	double old_err;
	double CG=0.01; //2e-7;				//収束判定(convergence)
	double CGtype=2;			//収束条件のタイプ　1:V/A 2:A
	while(ENDFLAG==OFF)
	{
		for(int i=1;i<=pn;i++){//初期化
			NUM[i]=0;
			for(int j=1;j<=wid[i];j++){
				G[i][j]=0;
				ROW[i][j]=0;
			}
		}
		for(int i=0;i<pn;i++) B[i]=0;//初期化
		loop_count++;

		for(int je=1;je<=nelm;je++)
		{   
			//辺－節点ﾃｰﾌﾞﾙ作成
			for(int i=1;i<=6;i++)
			{
				int iside=ELEM[je].sides[i];
				int ia=SIDE[iside].node[1];
				int ib=SIDE[iside].node[2];
				for(int j=1;j<=4;j++)
				{
					if(ELEM[je].node[j]==ia) table[i][1]=j;
					else if(ELEM[je].node[j]==ib) table[i][2]=j;
				}
			}////////////////
	
			for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
		
			double Xs=0;//要素の重心座標
			double Ys=0;
			double Zs=0;
			for(int j=1;j<=4;j++)
			{
				X[j]=NODE[N[j]].r[A_X];
				Y[j]=NODE[N[j]].r[A_Y];
				Z[j]=NODE[N[j]].r[A_Z];
				Xs+=X[j]*0.25;
				Ys+=Y[j]*0.25;
				Zs+=Z[j]*0.25;
			}
			////////////////////////////

			///ﾃﾞﾛｰﾆ分割の際に求めた体積は、ｽｰﾊﾟｰﾎﾞｯｸｽ内に移行して計算されているためもはや正しい値ではない。なのでここで求め直す。ただし１度求めたら２回目からは計算しなくていい
			if(loop_count==1) ELEM[je].volume=volume3D(NODE,N[1],N[2],N[3],N[4]);//体積の6倍であることに注意
	
			double delta6=ELEM[je].volume;//体積の6倍
		
			delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
	
			double delta=ELEM[je].volume/6;//本当の体積
	
			for(int i=1;i<=4;i++)
			{
				int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
				int m=j%4+1;
				int n=m%4+1;
	    
				b[i]=X[j]*(Y[n]*Z[m]-Y[m]*Z[n])+X[m]*(Y[j]*Z[n]-Y[n]*Z[j])+X[n]*(Y[m]*Z[j]-Y[j]*Z[m]);//iが偶数のとき
				c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
				d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
				e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
				if(i%2!=0)//iが奇数なら
				{
					b[i]*=-1;
					c[i]*=-1;
					d[i]*=-1;
					e[i]*=-1;
				}
			}
			/////////
	
			double v=v0;
			v/=RP[je]; //vは磁気低効率なのに注意
			double rp=RP[je];
/*			////要素ﾏﾄﾘｸｽ作成開始
			for(int i=1;i<=6;i++)
			{	
				int iside=ELEM[je].sides[i];//要素jeの辺番号
				if(SIDE[iside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
				{   
					int I=npp[iside]+1;///辺isideは行列のI番目
					int I1=SIDE[iside].node[1];//isideを構成する2点
					int I2=SIDE[iside].node[2];
					int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int k2=table[i][2];
					////解行列Ｂを作成
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
						int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[j][2];
					//	B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*v*delta6*delta6*delta6*A[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*delta6*delta6*delta6*A[jside]/rp;
					}
					if(ELEM[je].material==MAGNET)
					{
						B[I-1]+=1.0/18/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
						//B[I-1]+=1/18/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
					}
					else if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
					/////////////////////////////////////
					
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
					
						int m1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int m2=table[j][2];//
						
						if(SIDE[jside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
						{   
							int J=npp[jside]+1;///辺jsideは行列のJ番目
							int flag=0;
						
							int J1=SIDE[jside].node[1];//jsideを構成する2点
							int J2=SIDE[jside].node[2];
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(ELEM[je].material==MAGELAST)
									{
										double Uie=0;//p.106(4.88)の、uに関する∑の合計値
										double Uje=0;//p.106(4.88)の、kに関する∑の合計値
										for(int u=1;u<=6;u++)
										{
											int uside=ELEM[je].sides[u];
											int u1=table[u][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
											int u2=table[u][2];
											Uie+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*A[uside];
											Uje+=(d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])*A[uside];
										}
										double add=2/(18*18*delta*delta*delta*delta)*vB2[je]*Uje;
										G[I][h]+=add*Uie;
										//G[I][h]+=add*Uie*RP[je]*u0;
									}
					//				G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*v*delta6*delta6*delta6;//(三次元有限要素法 磁界解析技術の基礎 p.106(4.88)と同じ.)
									G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*delta6*delta6*delta6/rp;//(三次元有限要素法 磁界解析技術の基礎 p.106(4.88)と同じ.)
									flag=1;
								}
							}
							if(flag==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								if(ELEM[je].material==MAGELAST)
								{
									double Uie=0;//p.106(4.88)の、uに関する∑の合計値
									double Uje=0;//p.106(4.88)の、kに関する∑の合計値
									for(int u=1;u<=6;u++)
									{
										int uside=ELEM[je].sides[u];
										int u1=table[u][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
										int u2=table[u][2];
										Uie+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*A[uside];
										Uje+=(d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])*A[uside];
									}
									double add=2/(18*18*delta*delta*delta*delta)*vB2[je]*Uje;
									G[I][H]+=add*Uie;
									//G[I][H]+=add*Uie*RP[je]*u0;
								}
						//		G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*v*delta6*delta6*delta6;//(三次元有限要素法 磁界解析技術の基礎 p.106(4.88)と同じ.)					    
								G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*delta6*delta6*delta6/rp;//(三次元有限要素法 磁界解析技術の基礎 p.106(4.88)と同じ.)
								ROW[I][H]=J;
							}
						}
						////
						else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら いまは固定値＝０なので計算しないが、将来的にはどうにかすること
						{
							//int n=dn[jside];
							//B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*v*delta6*delta6*delta6*PHAT[n];
						}///////////
					}
				}
			}   	
		}
		///////////////////////*///確認
			for(int i=1;i<=6;i++)
		{	
			int iside=ELEM[je].sides[i];//要素jeの辺番号
			if(SIDE[iside].boundary_condition==0)//ディリクレでない。つまり未知なら
			{   
				int I=npp[iside]+1;///辺isideは行列のI番目
				int I1=SIDE[iside].node[1];//isideを構成する2点
				int I2=SIDE[iside].node[2];
				int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
				int k2=table[i][2];

	//			if(loop_count=1){
			    for(int j=1;j<=6;j++)
				{
					int jside=ELEM[je].sides[j];
					
					int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
					int u2=table[j][2];
						
					if(SIDE[jside].boundary_condition==0)//ディリクレでない。つまり未知なら
					{   
						int J=npp[jside]+1;///辺jsideは行列のJ番目
						int flag=0;
						
						int J1=SIDE[jside].node[1];//jsideを構成する2点
						int J2=SIDE[jside].node[2];
						for(int h=1;h<=NUM[I];h++)
						{
							if(J==ROW[I][h])//すでに同じJが格納されているなら
							{
							//	G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*v;
								G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
								flag=1;
							}
						}
						if(flag==0)
						{   
						    NUM[I]=NUM[I]+1;
						    int H=NUM[I];
						    
						   // G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*v;
							G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6/rp;
							ROW[I][H]=J;
						}
					}
					////
					else //jsideがディリクレ型境界節点なら
					{
					    int n=dn[jside];
					//	B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]*v;
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2.0/3.0*delta6*delta6*delta6*PHAT[n]/rp;
					}//////////
				}

				///B[I-1]を計算する（解行列）
				if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=u0*((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=u0*((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=u0*((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////
				else if(ELEM[je].material==MAGNET)
				{
					B[I-1]+=1.0/18.0/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
				}
//				}
//				else if(loop_count>1){
				////////////

	/*			////解行列Ｂを作成
					for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
						int u1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int u2=table[j][2];
					//	B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*v*delta6*delta6*delta6*A[jside];
						B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*delta6*delta6*delta6*A[jside]/rp;
					}
					if(ELEM[je].material==MAGNET)
					{
						B[I-1]+=1.0/18/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
						//B[I-1]+=1/18/delta*((d[k1]*e[k2]-e[k1]*d[k2])*Mx+(e[k1]*c[k2]-c[k1]*e[k2])*My+(c[k1]*d[k2]-d[k1]*c[k2])*Mz);
					}
					else if(ELEM[je].material==COIL)
				{
					j0x=current[A_X][je];
					j0y=current[A_Y][je];
					j0z=current[A_Z][je];
					B[I-1]+=((b[k1]*c[k2]-b[k2]*c[k1])+(d[k1]*c[k2]-c[k1]*d[k2])*Ys+(e[k1]*c[k2]-c[k1]*e[k2])*Zs)*j0x*delta6/6;
					B[I-1]+=((b[k1]*d[k2]-b[k2]*d[k1])+(c[k1]*d[k2]-d[k1]*c[k2])*Xs+(e[k1]*d[k2]-d[k1]*e[k2])*Zs)*j0y*delta6/6;
					B[I-1]+=((b[k1]*e[k2]-b[k2]*e[k1])+(c[k1]*e[k2]-e[k1]*c[k2])*Xs+(d[k1]*e[k2]-e[k1]*d[k2])*Ys)*j0z*delta6/6;
				}//////////*/

				for(int j=1;j<=6;j++)
					{
						int jside=ELEM[je].sides[j];
					
						int m1=table[j][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int m2=table[j][2];//
						
						if(SIDE[jside].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
						{   
							int J=npp[jside]+1;///辺jsideは行列のJ番目
							int flag=0;
						
							int J1=SIDE[jside].node[1];//jsideを構成する2点
							int J2=SIDE[jside].node[2];
							for(int h=1;h<=NUM[I];h++)
							{
								if(J==ROW[I][h])//すでに同じJが格納されているなら
								{
									if(ELEM[je].material==MAGELAST)
									{
										double Uie=0;//p.106(4.88)の、uに関する∑の合計値
										double Uje=0;//p.106(4.88)の、kに関する∑の合計値
										for(int u=1;u<=6;u++)
										{
											int uside=ELEM[je].sides[u];
											int u1=table[u][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
											int u2=table[u][2];
											Uie+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*A[uside];
											Uje+=(d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])*A[uside];
										}
										double add=2/(18*18*delta*delta*delta*delta)*vB2[je]*Uje;
										G[I][h]+=add*Uie;
										//G[I][h]+=add*Uie*RP[je]*u0;
									}
							//		G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*v*delta6*delta6*delta6;//(三次元有限要素法 磁界解析技術の基礎 p.106(4.88)と同じ.)
									G[I][h]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*delta6*delta6*delta6/rp;//(三次元有限要素法 磁界解析技術の基礎 p.106(4.88)と同じ.)
									flag=1;
								}
							}
							if(flag==0)
							{   
								NUM[I]=NUM[I]+1;
								int H=NUM[I];
								if(ELEM[je].material==MAGELAST)
								{
									double Uie=0;//p.106(4.88)の、uに関する∑の合計値
									double Uje=0;//p.106(4.88)の、kに関する∑の合計値
									for(int u=1;u<=6;u++)
									{
										int uside=ELEM[je].sides[u];
										int u1=table[u][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
										int u2=table[u][2];
										Uie+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*A[uside];
										Uje+=(d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])*A[uside];
									}
									double add=2/(18*18*delta*delta*delta*delta)*vB2[je]*Uje;
									G[I][H]+=add*Uie;
									//G[I][H]+=add*Uie*RP[je]*u0;
								}
						//		G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*v*delta6*delta6*delta6;//(三次元有限要素法 磁界解析技術の基礎 p.106(4.88)と同じ.)					    
								G[I][H]+=((d[k1]*e[k2]-e[k1]*d[k2])*(d[m1]*e[m2]-e[m1]*d[m2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[m1]*c[m2]-c[m1]*e[m2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[m1]*d[m2]-d[m1]*c[m2]))*2/3*delta6*delta6*delta6/rp;//(三次元有限要素法 磁界解析技術の基礎 p.106(4.88)と同じ.)
								ROW[I][H]=J;
							}
						}
						////
						else //jsideがﾃﾞｨﾘｸﾚ型境界節点なら いまは固定値＝０なので計算しないが、将来的にはどうにかすること
						{
							//int n=dn[jside];
							//B[I-1]-=((d[k1]*e[k2]-e[k1]*d[k2])*(d[u1]*e[u2]-e[u1]*d[u2])+(e[k1]*c[k2]-c[k1]*e[k2])*(e[u1]*c[u2]-c[u1]*e[u2])+(c[k1]*d[k2]-d[k1]*c[k2])*(c[u1]*d[u2]-d[u1]*c[u2]))*2/3*v*delta6*delta6*delta6*PHAT[n];
						}///////////
					}//
		//		}
			}
		}   	
    }
    ///////////////////////*/

		int number=0;//行列の非ゼロ要素数
		for(int i=1;i<=pn;i++) number+=NUM[i];
    
		///行列の実際の最大幅を求める
		if(loop_count==1)
		{
			int maxN=0;
			for(int i=1;i<=pn;i++) if(NUM[i]>maxN) maxN=NUM[i];
			cout<<"最大幅："<<maxN<<"/"<<mat_w<<endl;
			mat_w=maxN; //最大幅の書き換え
		}

		///このままではG[i][j]は[j]の小さい順に並んでいないので、GとROWを並び替える
		arrange_matrix(pn,NUM,ROW,G);

		//対称性チェック
		check_matrix_symmetry(pn,NUM,ROW,G);
    
		ofstream fout("matrix.dat");
		for(int i=1;i<=100;i++)
		{
			for(int j=1;j<=NUM[i];j++)
			{
				fout<<G[i][j]<<" ";
				//cout<<ROW[i][j]<<" ";
			}
			fout<<endl;
			//cout<<endl;
		}
		fout.close();/////////*/

		double *val = new double [number];
		int *ind = new int [number];//非ゼロ要素の列番号格納配列
		int *ptr = new int [pn+1];//各行の要素がvalの何番目からはじまるのかを格納
    
		/////////////////////val,ind ,ptrに値を格納
		int index=0;
		for(int n=0;n<pn;n++)
		{
			ptr[n]=index;
			for(int m=1;m<=NUM[n+1];m++)
			{
				val[index]=G[n+1][m];
				ind[index]=ROW[n+1][m]-1;
				index++;
			}
		}	    
		ptr[pn]=number;
		////////////////////*/

		///////////////////////行列計算開始
//		if(CON.get_FEMCG()==0) ICCG3D(val,ind,ptr,pn,ppn,B,A);
//		else
//		{
			double *XX=new double [pn];//行列の答え格納
			if(CON.get_FEMCG()==1) ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
			else if(CON.get_FEMCG()==2) parallel_ICCG3D2(CON,val,ind,ptr,pn,B,number,XX);
			else if(CON.get_FEMCG()==3) MRTR(CON,B,pn,XX,val,ind,ptr);
			else if(CON.get_FEMCG()==4) ICMRTR(CON,B,pn,XX,val,ind,ptr);
	
			for(int n=0;n<pn;n++)
			{
				int i=ppn[n];
				V[i]=XX[n];
			}
			delete [] XX;
//		}///////////////////////////*/
    
		//ﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙ更新
		err2=0;//誤差を初期化
		if(loop_count==1)
		{
			for(int i=1;i<=side_num;i++)
			{
				if(SIDE[i].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
				{
					double a=A[i];
					A[i]+=V[i];
					if(CGtype==1)
					{
						if(A[i]!=0) err2+=sqrt((V[i]/A[i])*(V[i]/A[i]));
					}
					else if(CGtype==2) err2+=sqrt((A[i]-a)*(A[i]-a));
				}
			}
		}
		else if(loop_count>1 && CGtype==1)
		{
			double *A2=new double [side_num+1];
			
			double alpha=2;//減速係数
			double err3=0;
			int breakflag=OFF;
			double tempA=alpha;
			double temperr=old_err;
			for(int n=1;n<=10;n++)
			{
				for(int i=1;i<=side_num;i++) A2[i]=A[i];
				err3=0;
				alpha*=0.5;
				for(int i=1;i<=side_num;i++)
				{
					if(SIDE[i].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
					{
						A2[i]+=alpha*V[i];
						if(A2[i]!=0) err3+=sqrt((V[i]/A2[i])*(V[i]/A2[i]));
					}
				}
				if(err3<temperr)
				{
					tempA=alpha;
					temperr=err3;
				}
			}///αがもとまった
			alpha=tempA;
			cout<<loop_count<<"回目のalpha="<<alpha;
			//そのαを使って実際にＡ[i]を更新
			for(int i=1;i<=side_num;i++)
			{
				if(SIDE[i].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
				{
					A[i]+=alpha*V[i];
					if(A[i]!=0) err2+=sqrt((V[i]/A[i])*(V[i]/A[i]));
				}
			}
			delete [] A2;
		}//////////*/
		else if(loop_count>1 && CGtype==2)
		{
			for(int i=1;i<=side_num;i++)
			{
				if(SIDE[i].boundary_condition==0)///ﾃﾞｨﾘｸﾚでない。つまり未知なら
				{
					double a=A[i];
					A[i]+=V[i];
				//	err2+=sqrt((A[i]-a)*(A[i]-a));
					err2+=sqrt(V[i]*V[i]);
				}
			}
		}//////////*/
		
		
		if(err2<CG || loop_count>5) ENDFLAG=ON;
		cout<<"  誤差="<<err2<<endl;
		old_err=err2;//前回のerrを保存
		//if(CGtype==2 && loop_count==1) CG=err2*1e-6;//誤差判定としてA[i]のみを使うときは、ﾓﾃﾞﾙにより最適な誤差判定が変わるので、このように定義する

		////////////////////////////新しい透磁率をもとめるための、磁性体要素の磁束をもとめる
		if(ENDFLAG==OFF)
		{
			for(int je=1;je<=nelm;je++)
			{   
				if(ELEM[je].material==MAGELAST)
				{
					//辺－節点ﾃｰﾌﾞﾙ作成
					for(int i=1;i<=6;i++)
					{
						int iside=ELEM[je].sides[i];
						int ia=SIDE[iside].node[1];
						int ib=SIDE[iside].node[2];
						for(int j=1;j<=4;j++)
						{
							if(ELEM[je].node[j]==ia) table[i][1]=j;
							else if(ELEM[je].node[j]==ib) table[i][2]=j;
						}
					}////////////////
			
					for(int j=1;j<=4;j++) N[j]=ELEM[je].node[j];
					
					double Xs=0;//要素の重心座標
					double Ys=0;
					double Zs=0;
					for(int j=1;j<=4;j++)
					{
						X[j]=NODE[N[j]].r[A_X];
						Y[j]=NODE[N[j]].r[A_Y];
						Z[j]=NODE[N[j]].r[A_Z];
						Xs+=X[j];
						Ys+=Y[j];
						Zs+=Z[j];
					}
					Xs/=4;Ys/=4;Zs/=4;
					////////////////////////////
		
					double delta6=ELEM[je].volume;//体積の6倍(正しい体積の値はすでにベクトルポテンシャルを求める際に計算している)
			
					delta6=1/delta6;//計算に必要なのは逆数なのでここで反転しておく
			
					double delta=ELEM[je].volume/6;//本当の体積
			
					for(int i=1;i<=4;i++)
					{
						int j=i%4+1;///i,j,m,nは循環する整数(三次元有限要素法 磁界解析技術の基礎 P14参照 ただし式1.39のxm-znはxm-xnの間違い) 
						int m=j%4+1;
						int n=m%4+1;
			    
						c[i]=(Y[j]*(Z[m]-Z[n])+Y[m]*(Z[n]-Z[j])+Y[n]*(Z[j]-Z[m]));//iが偶数のとき
						d[i]=(Z[j]*(X[m]-X[n])+Z[m]*(X[n]-X[j])+Z[n]*(X[j]-X[m]));//iが偶数のとき
						e[i]=(X[j]*(Y[m]-Y[n])+X[m]*(Y[n]-Y[j])+X[n]*(Y[j]-Y[m]));//iが偶数のとき
						if(i%2!=0)//iが奇数なら
						{
						    c[i]*=-1;
							d[i]*=-1;
							e[i]*=-1;
						}
					}
					/////////
		
					for(int D=0;D<3;D++) Be[D][je]=0;//初期化
		
					for(int i=1;i<=6;i++)
					{
						int s=ELEM[je].sides[i];
						int k1=table[i][1];//各点は要素の第何番節点か(三次元有限要素法 磁界解析技術の基礎 p.42(3.11)と同じ記号)
						int k2=table[i][2];
			
						Be[A_X][je]+=(d[k1]*e[k2]-e[k1]*d[k2])*A[s];
						Be[A_Y][je]+=(e[k1]*c[k2]-c[k1]*e[k2])*A[s];
						Be[A_Z][je]+=(c[k1]*d[k2]-d[k1]*c[k2])*A[s];
					}
			
					for(int D=0;D<3;D++) Be[D][je]*=delta6*delta6*2;
			
					double BB=0;//磁束密度の２乗
					for(int D=0;D<3;D++) BB+=Be[D][je]*Be[D][je];
//					cout<<"je="<<je<<", B="<<sqrt(BB)<<endl;
					/////アキマ補間によるv-B^2曲線から要素のRP[i]をもとめる
					int flagB=0;
					for(int n=2;n<=Nn+1;n++)
					{
						if(BB>=Bflux2[n] && BB<Bflux2[n+1])
						{
							double h=Bflux2[n+1]-Bflux2[n];
							double aj0=vm[n];
							double aj1=dvdB2[n];
							double aj2=(3*mm3[n]-2*dvdB2[n]-dvdB2[n+1])/h;
							double aj3=(dvdB2[n]+dvdB2[n+1]-2*mm3[n])/(h*h);
							double val=aj0+aj1*(BB-Bflux2[n])+aj2*(BB-Bflux2[n])*(BB-Bflux2[n])+aj3*(BB-Bflux2[n])*(BB-Bflux2[n])*(BB-Bflux2[n]);
							RP[je]=1/val/u0;
//							cout<<"je="<<je<<",RP="<<RP[je]<<"aj0="<<aj0<<"aj1="<<aj1<<"aj2="<<aj2<<"aj3="<<aj3<<"val="<<val<<endl;
							double bj0=dvdB2[n];
							double bj1=t4[n];
							double bj2=(3*mm4[n]-2*t4[n]-t4[n+1])/h;
							double bj3=(t4[n]+t4[n+1]-2*mm4[n])/(h*h);
							double val2=bj0+bj1*(BB-Bflux2[n])+bj2*(BB-Bflux2[n])*(BB-Bflux2[n])+bj3*(BB-Bflux2[n])*(BB-Bflux2[n])*(BB-Bflux2[n]);
							vB2[je]=val2;
							flagB=1;
						}
						
					}
					if(flagB==0) cout<<"磁束密度が想定領域外です B="<<sqrt(BB)<<endl;
					
				}
			}
		}
		delete [] val;
		delete [] ind;
		delete [] ptr;
	}
	///本来プログラムはV[i]を節点iのﾍﾞｸﾄﾙﾎﾟﾃﾝｼｬﾙとして計算するので、ここでV[i]=A[i]として値を渡しておく
	for(int i=1;i<=side_num;i++) V[i]=A[i];
	

    ///////////
	ofstream fp2("V.dat");
    for(int i=1;i<=node;i++)
    {
		if(NODE[i].boundary_condition==1) fp2<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
    }
    fp2.close();
	ofstream fp4("RP.dat");
	for(int i=1;i<=nelm;i++)
	{
		if(ELEM[i].material==MAGELAST)
		{
			double R[3]={0,0,0};
			for(int j=1;j<=4;j++) for(int D=0;D<3;D++) R[D]+=NODE[ELEM[i].node[j]].r[D]*0.25;
		
			if(R[1]<0.0015 && R[1]>0) 
			{
				fp4<<R[0]<<" "<<R[2]<<" "<<RP[i]<<endl;
			}
		}
	}
	fp4.close();

    for(int i=1;i<=pn;i++) delete [] G[i];
    delete [] G;
    for(int i=1;i<=pn;i++) delete [] ROW[i];
    delete [] ROW;
    delete [] NUM;
	delete [] wid;
    ///////////////////////*/
    delete [] dn;
    delete [] PHAT;
    delete [] npp;
    delete [] ppn;
    delete [] B;
	delete [] A;
    
    

	delete [] vB2;

	delete [] mm;
	delete [] vm;
	delete [] mm2;
	delete [] t2;
	delete [] Bflux2;
	delete [] mm3;
	delete [] dvdB2;
	delete [] mm4;
	delete [] t4;
}
//ICMRTR法//
void ICMRTR(mpsconfig &CON,double *B,int pn,double *X,double *val,int *ind,int *ptr)
{
	////////////////////////////////////////////////////////////////
	//val :ゼロ要素の値
	//ind:非ゼロ要素の列番号格納配列
	//ptr:各行の要素がvalの何番目からはじまるのかを格納
	//X[n]:解
	int count=0;
	double E=1;//誤差
	double BB=0;
	double rr=0;
		
	double alp;
	double beta;
	double zeta;
	double zeta_old;
	double eta;
	double nu;

	double Ar_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	double Ar_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	double y_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	double Ar_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)

	double w_r;
	double v_w;
	double y_w;
	double w_y;

	double *r= new double[pn];
	double *Ar = new double [pn];  //   _ 
	double *P = new double [pn];
	double *y = new double [pn];
	double *u = new double [pn];
	double *v = new double [pn];
	double *w = new double [pn];
	double *z = new double [pn];

	zeta=0;
	zeta_old=0;
	eta=0;
	nu=1;

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++)//初期値
	{
		X[n]=0;
		r[n]=B[n];
		P[n]=r[n];
		y[n]=0;
	}

	//////////前処理//////////
	double accel_re=CON.get_CGaccl();
	double accel;	//加速ファクタ
	accel= accel_re;

	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;

	double *val2=new double [num2];//係数行列を保存(対角成分を含む下三角行列) 非ゼロ要素だけを1次配列として保存
	int *ind2 = new int [num2];
	int *ptr2 = new int [pn+1];

	num2=0;
	double one=1;
	double sum=0;
	for(int k=0;k<pn;k++)
	{	
		ptr2[k]=num2;
	    for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
	    {
			if(ind[m]<=k)
			{
				val2[num2]=val[m];
				ind2[num2]=ind[m];
				if(ind[m]==k) val2[num2]*=accel;
				num2++;
			}
		}
	}
	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない

	int *NUM = new int [pn];			//列方向にみた、各列の要素数
	for(int k=0;k<pn;k++) NUM[k]=0;
	
	for(int k=0;k<pn;k++)
	{	
	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
	    {
			int J=ind2[m];
			NUM[J]=NUM[J]+1;
		}
	}
	double **VAL=new double *[pn];						//ゼロ要素の値 VAL[i][k]はi列のk番目の非ゼロ要素
	for(int i=0;i<pn;i++) VAL[i]=new double [NUM[i]];
	int **IND = new int *[pn];							//非ゼロ要素の行番号格納配列
	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
	
	/////////////////////////////////////iccg法
	double rLDLt_r;
	double *y2= new double [pn];
	double *LDLt_r= new double [pn];
//	double *L=new double[num2];//不完全コレスキー分解後の下三角行列L格納
	double *D1 = new double [pn];//D行列
	/////不完全コレスキー分解
//	Incomplete_Cholesky_Decomposition(CON,L,D1,val2,ptr2,ind2,pn,num2);//LとD1に値が書き込まれる
	for(int k=0;k<pn;k++)
	{
		for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非0要素
		{
			int i=ind2[m];//列番号
			if(i==0)
			{
				val2[m]=val2[m];
			}
			if(i>0 && i<k)
			{
				sum=0;
				for(int j=ptr2[k];j<m;j++)
				{
					for(int J=ptr2[i];J<ptr2[i+1];J++)
					{
						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
					}
				}
				val2[m]=val2[m]-sum;
			}
			if(i==k)
			{
				sum=0;
				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
				val2[m]=val2[m]-sum;
				D1[k]=one/val2[m];
			}
		}
	}
	//不完全コレスキー分解完了///////
//	delete [] val2;
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
		if(i==0) y2[0]=r[0]/val2[0]; //式（3.77） 
		else
		{
		    double sum=0;
			for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//式（3.78）
		    int m=ptr2[i+1]-1;
			y2[i]=(r[i]-sum)/val2[m];
		}
	}////y[i]がもとまった。
	for(int i=pn-1;i>=0;i--)
	{
	    double sum=0;
		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
	    LDLt_r[i]=y2[i]-D1[i]*sum;	
	}
	/////////////////*/

	for(int n=0;n<pn;n++) u[n]=LDLt_r[n];
	////////////////////////////////////////////////////////////////
	for(int n=0;n<pn;n++) BB+=abs(B[n])*abs(B[n]);
	ofstream ICEEe("ICMRTR.dat", ios::trunc);
	ICEEe.close();
	cout<<"ICMRTR法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON.get_omp_P()==OFF)//通常版
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		while(E>CON.get_MRTRep())// CON.get_MRTRep()//収束判定(convergence test)
		{
			if(count==pn) cout<<"count=pn"<<endl;
			ofstream ICEE("ICMRTR.dat", ios::app);
			ICEE<<count<<" "<<E<<endl;
			ICEE.close();
			count++;
			/////vの計算
			for(int n=0;n<pn;n++)
			{
				v[n]=0;
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					v[n]+=val[j]*u[ind[j]];
				}
			}
		
			/////Ar,Ar_c計算
			for(int n=0;n<pn;n++)
			{    
				Ar[n]=0.0;
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*u[ind[j]];
				}
			}
		
			//////wの計算
			//////////////////////////y[i]を求め、LDLt_r[i]を求める。
			for(int i=0;i<pn;i++)
			{
				if(i==0) y2[0]=v[0]/val2[0]; //式(3.77)
				else
				{
					sum=0;
					for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//式（3.78）
				    int m=ptr2[i+1]-1;
						y2[i]=(v[i]-sum)/val2[m];
				}
			}////y[i]がもとまった。
		
			for(int i=pn-1;i>=0;i--)
			{
				double sum=0;
				for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
				LDLt_r[i]=y2[i]-D1[i]*sum;	
			}
			/////////////////*/
			for(int n=0;n<pn;n++) w[n]=LDLt_r[n];

			w_r=0;
			v_w=0;
			y_w=0;
			w_y=0;
		
			//内積の計算
			for(int n=0;n<pn;n++)
			{
				w_r+=w[n]*r[n];
				v_w+=v[n]*w[n];
				y_w+=y[n]*w[n];
				w_y+=w[n]*y[n];
			}	

			zeta=nu*w_r/(nu*v_w-y_w*w_y);
			eta=-y_w*w_r/(nu*v_w-y_w*w_y);
			//////nuの計算//////////////
			nu=0.0;
			for(int n=0;n<pn;n++) nu+=zeta*r[n]*w[n];
			
			//////pの計算//////////////
			for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////Xの計算////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////yの計算////////////
			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*Ar[n];
	
			//////rの計算////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//////zの計算///////////
			for(int n=0;n<pn;n++) z[n]=eta*z[n]+zeta*w[n];

			//////uの計算///////////
			for(int n=0;n<pn;n++) u[n]-=z[n];

			/////////////////*/
			//誤差評価
			rr=0;
			for(int n=0;n<pn;n++) rr+=abs(r[n])*abs(r[n]);
			
			E=sqrt(rr/BB);
			c<<count<<" "<<E<<endl;
			if(count==30000) break;
	
			if(count>pn)
			{
				count=1;		//反復回数を１にして関数から脱出。反復回数が１だから、エラーだとわかる
				E=0;
			}
		}
	}
	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
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
	delete [] y2;

	delete [] r;
	delete [] Ar;
	delete [] P;
	delete [] y;
	delete [] u;
	delete [] v;
	delete [] w;
	delete [] z;
	
	/////////////////
	c.close();
}
//MRTR法//
void MRTR(mpsconfig &CON,double *B,int pn,double *X,double *val,int *ind,int *ptr)
{
	
	//unsigned int timeCG=GetTickCount();
	int count=0;
	double EP=CON.get_MRTRep();
	//complex<double> rr=(0,0);
	double E=1;//誤差
	double BB=0;
	double rr=0;
		
	double zeta;
	double zeta_old;
	double eta;
	double v;

	double Ar_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
	double Ar_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	double y_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	double Ar_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)

    double *r= new double[pn];
	double *Ar = new double [pn];  //   _ _
	
	double *P = new double [pn];
	double *y = new double [pn];

	zeta=0;
	zeta_old=0;
	eta=0.0;
	//v=complex<double>(1.0,0.0);//1回目の反復において、0でないことに注意
	v=1.0;

	ofstream c("convergence.dat");

	for(int n=0;n<pn;n++) //初期値
	{
		X[n]=0.0;			//近似解
		r[n]=B[n];			//右辺行列　Ax=b
		P[n]=r[n];			//rは残差
		y[n]=0.0;
	}

	for(int n=0;n<pn;n++) BB+=r[n]*r[n];
	//BB=sqrt(BB);
	//cout<<"||r0||2="<<sqrt(BB)<<endl;
	//cout<<"r[0]="<<r[0]<<" r[1]="<<r[1]<<endl;
	ofstream EEe("MRTR.dat", ios::trunc);
	EEe.close();
	cout<<"MRTR法スタート-----未知数="<<pn<<"  ---";
	unsigned int time=GetTickCount();
	if(CON.get_omp_P()==OFF)//通常版
	{
		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
		while(E>EP)// EP=CON->get_CGep();//収束判定(convergence test)
		{
			//if(count==pn) cout<<"count=pn"<<endl;
			count++;
			ofstream EE("MRTR.dat", ios::app);
			EE<<count<<" "<<E<<endl;
			EE.close();
			for(int n=0;n<pn;n++)//Ar,Ar_c計算
			{    
				Ar[n]=0.0;
				for(int j=ptr[n];j<ptr[n+1];j++)
				{
					Ar[n]+=val[j]*r[ind[j]];
				}
			}

			Ar_r=0.0;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
			Ar_Ar=0.0;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			y_Ar=0.0;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
			Ar_y=0.0;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
	
			//#pragma omp parallel for
			for(int n=0;n<pn;n++)  
			{
				Ar_r+=Ar[n]*r[n];
				Ar_Ar+=Ar[n]*Ar[n];
				y_Ar+=y[n]*Ar[n];
				Ar_y+=Ar[n]*y[n];
			}

			zeta=v*Ar_r/(v*Ar_Ar-y_Ar*Ar_y);
			eta=-y_Ar*Ar_r/(v*Ar_Ar-y_Ar*Ar_y);
		//	cout<<zeta<<" "<<eta<<endl;

			//////vの計算//////////////
			v=0.0;
			for(int n=0;n<pn;n++) v+=zeta*Ar[n]*r[n];

			//////pの計算//////////////
			for(int n=0;n<pn;n++) P[n]=r[n] + (zeta_old*eta/zeta)*P[n];
			zeta_old=zeta;

			///////Xの計算////////////
			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
	
			//////yの計算////////////
			for(int n=0;n<pn;n++)
			{
				y[n]=eta*y[n]+zeta*Ar[n];
				//cout<<y[n]<<endl;
			}

			//////rの計算////////////
			for(int n=0;n<pn;n++) r[n]-=y[n];

			//誤差評価
			rr=0;
			for(int n=0;n<pn;n++) rr+=r[n]*r[n];
			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
			
			E=sqrt(rr/BB);
			//cout<<E<<endl;
			c<<count<<" "<<E<<endl;
			//if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
	
			if(count>pn)
			{
				count=1;		//反復回数を１にして関数から脱出。反復回数が１だから、エラーだとわかる
				E=0;
			}
			if(count>=10000) break;
		}
	}
	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
	
	delete [] r;
	delete [] Ar;
	
	delete [] P;
	delete [] y;
	c.close();
}
////ic付きMRTR法
//void cs_ICMRTR(mpsconfig *CON,complex<double> *val,int *ind,int *ptr,int pn,complex<double> *B,int number,complex<double> *X)
//{
//	
//	//unsigned int timeCG=GetTickCount();
//	int count=0;
//	//complex<double> rr=(0,0);
//	double E=1;//誤差
//	double BB=0;
//	double rr=0;
//	
//	complex<double> alp;
//	complex<double> beta;
//	complex<double> zeta;
//	complex<double> zeta_old;
//	complex<double> eta;
//	complex<double> nu;
//
//	complex<double> cAr_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
//	complex<double> cAr_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//	complex<double> cy_Ar;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//	complex<double> cAr_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//
//	complex<double> cw_r;//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
//	complex<double> cv_w;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//	complex<double> cy_w;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//	complex<double> cw_y;//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//
//    complex<double> *r= new complex<double>[pn];
//	complex<double> *Ar = new complex<double> [pn];  //   _ _
//	complex<double> *cAr = new complex<double> [pn];//cAr:A r(複素共役をバーで表記)
//	complex<double> *P = new complex<double> [pn];
//	complex<double> *y = new complex<double> [pn];
//	complex<double> *u = new complex<double> [pn];
//	complex<double> *v = new complex<double> [pn];
//	complex<double> *w = new complex<double> [pn];
//	complex<double> *z = new complex<double> [pn];
//
//	zeta=complex<double>(0.0,0.0);
//	zeta_old=complex<double>(0.0,0.0);
//	eta=complex<double>(0.0,0.0);
//	nu=complex<double>(1.0,0.0);//1回目の反復において、0でないことに注意
//
//	ofstream c("convergence.dat");
//
//	for(int n=0;n<pn;n++) //初期値
//	{
//		X[n]=complex<double> (0.0,0.0);
//		r[n]=B[n];
//		P[n]=r[n];
//		y[n]=complex<double> (0.0,0.0);
//		u[n]=complex<double> (0.0,0.0);
//		v[n]=complex<double> (0.0,0.0);
//		w[n]=complex<double> (0.0,0.0);
//		z[n]=complex<double> (0.0,0.0);
//	}
//
//	//////前処理/////
//	double accel_re=CON->get_CGaccl();
//	complex<double> accel;//CON->get_CGaccl();//加速ファクタ
//	accel=complex<double> (accel_re,0);
//	double accel2=0.001;//複素シフト用加速ファクタ
//	
//	int num2=0;//対角成分を含む、下三角行列だけを考慮にいれた非ゼロ要素数
//	for(int k=0;k<pn;k++) for(int m=ptr[k];m<ptr[k+1];m++) if(ind[m]<=k) num2++;	
//	if(num2!=(number-pn)/2+pn) cout<<"ERROR"<<endl;
//	
//	complex<double> *val2=new complex<double> [num2];
//	int *ind2 = new int [num2];
//	int *ptr2 = new int [pn+1];
//
//	num2=0;
//
//	complex<double> one;
//	one=complex<double> (1.0,0);
//	complex<double> Im;
//	Im=complex<double> (0.0,1.0);
//	complex<double> sum;
//	sum=complex<double> (0.0,0.0);
//
//	for(int k=0;k<pn;k++)
//	{	
//		ptr2[k]=num2;
//	    for(int m=ptr[k];m<ptr[k+1];m++)///k行目の非０要素
//	    {
//			if(ind[m]<=k)
//			{
//				val2[num2]=val[m];
//				ind2[num2]=ind[m];
//				if(ind[m]==k) val2[num2]*=accel;//加速ﾌｧｸﾀ
//				//if(ind[m]==k) val2[num2]+=accel2*Im;//複素シフト
//				num2++;
//			}
//		}
//	}
//	ptr2[pn]=num2;//これをしておかないと、最後に(int m=ptr2[k];m<ptr2[k+1];m++)みたいなことができない
//
//	int *NUM = new int [pn];
//	for(int k=0;k<pn;k++) NUM[k]=0;
//	
//	for(int k=0;k<pn;k++)
//	{	
//	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
//	    {
//			int J=ind2[m];
//			NUM[J]=NUM[J]+1;
//		}
//	}
//	complex<double> **VAL=new complex<double> *[pn];//ゼロ要素の値
//	for(int i=0;i<pn;i++) VAL[i]=new complex<double> [NUM[i]];
//	int **IND = new int *[pn];//非ゼロ要素の行番号格納配列
//	for(int i=0;i<pn;i++) IND[i]=new int [NUM[i]];
//	
//	/////////////////////////////////////iccg法
//	complex<double> rLDLt_r;
//	complex<double> *y2=new complex<double> [pn];
//	complex<double> *LDLt_r= new complex<double> [pn];
//	complex<double> *D1 = new complex<double> [pn];//D行列
//	
//	/////不完全コレスキｰ分解
//	for(int k=0;k<pn;k++)
//	{	
//	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
//	    {
//	        int i=ind2[m];//列番号
//	        if(i==0)
//			{
//				val2[m]=val2[m];
//				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
//				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
//		    
//			}
//			if(i>0 && i<k)
//			{
//				sum=complex<double> (0,0);
//				
//				for(int j=ptr2[k];j<m;j++)
//				{	
//					for(int J=ptr2[i];J<ptr2[i+1];J++)//i行目のなかから列の一致するものを探している。少し手間か？
//					{
//						if(ind2[J]==ind2[j]) sum+=val2[j]*val2[J]*D1[ind2[j]];
//					}
//				}
//				val2[m]=val2[m]-sum;
//			}
//			if(i==k)
//			{
//				sum=complex<double> (0,0);
//				for(int j=ptr2[k];j<m;j++) sum+=val2[j]*val2[j]*D1[ind2[j]];
//				val2[m]=val2[m]-sum;
//				//if(val2[m]<0.0001 &&val2[m]>=0) val2[m]=0.0001;
//				//else if(val2[m]>-0.0001 &&val2[m]<=0) val2[m]=-0.0001;
//				//if(val2[m].real()<0.0001 &&val2[m].real()>=0) val2[m].real(0.0001);
//				//else if(val2[m].real()>-0.0001 &&val2[m].real()<=0) val2[m].real(-0.0001);
//				//if(val2[m].imag()<0.0001 &&val2[m].imag()>=0) val2[m].imag(0.0001);
//				//else if(val2[m].imag()>-0.0001 &&val2[m].imag()<=0) val2[m].imag(-0.0001);
//				
//				D1[k]=one/val2[m];
//				//if(val2[m]>0) cout<<"EE"<<endl;
//            }
//	    }
//	}    
//	///不完全ｺﾚｽｷｰ分解完了/////////*/
//
//	///列を基準にした配列に値を代入
//	for(int k=0;k<pn;k++) NUM[k]=0;
//	for(int k=0;k<pn;k++)
//	{	
//	    for(int m=ptr2[k];m<ptr2[k+1];m++)///k行目の非０要素
//	    {
//			int J=ind2[m];
//			VAL[J][NUM[J]]=val2[m];
//			IND[J][NUM[J]]=k;
//			NUM[J]=NUM[J]+1;
//		}
//	}////////*/
//
//	/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
//	for(int i=0;i<pn;i++)
//	{
//		if(i==0) y2[0]=r[0]/val2[0]; //式（3.77） 
//		else
//		{
//		    sum=complex<double> (0,0);
//		    /////////        
//		    for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//式（3.78）
//		    int m=ptr2[i+1]-1;
//		    y2[i]=(r[i]-sum)/val2[m];
//		}
//	}////y[i]がもとまった。
//	for(int i=pn-1;i>=0;i--)
//	{
//	    sum=complex<double> (0,0);
//		for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
//	    LDLt_r[i]=y2[i]-D1[i]*sum;	
//	}
//	/////////////////*///LDLtがM-1に相当？つまり、LDLt_rはu0と等しい
//	
//	for(int n=0;n<pn;n++) u[n]=LDLt_r[n];
//		/////////////////*/
//	//////////////////////
//
//	//for(int n=0;n<pn;n++) BB+=norm(B[n]);
//	for(int n=0;n<pn;n++) BB+=abs(B[n])*abs(B[n]);
//	//BB=sqrt(BB);
//	cout<<"||r0||2="<<sqrt(BB)<<endl;
//
//
//	 cout<<"cs_ICMRTR法スタート-----未知数="<<pn<<"  ---";
//	unsigned int time=GetTickCount();
//	if(CON->get_omp_P()==OFF)//通常版
//	{
//		//while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
//		while(E>CON->get_FEMCGep())// EP=CON->get_CGep();//収束判定(convergence test)
//		{
//			if(count==pn) cout<<"count=pn"<<endl;
//			count++;
//
//			/////vの計算
//			for(int n=0;n<pn;n++)
//			{    
//				v[n]=complex<double>(0.0,0.0);//vに相当
//				for(int j=ptr[n];j<ptr[n+1];j++)
//				{
//					v[n]+=val[j]*u[ind[j]];
//				}
//			}
//
//			/*///
//			for(int n=0;n<pn;n++)//Ar,Ar_c計算
//			{    
//				Ar[n]=complex<double>(0.0,0.0);//vに相当
//				cAr[n]=complex<double>(0.0,0.0);//conj(v)に相当
//				for(int j=ptr[n];j<ptr[n+1];j++)
//				{
//					Ar[n]+=val[j]*u[ind[j]];
//					cAr[n]+=conj(val[j]) * conj(u[ind[j]]);
//				}
//			}
//			///*/
//
//			////wの計算
//
//			/////////////////y[i]をもとめ、LDLt_r[i]をもとめる。
//			for(int i=0;i<pn;i++)
//			{
//				if(i==0) y2[0]=v[0]/val2[0]; //式（3.77） 
//				else
//				{
//					sum=complex<double> (0,0);
//					/////////        
//					for(int m=ptr2[i];m<ptr2[i+1]-1;m++) sum+=val2[m]*y2[ind2[m]];//式（3.78）
//					int m=ptr2[i+1]-1;
//					y2[i]=(v[i]-sum)/val2[m];
//				}
//			}////y[i]がもとまった。
//			for(int i=pn-1;i>=0;i--)
//			{
//				sum=complex<double> (0,0);
//				for(int h=1;h<NUM[i];h++) sum+=VAL[i][h]*LDLt_r[IND[i][h]];
//				LDLt_r[i]=y2[i]-D1[i]*sum;	
//			}
//			/////////////////*///LDLtがM-1に相当？つまり、LDLt_rはu0と等しい 同様に、LDL_vはw？
//	
//			for(int n=0;n<pn;n++) w[n]=LDLt_r[n];
//
//			cw_r=complex<double>(0.0,0.0);//cAr_r: (cAr,r)(内積)を指す。(u,v)=Σ((cu)*v)
//			cv_w=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//			cy_w=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//			cw_y=complex<double>(0.0,0.0);//cAr_Ar: (cAr,Ar)(内積)を指す。(u,v)=Σ((cu)*v)
//	
//			//内積の計算
//			for(int n=0;n<pn;n++)  
//			{
//				cw_r+=w[n]*r[n];
//				cv_w+=v[n]*w[n];
//				cy_w+=y[n]*w[n];
//				cw_y+=w[n]*y[n];
//			}
//
//			zeta=nu*cw_r/(nu*cv_w-cy_w*cw_y);
//			eta=-cy_w*cw_r/(nu*cv_w-cy_w*cw_y);
//
//			//////nuの計算//////////////
//			nu=complex<double>(0.0,0.0);
//			for(int n=0;n<pn;n++) nu+=zeta*r[n]*w[n];
//
//			//////pの計算//////////////
//			for(int n=0;n<pn;n++) P[n]=u[n] + (zeta_old*eta/zeta)*P[n];
//			zeta_old=zeta;
//
//			///////Xの計算////////////
//			for(int n=0;n<pn;n++) X[n]+=zeta*P[n];
//	
//			//////yの計算////////////
//			for(int n=0;n<pn;n++) y[n]=eta*y[n]+zeta*v[n];
//
//			//////rの計算////////////
//			for(int n=0;n<pn;n++) r[n]-=y[n];
//
//			//////zの計算////////////
//			for(int n=0;n<pn;n++) z[n]=eta*z[n]+zeta*w[n];
//
//			//////uの計算////////////
//			for(int n=0;n<pn;n++) u[n]-=z[n];
//
//			//誤差評価
//			rr=0;
//			//for(int n=0;n<pn;n++) rr+=norm(r[n]);
//			for(int n=0;n<pn;n++) rr+=abs(r[n])*abs(r[n]);
//			//for(int n=0;n<pn;n++) rr+=conj(r[n])*conj(r[n]);
//			
//			E=sqrt(rr/BB);
//			c<<count<<" "<<E<<endl;
//			if(count%1000==0) cout<<"反復回数="<<count<<" E="<<E<<endl;
//	
//			
//		}
//	}
//
//	cout<<"終了 反復回数="<<count<<" time="<<(GetTickCount()-time)*0.001<<"[sec]"<<endl;
//	
//	delete [] LDLt_r;
//	delete [] D1;
//
//	delete [] val2;
//	delete [] ind2;
//	delete [] ptr2;
//
//	for(int i=0;i<pn;i++)
//	{
//		delete [] VAL[i];
//		delete [] IND[i];
//	}
//	delete [] VAL;
//	delete [] IND;
//	delete [] NUM;
//	delete [] y2;
//
//	delete [] r;
//	delete [] Ar;
//	delete [] cAr;	
//	delete [] P;
//	delete [] y;
//	delete [] u;
//	delete [] v;
//	delete [] w;
//	delete [] z;
//
//	c.close();
//}

