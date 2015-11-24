#include"stdafx.h"	

#define FULL 1
#define HALF 2
#define HALFD 3
#define HALFD_shell 4

//using namespace std;

void writedata2(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double pax,double pay,double paz,double P,double h,int toBEM);
double get_volume(mpsconfig *CON);

//直線を分割するさいの最適な分割数と分割距離の算出関数 他所で使うためfunction.hに移動
void calc_N_and_L(double dis,double le,int *N,double *L);
//円周分割数計算関数
int calc_division_N_circle(double dis,double le);
//半径Rの円の外周
void set_circle_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R);
void set_circle_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R, double height);
//半径Rの円内部
void set_circle_in(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &type1,int *number,double le,double R,int edge_startID,int edge_lastID);
void set_circle_in_using_6_pieces(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int edge_startID,int edge_lastID);
void set_circle_in_using_6_pieces(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int edge_startID,int edge_lastID);
//半径Rの球作成関数
//長方形作成関数
void set_rectangular(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height);
void set_sphere(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag);
void set_sphere2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag,int *suf_num);
//長方形の辺作成関数
void set_rectangular_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V);
//長方形内部作成関数
void set_rectangular_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V,int edge_startID,int edge_lastID);
//円柱表面作成関数
void set_cylinder_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int circle_start_id,int circle_end_id,int top_flag);
//円柱内部設置関数
void set_cylinder_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int flag,int stert_id);
void set_cylinder_in(vector<double> &X,vector<double> &Y,vector<double> &Z,double erast_r,double erast_h,int *number,double le,double R,double height,int flag,int stert_id);
//ドーナツ作成
void set_doughnut2D(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_big,double R_smal,int edge_startID,int edge_lastID);
//FSWプローブ内部粒子セット関数
void set_hat_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID);
//箱作成関数
void set_box(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &surface,int *number,double le,double Width,double Height,double Depth);
//BOX内作成関数
void set_box_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height,double Depth,int BO_startID,int BO_lastID);
//るつぼ壁内作成関数（るつぼの肉厚部）
void set_crucible_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R_out,double height,double height_out,int flag,int fluid_number);

//分子動力学関数
void MD_2D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN);
void MD_3D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN,double r,double region[3][2]);

//物質合成関数
void make_fusion3D(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<double> &X2,vector<double> &Y2,vector<double> &Z2,vector<int> &surface2,int *number,double le);

void set_initial_placement_using_MD(mpsconfig *CON,int *particle_number)
{
	cout<<"初期粒子配置をMDにより最適化中--";
	unsigned int timeA=GetTickCount();
	int number=0;	//粒子数　最後にはparticle_numberに格納
	int model=CON->get_model_number();
	int dimension=CON->get_dimension();		//解析次元

	double val=0;
	double val2=0;

	double le=CON->get_distancebp();	//初期粒子間距離*0.8 le=0.001
	double A=sqrt(3.0)*0.5;				//よく使う係数(√3)/2
	double B=sqrt(2.0/3);				//よく使う係数√(2/3)
	

	vector<double> X;
	vector<double> Y;
	vector<double> Z;
	vector<int> type1;
	vector<int> surface;

	vector<double> X2;					//使い回し用
	vector<double> Y2;
	vector<double> Z2;
	vector<int> surface2;			//ONなら表面　OFFなら内部

	vector<double> X3;					//使い回し用 弾丸２用
	vector<double> Y3;
	vector<double> Z3;
	vector<int> surface3;			//ONなら表面　OFFなら内部

	vector<double> X4;					//使い回し用　弾丸２用
	vector<double> Y4;
	vector<double> Z4;

	vector<double> X5;					//使い回し用　弾丸1用
	vector<double> Y5;
	vector<double> Z5;
	vector<int> surface5;			//ONなら表面　OFFなら内部
	
	vector<double> X6;					//使い回し用　弾丸1用
	vector<double> Y6;
	vector<double> Z6;
	
	ofstream fq("initial_input.dat");

	cout<<model<<endl;
	//////////////////モデル1　球体//////////////////////////
	 if(model==1)
	{
		double R=CON->get_R1();			//作成する円の半径
		double Zg=R+5.2*le+10*le;					//球の中心高さ
		if(dimension==2)
		{
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
		}
		else if(dimension==3)
		{
			//円作成
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている

			//球作成
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,R,flag);
		}
		//初期情報書き込み
		//for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FRFLUID,1,0,Y[i]*40,0,0,0,0);
		double vx=0;
		double vy=0;
		double vz=0;
		if(dimension==3)
		{
			for(int i=0;i<number;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				double a=0.2;
		//		a=500*3;//クーロン分裂のとき
				a=200;	//遊び
				
				
				vx=-0.5*a*X[i];
				vy=-0.5*a*Y[i];
				vz=a*Z[i];
				 
			//	writedata2(fq,  i,X[i],Y[i],Z[i]+Zg,MAGELAST,1, OFF,val, vx,vy, vz,0,0,0,0, 0, 0);
				writedata2(fq,  i,X[i],Y[i],Z[i]+Zg,MAGELAST,1,1,val, 0,0, 0,0,0,0,0, 0, 1);
			}
			//for(int i=0;i<number;i++) if(fabs(Z[i])<0.001) writedata2(fq,  i,X[i],Y[i],Z[i]+Zg,FLUID,1, OFF,0, 0,0, 0, 0, 0, ON);

			//下壁
			//double Width=27*le;
			double Rw=CON->get_R1()*1.8;
			double Height=6*le;
			int numw=number;
			double mdis=le; //*0.5
			set_circle_edge(X,Y,Z,&number,mdis,Rw);//円外周　これがないと内部も作れない
			set_circle_in_using_6_pieces(X,Y,Z,&number,mdis,Rw,numw,number);//円内部    vector配列は参照渡ししている　//0をnum3に変更

			int end_id=number;	//円の粒子idを記憶
			////////

			int topw_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,mdis,Rw,Height,numw,end_id,topw_flag);//円柱表面座標作成

			set_cylinder_in(X,Y,Z,&number,mdis,Rw,Height,1,numw);//内部にパッキング　　　//最後にint num3(壁の最初の粒子数)を加える．MD_3Dで必要
			int ii=numw;
			for(int i=numw;i<number;i++){
				if(Z[i]>=Height/2){
					ii++;
				writedata2(fq,ii,X[i], Y[i], Z[i]+3*le-Height,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
				}
			}
			number=ii;
		}
		else if(dimension==2)
		{
			for(int i=0;i<number;i++)
			{
				double a=500*20;	
				//a*=3;
				vx=-a*X[i];
				vy=a*Y[i];
		
				writedata2(fq,  i,X[i],Y[i]+Zg,Z[i],ELASTIC,1, OFF,val,vx,vy, vz,0,0,0, 0, 0,0);
			}
		}
		//for(int i=0;i<number;i++) writedata2(fq,  i,X[i],Y[i],Z[i],FLUID,1, OFF,0, 0,Y[i]*40, 0, 0, 0, ON);
		cout<<"model1完了"<<endl;
	}
	////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////モデル2　円筒体/////////////////////////////////////////////
	else if(model==2)
	{
		double R=CON->get_fluidwidth()*le*0.8;//半径
		double height=(8*le*A)*2;//(6*le*A)*2;これは何？
		double vz=-2.0;
		if(dimension==3)
		{
			//円作成
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//円外周　これがないと内部も作れない
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
//			set_circle_edge(X,Y,Z,&number,le,R, height);//円外周　これがないと内部も作れない
//			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,height,0,number);//円内部    vector配列は参照渡ししている
			int circle_end_id=number;	//円の粒子idを記憶
			////////

			int top_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,le,R,height,circle_start_id,circle_end_id,top_flag);//円柱表面座標作成
			int face_n=number; //
			///////
			for(int s=0;s<face_n;s++){
				writedata2(fq,s,X[s],Z[s]-height/2,Y[s],TERMINAL1,1,1,0,0,0,vz,0,0,0,0,0,1);//粒子は,FACE
			}
			///////////////////////////////////////////////////////////////////////////////
			set_cylinder_in(X,Y,Z,&number,le,R,height,1,circle_start_id);//内部にパッキング

			int beforeNumber=number;
			
			//円柱表面+下壁の書き込み
			for(int i=face_n;i<number;i++) {
				 writedata2(fq,i,X[i],Z[i]-height/2,Y[i],TERMINAL1,1,0,0,0,0,vz,0,0,0,0,0,1);
			}
			//下壁
			//double Width=27*le;
			double Rw=CON->get_fluidwidth()*le*1.8;
			double Height=6*le;
			int num3=number;
			double mdis=le; //*0.5
			set_circle_edge(X,Y,Z,&number,mdis,Rw);//円外周　これがないと内部も作れない
			set_circle_in_using_6_pieces(X,Y,Z,&number,mdis,Rw,num3,number);//円内部    vector配列は参照渡ししている　//0をnum3に変更

			int end_id=number;	//円の粒子idを記憶
			////////

			int topw_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,mdis,Rw,Height,num3,end_id,topw_flag);//円柱表面座標作成

			set_cylinder_in(X,Y,Z,&number,mdis,Rw,Height,1,num3);//内部にパッキング　　　//最後にint num3(壁の最初の粒子数)を加える．MD_3Dで必要
			int ii=beforeNumber-1;
			for(int i=beforeNumber;i<number;i++){
				if(Z[i]>=Height/2){
					ii++;
				writedata2(fq,ii,X[i], Y[i], Z[i]-R-Height-4*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
				}
		}
			number=ii;
			//床上面の位置
			double ground=Z[beforeNumber]-R-Height-4*le;
			for(int i=beforeNumber;i<number;i++)
			{
				if((Z[i]-R-Height-4*le)>ground) ground=Z[i]-R-Height-4*le;
			}
			CON->set_ground_position(ground);
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////モデル　3　立方体///////////////////////////////////////////////////////////
	else if(model==3){
	if(dimension==2)
		{
			double origin[3]={0,0,0};//四角形の左下の点の座標
			double Width=CON->get_fluidwidth()*le;				//水平方向長さ
			double Height=CON->get_fluidwidth()*le;				//垂直方向長さ
			set_rectangular(X,Y,Z,&number,le,Width,Height);		//座標の原点は長方形の中心。よってあとで移動すること。
			
			//初期情報書き込み
			for(int i=0;i<number;i++) writedata2(fq,i,0.5*Width+origin[A_X]+X[i],Height+origin[A_Y]+Y[i],Z[i],MAGELAST,1,0,0,0,0,0,0,0,0,0,0,0);

		}
		else if(dimension==3)//立方体の箱
		{
			double Width=8*le;		
			double Height=10*le;	
			double Depth=10*le;
			set_box(X,Y,Z,surface,&number,le,Width,Height,Depth);//最後3つの引数は横、高さ、奥行き。デカルト座標の原点に対し、X正方向に横幅、Y正方向に奥行き、Z正方向に高さ
			for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],TERMINAL1,1, OFF,0, 0,0, 0,0,0,0, 0, 0, 0);
			for(int i=0;i<number;i++){
				X.push_back(X[i]);
				Y.push_back(Y[i]);
				Z.push_back(Z[i]);
			}
			for(int i=number;i<number*2;i++) writedata2(fq,i,X[i],Y[i]-Depth/2,Z[i]+Height+3*le,TERMINAL2,1, OFF,0, 0,0, -0.2, 0,0,0,0, 0, 0);
			number=number*2;
		
				//下壁
			//double Width=27*le;
			double Rw=CON->get_fluidwidth()*le*1.8;
			double Heightw=6*le;
			int num3=number;
			double mdis=le; //*0.5
			set_circle_edge(X,Y,Z,&number,mdis,Rw);//円外周　これがないと内部も作れない
			set_circle_in_using_6_pieces(X,Y,Z,&number,mdis,Rw,num3,number);//円内部    vector配列は参照渡ししている　//0をnum3に変更
			int end_id=number;	//円の粒子idを記憶
			////////

			int topw_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,mdis,Rw,Heightw,num3,end_id,topw_flag);//円柱表面座標作成

			set_cylinder_in(X,Y,Z,&number,mdis,Rw,Heightw,1,num3);//内部にパッキング　　　//最後にint num3(壁の最初の粒子数)を加える．MD_3Dで必要
			int ii=num3-1;
			for(int i=num3;i<number;i++){
				if(Z[i]>=Height/2){
					ii++;
				writedata2(fq,ii,X[i], Y[i], Z[i]-Heightw-4*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
				}
				
		}
			number=ii;
			cout<<"333"<<endl;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////モデル　4　試験片///////////////////////////////////////////////////////////////////
	else if(model==4){
		int length=25;	//[mm]引っ張り端子とは別35
		int width=11;//15
		int depth=11;
		int tan=0;//端子高さ粒子5
		double accel=0.005;
		le*=0.5;//*(11.0/12.0)
		
		////////////下端子////////////
		for(int i=0;i<tan;i++){
			for(int j=0;j<width;j++){
				for(int k=0;k<depth;k++){
					X.push_back(k*le);
					Y.push_back(j*le);
					Z.push_back(i*le);
					number++;
				}
			}
		}
		int low=number;
		//////////////////////////////
		/////////////MRE部分////////////
		for(int i=tan;i<length+tan;i++){
			for(int j=0;j<width;j++){
				for(int k=0;k<depth;k++){
					X.push_back(k*le);
					Y.push_back(j*le);
					Z.push_back(i*le);
					number++;
				}
			}
		}
		int mid=number;
		/////////////////////////////////
		//////////////上端子/////////////
		for(int i=length+tan;i<length+tan*2;i++){
			for(int j=0;j<width;j++){
				for(int k=0;k<depth;k++){
					X.push_back(k*le);
					Y.push_back(j*le);
					Z.push_back(i*le);
					number++;
				}
			}
		}
		////////////////////////////////
		
		for(int i=0;i<number;i++){ 
			if(i<low)writedata2(fq,i,X[i]-depth*le/2,Y[i]-width*le/2,Z[i]-(length+10)*le/2,WALL,0, OFF,0, 0, 0, -accel, 0,0,0,0, 0, 0);
			else if(i>=low && i<mid)writedata2(fq,i,X[i]-depth*le/2,Y[i]-width*le/2,Z[i]-(length+10)*le/2,MAGELAST,0, OFF,0, 0,0, 0,0,0,0, 0, 0, 0);
			else if(i>=mid)writedata2(fq,i,X[i]-depth*le/2,Y[i]-width*le/2,Z[i]-(length+10)*le/2,WALL,0, OFF,0, 0, 0, accel, 0,0,0,0, 0, 0);
		}
	}

	///////////////////////////////モデル　6　アクチュエータ////////////////////////////////////////////////////////////////
	else if(model==6)
	{
	double elast_r=0.0125;		//2.5cm
	double elast_h=0.03;		//3.0cm
	double mag_r=0.015;			//3.0cm
	double mag_h=0.04;			//4.0cm

	int start_ID=0;
	set_circle_edge(X,Y,Z,&number,le,mag_r);//円外周　これがないと内部も作れない
	set_circle_in_using_6_pieces(X,Y,Z,&number,le,mag_r,start_ID,number);//円内部    vector配列は参照渡ししている
	int circle_ID=number;			//円の粒子数　粒子IDはstart_IDからcircle_ID-1まで

	int top_flag=ON;		//円柱の上面も作成するからフラグをON
	set_cylinder_face(X,Y,Z,&number,le,mag_r,mag_h,start_ID,circle_ID,top_flag);//円柱表面座標作成
	int face_n=number;
	for(int s=0;s<face_n;s++){
		writedata2(fq,s,X[s],Y[s],Z[s]-(mag_h/2),MAGELAST,1,1,0,0,0,0,0,0,0,0,0,1);//粒子は,FACE
	}
	set_cylinder_in(X,Y,Z,&number,le,mag_r,mag_h,1,start_ID);//内部にパッキング
	int cylinder_ID=number;

	number=face_n;//粒子数のリセット
	for(int i=face_n;i<cylinder_ID;i++)
		{
			if((X[i]*X[i]+Y[i]*Y[i]<=pow(0.008,2) && Z[i]>0.012 && Z[i]<0.016+0.012)){
				if(X[i]*X[i]+Y[i]*Y[i]<=pow(0.005,2) && Z[i]>0.015 && Z[i]<0.010+0.015){;}//真ん中をくりぬく
				else{
					writedata2(fq,i,X[i],Y[i],Z[i]-(mag_h/2),WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
					number++;
				}
			}
			else if(X[i]*X[i]+Y[i]*Y[i]<=pow(elast_r,2) && Z[i]>0.005 && Z[i]<elast_h+0.0049){
				writedata2(fq,i,X[i],Y[i],Z[i]-(mag_h/2),ELASTIC,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はELASTIC
				number++;
			}
			else {
				writedata2(fq,i,X[i],Y[i],Z[i]-(mag_h/2),MAGELAST,1,0,0,0,0,0,0,0,0,0,0,1);//粒子はMAGELAST
				number++;
				}
		}

	//床上面の位置
			double ground=Z[number]-(elast_h/2)-10*le-2*le;
			for(int i=cylinder_ID;i<number;i++)
			{
				if((Z[i]-(elast_h/2)-10*le-2*le)>ground) ground=Z[i]-(elast_h/2)-10*le-2*le;
			}
			CON->set_ground_position(ground);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////モデル　7　アクチュエータ2////////////////////////////////////
	else if(model==7){
		le*=0.8;
		double r=0.018;	//0.036[m]
		double Wr=CON->get_fluidwidth()*le*1.5;//壁の半径
		double Wh=5*le;//壁の高さ
		double MRE_h=0.004;//[m]
		double silicone_h=0.006;
		double under=0;//-le*24
	////////////////////////////MRE///////////////////////////////
	int start_ID=0;
	set_circle_edge(X,Y,Z,&number,le,r);//円外周　これがないと内部も作れない
	set_circle_in_using_6_pieces(X,Y,Z,&number,le,r,start_ID,number);//円内部    vector配列は参照渡ししている
	int circle_ID=number;			//円の粒子数　粒子IDはstart_IDからcircle_ID-1まで

	int top_flag=ON;		//円柱の上面も作成するからフラグをON
	set_cylinder_face(X,Y,Z,&number,le,r,MRE_h,start_ID,circle_ID,top_flag);//円柱表面座標作成
	int surface=number;
	for(int i=start_ID;i<surface;i++){
		writedata2(fq,i,X[i],Y[i],Z[i]+under+silicone_h+1.5*le,MAGELAST,1,1,0,0,0,0,0,0,0,0,0,1);//粒子はMAGELAST
	}
	///内部パッキン
	set_cylinder_in(X,Y,Z,&number,le,r,MRE_h,1,start_ID);
	int cylinder_ID=number;

	for(int i=surface;i<cylinder_ID;i++){
		writedata2(fq,i,X[i],Y[i],Z[i]+under+silicone_h+1.5*le,MAGELAST,1,0,0,0,0,0,0,0,0,0,0,1);//粒子はMAGELAST
	}
	///////////////////////////////////////////////////////////////////

	///////////////////////////シリコーン///////////////////////////////
	int start2_ID=number;
	set_circle_edge(X,Y,Z,&number,le,r);//円外周　これがないと内部も作れない
	set_circle_in_using_6_pieces(X,Y,Z,&number,le,r,start2_ID,number);//円内部    vector配列は参照渡ししている
	int circle2_ID=number;			//円の粒子数　粒子IDはstart_IDからcircle_ID-1まで

			//円柱の上面も作成するからフラグをON
	set_cylinder_face(X,Y,Z,&number,le,r,silicone_h,start2_ID,circle2_ID,top_flag);//円柱表面座標作成
	set_cylinder_in(X,Y,Z,&number,le,r,silicone_h,1,start2_ID);//内部にパッキング
	int cylinder2_ID=number;

	for(int i=cylinder_ID;i<cylinder2_ID;i++){
		writedata2(fq,i,X[i],Y[i],Z[i]+under,ELASTIC,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はELASTIC
	}
	///////////////////////////////////////////////////////////
		//下壁
			int num3=number;
			double mdis=le; //*0.5
			set_circle_edge(X,Y,Z,&number,mdis,Wr);//円外周　これがないと内部も作れない
			set_circle_in_using_6_pieces(X,Y,Z,&number,mdis,Wr,num3,number);//円内部    vector配列は参照渡ししている　//0をnum3に変更
			int end_id=number;	//円の粒子idを記憶

			int topw_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,mdis,Wr,Wh,num3,end_id,topw_flag);//円柱表面座標作成
			set_cylinder_in(X,Y,Z,&number,mdis,Wr,Wh,1,num3);//内部にパッキング　　　//最後にint num3(壁の最初の粒子数)を加える．MD_3Dで必要
		
			for(int i=num3;i<number;i++){
				writedata2(fq,i,X[i], Y[i], Z[i]-Wh-1.5*le+under,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
		}
			
			//床上面の位置
			double ground=Z[num3]-Wh-1.5*le+under;
			for(int i=num3;i<number;i++)
			{
				if((Z[i]-Wh-1.5*le+under)>ground) ground=Z[i]-Wh-1.5*le+under;
			}
			CON->set_ground_position(ground);
			
	}
	//////////////////////////////////////////モデル　8　アクチュエータ3////////////////////////////////////
	else if(model==8){
		double r=0.018;	//0.036[m]
		double MRE_r=0.003;//[m]
		double MRE_h=0.003;//[m]
		double silicone_h=0.01;
		double under=-le*24;
	
	///////////////////////////シリコーン///////////////////////////////
	int start_ID=0;
	set_circle_edge(X,Y,Z,&number,le,r);//円外周　これがないと内部も作れない
	set_circle_in_using_6_pieces(X,Y,Z,&number,le,r,start_ID,number);//円内部    vector配列は参照渡ししている
	int circle_ID=number;			//円の粒子数　粒子IDはstart_IDからcircle_ID-1まで

	int top_flag=ON;		//円柱の上面も作成するからフラグをON		
	set_cylinder_face(X,Y,Z,&number,le,r,silicone_h,start_ID,circle_ID,top_flag);//円柱表面座標作成
	set_cylinder_in(X,Y,Z,&number,le,r,silicone_h,1,start_ID);//内部にパッキング
	int cylinder2_ID=number;

	for(int i=0;i<cylinder2_ID;i++){
		if(pow(X[i],2)+pow(Y[i]-MRE_r,2)<=pow(MRE_r,2) && (Z[i]>=silicone_h-MRE_h))
		{
			writedata2(fq,i,X[i],Y[i],Z[i]+under,MAGELAST,1,0,0,0,0,0,0,0,0,0,0,1);//粒子はELASTIC
		}
		else writedata2(fq,i,X[i],Y[i],Z[i]+under,ELASTIC,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はELASTIC
	}
	///////////////////////////////////////////////////////////
	//下壁
		//double Width=27*le;
		double Rw=CON->get_fluidwidth()*le*1.5;
		double Height=6*le;
		int num3=number;
		double mdis=le; //*0.5
		set_circle_edge(X,Y,Z,&number,mdis,Rw);//円外周　これがないと内部も作れない
		set_circle_in_using_6_pieces(X,Y,Z,&number,mdis,Rw,num3,number);//円内部    vector配列は参照渡ししている　//0をnum3に変更
		int end_id=number;	//円の粒子idを記憶

		int topw_flag=ON;		//円柱の上面も作成するからフラグをON
		set_cylinder_face(X,Y,Z,&number,mdis,Rw,Height,num3,end_id,topw_flag);//円柱表面座標作成
		set_cylinder_in(X,Y,Z,&number,mdis,Rw,Height,1,num3);//内部にパッキング　　　//最後にint num3(壁の最初の粒子数)を加える．MD_3Dで必要
		
		for(int i=num3;i<number;i++)
		{
			writedata2(fq,i,X[i], Y[i], Z[i]-Height-1.5*le+under,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
		}			
	}
	///////////////////////////モデル　9　薄膜//////////////////////
	else if(model==9){
		double Width=10*le;		
		double Height=5*le;	
		double Depth=20*le;

		set_box(X,Y,Z,surface,&number,le,Width,Height,Depth);//最後3つの引数は横、高さ、奥行き。デカルト座標の原点に対し、X正方向に横幅、Y正方向に奥行き、Z正方向に高さ
		for(int i=0;i<number;i++) writedata2(fq,i,X[i]-Width-Width/2-le,Y[i],Z[i]-Height/2,WALL,1, OFF,0, 0,0, 0,0,0,0, 0, 0, 0);

		for(int i=0;i<number;i++){
			X.push_back(X[i]);
			Y.push_back(Y[i]);
			Z.push_back(Z[i]);
		}
		for(int i=number;i<number*2;i++) writedata2(fq,i,X[i]+Width/2,Y[i],Z[i]-Height/2,WALL,1, OFF,0, 0,0,0,0,0, 0, 0, 0, 0);
		number=number*2;

		//弾性体
		int Wid=20;
		int Hei=1;
		int Dep=20;
		for(int wi=0;wi<Wid;wi++){
			for(int hi=0;hi<Hei;hi++){
				for(int di=0;di<Dep;di++){
			//		if(wi==5)writedata2(fq,number,(wi*le)-Wid*le/2,(di*le)-Dep*le/2,(hi*le)+Height/2+2*le,MAGELAST,1, OFF,0, 0,0, -0.1, 0, 0, 0);
					writedata2(fq,number,(wi*le)-Wid*le/2,(di*le),(hi*le)+Height/2+2*le,MAGELAST,1, OFF,0, 0,0, 0, 0,0,0,0, 0, 0);
					number++;
				}
			}
		}
		//
		//重し
		for(int i=0;i<Dep;i++){
			for(int j=0;j<2;j++){
			writedata2(fq,number,0.0,i*le,Height+(Hei+5+j)*le,WALL,1, OFF,0, 0, 0, -0.1,0,0,0, 0, 0, 0);
			number++;
			}
		}

	}
	/////////////////////////モデル　10 片持ち梁/////////////////////
	else if(model==10){
		int Length=(int)(0.1/CON->get_distancebp())+1; //長さ100mm/粒子の大きさ
		int Height=(int)(0.01/CON->get_distancebp());	
		int Width=(int)(0.01/CON->get_distancebp());
		CON->Set_length(Length-1);
		for(int k=0;k<Width;k++){
			for(int j=0;j<Height;j++){
				for(int i=0;i<Length;i++){
					if(i==0)writedata2(fq,number,i*le,k*le,j*le,WALL,1, OFF,0, 0, 0, 0,0,0, 0, 0, 0, 0);
					else if(i==(Length-1) && j==Height)writedata2(fq,number,i*le,k*le,j*le,MAGELAST,1, OFF,0, 0, 0, 0,0,0,0, 0, 0, 0);
					else writedata2(fq,number,i*le,k*le,j*le,MAGELAST,1, OFF,0, 0, 0, 0, 0,0,0,0, 0, 0);
				number++;
				}
			}
		}
	}
	///////////////////////////モデル　11　ミッキー//////////////////////
	else if(model==11){
		double ear_R=CON->get_R1()/2; //耳の半径
		double spac=CON->get_R1()*2;	//耳同士の間隔
		double ear_hi=CON->get_R1(); //耳の高さ
		double face_R=CON->get_R1(); //顔の半径
		cout<<"耳の作成"<<endl;

		int num=0;	//粒子数
	 //左耳作成
			set_circle_edge(X,Y,Z,&number,le,ear_R);//円外周
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,ear_R,0,number);//円内部    vector配列は参照渡ししている

			//球作成
			int flag=FULL;
			set_sphere(X,Y,Z,&number,le,ear_R,flag);

			for(int i=0;i<X.size();i++){
				X[i]=X[i]-spac/2;
				Z[i]=Z[i]+ear_hi+10*le+CON->get_R1();
			}
			
			for(int i=0;i<X.size();i++){
//				if(sqrt(pow(X[i],2)+pow(Y[i],2)+pow(Z[i],2)) > face_R+0.1*le)
				{
				 X2.push_back(X[i]);
				 Y2.push_back(Y[i]);
				 Z2.push_back(Z[i]);
				}
			}
			X.clear();
			Y.clear();
			Z.clear();
			num+=(int)X2.size();
		//右耳作成
			int number2=0;
			set_circle_edge(X,Y,Z,&number2,le,ear_R);//円外周
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number2,le,ear_R,0,number2);//円内部    vector配列は参照渡ししている

			//球作成
			set_sphere(X,Y,Z,&number2,le,ear_R,flag);

			
			for(int i=0;i<X2.size();i++){
				X[i]=X[i]+spac/2;
				Z[i]=Z[i]+ear_hi+10*le+CON->get_R1();
			}
			for(int i=0;i<X.size();i++){
//				if(sqrt(pow(X[i],2)+pow(Y[i],2)+pow(Z[i],2)) > face_R+0.1*le)
				{
				 X3.push_back(X[i]);
				 Y3.push_back(Y[i]);
				 Z3.push_back(Z[i]);
				}
			}
			num+=(int)X3.size();
			cout<<"顔の作成"<<endl;
			 //顔作成
			int number3=0;
			set_circle_edge(X4,Y4,Z4,&number3,le,face_R);//円外周
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X4,Y4,Z4,&number3,le,face_R,0,number3);//円内部    vector配列は参照渡ししている.

			//球作成
			set_sphere(X4,Y4,Z4,&number3,le,face_R,flag);

			for(int i=0;i<X4.size();i++){
//				if(sqrt(pow(spac/2-X4[i],2)+pow(Y4[i],2)+pow(ear_hi+10*le+CON->get_re()-Z4[i],2)) > face_R-0.2*le && sqrt(pow(spac/2+X4[i],2)+pow(Y4[i],2)+pow(ear_hi+10*le+CON->get_re()-Z4[i],2)) > face_R-0.2*le)
				surface3.push_back(0);
			}
			num+=(int)X4.size();
			//合成
//			make_fusion3D(X,Y,Z,X3,Y3,Z3,surface3,&num,le);
//			make_fusion3D(X,Y,Z,X2,Y2,Z2,surface2,&num,le);
			//出力
			for(int i=0;i<X2.size();i++){
				writedata2(fq,i,X2[i],Y2[i],Z2[i]+5*le,ELASTIC,1, 0,0, 0, 0, 0,0,0,0, 0, 0, 0);
			}
			for(int i=0;i<X3.size();i++){
				writedata2(fq,i+X2.size(),X3[i],Y3[i],Z3[i]+5*le,ELASTIC,1, 0,0, 0, 0, 0,0,0,0, 0, 0,0);
			}
			for(int i=0;i<X4.size();i++){
				writedata2(fq,i+X2.size()+X3.size(),X4[i],Y4[i],Z4[i]+10*le+CON->get_R1()+5*le,MAGELAST,1, 0,0, 0, 0, 0,0,0,0, 0, 0,1);
			}
			
			number=num;//*/
			cout<<"床の作成"<<endl;
			//下壁
			double Wr=20*le;//壁の半径
			double Wh=6*le;//壁の高さ
			int num3=num;
			int number4=0;
/*			double mdis=le; //*0.5
			set_circle_edge(X5,Y5,Z5,&number4,mdis,Wr);//円外周　これがないと内部も作れない
			set_circle_in_using_6_pieces(X5,Y5,Z5,&number4,mdis,Wr,0,number4);//円内部    vector配列は参照渡ししている　//0をnum3に変更
			int end_id=number;	//円の粒子idを記憶

			int topw_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X5,Y5,Z5,&number4,mdis,Wr,Wh,0,end_id,topw_flag);//円柱表面座標作成
			set_cylinder_in(X5,Y5,Z5,&number4,mdis,Wr,Wh,1,0);//内部にパッキング　　　//最後にint num3(壁の最初の粒子数)を加える．MD_3Dで必要
		int ii=0;
			for(int i=0;i<number4;i++){
				if(Z5[i]>=Wh/2){
				writedata2(fq,ii+num3,X5[i], Y5[i], Z5[i]-Wh-8*le-face_R/2,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
				ii++;
				}
		}
			number+=ii;//*/
			int ii=0;
			for(int i=0;i<30;i++){
				for(int j=0;j<30;j++){
					for(int k=0;k<3;k++){
						writedata2(fq,ii+num3,i*le-(15*le), j*le-(15*le), k*le+6*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
						ii++;
					}
				}
			}
			number+=ii;
	}

	//////////////////モデル13　球体2つ//////////////////////////
	else if(model==15)
	{
		double R=CON->get_R1();			//作成する円の半径
		double Zg=R*3;					//球の中心高さ
		double Xg=3*CON->get_R1();//球の座標
		int dimension=3;
		int number2=0;
		int flag=FULL;
		
		//初期情報書き込み
		//for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FRFLUID,1,0,Y[i]*40,0,0,0,0);
		double vx=0;
		double vy=0;
		double vz=0;

			//円作成
			set_circle_edge(X,Y,Z,&number,le,R);//円外周
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
			//球作成
			set_sphere(X,Y,Z,&number,le,R,flag);

			for(int i=0;i<number;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				double a=0.2;
		//		a=500*3;//クーロン分裂のとき
				a=200;	//遊び
								
				vx=-0.5*a*X[i];
				vy=-0.5*a*Y[i];
				vz=a*Z[i];
				 
			//	writedata2(fq,  i,X[i],Y[i],Z[i]+Zg,MAGELAST,1, OFF,val, vx,vy, vz,0,0,0,0, 0, 0);
				writedata2(fq,  i,X[i]+Xg,Y[i],Z[i]+Zg+2*R,MAGELAST,1,1,val, 0,0, 0,0,0,0,0, 0, 1);
			}
			cout<<"球１完成\n";

			//円作成
			set_circle_edge(X2,Y2,Z2,&number2,le,R);//円外周
			//set_circle_in(X,Y,Z,type1,&number,le,R,0,number);////円内部  vector配列は参照渡ししている
			set_circle_in_using_6_pieces(X2,Y2,Z2,&number2,le,R,0,number2);//円内部    vector配列は参照渡ししている
			//球作成
			set_sphere(X2,Y2,Z2,&number2,le,R,flag);

			for(int i=0;i<number2;i++)
			{
				double r=sqrt(X2[i]*X2[i]+Y2[i]*Y2[i]);
				double a=0.2;
		//		a=500*3;//クーロン分裂のとき
				a=200;	//遊び
				
				
				vx=-0.5*a*X2[i];
				vy=-0.5*a*Y2[i];
				vz=a*Z2[i];
				 
			//	writedata2(fq,  i,X[i],Y[i],Z[i]+Zg,MAGELAST,1, OFF,val, vx,vy, vz,0,0,0,0, 0, 0);
				writedata2(fq,  i+number,X2[i]-Xg,Y2[i],Z2[i]+Zg,ELASTIC,1,1,val, 0,0, 0,0,0,0,0, 0, 0);
			}
			cout<<"球２完成\n";

			number+=number2;

			//for(int i=0;i<number;i++) if(fabs(Z[i])<0.001) writedata2(fq,  i,X[i],Y[i],Z[i]+Zg,FLUID,1, OFF,0, 0,0, 0, 0, 0, ON);

			//下壁
			//double Width=27*le;
			cout<<"床の作成"<<endl;
			//下壁
			double Wr=60*le;//壁の半径
			double Wh=6*le;//壁の高さ
			int ii=0;
			for(int i=0;i<60;i++){
				for(int j=0;j<60;j++){
					for(int k=0;k<3;k++){
						writedata2(fq,ii+number,i*le-(30*le), j*le-(30*le), k*le+6*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
						ii++;
					}
				}
			}
			number+=ii;
			cout<<"model13完成\n";
		}

	 //////////////モデル16　層/////////////////////////////////////////////////////////////////////////
	 else if(model==16)
	 {
		 double h=5;
		 double lm=27;
		 double ls=lm+15;
		 double lw=ls+10;
		 double d=0;
		 double f=0;

		 int num=0;
		 int num2=0;
		 int num3=0;

		 for(int i=0;i<ls;i++)
		 {
			 double ii=abs(i-ls/2);
			 for(int j=0;j<ls;j++)
			 {
				 double jj=abs(j-ls/2);
				 for(int k=0;k<lm;k++)
				 {
					 if(ii<=lm/2 && jj<=h/2)
					 {
						num2++;
						X2.push_back(i*le);
						Y2.push_back(j*le);
						Z2.push_back(k*le);
					 }
					 else
					 {
						num++;
						X.push_back(i*le);
						Y.push_back(j*le);
						Z.push_back(k*le);
					 }
				 }
			 }
		 }

		 for(int i=0;i<num;i++)
			 writedata2(fq,num,X[i]-ls*le/2,Y[i]-ls*le/2,Z[i]+5*le,ELASTIC,1,1,val,0,0,0,0,0,0,0,0,0);
		 for(int i=0;i<num2;i++)
			 writedata2(fq,num2+num,X2[i]-ls*le/2,Y2[i]-ls*le/2,Z2[i]+5*le,MAGELAST,1,1,val,0,0,0,0,0,0,0,0,1);

		 number+=num+num2;
		 cout<<"モデル完成\n";

		 for(int i=0;i<lw;i++)
		 {
			for(int j=0;j<lw;j++)
			{
				for(int k=0;k<6;k++)
				{
					if(i%2==0)
						d=0.5;
					else
						d=-0.5;
					if(j%2==0)
						f=1;
					else
						f=-1;

					num3++;
					writedata2(fq,num3+number,(i-lw/2+d*f)*le,(j-lw/2)*le,(k-6)*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);
				}
			}
		 }
		 number+=num3;
		 cout<<"壁完成\n";
	 }


	////////////////////////モデル17　層２////////////////////////////////////////////////////////////////
	 else if(model==17)
	 {
		 double h=5;
		 double lm=27;
		 double ls=lm+15;
		 double lw=ls+10;
		 double d=0;
		 double f=0;

		 int num=0;
		 int num2=0;
		 int num3=0;

		 for(int i=0;i<ls;i++)
		 {
			 double ii=abs(i-ls/2);
			 for(int j=0;j<ls;j++)
			 {
				 double jj=abs(j-ls/2);
				 for(int k=0;k<lm;k++)
				 {
					 if((ii<=lm/2+h && jj<=h/2) || (ii<=lm/2+h && ii>=lm/2 && jj<=lm/2))
					 {
						num2++;
						X2.push_back(i*le);
						Y2.push_back(j*le);
						Z2.push_back(k*le);
					 }
					 else
					 {
						num++;
						X.push_back(i*le);
						Y.push_back(j*le);
						Z.push_back(k*le);
					 }
				 }
			 }
		 }
		 for(int i=0;i<num;i++)writedata2(fq,num,X[i]-ls*le/2,Y[i]-ls*le/2,Z[i]+3*le,ELASTIC,1,1,val,0,0,0,0,0,0,0,0,0);			 
		 for(int i=0;i<num2;i++)writedata2(fq,num2+num,X2[i]-ls*le/2,Y2[i]-ls*le/2,Z2[i]+3*le,MAGELAST,1,1,val,0,0,0,0,0,0,0,0,1);
		 number+=num+num2;
		 cout<<"モデル完成\n";

		 for(int i=0;i<lw;i++)
		 {
			for(int j=0;j<lw;j++)
			{
				for(int k=0;k<6;k++)
				{
					if(i%2==0)
						d=0.5;
					else
						d=-0.5;
					if(j%2==0)
						f=1;
					else
						f=-1;

					num3++;
					writedata2(fq,num3+number,(i-lw/2+d*f)*le,(j-lw/2)*le,(k-6)*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);
				}
			}
		 }
		 number+=num3;
		 cout<<"壁完成\n";

	 }
	 ////////////////////////////////モデル18　実機////////////////////////////////////////////////
	 else if(model==18)
	 {
		 double h=2;
		 double lm=8.7;
		 double ls=17.3;
		 double flag_model=1;

		 int num=0;
		 int num2=0;

		 for(int i=0;i<ls;i++)
		 {
			 double ii=fabs(i-ls/2);
			 for(int j=0;j<ls;j++)
			 {
				 double jj=fabs(j-ls/2);
				 for(int k=0;k<ls;k++)
				 {
					 double kk=fabs(k-ls/2);

					 if(flag_model==0)
					 {
						 if(ii<=lm/2 && kk<=h/2)
						 {
							num2++;
							X2.push_back(i*le);
							Y2.push_back(j*le);
							Z2.push_back(k*le);
						 }
						 else
						 {
							num++;
							X.push_back(i*le);
							Y.push_back(j*le);
							Z.push_back(k*le);
						 }
					 }

					 if(flag_model==1)
					 {
						 if(ii<=h/2 && kk<=lm/2)
						 {
							num2++;
							X2.push_back(i*le);
							Y2.push_back(j*le);
							Z2.push_back(k*le);
						 }
						 else
						 {
							num++;
							X.push_back(i*le);
							Y.push_back(j*le);
							Z.push_back(k*le);
						 }
					 }

					 if(flag_model==2)
					 {
						 if(ii<=lm/2 && jj<=h/2)
						 {
							num2++;
							X2.push_back(i*le);
							Y2.push_back(j*le);
							Z2.push_back(k*le);
						 }
						 else
						 {
							num++;
							X.push_back(i*le);
							Y.push_back(j*le);
							Z.push_back(k*le);
						 }
					 }
				 }
			 }
		 }

		 for(int i=0;i<num;i++)	writedata2(fq,num,X[i]-ls*le/2,Y[i]-ls*le/2,Z[i]+2*le,ELASTIC,1,1,val,0,0,0,0,0,0,0,0,0);
		 for(int i=0;i<num2;i++)	writedata2(fq,num2+num,X2[i]-ls*le/2,Y2[i]-ls*le/2,Z2[i]+2*le,MAGELAST,1,1,val,0,0,0,0,0,0,0,0,1);			 
		 cout<<"モデル完成\n";
 		 number=num+num2;
		 ///////////////////////下壁作成///////////////////
		 double Rw=14*le;
		 double hw=3*le;
		 int number2=0;
		 int ii=0;
		 set_circle_edge(X3,Y3,Z3,&number2,le,Rw);
		 set_circle_in_using_6_pieces(X3,Y3,Z3,&number2,le,Rw,0,number2);

		 int end_id=number2;
		 int topw_flag=ON;
		 set_cylinder_face(X3,Y3,Z3,&number2,le,Rw,hw,0,end_id,topw_flag);
		 set_cylinder_in(X3,Y3,Z3,&number2,le,Rw,hw,1,0);

		 for(int i=0;i<number2;i++)
		 {
			 if(Z3[i]<le)
			 {
				 ii++;
				 writedata2(fq,ii+number,X3[i],Y3[i],Z3[i]-le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);
			 }
		 }
		 cout<<"壁完成\n";
		 number+=ii;
	 }
	 ///////////////////////////////////////////////////////////////////////////////////////////////////
	 
	 /////////////////////////////////////////////モデル19　表情筋/////////////////////////////////////////////////////
	 else if(model==19)
	 {
		 double Ra=20.0;
		 double Rb=10.0;
		 double hight=5.0;
		 int number2=0;

		 for(int i=0;i<Ra;i++)
		 {
			 for(int j=0;j<Rb;j++)
			 {
				 for(int k=0;k<hight;k++)
				 {
					 X.push_back(i*le);
					 Y.push_back(j*le);
					 Z.push_back(k*le);
					 number++;
				 }
			 }
		 }
		 for(int i=0;i<number;i++)	 writedata2(fq,i,X[i]-(Ra*le+5*le),Y[i]-Rb*le/2,Z[i]+4*le,MAGELAST,1,1,0,0,0,0,0,0,0,0,0,1);
		 cout<<"MRE完成\n";

		 for(int i=0;i<Ra*3;i++)
		 {
			 for(int j=0;j<Rb*1.5;j++)
			 {
				 for(int k=0;k<3;k++)
				 {
					 number2++;
					 writedata2(fq,number+number2,i*le-Ra*3/2*le,j*le-Rb*1.5/2*le,k*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);
				 }
			 }
		 }
		 cout<<"壁完成\n";
		 number+=number2;
	 }

	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	 	 /////////////////////////////////////////////モデル20　表情筋2/////////////////////////////////////////////////////
	 else if(model==20)
	 {
		 double Ra=20.0;
		 double Rb=10.0;
		 double hight=5.0;
		 int number2=0;
		 int number3=0;
		 double f=0;
		 double a=0;

		 for(int i=0;i<Ra;i++)
		 {
			 for(int j=0;j<Rb;j++)
			 {
				 for(int k=0;k<hight;k++)
				 {
					 double f=k-0.75*(i-(Ra/2-5));
					 double a=hight/2-1;

					 if (fabs(f)<=fabs(a))
					 {
						 X2.push_back(i*le);
						 Y2.push_back(j*le);
						 Z2.push_back(k*le);
						 number2++;
					 }
					 else
					 {
					 X.push_back(i*le);
					 Y.push_back(j*le);
					 Z.push_back(k*le);
					 number++;
					 }
				 }
			 }
		 }
		 for(int i=0;i<number;i++)	writedata2(fq,i,X[i]-(Ra*le+5*le),Y[i]-Rb*le/2,Z[i]+4*le,ELASTIC,1,1,0,0,0,0,0,0,0,0,0,0);
		  for(int i=0;i<number2;i++)	 writedata2(fq,i+number,X2[i]-(Ra*le+5*le),Y2[i]-Rb*le/2,Z2[i]+4*le,MAGELAST,1,1,0,0,0,0,0,0,0,0,0,1);
		  number+=number2;

		 cout<<"model完成\n";
		 for(int i=0;i<Ra*3;i++)
		 {
			 for(int j=0;j<Rb*1.5;j++)
			 {
				 for(int k=0;k<3;k++)
				 {
					 number3++;
					 writedata2(fq,number+number3,i*le-Ra*3/2*le,j*le-Rb*1.5/2*le,k*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);
				 }
			 }
		 }
		 cout<<"壁完成\n";
		 number+=number3;
	 }
	 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	 ///////////////////////////////////////////モデル21　超弾性体///////////////////////////////////////////////////////////
	 else if(model==21)	//越塚先生先行研究の角柱
	 {
		 double height=18;
		 double base=3;

		 for(int i=0;i<base;i++)
		 {
			 for(int j=0;j<base;j++)
			 {
				 for(int k=0;k<height;k++)
				 {
					 writedata2(fq,number,(i-(base-1)/2)*le,(j-(base-1)/2)*le,(k-(height-1)/2)*le,HYPERELAST,1,ON,0,0,0,0,0,0,0,0,0,1);
					 number++;
				 }
			 }
		 }
		 cout<<"model完成\n";
	 }

	 ///////////////////////////////////////
	 else if(model==22)	//越塚先生先行研究の円筒型
	 {
		 double height=3;
		 double r1=1;
		 double r2=2;
		 double theta;
		 double x=0,y=0,z=0;

		 int n1=12;
		 int n2=24;
		 number=0;

		 for(int h=0;h<height;h++)
		 {
			 z=h-(height-1)/2;
			 X.push_back(0);
			 Y.push_back(0);
			 Z.push_back(z);
			 number++;
			
			 theta=0;
			 for(int n=0;n<n1;n++)
			 {
				theta=n*(2*PI)/n1;
				x=r1*cos(theta);
				y=r1*sin(theta);
				
				X.push_back(x);
				Y.push_back(y);
				Z.push_back(z);
				number++;
			 }

			 theta=0;
			 for(int n=0;n<n2;n++)
			 {
				 theta=n*(2*PI)/n2;
				 x=r2*cos(theta);
				 y=r2*sin(theta);

				 X.push_back(x);
				 Y.push_back(y);
				 Z.push_back(z);
				 number++;
			 }
		 }

		 for(int i=0;i<number;i++)
		 {
			 writedata2(fq,i,X[i]*le,Y[i]*le,Z[i]*le,HYPERELAST,1,ON,0,0,0,0,0,0,0,0,0,0);
		 }
		 cout<<"model完成\n";
	 }
	 ///////////////////////////////////////
	 /////////////////////////////////ただの立方体////////////////////////////////////
	 else if(model==23)
	 {
		int height=3;
		int base=3;
		vector<int>suf;
		vector<int>num;
		//writedata2内決め打ち有

		double le2=le*sqrt(2.0);

		 for(int k=0;k<height*2-1;k++)
		 {
			 if(k%2==0)
			 {
				 for(int i=0;i<base*2-1;i++)
				 {
					 if(i%2==0)
					 {
						 for(int j=0;j<base;j++)
						 {
							 X.push_back(i*0.5);
							 Y.push_back(j);
							 Z.push_back(k*0.5);
							if(i==0||j==0||k==0||i==base-1||j==base-1||k==height*2-2) suf.push_back(1);
							else suf.push_back(0);
							number++;
						 }
					 }
					 else 
					 {
						 for(int j=0;j<base-1;j++)
						 {
							X.push_back(i*0.5);
							Y.push_back(j+0.5);
							Z.push_back(k*0.5);
							if(i==0||j==0||k==0||i==base-1||j==base-1||k==height*2-2) suf.push_back(1);
							else suf.push_back(0);
							number++;
						 }
					 }
				 }
			 }
			 else
			{
				 for(int i=0;i<base*2-1;i++)
				 {
					 if(i%2==1)
					 {
						 for(int j=0;j<base;j++)
						 {
							 X.push_back(i*0.5);
							 Y.push_back(j);
							 Z.push_back(k*0.5);
							if(i==0||j==0||k==0||i==base-1||j==base-1||k==height*2-2) suf.push_back(1);
							else suf.push_back(0);
							number++;
						 }
					 }
					 else 
					 {
						 for(int j=0;j<base-1;j++)
						 {
							X.push_back(i*0.5);
							Y.push_back(j+0.5);
							Z.push_back(k*0.5);
							if(i==0||j==0||k==0||i==base-1||j==base-1||k==height*2-2) suf.push_back(1);
							else suf.push_back(0);
							number++;
						 }
					 }
				}
			 }
		 }
		 for(int i=0;i<number;i++)	writedata2(fq,i,(X[i]-1.0)*le2,(Y[i]-1.0)*le2,(Z[i]+2.0)*le2,HYPERELAST,1,suf[i],0,0,0,0,0,0,0,0,0,ON);

		 cout<<"超弾性体完成\n";
		 
		 int number2=0;
	 	int number3=0;
		 int w_base=4;
		 vector<int> w_suf;
		 for(int k=0;k<3*2-1;k++)
		 {
			 if(k%2==0)
			 {
				 for(int i=0;i<w_base*2-1;i++)
				 {
					 if(i%2==0)
					 {
						 for(int j=0;j<w_base;j++)
						 {
							 X2.push_back(i*0.5);
							 Y2.push_back(j);
							 Z2.push_back(k*0.5);
							number2++;
							if(i==0||j==0||k==0||i==w_base*2-2||j==w_base||k==3*2-2) w_suf.push_back(1);
							else w_suf.push_back(0);
						 }
					 }
					 else 
					 {
						 for(int j=0;j<w_base-1;j++)
						 {
							X2.push_back(i*0.5);
							Y2.push_back(j+0.5);
							Z2.push_back(k*0.5);
							number2++;
							if(i==0||j==0||k==0||i==w_base*2-2||j==w_base-2||k==3*2-2) w_suf.push_back(1);
							else w_suf.push_back(0);
						 }
					 }
				 }
			 }
			 else
			{
				 for(int i=0;i<w_base*2-1;i++)
				 {
					 if(i%2==1)
					 {
						 for(int j=0;j<w_base;j++)
						 {
							 X2.push_back(i*0.5);
							 Y2.push_back(j);
							 Z2.push_back(k*0.5);
							number2++;
							if(i==0||j==0||k==0||i==w_base*2-2||j==w_base-1||k==3*2-2) w_suf.push_back(1);
							else w_suf.push_back(0);
						 }
					 }
					 else 
					 {
						 for(int j=0;j<w_base-1;j++)
						 {
							X2.push_back(i*0.5);
							Y2.push_back(j+0.5);
							Z2.push_back(k*0.5);
							number2++;
							if(i==0||j==0||k==0||i==w_base*2-2||j==w_base-2||k==3*2-2) w_suf.push_back(1);
							else w_suf.push_back(0);
						 }
					 }
				 }
			 }
		 }
		 for(int i=0;i<number2;i++)		 writedata2(fq,i+number,(X2[i]-1.5)*le2,(Y2[i]-1.5)*le2,(Z2[i]-1.0)*le2,WALL,1,w_suf[i],0,0,0,0,0,0,0,0,0,0);
		 cout<<"number2"<<number2<<endl;
		number+=number2;/*
		for(int k=0;k<3*2-1;k++)
		 {
			if(k%2==0)
			{
				for(int i=0;i<w_base;i++)
				{
					for(int j=0;j<w_base;j++)
					{
						X3.push_back(i);
						Y3.push_back(j);
						Z3.push_back(k*0.5);
						number3++;
					}
				}
			}
		else
		{
			for(int j=0;j<w_base*2-1;j++)
			{
				if(j%2==0)
				{
					for(int i=0;i<w_base;i++)
					{
						if(i+0.5<w_base-1)
						{
							X3.push_back(i+0.5);
							Y3.push_back(j*0.5);
							Z3.push_back(k*0.5);
							number3++;
						}
					}
				}
				else
				{
					for(int i=0;i<w_base;i++)
					{
							X3.push_back(i);
							Y3.push_back(j*0.5);
							Z3.push_back(k*0.5);
							number3++;
						}
					}
				}
			}
		}
 		 for(int i=0;i<number3;i++)		 writedata2(fq,i+number2+number,(X3[i]-2.5)*le,(Y3[i]-2.5)*le,(Z3[i]-1.0+6.0)*le,WALL,1,0,0,0,0,-10.0*le,0,0,0,0,0,0);
		 cout<<"number3"<<number3<<endl;*/
		 		 
		 //number+=number2+number3;
		 cout<<"壁完成\n";
		 cout<<"model完成\n";
	 }
	 /////////////////////////
	 ///////////////////////円筒型磁性エラストマー課題用
	else if(model==24)	
	{
		double R=CON->get_fluidwidth()*le*0.8/2;//半径
		double height=(8*le*A);//(6*le*A)*2;これは何？
		double vz=-2.0;
		if(dimension==3)
		{
			//円作成
			int circle_start_id=0; 
			set_circle_edge(X,Y,Z,&number,le,R);//円外周　これがないと内部も作れない
			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,0,number);//円内部    vector配列は参照渡ししている
//			set_circle_edge(X,Y,Z,&number,le,R, height);//円外周　これがないと内部も作れない
//			set_circle_in_using_6_pieces(X,Y,Z,&number,le,R,height,0,number);//円内部    vector配列は参照渡ししている
			int circle_end_id=number;	//円の粒子idを記憶
			////////

			int top_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,le,R,height,circle_start_id,circle_end_id,top_flag);//円柱表面座標作成
			int face_n=number; //
			///////
			for(int s=0;s<face_n;s++)writedata2(fq,s,X[s],Y[s],Z[s]-height/2,HYPERELAST,1,1,0,0,0,vz,0,0,0,0,0,1);//粒子は,FACE
			///////////////////////////////////////////////////////////////////////////////
			set_cylinder_in(X,Y,Z,&number,le,R,height,1,circle_start_id);//内部にパッキング

			int beforeNumber=number;
			
			//円柱表面+下壁の書き込み
			for(int i=face_n;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i]-height/2,HYPERELAST,1,0,0,0,0,vz,0,0,0,0,0,1);
			//下壁
			//double Width=27*le;
			double Rw=CON->get_fluidwidth()*le*1.8/2;
			double Height=6*le;
			int num3=number;
			double mdis=le; //*0.5
			set_circle_edge(X,Y,Z,&number,mdis,Rw);//円外周　これがないと内部も作れない
			set_circle_in_using_6_pieces(X,Y,Z,&number,mdis,Rw,num3,number);//円内部    vector配列は参照渡ししている　//0をnum3に変更

			int end_id=number;	//円の粒子idを記憶
			////////

			int topw_flag=ON;		//円柱の上面も作成するからフラグをON
			set_cylinder_face(X,Y,Z,&number,mdis,Rw,Height,num3,end_id,topw_flag);//円柱表面座標作成
			set_cylinder_in(X,Y,Z,&number,mdis,Rw,Height,1,num3);//内部にパッキング　　　//最後にint num3(壁の最初の粒子数)を加える．MD_3Dで必要
			int ii=beforeNumber-1;
			for(int i=beforeNumber;i<number;i++)
			{
				if(Z[i]>=Height/2)
				{
					ii++;
					writedata2(fq,ii,X[i], Y[i], Z[i]-1/2*Height-height-5*le,WALL,1,0,0,0,0,0,0,0,0,0,0,0);//粒子はWALL
				}
			}
			number=ii;
			//床上面の位置
			double ground=Z[beforeNumber]-1/2*Height-height-5*le;
			for(int i=beforeNumber;i<number;i++)
			{
				if((Z[i]-1/2*Height-height-5*le)>ground) ground=Z[i]-1/2*Height-height-5*le;
			}
			CON->set_ground_position(ground);
		}
	}


	 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////例外処理///////////////////////////////////////
	else cout<<"modelエラー"<<endl;
	////////////////////////////////////////////////////////////////////////////////////////////////
	fq.close();
		
	*particle_number=number;
	ofstream fn("particle_number.dat");
	fn<<number<<endl;
	fn.close();
	cout<<"ok time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;
}

void writedata2(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double pax,double pay,double paz,double P,double h,int toBEM)
{	
	//ファイル出力
	fp<<id<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<pax<<"\t"<<pay<<"\t"<<paz<<"\t"<<P<<"\t"<<h<<"\t"<<val<<"\t"<<type<<"\t"<<materialID<<"\t"<<surface<<"\t"<<toBEM<<"\t"<<endl;
}

/*void writedata(ofstream &fp, int number, double x, double y,double z, int type,double vx,double vy,double vz,double P,double h,int toFEM)
{	
	//fprintf( fp, "%5.10f\t", 0); とかいう数字の放り込みはだめ？
	
	double angle=0;//回転角
	double angle2=0;
	double angle3=0;
	double angle_s=1;
	double anglar_u=0;//角速度
	double anglar_u2=0;//角速度
	double anglar_u3=0;//角速度
	
	fp<<number<<"\t";
	fp<<x<<"\t";
	fp<<y<<"\t";
	fp<<z<<"\t";
	fp<<vx<<"\t";					//速度x成分
	fp<<vy<<"\t";					//速度y成分
	fp<<vz<<"\t";					//速度z成分
	fp<<P<<"\t";					//圧力
	fp<<h<<"\t";					//エンタルピー
	fp<<angle<<"\t";					//回転角
	fp<<type<<"\t";
	fp<<toFEM<<endl;
}*/


//半径Rの円の外周
void set_circle_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R)
{
	//対称性を考えた時、円周を記述する粒子は偶数個でなければならない。

	int N=calc_division_N_circle(2*PI*R,le);//円周の分割数
	
	double L=2*PI*R/N;				//粒子間距離
	double theta=L/R;				//粒子を配置する角度

	for(int n=0;n<N;n++)
	{
		X.push_back(R*cos(theta*n));
		Y.push_back(R*sin(theta*n));
		Z.push_back(0);
	}
	*number=*number+N;
}

//半径Rの円の外周+z方向ずらす（オーバーロード）
void set_circle_edge(vector<double> &X, vector<double> &Y, vector<double> &Z, int *number, double le,double R, double height)
{
	//対称性を考えた時、円周を記述する粒子は偶数個でなければならない。

	int N=calc_division_N_circle(2*PI*R,le);//円周の分割数
	
	double L=2*PI*R/N;				//粒子間距離
	double theta=L/R;				//粒子を配置する角度

	for(int n=0;n<N;n++)
	{
		X.push_back(R*cos(theta*n));
		Y.push_back(R*sin(theta*n));
		Z.push_back(-height/2);
	}
	*number=*number+N;
}

//円周分割数計算関数
int calc_division_N_circle(double dis,double le)
{
	//対称性を考えた時、円周を記述する粒子は偶数個でなければならない。だから他の辺分割数とは扱いが少し特殊
	//dis:分割する距離(円周)
	double temp_num=dis/le;		//円外周に設置する『仮の』粒子数。ただし外周がうまくleで割り切れるとは限らない

	int N1=(int)(temp_num/2);
	N1*=2;							
	int N2=N1+2;					//temp_numはN1とN2の間にある。ここでN1,N2は偶数

	double dif1=temp_num-N1;		//各Nとの差
	double dif2=N2-temp_num;
	int N=N1;						//周方向分割数
	if(dif2<dif1) N=N2;				//差の小さい方をNとして採用する。

	return N;
}

//半径Rの円内部
void set_circle_in(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &type1,int *number,double le,double R,int edge_startID,int edge_lastID)
{
	//edge_startIDから(edge_lastID-1)までの粒子が、円の外周を構成する粒子に該当する
	//vector型配列は参照渡ししている。vector<double> *Xではなくvector<double> &Xであることに注意。これで各配列は通常通りに仕様可能。アロー演算子もいらない
	//参照渡しでなく通常のやりかたでもいいけど、その場合、例えばa=X[5]と書いても利用できない。a=(*X)[5]などとしなければならない

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(R/le)+1;  //円を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //円を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-le*0.5;				//少し小さめの半径を設定
	
	//初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii*ii+jj*jj<R2*R2)
			{
				X.push_back(ii);
				Y.push_back(jj);
				Z.push_back(0);
				newN++;
			}
		}
	}

	//分子動力学により位置を最適化
	MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);

	*number=*number+newN;
}

//半径Rの円内部
void set_circle_in_using_6_pieces(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int edge_startID,int edge_lastID)
{
	//set_circle_in()と違い、60度分だけ計算し、それを6つｺﾋﾟｰして円を構成する。時間短縮と配置の均等化が目的
	//edge_startIDから(edge_lastID-1)までの粒子が、円の外周を構成する粒子に該当する
	//vector型配列は参照渡ししている。vector<double> *Xではなくvector<double> &Xであることに注意。これで各配列は通常通りに仕様可能。アロー演算子もいらない
	//参照渡しでなく通常のやりかたでもいいけど、その場合、例えばa=X[5]と書いても利用できない。a=(*X)[5]などとしなければならない

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(R/le)+1;  //円を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //円を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-le*0.5;				//少し小さめの半径を設定


	//////////////////////円を6つにわけるための直線6つを生成

	double temp_R_num=R/le;			//半径方向に設置する『仮の』粒子数
	int    R_num=(int)temp_R_num;	//真の粒子数　とりあえず仮の粒子数の整数番とする。ここで、temp_R_num>R_numが成立している。
	double difference=temp_R_num-R_num;	//仮の数と真の数の差
	
	//仮の数と真の数の差が0.5までなら、真の数はNとする。0.5以上ならN+1とする
	if(difference>0.5) R_num++;
	double L=R/R_num;					//粒子間距離

	X.push_back(0);						//中心粒子を追加
	Y.push_back(0);
	Z.push_back(0);
	newN++;

	for(int k=0;k<6;k++)				//6つの直線のloop
	{
		double theta=PI/3*k;			//直線の角度
		for(int n=1;n<R_num;n++)		//中心粒子と最外周粒子はもうあるから、ここでのloopはそれをカウントしない
		{
			double r=L*n;				//中心からの距離
			X.push_back(r*cos(theta));
			Y.push_back(r*sin(theta));
			Z.push_back(0);
			newN++;
		}
	}
	*number=*number+newN;
	newN=0;
	beforeN=*number;
	edge_lastID=beforeN;
	/////////////////直線6つを生成完了

	//内部初期位置 ただし最初の1ピースのみ
	for(int i=0;i<=half_WX;i++)
	{
		for(int j=1;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii*ii+jj*jj<R2*R2)
			{
				if(jj<sqrt(3.0)*ii-le)		//最初のピースの斜め線より低い領域に設置。ただし直線ぎりぎりはまずいので、保険で-leしている
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(0);
					newN++;
				}
			}
		}
	}////////////////////////

	//分子動力学により位置を最適化
	//MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	///内部粒子を周方向に6つｺﾋﾟｰ
	for(int angle=1;angle<6;angle++)
	{
		double theta=PI/3*angle;//回転する角度
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;
			double x=cos(theta)*X[i]-sin(theta)*Y[i];//回転後の座標
			double y=sin(theta)*X[i]+cos(theta)*Y[i];

			X.push_back(x);
			Y.push_back(y);
			Z.push_back(0);
			//newN++;
		}
	}///////////////////*/

	*number=*number+newN*6;//newNはひとつのピース内の粒子数を表しているからここでは6倍
}

//半径Rの円内部(高さ方向オーバーロード)
void set_circle_in_using_6_pieces(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int edge_startID,int edge_lastID)
{
	int newN=0;					//この関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(R/le)+1;  //円を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //円を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-le*0.5;				//少し小さめの半径を設定


	//////////////////////円を6つにわけるための直線6つを生成

	double temp_R_num=R/le;			//半径方向に設置する『仮の』粒子数
	int    R_num=(int)temp_R_num;	//真の粒子数　とりあえず仮の粒子数の整数番とする。ここで、temp_R_num>R_numが成立している。
	double difference=temp_R_num-R_num;	//仮の数と真の数の差
	
	//仮の数と真の数の差が0.5までなら、真の数はNとする。0.5以上ならN+1とする
	if(difference>0.5) R_num++;
	double L=R/R_num;					//粒子間距離

	X.push_back(0);						//中心粒子を追加
	Y.push_back(0);
	Z.push_back(-height/2);
	newN++;

	for(int k=0;k<6;k++)				//6つの直線のloop
	{
		double theta=PI/3*k;			//直線の角度
		for(int n=1;n<R_num;n++)		//中心粒子と最外周粒子はもうあるから、ここでのloopはそれをカウントしない
		{
			double r=L*n;				//中心からの距離
			X.push_back(r*cos(theta));
			Y.push_back(r*sin(theta));
			Z.push_back(-height/2);
			newN++;
		}
	}
	*number=*number+newN;
	newN=0;
	beforeN=*number;
	edge_lastID=beforeN;
	/////////////////直線6つを生成完了

	//内部初期位置 ただし最初の1ピースのみ
	for(int i=0;i<=half_WX;i++)
	{
		for(int j=1;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii*ii+jj*jj<R2*R2)
			{
				if(jj<sqrt(3.0)*ii-le)		//最初のピースの斜め線より低い領域に設置。ただし直線ぎりぎりはまずいので、保険で-leしている
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(-height/2);
					newN++;
				}
			}
		}
	}////////////////////////

	//分子動力学により位置を最適化
	//MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	///内部粒子を周方向に6つコピー
	for(int angle=1;angle<6;angle++)
	{
		double theta=PI/3*angle;//回転する角度
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;
			double x=cos(theta)*X[i]-sin(theta)*Y[i];//回転後の座標
			double y=sin(theta)*X[i]+cos(theta)*Y[i];

			X.push_back(x);
			Y.push_back(y);
			Z.push_back(-height/2);
			//newN++;
		}
	}///////////////////*/

	*number=*number+newN*6;//newNはひとつのピース内の粒子数を表しているからここでは6倍
}

//半径Rの球作成関数
void set_sphere(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag)
{
	//まずは半球を作る。そのためには半球表面を作成する必要がある。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数
	int beforeN2=*number;		//関数呼び出し時の粒子数。この関数の最後まで記憶しておく

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int half_WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int half_WZ=(int)(R/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-0.5*le;				//少し小さめの半径を設定

	///////////半球表面
	int Nt;						//球表面の、θ方向の分割数
	double Lt;					//球表面の、θ方向の分割距離
	calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//半円の分割数は偶数・奇数どちらでもよい
	double d_theta=Lt/R;		//弧の長さがLtになる角度

	for(int k=0;k<Nt;k++)//loopはk<Ntで終わらせる。Ntに該当するところはすでに設置済み
	{
		double THETA=k*d_theta;	//θ
		double r=R*sin(THETA);	//その高さにおける円の半径
		double round=2*PI*r;//その高さにおける円周

		int Nf=calc_division_N_circle(round,le);//球表面の、θ方向の分割数
		double Lf=round/Nf;						//球表面の、θ方向の分割距離
		double d_fai=Lf/r;						//弧の長さがLfになる角度
		
		for(int i=0;i<Nf;i++)
		{
			double fai=d_fai*i;
			if(Nt%2==0)
			{
				if(k%2!=0) fai+=0.5*d_fai;//Ntが偶数なら、作成済みの円と接するときは奇数番目。よって奇数をずらす
			}
			else
			{
				if(k%2==0) fai+=0.5*d_fai;//Ntが奇数なら、作成済みの円と接するときは偶数番目。よって奇数をずらす
			}
			double x=r*cos(fai);
			double y=r*sin(fai);
			double z=R*cos(THETA);
			X.push_back(x);
			Y.push_back(y);
			Z.push_back(z);
			
			newN++;
		}
	}
	if(Nt%2!=0)//Ntが奇数のときは、頂上に粒子が置かれなければならない。しかし上のloopはそれが不可能。よってここで追加
	{
		X.push_back(0);
		Y.push_back(0);
		Z.push_back(R);
		newN++;
	}
	//////////////////////////////////////

	*number=*number+newN;

	if(flag!=HALFD_shell)
	{
	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=*number;		//この関数呼び出し時における粒子数
	
	//内部初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(ii*ii+jj*jj+kk*kk<R2*R2*0.6*0.6)//半径*0.7内の粒子は位置固定
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;

	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=*number;		//この関数呼び出し時における粒子数
	///////////////////////*/

	

	//内部初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				//if(ii*ii+jj*jj+kk*kk<R2*R2)
				if(ii*ii+jj*jj+kk*kk<R2*R2 && ii*ii+jj*jj+kk*kk>=R2*R2*0.6*0.6)
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}///////////////////////*/
	

	//分子動力学により位置を最適化
	double r=1.50;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-0.2*R; rigion[A_Z][1]=1.2*R;

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,rigion);

	*number=*number+newN;
	}

	///上半球を下半球へｺﾋﾟｰし球を完成
	if(flag==FULL)
	{
		newN=0;					//新しく追加する粒子数
		beforeN=*number;		//この時における粒子数

		for(int k=0;k<beforeN;k++)
		{
			if(Z[k]>0.4*le)
			{
				X.push_back(X[k]);
				Y.push_back(Y[k]);
				Z.push_back(-Z[k]);
				newN++;
			}
		}
		*number=*number+newN;
	}///////////////////////////////
	if(flag==HALFD || flag==HALFD_shell)		//下半球がほしいときに、つくった上半球を上下反転させる
	{
		for(int k=beforeN2;k<*number;k++) Z[k]*=-1;
	}
	
}

void set_sphere2(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,int flag,int *suf_num)
{
	//まずは半球を作る。そのためには半球表面を作成する必要がある。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数
	int beforeN2=*number;		//関数呼び出し時の粒子数。この関数の最後まで記憶しておく

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int half_WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int half_WZ=(int)(R/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-0.5*le;				//少し小さめの半径を設定

	///////////半球表面
	int Nt;						//球表面の、θ方向の分割数
	double Lt;					//球表面の、θ方向の分割距離
	calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//半円の分割数は偶数・奇数どちらでもよい
	double d_theta=Lt/R;		//弧の長さがLtになる角度

	for(int k=0;k<Nt;k++)//loopはk<Ntで終わらせる。Ntに該当するところはすでに設置済み
	{
		double THETA=k*d_theta;	//θ
		double r=R*sin(THETA);	//その高さにおける円の半径
		double round=2*PI*r;//その高さにおける円周

		int Nf=calc_division_N_circle(round,le);//球表面の、θ方向の分割数
		double Lf=round/Nf;						//球表面の、θ方向の分割距離
		double d_fai=Lf/r;						//弧の長さがLfになる角度
		
		for(int i=0;i<Nf;i++)
		{
			double fai=d_fai*i;
			if(Nt%2==0)
			{
				if(k%2!=0) fai+=0.5*d_fai;//Ntが偶数なら、作成済みの円と接するときは奇数番目。よって奇数をずらす
			}
			else
			{
				if(k%2==0) fai+=0.5*d_fai;//Ntが奇数なら、作成済みの円と接するときは偶数番目。よって奇数をずらす
			}
			double x=r*cos(fai);
			double y=r*sin(fai);
			double z=R*cos(THETA);
			X.push_back(x);
			Y.push_back(y);
			Z.push_back(z);
			
			newN++;
		}
	}
	if(Nt%2!=0)//Ntが奇数のときは、頂上に粒子が置かれなければならない。しかし上のloopはそれが不可能。よってここで追加
	{
		X.push_back(0);
		Y.push_back(0);
		Z.push_back(R);
		newN++;
	}
	//////////////////////////////////////

	*number=*number+newN;
	*suf_num=*suf_num+newN;

	if(flag!=HALFD_shell)
	{
	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=*number;		//この関数呼び出し時における粒子数
	
	//内部初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(ii*ii+jj*jj+kk*kk<R2*R2*0.6*0.6)//半径*0.7内の粒子は位置固定
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;

	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=*number;		//この関数呼び出し時における粒子数
	///////////////////////*/

	

	//内部初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			for(int k=1;k<half_WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				//if(ii*ii+jj*jj+kk*kk<R2*R2)
				if(ii*ii+jj*jj+kk*kk<R2*R2 && ii*ii+jj*jj+kk*kk>=R2*R2*0.6*0.6)
				{
					X.push_back(ii);
					Y.push_back(jj);
					Z.push_back(kk);
					newN++;
				}
			}
		}
	}///////////////////////*/
	

	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-0.2*R; rigion[A_Z][1]=1.2*R;

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,rigion);

	*number=*number+newN;
	}

	///上半球を下半球へｺﾋﾟｰし球を完成
	if(flag==FULL)
	{
		newN=0;					//新しく追加する粒子数
		beforeN=*number;		//この時における粒子数

		for(int k=0;k<beforeN;k++)
		{
			if(Z[k]>0.4*le)
			{
				X.push_back(X[k]);
				Y.push_back(Y[k]);
				Z.push_back(-Z[k]);
				newN++;
			}
		}
		*number=*number+newN;
	}///////////////////////////////
	if(flag==HALFD || flag==HALFD_shell)		//下半球がほしいときに、つくった上半球を上下反転させる
	{
		for(int k=beforeN2;k<*number;k++) Z[k]*=-1;
	}
	
}


//直線を分割するさいの最適な分割数と分割距離の算出関数
void calc_N_and_L(double dis,double le,int *N,double *L)
{
	double temp_N=dis/le;			//仮の分割数。leで割り切れたら一番いいけど、そうもいかないときがある
	int Ns=(int) temp_N;				//真の分割数
	double difference=temp_N-Ns;		//仮と真の差
	if(difference>0.5) Ns++;
	*L=dis/Ns;			//粒子の距離
	*N=Ns;
}

//分子動力学関数
void MD_2D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN)
{
	//分子動力学によりnewN個の粒子の位置を最適化　IDがBstartIDからBendIDまでのは境界粒子なので動かさない

	double region[2][2];	//解析領域

	/////////////////////解析領域の決定
	region[A_X][0]=100; region[A_X][1]=-100;
	region[A_Y][0]=100; region[A_Y][1]=-100;
	for(int i=BstartID;i<BendID;i++)
	{
		if(X[i]<region[A_X][0]) region[A_X][0]=X[i];
		else if(X[i]>region[A_X][1]) region[A_X][1]=X[i];

		if(Y[i]<region[A_Y][0]) region[A_Y][0]=Y[i];
		else if(Y[i]>region[A_Y][1]) region[A_Y][1]=Y[i];
	}
	for(int D=0;D<2;D++)
	{
		region[D][0]-=5*le;	//少し領域を広めにとる
		region[D][1]+=5*le;
	}//////////////////////////

	//パラメータ
	double k0=1;
	double r=1.5;
	double dt=0.001;
	
	//力はax^3+bx^2+dの式を採用。文献[Bubble Mesh Automated Triangular Meshing of Non-Manifold Geometry by Sphere Packing]を参照
	double a=(r+1)/(2*r*r-r-1)*k0/(le*le);
	double b=-0.5*k0/le-1.5*a*le;
	double d=-a*le*le*le-b*le*le;
	/////////////
	int lastN=beforeN+newN;

	vector<double> Fx(newN);	//各粒子に働くX方向力
	vector<double> Fy(newN);	//各粒子に働くY方向力
	vector<double> Ax(newN,0);	//X方向加速度
	vector<double> Ay(newN,0);	//Y方向加速度
	vector<double> U(newN,0);	//X方向速度
	vector<double> V(newN,0);	//Y方向速度
	vector<double> visX(newN);	//X方向粘性係数
	vector<double> visY(newN);	//Y方向粘性係数

	//計算の高速化のために格子を形成 解析幅がr*leで割り切れるとは限らないので、はみ出したところは切り捨て。なので各軸とも正の方向には余裕を持つこと
	double grid_width=le*((int)(r+1));								//格子の幅。rを含む整数*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X方向の格子の個数
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//各内部粒子を含む格子番号
	vector<int> *MESH=new vector<int>[plane_SIZE];				//各メッシュに格納される粒子ID格納

	for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);	//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);	//Y方向に何個目の格子か
		int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
	{
		int i=beforeN+k;
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
		MESH[number].push_back(i);
		index[k]=number;
	}//////////////////////////////////////////

	
	//計算開始
	for(int t=0;t<100;t++)
	{
		if(t%10==0 &&t>0)//MESHを作り直す
		{
			//まずはMESHを一度破壊する。
			for(int n=0;n<plane_SIZE;n++)
			{
				size_t size=MESH[n].size();
				for(int k=0;k<size;k++) MESH[n].pop_back();
			}
			
			for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
			{
				int xn=(int)((X[i]-region[A_X][0])/grid_width);	//X方向に何個目の格子か 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);	//Y方向に何個目の格子か
				int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
			{
				int i=beforeN+k;
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
				int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
				MESH[number].push_back(i);
				index[k]=number;
			}
		}////////////

		for(int k=0;k<newN;k++)
		{
			Fx[k]=0; Fy[k]=0;					//初期化
			int i=beforeN+k;					//対応する粒子番号
			double kx=0;						//X方向バネ係数
			double ky=0;
			int G_id=index[k];				//格納する格子番号
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					int M_id=II+JJ;
					for(int L=0;L<MESH[M_id].size();L++)
					{
						int j=MESH[M_id][L];
						double x=X[j]-X[i];
						double y=Y[j]-Y[i];
						double dis=sqrt(x*x+y*y);
						if(dis<r*le && dis!=0)			//このloopは自分自身も通過するから、dis!=0は必要
						{
							double F=a*dis*dis*dis+b*dis*dis+d;
							Fx[k]-=F*x/dis;					//Fの値が正のときは斥力なので、-=にする
							Fy[k]-=F*y/dis;
							double K=3*a*dis*dis+2*b*dis;//バネ係数　力の式の微分に相当
							K=sqrt(K*K);					//正の値が欲しい。だから負のときに備えて正に変換
							kx+=K*x*x/(dis*dis);			//kを各方向に分配。ここで、常に正の量が分配されるようにx*x/(dis*dis)となっている
							ky+=K*y*y/(dis*dis);
						}
					}
				}
			}
			visX[k]=1.414*sqrt(kx);//このように各軸方向の粘性係数を決める。文献「物理モデルによる自動メッシュ分割」P6参照。ただし質量は1としている。
			visY[k]=1.414*sqrt(ky);
			Ax[k]=(Fx[k]-visX[k]*U[k]);
			Ay[k]=(Fy[k]-visY[k]*V[k]);
		}//各粒子の加速度が求まった。
		
		if(t==0)	//最初のｽﾃｯﾌﾟ時にdtを決定
		{
			double MaxAccel=0;
			for(int k=0;k<newN;k++)
			{
				double accel2=Ax[k]*Ax[k]+Ay[k]*Ay[k];
				if(accel2>MaxAccel) MaxAccel=accel2;
			}
			MaxAccel=sqrt(MaxAccel);//最大加速度が求まった
			dt=sqrt(0.02*le/MaxAccel);
		}

		for(int k=0;k<newN;k++)//速度と位置の更新
		{
			int i=beforeN+k;
			double u=U[k];
			double v=V[k];
			U[k]+=dt*Ax[k];
			V[k]+=dt*Ay[k];
			X[i]+=dt*(U[k]+u)*0.5;
			Y[i]+=dt*(V[k]+v)*0.5;
		}

		//再近接距離がle以下の場合はこれを修正
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//対応する粒子番号
			int G_id=index[k];				//格納する格子番号
			double mindis=le;
			int J=k;						//最近接距離の相手粒子
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					int M_id=II+JJ;
					for(int L=0;L<MESH[M_id].size();L++)
					{
						int j=MESH[M_id][L];
						double x=X[j]-X[i];
						double y=Y[j]-Y[i];
						double dis=sqrt(x*x+y*y);
						if(dis<mindis && i!=j)
						{
							mindis=dis;
							J=j;
						}
					}
				}
			}
			if(J!=i && J<beforeN)//leより近接している相手が境界粒子なら
			{
				double L=le-mindis;//開くべき距離
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;			
			}
			else if(J!=i && J>=beforeN)//leより近接している相手が内部粒子なら
			{
				double L=0.5*(le-mindis);//開くべき距離
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;
				X[J]+=dX/mindis*L;
				Y[J]+=dY/mindis*L;
			}
		}//////////*/
	}/////MD終了

	delete [] index;
	delete [] MESH;
}

//分子動力学関数
void MD_3D(vector<double> &X,vector<double> &Y,vector<double> &Z,double le,int BstartID,int BendID,int beforeN,int newN,double r,double region[3][2])
{
	//分子動力学によりnewN個の粒子の位置を最適化　IDがBstartIDからBendIDまでのは境界粒子なので動かさない
	double k0=1;
	double dt=0.001;
	
	//力はax^3+bx^2+dの式を採用。文献[Bubble Mesh Automated Triangular Meshing of Non-Manifold Geometry by Sphere Packing]を参照
	double a=(r+1)/(2*r*r-r-1)*k0/(le*le);
	double b=-0.5*k0/le-1.5*a*le;
	double d=-a*le*le*le-b*le*le;
	/////////////
	int lastN=beforeN+newN;

	//cout<<"F="<<a*le*le*le+b*le*le+d<<" "<<a*1.5*le*1.5*le*1.5*le+b*1.5*le*1.5*le+d<<endl;

	vector<double> Fx(newN);	//各粒子に働くX方向力
	vector<double> Fy(newN);	//各粒子に働くY方向力
	vector<double> Fz(newN);	//各粒子に働くZ方向力
	vector<double> Ax(newN,0);	//X方向加速度
	vector<double> Ay(newN,0);	//Y方向加速度
	vector<double> Az(newN,0);	//Z方向加速度
	vector<double> U(newN,0);	//X方向速度
	vector<double> V(newN,0);	//Y方向速度
	vector<double> W(newN,0);	//Z方向速度
	vector<double> visX(newN);	//X方向粘性係数
	vector<double> visY(newN);	//Y方向粘性係数
	vector<double> visZ(newN);	//Y方向粘性係数
	vector<double> KX(newN);	//X方向バネ係数
	vector<double> KY(newN);	//Y方向バネ係数
	vector<double> KZ(newN);	//Y方向バネ係数

	//計算の高速化のために格子を形成 解析幅がr*leで割り切れるとは限らないので、はみ出したところは切り捨て。なので各軸とも正の方向には余裕を持つこと
	double grid_width=le*((int)(r));								//格子の幅。rを含む整数*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X方向の格子の個数
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int grid_sizeZ=(int)((region[A_Z][1]-region[A_Z][0])/grid_width);
	int grid_SIZE=grid_sizeX*grid_sizeY*grid_sizeZ;
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//各内部粒子を含む格子番号
	vector<int> *MESH=new vector<int>[grid_SIZE];				//各メッシュに格納される粒子ID格納
//	cout<<"grid="<<grid_SIZE<<endl;
	for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
//		cout<<number<<endl;
//		if(X[i]<region[A_X][0] || X[i]>region[A_X][1])cout<<"Xsize="<<X[i]<<endl;
//		if(Y[i]<region[A_Y][0] || Y[i]>region[A_Y][1])cout<<"Ysize="<<Y[i]<<endl;
//		if(Z[i]<region[A_Z][0] || Z[i]>region[A_Z][1])cout<<"Zsize="<<Z[i]<<endl;
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
	{
		int i=beforeN+k;
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
		MESH[number].push_back(i);
		index[k]=number;
	}

	//計算開始
	for(int t=0;t<100;t++)
	{
		if(t%10==0 && t>0)
		{
			//MESHを一度破壊する。
			for(int n=0;n<grid_SIZE;n++)
			{
				size_t size=MESH[n].size();
				for(int k=0;k<size;k++) MESH[n].pop_back();
			}
			
			for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
			{
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
				int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
			{
				int i=beforeN+k;
				int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
				int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
				MESH[number].push_back(i);
				index[k]=number;
			}
		}

		for(int k=0;k<newN;k++)
		{
			Fx[k]=0; Fy[k]=0, Fz[k]=0;			//初期化
			KX[k]=0;KY[k]=0; KZ[k]=0;			//バネ係数
		}

		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//対応する粒子番号
			int G_id=index[k];				//格納する格子番号
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
					{
						int M_id=II+JJ+KK;
						for(int L=0;L<MESH[M_id].size();L++)
						{
							int j=MESH[M_id][L];
							if(j>=beforeN && j>i)	//同じ領域内でかつiより大きな番号なら
							{
								int J=j-beforeN;	//newN内での番号
								double x=X[j]-X[i];
								double y=Y[j]-Y[i];
								double z=Z[j]-Z[i];
								double dis=sqrt(x*x+y*y+z*z);
								if(dis<r*le)			//このloopは自分自身も通過するから、dis!=0は必要
								{
									double F=a*dis*dis*dis+b*dis*dis+d;
									Fx[k]-=F*x/dis;					//Fの値が正のときは斥力なので、-=にする
									Fy[k]-=F*y/dis;
									Fz[k]-=F*z/dis;
									Fx[J]+=F*x/dis;					//相手粒子の力もここで計算。符号は反転させる
									Fy[J]+=F*y/dis;
									Fz[J]+=F*z/dis;
									double K=3*a*dis*dis+2*b*dis;//バネ係数　力の式の微分に相当
									K=sqrt(K*K);					//正の値が欲しい。だから負のときに備えて正に変換
									KX[k]+=K*x*x/(dis*dis);			//kを各方向に分配。ここで、常に正の量が分配されるようにx*x/(dis*dis)となっている
									KY[k]+=K*y*y/(dis*dis);
									KZ[k]+=K*z*z/(dis*dis);
									KX[J]+=K*x*x/(dis*dis);			//kを相手粒子にも分配
									KY[J]+=K*y*y/(dis*dis);
									KZ[J]+=K*z*z/(dis*dis);
								}
							}
							if(j<BendID && j>=BstartID)
							{
								double x=X[j]-X[i];
								double y=Y[j]-Y[i];
								double z=Z[j]-Z[i];
								double dis=sqrt(x*x+y*y+z*z);
								if(dis<r*le && dis>0)			//このloopは自分自身は通過しない、dis!=0は不要
								{
									double F=a*dis*dis*dis+b*dis*dis+d;
									Fx[k]-=F*x/dis;					//Fの値が正のときは斥力なので、-=にする
									Fy[k]-=F*y/dis;
									Fz[k]-=F*z/dis;
									double K=3*a*dis*dis+2*b*dis;//バネ係数　力の式の微分に相当
									K=sqrt(K*K);					//正の値が欲しい。だから負のときに備えて正に変換
									KX[k]+=K*x*x/(dis*dis);			//kを各方向に分配。ここで、常に正の量が分配されるようにx*x/(dis*dis)となっている
									KY[k]+=K*y*y/(dis*dis);
									KZ[k]+=K*z*z/(dis*dis);
								}
							}
						}
					}
				}
			}
			//visX[k]=1.414*sqrt(KX[k]);//このように各軸方向の粘性係数を決める。文献「物理モデルによる自動メッシュ分割」P6参照。ただし質量は1としている。
			//visY[k]=1.414*sqrt(KY[k]);
			//visZ[k]=1.414*sqrt(KZ[k]);
			visX[k]=1.414*sqrt(KX[k]);//このように各軸方向の粘性係数を決める。文献「物理モデルによる自動メッシュ分割」P6参照。ただし質量は1としている。
			visY[k]=1.414*sqrt(KY[k]);
			visZ[k]=1.414*sqrt(KZ[k]);
			Ax[k]=(Fx[k]-visX[k]*U[k]);
			Ay[k]=(Fy[k]-visY[k]*V[k]);
			Az[k]=(Fz[k]-visZ[k]*W[k]);
		}//各粒子の加速度が求まった。
		
		if(t==0)	//最初のｽﾃｯﾌﾟ時にdtを決定
		{
			double MaxAccel=0;
			for(int k=0;k<newN;k++)
			{
				double accel2=Ax[k]*Ax[k]+Ay[k]*Ay[k]+Az[k]*Az[k];
				if(accel2>MaxAccel) MaxAccel=accel2;
			}
			MaxAccel=sqrt(MaxAccel);//最大加速度が求まった
			if(MaxAccel!=0)
			{
				dt=sqrt(0.02*le/MaxAccel);
			}
		}

		for(int k=0;k<newN;k++)//速度と位置の更新
		{
			int i=beforeN+k;
			double u=U[k];
			double v=V[k];
			double w=W[k];
			U[k]+=dt*Ax[k];
			V[k]+=dt*Ay[k];
			W[k]+=dt*Az[k];
			X[i]+=dt*(U[k]+u)*0.5;
			Y[i]+=dt*(V[k]+v)*0.5;
			Z[i]+=dt*(W[k]+w)*0.5;
		}

		//再近接距離がle以下の場合はこれを修正
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;					//対応する粒子番号
			int G_id=index[k];				//格納する格子番号
			double mindis=le;
			int J=k;						//最近接距離の相手粒子
			for(int II=G_id-1;II<=G_id+1;II++)
			{       
				for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
				{
					for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
					{
						int M_id=II+JJ+KK;
						for(int L=0;L<MESH[M_id].size();L++)
						{
							int j=MESH[M_id][L];
							double x=X[j]-X[i];
							double y=Y[j]-Y[i];
							double z=Z[j]-Z[i];
							double dis=sqrt(x*x+y*y+z*z);
							if(dis<mindis && i!=j)
							{
								mindis=dis;
								J=j;
							}
						}
					}
				}
			}
			if(J!=i && J<beforeN)//leより近接している相手が境界粒子なら
			{
				double L=le-mindis;//開くべき距離
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				double dZ=Z[J]-Z[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;	
				Z[i]-=dZ/mindis*L;	
			}
			else if(J!=i && J>=beforeN)//leより近接している相手が内部粒子なら
			{
				double L=0.5*(le-mindis);//開くべき距離
				double dX=X[J]-X[i];
				double dY=Y[J]-Y[i];
				double dZ=Z[J]-Z[i];
				X[i]-=dX/mindis*L;
				Y[i]-=dY/mindis*L;
				Z[i]-=dZ/mindis*L;
				X[J]+=dX/mindis*L;
				Y[J]+=dY/mindis*L;
				Z[J]+=dZ/mindis*L;
			}
		}//////////*/
	}/////MD終了

	delete [] index;
	delete [] MESH;
}


//長方形の辺作成関数
void set_rectangular_edge(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V)
{
	//水平方向長さLen_H,垂直長さLen_Vの長方形の辺を作成する
	
	int newN=0;					//この関数で追加される粒子数
	
	////////////////////まず水平方向作成

	int Nh;				//水平辺の分割数
	double dL_H;		//水平辺の分割長さ
	calc_N_and_L(Len_H,le,&Nh,&dL_H);

	for(int n=0;n<=Nh;n++)//底辺 ここではY,Z座標はゼロとしておく。
	{
		X.push_back(-Len_H*0.5+dL_H*n);
		Y.push_back(-Len_V*0.5);
		Z.push_back(0);
		newN++;
	}
	for(int n=0;n<=Nh;n++)//上辺
	{
		X.push_back(-Len_H*0.5+dL_H*n);
		Y.push_back(Len_V*0.5);
		Z.push_back(0);
		newN++;
	}////////////////////////////////////////////

	////////////////////次に垂直方向作成

	int Nv;				//水直辺の分割数
	double dL_V;		//水直辺の分割長さ
	calc_N_and_L(Len_V,le,&Nv,&dL_V);

	for(int n=1;n<Nv;n++)//左辺 ここではX座標はゼロとしておく。またn=0に該当する点はすでに設置済み
	{
		X.push_back(-Len_H*0.5);
		Y.push_back(-Len_V*0.5+n*dL_V);
		Z.push_back(0);
		newN++;
	}
	for(int n=1;n<Nv;n++)//右辺 
	{
		X.push_back(Len_H*0.5);
		Y.push_back(-Len_V*0.5+n*dL_V);
		Z.push_back(0);
		newN++;
	}////////////////////////////////////////////

	*number=*number+newN;
}

//長方形作成関数
void set_rectangular_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Len_H,double Len_V,int edge_startID,int edge_lastID)
{
	//水平方向長さLen_H,垂直長さLen_Vの長方形の辺を作成する
	
	int newN=0;					//この関数で追加される粒子数
	int beforeN=*number;			//この関数呼び出し時にすでに存在している粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(Len_H/le)+1;  //長方形を十分含む幅
	int half_WY=(int)(Len_V/(le*A))+1;
	double gap=0.4*le;				//辺ぎりぎりに内部粒子を配置しないよう、隙間を設ける
	
	//初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WX;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii>-Len_H*0.5+gap && ii<Len_H*0.5-gap)
			{
				if(jj>-Len_V*0.5+gap && jj<Len_V*0.5-gap)
				{
					if(ii*ii<Len_H*0.25*Len_H*0.25 && jj*jj<Len_V*0.25*Len_V*0.25)//十分内部は配置を決め打ち
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(0);
						newN++;
					}
				}
			}
		}
	}
	*number=*number+newN;
	newN=0;					
	beforeN=*number;
	edge_lastID=*number;//境界IDを変更　ここまでに設置された粒子が境界粒子になる。

	//初期位置
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WX;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii>-Len_H*0.5+gap && ii<Len_H*0.5-gap)
			{
				if(jj>-Len_V*0.5+gap && jj<Len_V*0.5-gap)
				{
					if(ii*ii>=Len_H*0.25*Len_H*0.25 || jj*jj>=Len_V*0.25*Len_V*0.25)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(0);
						newN++;
					}
				}
			}
		}
	}

	//分子動力学により位置を最適化
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	*number=*number+newN;
}

//円柱表面作成関数
void set_cylinder_face(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int circle_start_id,int circle_end_id,int top_flag)
{
	//半径R,高さheightの円柱の面を作成する。ただしこの関数呼び出し時において、すでに下面の円(Z=0)は作成済みとする
	//top_flag=ONなら円柱上面を作成する。OFFならしないが、側面だけは作成する。
	int beforeN=*number;
	int newN=0;

	int Nv;				//水直の分割数
	double dL_V;		//水直の分割長さ
	double A=sqrt(3.0)/2;		//よく使う係数
	calc_N_and_L(height,le*A,&Nv,&dL_V);

	int Nr=calc_division_N_circle(2*PI*R,le);//円周の分割数
	double Lr=2*PI*R/Nr;				//円周分割距離

	double gap=0.4*le;				//辺ぎりぎりに内部粒子を配置しないよう、隙間を設ける

	///////////////////////////////////側面
	for(int i=0;i<Nr;i++)
	{
		for(int j=1;j<Nv;j++)//j=0,j=Nvは下面、上面に該当するのでここではぬかす
		{
			double jj=j*dL_V;
			double ii=i*Lr;
			if(j%2!=0) ii+=0.5*Lr;//jが奇数ならiiを0.5格子だけずらす
			if(ii<2*PI*R-gap)
			{
				if(jj<height-gap)	
				{
					double theta=2*PI*(ii/(2*PI*R));
					X.push_back(R*cos(theta));
					Y.push_back(R*sin(theta));
					Z.push_back(jj);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;
	////////////////////////


	//上面作成（下面のコピー）ただしNvが奇数なら上面は下面と半格子ずれなければならない
	beforeN=*number;
	newN=0;
	if(Nv%2==0)	//偶数ならそのままｺﾋﾟｰ
	{
		if(top_flag==ON)
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				X.push_back(X[i]);
				Y.push_back(Y[i]);
				Z.push_back(height);
				newN++;
			}
		}
		else if(top_flag==HALF)//内壁部より内側はなし
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R && r<R+4*le-0.1*le)//外周のみ作成
				{
				X.push_back(X[i]);
				Y.push_back(Y[i]);
				Z.push_back(height);
				newN++;
				}
			}
		}
		else if(top_flag==OFF)//上面が必要ないなら
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R-0.1*le)//外周のみ作成
				{
					X.push_back(X[i]);
					Y.push_back(Y[i]);
					Z.push_back(height);
					newN++;
				}
			}
		}
	}
	else
	{
		double d_theta=0.5*Lr/R;//この微小角度だけ回転させる。
		if(top_flag==ON)
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//回転後の座標
				double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
				X.push_back(X2);
				Y.push_back(Y2);
				Z.push_back(height);
				newN++;
			}
		}
		else if(top_flag==HALF)//上面が必要ないなら
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				//if(r>R-0.1*le)//外周のみ作成
				if(r>R && r<R+4*le-0.1*le)//外周のみ作成
				{
					double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//回転後の座標
					double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
					X.push_back(X2);
					Y.push_back(Y2);
					Z.push_back(height);
					newN++;
				}
			}
		}
		else if(top_flag==OFF)//上面が必要ないなら
		{
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double r=sqrt(X[i]*X[i]+Y[i]*Y[i]);
				if(r>R-0.1*le)//外周のみ作成
				{
					double X2=X[i]*cos(d_theta)-Y[i]*sin(d_theta);//回転後の座標
					double Y2=X[i]*sin(d_theta)+Y[i]*cos(d_theta);
					X.push_back(X2);
					Y.push_back(Y2);
					Z.push_back(height);
					newN++;
				}
			}
		}
	}
	*number=*number+newN;
	
	/////////////////////////////

}

//円柱内部設置関数
void set_cylinder_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double height,int flag,int stert_id)
{
	//半径R,高さheightの円柱内部を作成する。この関数呼び出し時に0<=i<numberの粒子で円柱表面が形成されているとする。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int WZ=(int)(height/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-0.3*le;				//少し小さめの半径を設定
	double height2=height-0.3*le;				//少し小さめの高さを設定

	//内部初期位置
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(ii*ii+jj*jj<R2*R2)
				{
					if(kk<height2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/

	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,stert_id,beforeN,beforeN,newN,r,rigion);  //0をstert_idに変更

	*number=*number+newN;

}
void set_cylinder_in(vector<double> &X,vector<double> &Y,vector<double> &Z,double erast_r,double erast_h,int *number,double le,double R,double height,int flag,int stert_id)
{
	//半径R,高さheightの円柱内部を作成する。この関数呼び出し時に0<=i<numberの粒子で円柱表面が形成されているとする。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int WZ=(int)(height/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-0.3*le;				//少し小さめの半径を設定
	double height2=height-0.3*le;				//少し小さめの高さを設定

	//内部初期位置
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(ii*ii+jj*jj<R2*R2)
				{
					if(kk<height2)
					{
						if(!(ii*ii+jj*jj<erast_r*erast_r && (kk<erast_h+0.005 && kk>0.005))){
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
						}
					}
				}
			}
		}
	}///////////////////////*/

	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R; rigion[A_X][1]=1.2*R;
	rigion[A_Y][0]=-1.2*R; rigion[A_Y][1]=1.2*R;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,stert_id,beforeN,beforeN,newN,r,rigion);  //0をstert_idに変更

	*number=*number+newN;

}
void set_crucible_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R,double R_out,double height,double height_out,int flag,int fluid_number)
{
	//半径R,高さheightの円柱内部を作成する。この関数呼び出し時に0<=i<numberの粒子で円柱表面が形成されているとする。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int WX=(int)(R_out/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int WY=(int)(R_out/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int WZ=(int)(height/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R+0.3*le;				//少し小さめの半径を設定
	double R2_out=R_out-0.3*le;				//少し小さめの半径を設定
	double height2=height-0.3*le;				//少し小さめの高さを設定

	//内部初期位置
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=-WZ;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk>=0 && kk<height2)//内壁円柱部
				{
					if(ii*ii+jj*jj<R2_out*R2_out &&ii*ii+jj*jj>R2*R2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				else if(kk<0 && kk>-R2_out)//内壁円柱部
				{
					if(ii*ii+jj*jj+kk*kk<R2_out*R2_out &&ii*ii+jj*jj+kk*kk>R2*R2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				
			}
		}
	}///////////////////////*/

	///分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R_out; rigion[A_X][1]=1.2*R_out;
	rigion[A_Y][0]=-1.2*R_out; rigion[A_Y][1]=1.2*R_out;
	rigion[A_Z][0]=-1.2*R_out;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,fluid_number,beforeN,beforeN,newN,r,rigion);
	
	*number=*number+newN;

}

//ドーナツ作成
void set_doughnut2D(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_big,double R_smal,int edge_startID,int edge_lastID)
{
	//すでに作成された円周を用いて、それより大きな半径のドーナツを作成する。
	//edge_startIDから(edge_lastID-1)までの粒子が、すでに作成されている円の外周を構成する粒子に該当する
	//vector型配列は参照渡ししている。vector<double> *Xではなくvector<double> &Xであることに注意。これで各配列は通常通りに仕様可能。アロー演算子もいらない
	//参照渡しでなく通常のやりかたでもいいけど、その場合、例えばa=X[5]と書いても利用できない。a=(*X)[5]などとしなければならない

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;		//よく使う係数

	int half_WX=(int)(R_big/le)+1;  //円を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R_big/(le*A))+1;  //円を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R_big-le*0.4;				//少し小さめの半径を設定
	double R1=R_smal+le*0.4;				//少し大きめの半径を設定


	//もうひとつの円周を作成 粒子数は内部で増加することに注意
	set_circle_edge(X,Y,Z,number,le,R_big);

	edge_lastID=*number;//MD2D()における境界粒子はedge_startID<=i<edge_lastID
	beforeN=*number;

	//内部初期位置 ただし最初の1ピースのみ
	for(int i=-half_WX;i<=half_WX;i++)
	{
		for(int j=-half_WY;j<=half_WY;j++)
		{
			double jj=j*le*A;
			double ii=i*le;
			if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
			if(ii*ii+jj*jj<R2*R2 && ii*ii+jj*jj>R1*R1)
			{
				X.push_back(ii);
				Y.push_back(jj);
				Z.push_back(0);
				newN++;
			}
		}
	}////////////////////////

	//分子動力学により位置を最適化
	MD_2D(X,Y,Z,le,edge_startID,edge_lastID,beforeN,newN);

	*number=*number+newN;
}

//FSWプローブ内部粒子セット関数
void set_hat_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double R_smal,double R_big,double H_hat,double H_flange,int bound_startID,int bound_endID)
{
	//帽子型の物質内部を作成する。例えばFSWのツール形状。
	
	//帽子の頭に該当する半径をR_smal,つばに該当する半径をR_big,つばの幅をH_flange,頭の幅をH_hat
	//ここで作成する帽子型の姿勢は、FSWのツールと同じで、頭を下にしてつばが上。
	//プローブ底辺のZ=0とする

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	double height=H_hat+H_flange;

	int WX=(int)(R_big/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int WY=(int)(R_big/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int WZ=(int)(height/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double gap=0.4*le;					//隙間
	double R_big2=R_big-gap;				//少し小さめの半径を設定
	double R_smal2=R_smal-gap;				//少し小さめの半径を設定
	double height2=height-0.3*le;				//少し小さめの高さを設定
	

	//内部初期位置
	for(int i=-WX;i<=WX;i++)
	{
		for(int j=-WY;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk<=H_hat+gap)
				{
					if(ii*ii+jj*jj<R_smal2*R_smal2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
				else if(kk<height-gap)
				{
					if(ii*ii+jj*jj<R_big2*R_big2)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}///////////////////////*/

	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-1.2*R_big; rigion[A_X][1]=1.2*R_big;
	rigion[A_Y][0]=-1.2*R_big; rigion[A_Y][1]=1.2*R_big;
	rigion[A_Z][0]=-5*le;  rigion[A_Z][1]=height*1.5;

	MD_3D(X,Y,Z,le,bound_startID,bound_endID,beforeN,newN,r,rigion);

	*number=*number+newN;

}

//箱作成関数
void set_box(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<int> &surface,int *number,double le,double Width,double Height,double Depth)
{
	//横(Width)×高さ(Height)×奥行き(Depth)の箱を作成する。デカルト座標の原点に対し、X正方向に横幅、Y正方向に奥行き、Z正方向に高さ

	int BOX_startID=*number;	//この関数開始時の粒子数
	int beforeN=*number;
	int newN=0;
	double A=sqrt(3.0)/2;		//よく使う係数

	vector<double> X2;					//使い回し用
	vector<double> Y2;
	vector<double> Z2;
	int number2;						//使いまわしよう粒子数

	int Ndepth;					//奥行き方向の分割数
	int Nwidth;					//奥行き方向の分割数
	int Nheight;				//奥行き方向の分割数
	double dL_depth;			//行き方向の分割長さ
	double dL_width;			//行き方向の分割長さ
	double dL_height;			//行き方向の分割長さ

	double gap=0.4*le;
	
	calc_N_and_L(Depth,le,&Ndepth,&dL_depth);//各方向の分割数とその長さが求まる
	calc_N_and_L(Width,le,&Nwidth,&dL_width);
	calc_N_and_L(Height,le,&Nheight,&dL_height);

	//XY平面(上下面)作成//////長方形作成関数を使うのではなく、ここできちんと作成する。理由は、ひとつの長方形の辺を別の長方形が使用するから

	set_rectangular(X,Y,Z,number,le,Width,Depth);		//座標の原点は長方形の中心。よってあとで移動すること。
	for(int i=beforeN;i<*number;i++)					//重心移動
	{
		X[i]+=Width*0.5;
		Y[i]+=Depth*0.5;
	}

	for(int i=beforeN;i<*number;i++)					//上面にｺﾋﾟｰ
	{
		X.push_back(X[i]);
		Y.push_back(Y[i]);
		Z.push_back(Height);
		newN++;
	}
	*number=*number+newN;
	/////////////////////////////////////

	//XZ平面(正面・背面)作成
	number2=0;
	set_rectangular(X2,Y2,Z2,&number2,le,Width,Height);		//座標の原点は長方形の中心。よってあとで移動すること。
	for(int i=0;i<number2;i++)								//重心移動
	{
		X2[i]+=Width*0.5;
		Y2[i]+=Height*0.5;		//set_rectangular()は2D用なので、Zは値がゼロに注意
	}
		//粒子追加
	newN=0;
	beforeN=*number;
	for(int i=0;i<number2;i++)								//正式な座標に移動し、正式配列にｺﾋﾟｰ
	{
		if(Y2[i]>gap && Y2[i]<Height-gap)					//Y2=0,Y2=Heightの粒子はすでに作成済みなので省く
		{
			X.push_back(X2[i]);
			Y.push_back(0.0);			//Y=0すなわち正面
			Z.push_back(Y2[i]);			//Y2から値をもらうことに注意
			newN++;
		}
	}
	*number=*number+newN;				//正面粒子配置完了

	newN=0;
	for(int i=beforeN;i<*number;i++)
	{
		X.push_back(X[i]);
		Y.push_back(Depth);			//Y=Depthすなわち背面
		Z.push_back(Z[i]);
		newN++;
	}
	*number=*number+newN;				//背面粒子配置完了
	////////////////////////////////////////////////////////////


	//YZ平面(側面)作成
	newN=0;
	size_t vector_size=X2.size();
	for(int k=0;k<vector_size;k++)//X2,Y2,Z2をいったん消去
	{
		X2.pop_back();
		Y2.pop_back();
		Z2.pop_back();
	}

	number2=0;
	set_rectangular(X2,Y2,Z2,&number2,le,Depth,Height);		//座標の原点は長方形の中心。よってあとで移動すること。
	for(int i=0;i<number2;i++)								//重心移動
	{
		X2[i]+=Depth*0.5;
		Y2[i]+=Height*0.5;									//set_rectangular()は2D用なので、Zは値がゼロに注意
	}
		//粒子追加
	newN=0;
	beforeN=*number;
	for(int i=0;i<number2;i++)								//正式な座標に移動し、正式配列にｺﾋﾟｰ
	{
		if(X2[i]>gap && X2[i]<Depth-gap)					//辺粒子はすでに4辺とも作成ずみなので省く
		{
			if(Y2[i]>gap && Y2[i]<Height-gap)				
			{
				X.push_back(0.0);			//Z=0すなわち左面
				Y.push_back(X2[i]);			//X2から値をもらうことに注意
				Z.push_back(Y2[i]);			//Y2から値をもらうことに注意
				newN++;
			}
		}
	}
	*number=*number+newN;				//左面粒子配置完了

	newN=0;
	for(int i=beforeN;i<*number;i++)
	{
		X.push_back(Width);
		Y.push_back(Y[i]);			//Y=Depthすなわち背面
		Z.push_back(Z[i]);
		newN++;
	}
	*number=*number+newN;				//背面粒子配置完了
	////////////////////////////////////////////////////////////
	
	for(int i= BOX_startID;i<*number;i++)  surface.push_back(ON);//ここまでは表面
	beforeN=*number;

	set_box_in(X,Y,Z,number,le,Width,Height,Depth,BOX_startID,*number);
	
	for(int i=beforeN;i<*number;i++) surface.push_back(OFF);//内部
}

//長方形作成関数
void set_rectangular(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height)
{
	int BO_startID=*number;
	set_rectangular_edge(X,Y,Z,number,le,Width,Height);	//座標の原点は長方形の中心。よってあとで移動すること。
	int BO_lastID=*number;
	set_rectangular_in(X,Y,Z,number,le,Width,Height,BO_startID,BO_lastID);
}

//BOX内作成関数
void set_box_in(vector<double> &X,vector<double> &Y,vector<double> &Z,int *number,double le,double Width,double Height,double Depth,int BO_startID,int BO_lastID)
{
	//横(Width)×高さ(Height)×奥行き(Depth)の箱を作成する。デカルト座標の原点に対し、X正方向に横幅、Y正方向に奥行き、Z正方向に高さ

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=*number;		//この関数呼び出し時における粒子数
	
	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数

	int WX=(int)(Width/le)+1;		
	int WY=(int)(Depth/(le*A))+1; 
	int WZ=(int)(Height/(le*B))+1;  
	double gap=0.4*le;					//隙間
	
	//内部固定初期位置
	for(int i=1;i<=WX;i++)
	{
		for(int j=1;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk<Height*0.75 && jj<Depth*0.75 && ii<Width*0.75)
				{
					if(kk<Height*0.25 && jj<Depth*0.25 && ii<Width*0.25)
					{
						X.push_back(ii);//十分内部は決め打ち
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}
	*number=*number+newN;
	newN=0;					//この関数で新しく追加する粒子数
	beforeN=*number;
	///////////////////////*/

	//内部流動初期位置
	for(int i=1;i<=WX;i++)
	{
		for(int j=1;j<=WY;j++)
		{
			for(int k=1;k<=WZ;k++)
			{
				double ii=i*le;
				double jj=j*le*A;
				double kk=k*le*B;
				if(j%2!=0) ii+=0.5*le;//jが奇数ならiiを0.5格子だけずらす
				if(k%2!=0) {ii+=0.5*le; jj+=sqrt(3.0)/6*le;}//kが奇数ならiiとjjをずらす
				if(kk<Height-gap && jj<Depth-gap && ii<Width-gap)
				{
					if(kk<=Height*0.25 || jj<=Depth*0.25 || ii<=Width*0.25 || kk>=Height*0.75 || jj>=Depth*0.75 || ii>=Width*0.75)
//					if(kk<=Height*0.2 || jj<=Depth*0.2 || ii<=Width*0.2 || kk>=Height*0.8 || jj>=Depth*0.8 || ii>=Width*0.8)
					{
						X.push_back(ii);
						Y.push_back(jj);
						Z.push_back(kk);
						newN++;
					}
				}
			}
		}
	}


	//分子動力学により位置を最適化
	double r=1.5;
	double rigion[3][2];	//解析領域
	rigion[A_X][0]=-5*le; rigion[A_X][1]=1.5*Width;
	rigion[A_Y][0]=-5*le; rigion[A_Y][1]=1.5*Depth;
	rigion[A_Z][0]=-5*le; rigion[A_Z][1]=Height*1.5;

	MD_3D(X,Y,Z,le,BO_startID,BO_lastID,beforeN,newN,r,rigion);

	*number=*number+newN;

}

//物質合成関数
void make_fusion3D(vector<double> &X,vector<double> &Y,vector<double> &Z,vector<double> &X2,vector<double> &Y2,vector<double> &Z2,vector<int> &surface2,int *number,double le)
{
	//存在が優先される粒子の座標がX,Y,Zに、存在が優先されない粒子の座標がX2,Y2,Z2に格納されている。
	//また、存在が優先されない粒子が、表面粒子か内部粒子かがsurface2[i]に格納されている。表面粒子は分子動力学で動かさないようにしないといけない

	size_t pri_num=X.size();			//優先物体の構成粒子数
	size_t neg_num=X2.size();			//消される物体の構成粒子数

	//cout<<X2.size()<<" "<<surface2.size()<<endl;

	double r=1.5;
	double region[3][2];	//解析領域
	int *flag=new int[neg_num];			//flagがONなら存在を許される。OFFなら消される
	int newN;
	int beforeN=*number;

	////////////////////解析領域の決定
	region[A_X][0]=100; region[A_X][1]=-100;
	region[A_Y][0]=100; region[A_Y][1]=-100;
	region[A_Z][0]=100; region[A_Z][1]=-100;
	for(int i=0;i<pri_num;i++)
	{
		if(X[i]<region[A_X][0]) region[A_X][0]=X[i];
		else if(X[i]>region[A_X][1]) region[A_X][1]=X[i];
		if(Y[i]<region[A_Y][0]) region[A_Y][0]=Y[i];
		else if(Y[i]>region[A_Y][1]) region[A_Y][1]=Y[i];
		if(Z[i]<region[A_Z][0]) region[A_Z][0]=Z[i];
		else if(Z[i]>region[A_Z][1]) region[A_Z][1]=Z[i];
	}
	for(int i=0;i<neg_num;i++)
	{
		if(X2[i]<region[A_X][0]) region[A_X][0]=X2[i];
		else if(X2[i]>region[A_X][1]) region[A_X][1]=X2[i];
		if(Y2[i]<region[A_Y][0]) region[A_Y][0]=Y2[i];
		else if(Y2[i]>region[A_Y][1]) region[A_Y][1]=Y2[i];
		if(Z2[i]<region[A_Z][0]) region[A_Z][0]=Z2[i];
		else if(Z2[i]>region[A_Z][1]) region[A_Z][1]=Z2[i];
	}

	for(int D=0;D<3;D++)
	{
		region[D][0]-=5*le;//保険の意味で少し広めにとる
		region[D][1]+=5*le;
	}
	//////////////////////////////

	//計算の高速化のために格子を形成 解析幅がr*leで割り切れるとは限らないので、はみ出したところは切り捨て。なので各軸とも正の方向には余裕を持つこと
	double grid_width=le*((int)(r+1));								//格子の幅。rを含む整数*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X方向の格子の個数
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int grid_sizeZ=(int)((region[A_Z][1]-region[A_Z][0])/grid_width);
	int grid_SIZE=grid_sizeX*grid_sizeY*grid_sizeZ;
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[neg_num];									//消される物体の構成粒子を含む格子番号
	vector<int> *MESH_pri=new vector<int>[grid_SIZE];				//各メッシュに格納される優先粒子ID格納
	//vector<int> *MESH_neg=new vector<int>[grid_SIZE];				//各メッシュに格納される非優先粒子ID格納

	for(int i=0;i<pri_num;i++)	//まずは優先粒子を格子に格納
	{
		int xn=(int)((X[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((Z[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
		MESH_pri[number].push_back(i);
	}
	for(int i=0;i<neg_num;i++)	//次に非優先粒子を格子に格納
	{
		int xn=(int)((X2[i]-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((Y2[i]-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((Z2[i]-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
		//MESH_neg[number].push_back(i);
		index[i]=number;
	}

	//flag[i]計算開始
	for(int i=0;i<neg_num;i++)
	{
		int G_id=index[i];				//格納する格子番号
		flag[i]=ON;
		for(int II=G_id-1;II<=G_id+1;II++)
		{       
			for(int JJ=-1*grid_sizeX;JJ<=grid_sizeX;JJ+=grid_sizeX)
			{
				for(int KK=-1*plane_SIZE;KK<=plane_SIZE;KK+=plane_SIZE)
				{
					int M_id=II+JJ+KK;
					for(int L=0;L<MESH_pri[M_id].size();L++)//近隣の優先粒子を探索
					{
						int j=MESH_pri[M_id][L];
						
						double x=X[j]-X2[i];
						double y=Y[j]-Y2[i];
						double z=Z[j]-Z2[i];
						double dis=sqrt(x*x+y*y+z*z);
						if(dis<0.7*le)
						{
							flag[i]=OFF;	
						}
					}
				}
			}
		}
	}///////////////////

	//flag[i]=ONの粒子のみX,Y,Zに追加
	newN=0;
	for(int i=0;i<neg_num;i++)
	{
		if(flag[i]==ON && surface2[i]==ON)//flag=ONかつ表面粒子
		{
			X.push_back(X2[i]);
			Y.push_back(Y2[i]);
			Z.push_back(Z2[i]);
			newN++;
		}
	}
	*number=*number+newN;//分子動力学ではここまでの粒子が固定粒子扱い
	beforeN=*number;

	newN=0;
	for(int i=0;i<neg_num;i++)
	{
		if(flag[i]==ON && surface2[i]==OFF)//flag=ONかつ内部粒子
		{
			X.push_back(X2[i]);
			Y.push_back(Y2[i]);
			Z.push_back(Z2[i]);
			newN++;
		}
	}

	MD_3D(X,Y,Z,le,0,beforeN,beforeN,newN,r,region);

	*number=*number+newN;

	delete [] index;
	delete [] MESH_pri;
	//delete [] MESH_neg;
	delete [] flag;
	
}

