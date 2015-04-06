#include "stdafx.h"	

#define  A_R 0
#define  A_t 1 //θ
		
using namespace std;

void make_cube_region(mpsconfig &CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2]);//直方体解析領域作成関数
void make_cylinder_region(mpsconfig &CON,vector<point3D> &NODE,int *node, int *divN,double* regionR,double* regionZ);//円筒解析領域作成関数

void MPSTOFEM3D_MRE(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number);//MRE
void MPSTOFEM3D_droplet(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number);//液滴
void MPSTOFEM3D_nanoe(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number);

//MPS_TO_FEM3Dmain関数
void MPS_TO_FEM3D(mpsconfig &CON, int *node_num, vector<point3D> &NODE, vector<mpselastic> &PART, int fluid_number, int particle_number)
{
	int model=CON.get_model_number();

	if(model==3) MPSTOFEM3D_droplet(CON, node_num, NODE,PART, fluid_number, particle_number);//液滴

	if(model>=5 || model<=8 || model==1 || model==11) MPSTOFEM3D_MRE(CON, node_num, NODE, PART, fluid_number, particle_number); 

	if(model==14) MPSTOFEM3D_nanoe(CON,node_num,NODE,PART, fluid_number, particle_number);//液滴

	//ofstream fp("node_c.dat");
	//for(int i=1;i<=*node_num;i++) fp<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	//fp.close();
	cout<<"MPSTOFEM3D完了 節点数="<<*node_num<<endl;
}

void MPSTOFEM3D_MRE(mpsconfig &CON, int *node_num, vector<point3D> &NODE, vector<mpselastic> &PART, int fluid_number, int particle_number)
{
	//粒子位置座標を節点要素座標として変換する。
	//エラストマー→空気領域への強制節点・境界条件設定→磁石の順番
	int num=0;//節点数
	double le=CON.get_distancebp();
    double err=1e-10; //1e-10;
	
	point3D NODE0;
	NODE.clear();
	NODE.push_back(NODE0);//NODEは節点番号1からスタートするから、ここでひとつ確保しておく

    for(unsigned i=0;i<PART.size();i++)
    {
		if(PART[i].type!=WALL){
		//if(PART[i].surface==ON || i%4==0)
		if(PART[i].toFEM==ON)//解析中、FEMに渡す粒子番号を統一したい 2012/05/22
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;	//未知数
			NODE[num].material=PART[i].type;		//FLUID;
			NODE[num].particleID=i;			//節点iに対応する粒子番号はi
			NODE[num].remesh=ON;			//リメッシュON
		}
/*		else{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;	//未知数
			NODE[num].material=AIR;		//FLUID;
			NODE[num].particleID=i;			//節点iに対応する粒子番号はi
			NODE[num].remesh=ON;			//リメッシュON
		}//*/
		}
    }////////////*/

/*	////表面粒子判定．弾性体では不要？
	for(int i=fluid_number;i<particle_number;i++)
    {
		if(PART[i].type==WALL && PART[i].surface==OFF)//流体と接するINWALL
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=ELASTIC;
			NODE[num].particleID=-1;		//節点iに対応する粒子は存在しない
			NODE[num].remesh=ON;			//リメッシュON
		}
    }*/

	//円筒空気領域作成
	if(CON.get_region_shape()==1)
	{
		int divN[3];
		if(CON.get_model_number()==23)	//FEMで追加する節点数を減らす15/2/4
		{
			divN[A_R]=20;
			divN[A_t]=20;
			divN[A_Z]=20;

		}
		else
		{
			divN[A_R]=20;
			divN[A_t]=20;
			divN[A_Z]=20;
		}
		double regionR[2]={0.0, CON.get_RU()};
		double regionZ[2]={CON.get_ZD(), CON.get_ZU()};
		
		int previous_points=num;//現時点での節点数
		
		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//解析領域は円筒

		//次に、いま作成した解析領域の少し内側に節点を追加する。これはFINE3Dで境界条件面によけいな節点が追加されることを防ぐためである。
		double dR=(regionR[1]-regionR[0])/divN[A_R];				//先ほど作成した解析領域のメッシュ長さdR
		double dZ=(regionZ[1]-regionZ[0])/divN[A_Z];				//先ほど作成した解析領域のメッシュ長さdZ
		//以下のように定義した寸法の箱領域を作成すれば、解析領域と新しい箱との隙間には良質な正四面体に近いメッシュが作成される

		regionR[0]=0.0;
		regionR[1]=CON.get_RU()-dR;
		regionZ[0]=CON.get_ZD()+dZ; 
		regionZ[1]=CON.get_ZU()-dZ;

		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//内側の円筒作成

		//固定境界の設定
		for(int i=previous_points+1;i<=num;i++)
		{
			double X=NODE[i].r[A_X];
			double Y=NODE[i].r[A_Y];
			double Z=NODE[i].r[A_Z];
			double R=sqrt(X*X+Y*Y);

			if(Z>CON.get_ZU()-err) NODE[i].boundary_condition=1;
			else if(Z<CON.get_ZD()+err) NODE[i].boundary_condition=1;
			else if(R>CON.get_RU()-err) NODE[i].boundary_condition=1;
			else NODE[i].boundary_condition=0;
			NODE[i].remesh=OFF;
		}

		//remesh領域作成・・・これifで分けるべきでは？
		regionR[0]=0.0;
		regionR[1]=CON.get_RU()*0.7;
		regionZ[0]=CON.get_ZD()*0.7;
		regionZ[1]=CON.get_ZU()*0.7;
		previous_points=num;//現時点での節点数

		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//解析領域は円筒

		for(int i=previous_points+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//境界条件
			NODE[i].remesh=ON;
		}

	}else if(CON.get_region_shape()==0){
		cout<<"直方体解析領域は未解決"<<endl;
		cout<<"計算終了"<<endl;
		exit(1);
	}

	//磁石
	point3D NODE01;
	int divN[3];
	if(CON.get_model_number()==23)	//節点数を減らす15/2/4
	{
		divN[A_R]=20;//半径方向分割数
		divN[A_t]=20;//角度方向分割数（正何角形で近似するか）
		divN[A_Z]=20;//高さ方向分割数
	}
	else
	{
		divN[A_R]=20;//半径方向分割数
		divN[A_t]=20;//角度方向分割数（正何角形で近似するか）
		divN[A_Z]=20;//高さ方向分割数
	}
	double Rmin=0.0;
	double Rmax=CON.get_magnet_r();

	double Zmin=CON.get_magnet_Z()-0.5*CON.get_magnet_H();
	double Zmax=CON.get_magnet_Z()+0.5*CON.get_magnet_H();
	double divL[3];
	divL[A_R]=(Rmax-Rmin)/divN[A_R];
	divL[A_t]=(2*PI)/divN[A_t];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];
	double theta2=CON.get_magnet_angle()*PI*2/360;//磁石の回転角度
	if(CON.get_EM_calc_type()==2){	//永久磁石静磁場
	for(int n=0;n<=divN[A_R];n++)
	{
		if(n==0)//中心
		{
			for(int k=0;k<=divN[A_Z];k++)
			{
				num++;
				NODE.push_back(NODE01);
				NODE[num].r[A_X]=0;
				NODE[num].r[A_Y]=0;
				NODE[num].r[A_Z]=Zmin+divL[A_Z]*k;
				NODE[num].material=MAGNET;
				NODE[num].particleID=-1;					//対応する粒子が存在しないから-1を格納
				NODE[num].boundary_condition=0;			//境界条件
				NODE[num].remesh=OFF;

				if(fabs(theta2)>0.0175)//doubleを0と比較してはいけない（一致するわけがない）
				{
					double XX=NODE[num].r[A_X];
					double ZZ=NODE[num].r[A_Z]-CON.get_magnet_Z();
					double newX=XX*cos(theta2)-ZZ*sin(theta2);
					double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
					NODE[num].r[A_X]=newX;
					NODE[num].r[A_Z]=newZ;
				}
			}
		}
		else
		{
			for(int m=0;m<divN[A_t];m++)
			{
				for(int k=0;k<=divN[A_Z];k++)
				{
					num++;
					NODE.push_back(NODE01);
					double r=divL[A_R]*n;
					double theta=divL[A_t]*m;
					NODE[num].r[A_X]=r*cos(theta);
					NODE[num].r[A_Y]=r*sin(theta);
					NODE[num].r[A_Z]=Zmin+divL[A_Z]*k;
					NODE[num].material=MAGNET;
					NODE[num].particleID=-1;					//対応する粒子が存在しないから-1を格納
					NODE[num].boundary_condition=0;			//境界条件
					NODE[num].remesh=OFF;

					if(fabs(theta2)>0.0175)//磁石が傾いている場合(1degより大きい場合)・・・doubleを0と比較してはいけない
					//if(theta2!=0)これはダメ！！
					{
						double XX=NODE[num].r[A_X];
						double ZZ=NODE[num].r[A_Z]-CON.get_magnet_Z();
						double newX=XX*cos(theta2)-ZZ*sin(theta2);
						double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
						NODE[num].r[A_X]=newX;
						NODE[num].r[A_Z]=newZ;
					}
				}
			}
		}
	}
	}
	else if(CON.get_EM_calc_type()==3){	//電磁石　動磁場
	for(int n=0;n<=divN[A_R];n++)
	{
		if(n==0)//中心
		{
			for(int k=0;k<=divN[A_Z];k++)
			{
				num++;
				NODE.push_back(NODE01);
				NODE[num].r[A_X]=0;
				NODE[num].r[A_Y]=0;
				NODE[num].r[A_Z]=Zmin+divL[A_Z]*k;
				NODE[num].material=IRON;
				NODE[num].particleID=-1;					//対応する粒子が存在しないから-1を格納
				NODE[num].boundary_condition=0;			//境界条件
				NODE[num].remesh=OFF;

				if(fabs(theta2)>0.0175)//doubleを0と比較してはいけない（一致するわけがない）
				{
					double XX=NODE[num].r[A_X];
					double ZZ=NODE[num].r[A_Z]-CON.get_magnet_Z();
					double newX=XX*cos(theta2)-ZZ*sin(theta2);
					double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
					NODE[num].r[A_X]=newX;
					NODE[num].r[A_Z]=newZ;
				}
			}
		}
		else
		{
			for(int m=0;m<divN[A_t];m++)
			{
				for(int k=0;k<=divN[A_Z];k++)
				{
					num++;
					NODE.push_back(NODE01);
					double r=divL[A_R]*n;
					double theta=divL[A_t]*m;
					NODE[num].r[A_X]=r*cos(theta);
					NODE[num].r[A_Y]=r*sin(theta);
					NODE[num].r[A_Z]=Zmin+divL[A_Z]*k;
					if(n<6){
					NODE[num].material=IRON;	
					NODE[num].boundary_condition=0;			//境界条件
					}
					else {
					NODE[num].material=COIL;
					 if((n==divN[A_R]) && (k%3==0))NODE[num].boundary_condition=23;//外周側面
					else if((n==divN[A_R]) && (k%3==1))NODE[num].boundary_condition=22;
					else if((n==divN[A_R]) && (k%3==2))NODE[num].boundary_condition=21;
		/*			else if((n==6) && (k%3==0)) NODE[num].boundary_condition=21;	//鉄心との境界面
					else if((n==6) && (k%3==1))NODE[num].boundary_condition=22;
					else if((n==6) && (k%3==2))NODE[num].boundary_condition=23;  */
					else NODE[num].boundary_condition=0;			//境界条件
					}
					NODE[num].particleID=-1;					//対応する粒子が存在しないから-1を格納
					NODE[num].remesh=OFF;

					if(fabs(theta2)>0.0175)//磁石が傾いている場合(1degより大きい場合)・・・doubleを0と比較してはいけない
					//if(theta2!=0)これはダメ！！
					{
						double XX=NODE[num].r[A_X];
						double ZZ=NODE[num].r[A_Z]-CON.get_magnet_Z();
						double newX=XX*cos(theta2)-ZZ*sin(theta2);
						double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
						NODE[num].r[A_X]=newX;
						NODE[num].r[A_Z]=newZ;
					}
				}
			}
		}
	}
	}
	//磁石の外側に空気層設置
	//全体の境界と被らないように気をつけること
	double regionR[2];
	double regionZ[2];
	int number_of_layers=1+CON.get_magnetic_layers();
	divN[2]=22;
	regionR[0]=0.0;//変更なし
	
	for(int layers=1;layers<=number_of_layers;layers++)
	{
		regionR[1]=Rmax+(divL[A_R])*layers;
		regionZ[0]=Zmin-(divL[A_Z])*layers/2;
		regionZ[1]=Zmax+(divL[A_Z])*layers/2;
		int before=num;//現時点での節点数

		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//解析領域は円筒
		
		for(int i=before+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//境界条件
			NODE[i].remesh=OFF;
		}

		if(fabs(theta2)>0.0175)//磁石が傾いている場合(1degより大きい場合)・・・doubleを0と比較してはいけない
		{
			for(int i=before +1;i<=num;i++)
			{
				double XX=NODE[i].r[A_X];
				double ZZ=NODE[i].r[A_Z]-CON.get_magnet_Z();
				double newX=XX*cos(theta2)-ZZ*sin(theta2);
				double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
				NODE[i].r[A_X]=newX;
				NODE[i].r[A_Z]=newZ;
			}
		}
	}

	/*//磁石の外側に空気層設置
	double regionR[2]={0,Rmax+divL[A_R]};//{0,Rmax+divL[A_R]/2};
	double regionZ[2]={Zmin-divL[A_Z],Zmax+divL[A_Z]};//{Zmin-divL[A_Z]/2,Zmax+divL[A_Z]/2};
	int count=num;//現時点での節点数
		
	make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//解析領域は円筒

	for(int i=count+1;i<=num;i++)
	{
		NODE[i].boundary_condition=0;			//境界条件
		NODE[i].remesh=OFF;
	}

	if(theta2!=0)
	{
		for(int i=count+1;i<=num;i++)
		{
			double XX=NODE[i].r[A_X];
			double ZZ=NODE[i].r[A_Z]-CON.get_magnet_Z();
			double newX=XX*cos(theta2)-ZZ*sin(theta2);
			double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
			NODE[i].r[A_X]=newX;
			NODE[i].r[A_Z]=newZ;
		}
	}*/

	*node_num=num;
}

//液滴
void MPSTOFEM3D_droplet(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number)
{
	int num=0;						//節点数
	double le=CON.get_distancebp();
    double err=1e-10;
	
	point3D NODE0;
	NODE.clear();
	NODE.push_back(NODE0);			//NODEは節点番号1からスタートするから、ここでひとつ確保しておく

	////流体粒子出力
    for(int i=0;i<fluid_number;i++)
    {
		//if(PART[i].surface==ON || i%4==0)
		//if(PART[i].toFEM==ON)//解析中、FEMに渡す粒子番号を統一したい
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=FLUID;
			NODE[num].particleID=i;				//節点iに対応する粒子番号はi
			NODE[num].remesh=ON;			//リメッシュON
		}
    }////////////*/

	/*///////表面粒子の法線ﾍﾞｸﾄﾙ
    double *direct[DIMENTION];
    for(int D=0;D<DIMENTION;D++) direct[D]=new double [fluid_number];
    
  //  double inpoint[3];	//FEMに送る内部節点座標(X,Y,Z)
	int count=0;			
	double Le=CON.get_distancebp()*0.5;
	for(int i=0;i<fluid_number;i++)
	{
		if(PART[i].surface==ON)
		{
			count++;
			//if(count%2==0)//間引き
			{
				direct_f(CON,PART,i,direct);
				for(int D=0;D<DIMENTION;D++) direct[D][i]*=-1;//外向き法線ﾍﾞｸﾄﾙが欲しいから符合をかえる
				
				for(int n=1;n<=4;n+=1)
				{   
					num++;
					NODE.push_back(NODE0);
					for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D]+direct[D][i]*Le*n;
					NODE[num].boundary_condition=0;
					NODE[num].material=AIR;
					NODE[num].particleID=-1;				
					NODE[num].remesh=ON;			//リメッシュON
				}
			}
		}
	}
	for(int D=0;D<DIMENTION;D++) delete [] direct[D];
	//////*/

	if(CON.get_region_shape()==0)
	{
		int divN[3];
		divN[A_X]=10;
		divN[A_Y]=10;
		divN[A_Z]=10;
		double regionX[2]={CON.get_XL(),CON.get_XR()};
		double regionY[2]={CON.get_YD(),CON.get_YU()};
		double regionZ[2]={CON.get_ZD(),CON.get_ZU()};
		int count=num;//現時点での節点数
		
		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//解析領域は直方体

		//次に、いま作成した解析領域の少し内側に節点を追加する。これはFINE3Dで境界条件面によけいな節点が追加されることを防ぐためである。
		double dX=(regionX[1]-regionX[0])/divN[A_X];				//先ほど作成した解析領域のメッシュ長さdX
		double dY=(regionY[1]-regionY[0])/divN[A_Y];				//先ほど作成した解析領域のメッシュ長さdY
		double dZ=(regionZ[1]-regionZ[0])/divN[A_Z];				//先ほど作成した解析領域のメッシュ長さdZ
		//以下のように定義した寸法の箱領域を作成すれば、解析領域と新しい箱との隙間には良質な正四面体に近いメッシュが作成される
		regionX[0]=CON.get_XL()+dX*sqrt(3.0)*0.5; regionX[1]=CON.get_XR()-dX*sqrt(3.0)*0.5;
		regionY[0]=CON.get_YD()+dY*sqrt(3.0)*0.5; regionY[1]=CON.get_YU()-dY*sqrt(3.0)*0.5;
		regionZ[0]=CON.get_ZD()+dZ*sqrt(3.0)*0.5; regionZ[1]=CON.get_ZU()-dZ*sqrt(3.0)*0.5;

		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//内側の箱作成

		for(int i=count+1;i<=num;i++)				//境界条件
		{
			double Z=NODE[i].r[A_Z];
			if(Z>CON.get_ZU()-err) NODE[i].boundary_condition=2;
			else if(Z<CON.get_ZD()+err) NODE[i].boundary_condition=1;
			else NODE[i].boundary_condition=0;
			NODE[i].remesh=OFF;
		}

		//remesh領域作成
		regionX[0]=CON.get_XL()*0.5; regionX[1]=CON.get_XR()*0.5;
		regionY[0]=CON.get_YD()*0.5; regionY[1]=CON.get_YU()*0.5;
		regionZ[0]=CON.get_ZD()*0.7; regionZ[1]=CON.get_ZU()*0.7;
		count=num;//現時点での節点数
		make_cube_region(CON,NODE,&num,divN,regionX,regionY,regionZ);//解析領域は直方体

		for(int i=count+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//境界条件
			NODE[i].remesh=ON;
		}
		
	}
	else if(CON.get_region_shape()==1)
	{
		int divN[3];
		divN[A_R]=10;
		divN[A_t]=30;
		divN[A_Z]=10;
		double regionR[2]={0,CON.get_RU()};
		double regionZ[2]={CON.get_ZD(),CON.get_ZU()};
		int count=num;//現時点での節点数
		
		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//解析領域は円筒

		//次に、いま作成した解析領域の少し内側に節点を追加する。これはFINE3Dで境界条件面によけいな節点が追加されることを防ぐためである。
		double dR=(regionR[1]-regionR[0])/divN[A_R];				//先ほど作成した解析領域のメッシュ長さdR
		double dZ=(regionZ[1]-regionZ[0])/divN[A_Z];				//先ほど作成した解析領域のメッシュ長さdZ
		//以下のように定義した寸法の箱領域を作成すれば、解析領域と新しい箱との隙間には良質な正四面体に近いメッシュが作成される
		regionR[0]=0; regionR[1]=CON.get_RU()-dR;
		regionZ[0]=CON.get_ZD()+dZ; regionZ[1]=CON.get_ZU()-dZ;

		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//内側の円筒作成

		for(int i=count+1;i<=num;i++)				//境界条件
		{
			double X=NODE[i].r[A_X];
			double Y=NODE[i].r[A_Y];
			double Z=NODE[i].r[A_Z];
			double R=sqrt(X*X+Y*Y);
			if(Z>CON.get_ZU()-err) NODE[i].boundary_condition=1;
			else if(Z<CON.get_ZD()+err) NODE[i].boundary_condition=1;
			else if(R>CON.get_RU()-err) NODE[i].boundary_condition=1;
			else NODE[i].boundary_condition=0;
			NODE[i].remesh=OFF;
		}

		//remesh領域作成
		regionR[0]=0; regionR[1]=CON.get_RU()*0.7;
		regionZ[0]=CON.get_ZD()*0.7; regionZ[1]=CON.get_ZU()*0.7;
		count=num;//現時点での節点数
		make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//解析領域は円筒

		for(int i=count+1;i<=num;i++)
		{
			NODE[i].boundary_condition=0;			//境界条件
			NODE[i].remesh=ON;
		}
	}///*/

	

	

	*node_num=num;

}

///3D用静電霧化
void MPSTOFEM3D_nanoe(mpsconfig &CON,int *node_num,vector<point3D> &NODE,vector<mpselastic> &PART, int fluid_number, int particle_number)
{
    double le=CON.get_distancebp();	//初期粒子間距離
    double P_R=0.0005;					//電極半径
//    double dX,dY,dZ;					//各軸方向の分割幅
//    int Nx,Ny,Nz;						//各軸方向の分割数
    double ex_r=0.0006;					//流体がくると想定してる領域半径
//    double hight;						//流体がくると想定してる高さ
	int TOUCH=0;						//接触
	int UNTOUCH=1;						//非接触
	double err=1e-10;
	double beta=CON.get_beta();
	int divN[3];

	int num=0;						//節点数
	
	point3D NODE0;
	NODE.clear();
	NODE.push_back(NODE0);			//NODEは節点番号1からスタートするから、ここでひとつ確保しておく

	////流体粒子出力
    for(int i=0;i<fluid_number;i++)
    {
       // if(PART[i].surface==ON || i%4==0)
		//if(PART[i].toFEM==ON)//解析中、FEMに渡す粒子番号を統一したい
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=FLUID;
			NODE[num].particleID=i;				//節点iに対応する粒子番号はi
			NODE[num].remesh=ON;			//リメッシュON
		}
    }////////////*/

    ///電極粒子出力
	double minZ=0;		//出力する壁粒子の最少高さ
    for(int i=fluid_number;i<particle_number;i++)
    {   
       // if(PART[i].type==BDWALL || PART[i].type==INWALL || i%5==0)
		{
			num++;
			NODE.push_back(NODE0);
			for(int D=0;D<3;D++) NODE[num].r[D]=PART[i].r[D];
			NODE[num].boundary_condition=0;
			NODE[num].material=ELECTRODE;
			NODE[num].particleID=i;				//節点iに対応する粒子番号はi
			NODE[num].remesh=ON;			//リメッシュON*/
			if(minZ>PART[i].r[A_Z]) minZ=PART[i].r[A_Z];
        }
    }////////////////*/

	//下電極作成
	divN[A_R]=10;
	divN[A_t]=20;
	divN[A_Z]=200;
	double regionR[2]={0,P_R+le};
	double regionZ[2]={CON.get_ZD(),minZ-3*le};		//作成する円柱の高さ
	int count=num;//現時点での節点数
		
	make_cylinder_region(CON,NODE,&num,divN,regionR,regionZ);//解析領域は円筒

	for(int i=count+1;i<=num;i++)				//境界条件
	{
		NODE[i].material=ELECTRODE;
		NODE[i].boundary_condition=1;
		NODE[i].remesh=OFF;
	}
    
    
   /* double WX=CON.get_XR()-CON.get_XL();///解析領域X方向幅
    double WY=CON.get_YU()-CON.get_YD();///解析領域Y方向幅
    double WZ=CON.get_ZU()-CON.get_ZD();///解析領域Z方向幅
    
    /////解析領域上 (流体はこないと想定している領域)
    
	//空気領域1
    Nx=40;//120;//15;
    Ny=40;//120;//15;
    Nz=75*2;//60;//75;//15;
	double Xmin1=CON.get_XL()*3/20; double Xmax1=CON.get_XR()*3/20;//空気領域1の最小・最大X座標
	double Ymin1=CON.get_YD()*3/20; double Ymax1=CON.get_YU()*3/20;//空気領域1の最小・最大Y座標
	double Zmin1=CON.get_ZD(); double Zmax1=CON.get_ZU()*0.24;//空気領域1の最小・最大Z座標
	dX=(Xmax1-Xmin1)/Nx;
    dY=(Ymax1-Ymin1)/Ny;
    dZ=(Zmax1-Zmin1)/Nz;

	
    for(double X=Xmin1;X<=Xmax1+err;X+=dX)
    {   
        for(double Y=Ymin1;Y<=Ymax1+err;Y+=dY)
		{   
			for(double Z=Zmin1;Z<=Zmax1+err;Z+=dZ)
			{
				if(X*X+Y*Y>ex_r*ex_r || Z>=hight)
				{
					num++;
					int BC=0;///境界条件
					if(Z>=CON.get_ZU()-err) BC=2;
					if(Z<CON.get_ZD()+err) BC=1;
					fout<<X<<" "<<Y<<" "<<Z<<" "<<AIR<<" "<<BC<<endl;//最後の数字は境界条件
					fprintf(fp,"%lf %lf %lf %d %d\n",X,Y,Z,AIR,BC);//最後の数字は境界条件
				}
			} 
		}
    }////////*/

	/*/空気領域2
	Nx=120/2;//15;
    Ny=120/2;//15;
    Nz=75/2;//15;
    dX=WX/Nx;
    dY=WY/Ny;
    dZ=WZ/Nz;
	
    for(double X=CON.get_XL();X<=CON.get_XR()+err;X+=dX)
    {   
        for(double Y=CON.get_YD();Y<=CON.get_YU()+err;Y+=dY)
		{   
			for(double Z=CON.get_ZD();Z<=CON.get_ZU()+err;Z+=dZ)
			{
				int flag=1;//flag=1なら空気領域1　2なら空気領域2
				if(X>Xmax1+le || X<Xmin1-le) flag=2;
				if(Y>Ymax1+le || Y<Ymin1-le) flag=2;
				if(Z>Zmax1+le || Z<Zmin1-le) flag=2;
				if(flag==2)
				{
					num++;
					int BC=0;///境界条件
					if(Z>=CON.get_ZU()-err) BC=2;
					if(Z<CON.get_ZD()+err) BC=1;
					fout<<X<<" "<<Y<<" "<<Z<<" "<<AIR<<" "<<BC<<endl;//最後の数字は境界条件
					fprintf(fp,"%lf %lf %lf %d %d\n",X,Y,Z,AIR,BC);//最後の数字は境界条件
				}
			} 
		}
    }////////*/

	*node_num=num;
   
}


//直方体解析領域作成関数
void make_cube_region(mpsconfig &CON,vector<point3D> &NODE,int *node, int *divN,double regionX[2],double regionY[2],double regionZ[2])
{
	int node_num=*node;
	//divN[3];								//各辺の分割数
	
	double Xmin=regionX[0];				//解析領域
	double Xmax=regionX[1];
	double Ymin=regionY[0];
	double Ymax=regionY[1];
	double Zmin=regionZ[0];
	double Zmax=regionZ[1];

	double divL[3];								//分割幅
	divL[A_X]=(Xmax-Xmin)/divN[A_X];
	divL[A_Y]=(Ymax-Ymin)/divN[A_Y];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];

	point3D NODE01;

	//底面
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin;					//解析領域の底面
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}

	//上面
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=0;m<=divN[A_Y];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin+m*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmax;					//解析領域の上面
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}////

	
	//側面Y
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymin;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];		
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	for(int n=0;n<=divN[A_X];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin+n*divL[A_X];
			NODE[node_num].r[A_Y]=Ymax;
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}

	//側面
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmin;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	for(int n=1;n<divN[A_Y];n++)
	{
		for(int m=1;m<divN[A_Z];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			NODE[node_num].r[A_X]=Xmax;
			NODE[node_num].r[A_Y]=Ymin+n*divL[A_Y];
			NODE[node_num].r[A_Z]=Zmin+m*divL[A_Z];	
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	*node=node_num;
}

//円筒解析領域作成関数
void make_cylinder_region(mpsconfig &CON,vector<point3D> &NODE,int *node, int *divN, double* regionR, double* regionZ)
{
	int node_num=*node;
	//divN[3];					//各辺の分割数
	
	double Rmin=0;				//解析領域
	double Rmax=regionR[1];
	double Zmin=regionZ[0];
	double Zmax=regionZ[1];

	double divL[3];//分割幅
	divL[A_R]=(Rmax-Rmin)/divN[A_R];
	divL[A_t]=(2*PI)/divN[A_t];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];

	point3D NODE01;

	/////////////////////底面
	node_num++;					//中心点
	NODE.push_back(NODE01);
	NODE[node_num].r[A_X]=0;
	NODE[node_num].r[A_Y]=0;
	NODE[node_num].r[A_Z]=Zmin;					//解析領域の底面
	NODE[node_num].material=AIR;
	NODE[node_num].particleID=-1;

	for(int n=1;n<=divN[A_R];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			double r=divL[A_R]*n;
			double theta=divL[A_t]*m;
			NODE[node_num].r[A_X]=r*cos(theta);
			NODE[node_num].r[A_Y]=r*sin(theta);
			NODE[node_num].r[A_Z]=Zmin;					//解析領域の底面
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	//////////////////////////////////上面
	node_num++;					//中心点
	NODE.push_back(NODE01);
	NODE[node_num].r[A_X]=0;
	NODE[node_num].r[A_Y]=0;
	NODE[node_num].r[A_Z]=Zmax;					//解析領域の底面
	NODE[node_num].material=AIR;
	NODE[node_num].particleID=-1;	
	for(int n=1;n<=divN[A_R];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			double r=divL[A_R]*n;
			double theta=divL[A_t]*m;
			NODE[node_num].r[A_X]=r*cos(theta);
			NODE[node_num].r[A_Y]=r*sin(theta);
			NODE[node_num].r[A_Z]=Zmax;					//解析領域の底面
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}

	//側面
	double RR=divL[A_R]*divN[A_R];
	for(int n=1;n<divN[A_Z];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			node_num++;
			NODE.push_back(NODE01);
			double theta=divL[A_t]*m;
			NODE[node_num].r[A_X]=RR*cos(theta);
			NODE[node_num].r[A_Y]=RR*sin(theta);;
			NODE[node_num].r[A_Z]=Zmin+n*divL[A_Z];		
			NODE[node_num].material=AIR;
			NODE[node_num].particleID=-1;					//対応する粒子が存在しないから-1を格納
		}
	}
	*node=node_num;
}