#include "stdafx.h"
#include "Model.h"

#define FULL 1
#define HALF 2
#define HALFD 3
#define HALFD_shell 4

Model :: Model()
	:fq("initial_input.dat")
{
	
}

Model :: ~Model()
{
	fq.close();
}

////////////////////////モデルセット//////////////////////
int Model :: Model_set(){
	////初期化されたmpsconfigクラスのデータ
	mpsconfig CON;		
	////////////////////////////////////
	int number;
	if(CON.get_model_number()==0);
	else if(CON.get_model_number()==1 || CON.get_model_number()==11) number=Model :: Set_sphere_model() ;
	else if(CON.get_model_number()==2) number=Model :: Set_cylinder_model() ;
	else if(CON.get_model_number()==3) number=Model :: Set_cube_model() ;
	else if(CON.get_model_number()==4) number=Model :: Set_tensiontest_model();
	else if(CON.get_model_number()==5) number=Model :: Set_actuator_model();
	else if(CON.get_model_number()==6) number=Model :: Set_benchmark_model();
	else {cout<<"non-existent Model number"<<endl;
	exit(0);}
	return (number);
}

int Model :: Set_sphere_model()
{
	////初期化されたmpsconfigクラスのデータ
	mpsconfig CON;		
	////////////////////////////////////
	double le=CON.get_distancebp();
	int initial_ID=0;
	int number=0;
	double R=CON.get_fluidwidth()*le*0.5;			//作成する円の半径
	double Zg=CON.get_height();					//球の中心高さ
	

			//球作成
			int flag=FULL;
			number=Model :: Set_sphere_function(initial_ID,  le,  R);
		
		//初期情報書き込み
		//for(int i=0;i<number;i++) writedata2(fq,i,X[i],Y[i],Z[i],FRFLUID,1,0,Y[i]*40,0,0,0,0);
		double vx=0;
		double vy=0;
		double vz=0;
		
			for(int i=0;i<number;i++)
			{
				double r=sqrt(PART[i].Get_X()*PART[i].Get_X()+PART[i].Get_Y()*PART[i].Get_Y());
				double a=0.2;
		//		a=500*3;//クーロン分裂のとき
				a=200;	//遊び
				
				
				vx=-0.5*a*PART[i].Get_X();
				vy=-0.5*a*PART[i].Get_Y();
				vz=a*PART[i].Get_Z();
				 
				writedata(fq,  i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()+Zg,ELASTIC,1, OFF,0, vx,vy, vz,0, 0, 0);
		}
			return (number);
}

int Model :: Set_cylinder_model()
{
	return (0);
}

int Model :: Set_cube_model()
{
	mpsconfig CON;					//初期化されたmpsconfigクラスのデータ
	double le=CON.get_distancebp();
	double z=10.0;
	double y=25.0;
	double x=25.0;
	double zw=3.0;
	double yw=31.0;
	double xw=31.0;
	int num=0;

	for(int i=0;i<x;i++){
		for(int j=0;j<y;j++){
			for(int k=0;k<z;k++){
				PART.push_back(PART0);
				PART[num].Set(i*le, j*le, k*le, le);
				num++;
			}
		}
	}
	int elasnum=num;  //1～の数 壁粒子IDの1番目

/*	for(int i=0;i<xw;i++){
		for(int j=0;j<yw;j++){
			for(int k=0;k<zw;k++){
				PART.push_back(PART0);
				PART[num].Set(i*le, j*le, k*le, le);
				num++;
			}
		}
	}*/
	for(int i=0;i<elasnum;i++) writedata(fq,i,PART[i].Get_X()-((x-1)*le)/2.0,PART[i].Get_Y()-((y-1)*le)/2.0,PART[i].Get_Z(),MAGELAST,1,0,0,0,0,0,0,0,1);//粒子はMAGELAST
//	for(int i=elasnum;i<num;i++) writedata(fq,i,PART[i].Get_X()-(xw-1)*le/2.0,PART[i].Get_Y()-(yw-1)*le/2.0,PART[i].Get_Z()-(2.0+zw)*le,WALL,1,0,0,0,0,0,0,0,0);//粒子はWALL
	return num;
}
int Model :: Set_tensiontest_model()
{
	mpsconfig CON;					//初期化されたmpsconfigクラスのデータ
	int length=25;	//[mm]引っ張り端子とは別35
	int width=11;//15
	int depth=11;
	int tan=0;//端子高さ粒子5
	int number=0;//初期粒子数
	double accel=0.005;
	double le=CON.get_distancebp();
	le*=0.5;//*(11.0/12.0)

	int low_ID=(((width*depth)-1)/2)+(width*depth*(8-1));
	int left_ID=(width*depth)*(length-1)/2+(width*2)+(width-1)/2;
	int right_ID=(width*depth)*(length-1)/2+(width*(width-3))+(width-1)/2;
	int top_ID=(((width*depth)-1)/2)+(width*depth*(18-1));
	Model :: Set_point(low_ID,top_ID,left_ID,right_ID);
	cout<<low_ID<<","<<top_ID<<","<<left_ID<<","<<right_ID<<endl;
	////////////下端子////////////
	for(int i=0;i<tan;i++){
		for(int j=0;j<width;j++){
			for(int k=0;k<depth;k++){
				PART.push_back(PART0);
				PART[number].Set(k*le, j*le, i*le, le);
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
				PART.push_back(PART0);
				PART[number].Set(k*le, j*le, i*le, le);
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
				PART.push_back(PART0);
				PART[number].Set(k*le, j*le, i*le, le);
				number++;
			}
		}
	}
	////////////////////////////////
		
	for(int i=0;i<number;i++){ 
		if(i<low)writedata(fq,i,PART[i].Get_X()-depth*le/2,PART[i].Get_Y()-width*le/2,PART[i].Get_Z()-(length+10)*le/2,WALL,0, OFF,0, 0, 0, -accel, 0, 0, 0);
		else if(i>=low && i<mid)writedata(fq,i,PART[i].Get_X()-depth*le/2,PART[i].Get_Y()-width*le/2,PART[i].Get_Z()-(length+10)*le/2,MAGELAST,0, OFF,0, 0,0, 0, 0, 0, 0);
		else if(i>=mid)writedata(fq,i,PART[i].Get_X()-depth*le/2,PART[i].Get_Y()-width*le/2,PART[i].Get_Z()-(length+10)*le/2,WALL,0, OFF,0, 0, 0, accel, 0, 0, 0);
	}
	return number;
}

///////////////アクチュエータモデル/////////////////
int Model :: Set_actuator_model()	
{
	//////アクチュエータの寸法///////
	double elast_R=0.012;	//半径
	double elast_H=0.03;
	double MRE_R=0.015;
	double MRE_H=0.04;
	////////////////////////////////

	////初期化されたmpsconfigクラスのデータ
	mpsconfig CON;		
	////////////////////////////////////

	//////粒子の直径////
	double le=CON.get_distancebp();
	////////////////////

	//////粒子数格納//////
	int total_number;	//この関数での全粒子数
	int start_ID=0;		//初期粒子数
	/////////////////////

	//////円筒形作成/////
	total_number=Model :: Set_cylinder(start_ID, le, MRE_R, MRE_H);
	////////////////////

	for(int i=start_ID;i<total_number;i++)
		{
			if(pow(PART[i].Get_X(),2)+pow(PART[i].Get_Y(),2)<=pow(0.01,2) && PART[i].Get_Z()>0.01 && PART[i].Get_Z()<0.02+0.01){
					writedata(fq,i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()-(MRE_H/2),WALL,1,0,0,0,0,0,0,0,0);//粒子はWALL
			}
			else if(pow(PART[i].Get_X(),2)+pow(PART[i].Get_Y(),2)<=pow(elast_R,2) && PART[i].Get_Z()>0.005 && PART[i].Get_Z()<elast_H+0.0045){
				writedata(fq,i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()-(MRE_H/2),ELASTIC,1,0,0,0,0,0,0,0,0);//粒子はELASTIC
			}
			else writedata(fq,i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()-(MRE_H/2),MAGELAST,1,0,0,0,0,0,0,0,1);//粒子はMAGELAST
			
		}

	return total_number;

}

int Model :: Set_benchmark_model()
{
	double bord_r=0.015;
	double bord_h=0.005;
	double elast_r=0.012;
	double elast_h=0.010;
	int total_number;
	int start_ID=0;
	mpsconfig CON;					//初期化されたmpsconfigクラスのデータ
	double le=CON.get_distancebp();

	total_number=Model :: Set_cylinder(start_ID, le, bord_r, bord_h);
	for(int i=start_ID;i<total_number;i++) writedata(fq,i,PART[i].Get_X(),PART[i].Get_Y(),PART[i].Get_Z()-(elast_h/2+bord_h),WALL,1,0,0,0,0,0,0,0,0);//粒子はWALL

	start_ID=total_number;
	total_number+=Model :: Set_cylinder(start_ID, le, elast_r, elast_h);
	for(int j=start_ID;j<total_number;j++) writedata(fq,j,PART[j].Get_X(),PART[j].Get_Y(),PART[j].Get_Z()-(elast_h/2),ELASTIC,1,0,0,0,0,0,0,0,0);//粒子はELASTIC

/*	start_ID=total_number;
	total_number+=Model :: Set_cylinder(start_ID, le, bord_r, bord_h);
	for(int k=start_ID;k<total_number;k++) writedata(fq,k,PART[k].Get_X(),PART[k].Get_Y(),PART[k].Get_Z()+(elast_h/2),WALL,1,0,0,0,0,0,0,0,0);//粒子はWALL
	*/
	return total_number;
}

///////////////粒子の情報を外部ファイルに出力/////////////////
void Model :: writedata(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM)
{
	//ファイル出力
	fp<<id<<"\t"<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<"\t"<<P<<"\t"<<h<<"\t"<<val<<"\t"<<type<<"\t"<<materialID<<"\t"<<surface<<"\t"<<toBEM<<"\t"<<endl;
}

int Model :: Set_sphere_function(int initial_ID, double le, double R)
{
	int number=0;
	int flag=FULL;
	cout<<"Set_circle_edge開始"<<endl;
	number+=Set_circle_edge(initial_ID, le, R);
	cout<<"Set_circle_in_using_6_pieces開始"<<endl;
	number+=Set_circle_in_using_6_pieces(number,le, R, initial_ID);
	cout<<"Set_sphere開始"<<endl;
	number+=Set_sphere(number,le, R, flag);
	return number;
}

int Model :: Set_cylinder(int initial_ID, double le, double R, double height)	//円筒作成関数
{
	int number=0;
	cout<<"Set_circle_edge開始"<<endl;
	number+=Set_circle_edge(initial_ID, le, R);
	cout<<"Set_circle_in開始"<<endl;
	number+=Set_circle_in(number, le, R, initial_ID);
	cout<<"Set_cylinder_face開始"<<endl;
	number+=Set_cylinder_face(number, le, R, height, initial_ID);
	cout<<"Set_cylinder_in開始"<<endl;
	number+=Set_cylinder_in(number, le, R, height, initial_ID);
	cout<<"終了"<<endl;
	return number;
}

int Model :: Set_circle_edge(int initial_number, double le, double R)	//円外周関数
{
	int N=(int)((2*PI*R)/le);//円周の分割数
	
	double L=2*PI*R/N;				//粒子間距離
	double theta=L/R;				//粒子を配置する角度

	for(int n=initial_number;n<N;n++)
	{
		PART.push_back(PART0);
		PART[n].Set(R*cos(theta*n), R*sin(theta*n), 0, le);
	}
	return N;
}

int Model :: Set_circle_in(int number, double le, double R, int initial_number)	//円内部充填関数
{
	//set_circle_in()と違い、60度分だけ計算し、それを6つｺﾋﾟｰして円を構成する。時間短縮と配置の均等化が目的
	//edge_startIDから(edge_lastID-1)までの粒子が、円の外周を構成する粒子に該当する
	//vector型配列は参照渡ししている。vector<double> *Xではなくvector<double> &Xであることに注意。これで各配列は通常通りに仕様可能。アロー演算子もいらない
	//参照渡しでなく通常のやりかたでもいいけど、その場合、例えばa=X[5]と書いても利用できない。a=(*X)[5]などとしなければならない

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=number;		//この関数呼び出し時における粒子数
	

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

	PART.push_back(PART0);
	PART[beforeN].Set(0, 0, 0, le);						//中心粒子を追加
	
	newN++;
	for(int k=0;k<6;k++)				//6つの直線のloop
	{
		double theta=PI/3*k;			//直線の角度
		for(int n=1;n<R_num;n++)		//中心粒子と最外周粒子はもうあるから、ここでのloopはそれをカウントしない
		{
			double r=L*n;				//中心からの距離
			PART.push_back(PART0);
			PART[beforeN+newN].Set(r*cos(theta), r*sin(theta), 0, le);
			newN++;
		}
	}	
	beforeN+=newN;	//ここまでの粒子数
	newN=0;			//これから増える粒子数
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
					PART.push_back(PART0);
					PART[beforeN+newN].Set(ii, jj, 0, le);	
					newN++;
				}
			}
		}
	}////////////////////////
	//分子動力学により位置を最適化
	//MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);
	MD_2D(le,initial_number,beforeN,newN);
	int M_N=newN;
	///内部粒子を周方向に6つｺﾋﾟｰ
	for(int angle=1;angle<6;angle++)
	{
		double theta=PI/3*angle;//回転する角度
		for(int k=0;k<M_N;k++)
		{
			int i=beforeN+k;
			double x=cos(theta)*PART[i].Get_X()-sin(theta)*PART[i].Get_Y();//回転後の座標
			double y=sin(theta)*PART[i].Get_X()+cos(theta)*PART[i].Get_Y();

			PART.push_back(PART0);
			PART[beforeN+newN].Set(x, y, 0, le);
			newN++;
		}
	}///////////////////*/
	return (beforeN+newN-number);	//固定枠粒子数＋流動粒子数
}

int Model:: Set_circle_in_using_6_pieces(int number,double le,double R,int  initial_ID)
{
	//set_circle_in()と違い、60度分だけ計算し、それを6つｺﾋﾟｰして円を構成する。時間短縮と配置の均等化が目的
	//edge_startIDから(edge_lastID-1)までの粒子が、円の外周を構成する粒子に該当する
	//vector型配列は参照渡ししている。vector<double> *Xではなくvector<double> &Xであることに注意。これで各配列は通常通りに仕様可能。アロー演算子もいらない
	//参照渡しでなく通常のやりかたでもいいけど、その場合、例えばa=X[5]と書いても利用できない。a=(*X)[5]などとしなければならない

	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=number;		//この関数呼び出し時における粒子数

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
	PART.push_back(PART0);
	PART[newN].Set(0,0,0,le);	//中心粒子を追加			
	newN++;

	for(int k=0;k<6;k++)				//6つの直線のloop
	{
		double theta=PI/3*k;			//直線の角度
		for(int n=1;n<R_num;n++)		//中心粒子と最外周粒子はもうあるから、ここでのloopはそれをカウントしない
		{
			double r=L*n;				//中心からの距離
			PART.push_back(PART0);
			PART[newN].Set(r*cos(theta),r*sin(theta),0,le);	//中心粒子を追加			
			newN++;
		}
	}
	number=number+newN;
	newN=0;
	beforeN=number;
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
					PART.push_back(PART0);
					PART[newN].Set(ii,jj,0,le);	//中心粒子を追加			
					newN++;
				}
			}
		}
	}////////////////////////

	//分子動力学により位置を最適化
	//MD_2D(X,Y,Z,le,0,beforeN,beforeN,newN);
	MD_2D(le, initial_ID,beforeN,newN);

	///内部粒子を周方向に6つｺﾋﾟｰ
	for(int angle=1;angle<6;angle++)
	{
		double theta=PI/3*angle;//回転する角度
		for(int k=0;k<newN;k++)
		{
			int i=beforeN+k;
			double x=cos(theta)*PART[i].Get_X()-sin(theta)*PART[i].Get_Y();//回転後の座標
			double y=sin(theta)*PART[i].Get_X()+cos(theta)*PART[i].Get_Y();

			PART.push_back(PART0);
			PART[newN].Set(x,y,0,le);	//中心粒子を追加			
			newN++;
			}
	}///////////////////*/

	return (number+newN*6);//newNはひとつのピース内の粒子数を表しているからここでは6倍
}


//半径Rの球作成関数
int Model :: Set_sphere(int number,double le,double R,int flag)
{
	//まずは半球を作る。そのためには半球表面を作成する必要がある。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=number;		//この関数呼び出し時における粒子数
	int beforeN2=number;		//関数呼び出し時の粒子数。この関数の最後まで記憶しておく

	double A=sqrt(3.0)/2;				//よく使う係数
	double B=sqrt(2.0/3);						////よく使う係数
	int half_WX=(int)(R/le)+1;		 //球を十分含む四角形を想定する。その四角形の幅*0.5
	int half_WY=(int)(R/(le*A))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	int half_WZ=(int)(R/(le*B))+1;  //球を十分含む正四角形を想定する。その四角形の幅*0.5
	double R2=R-0.5*le;				//少し小さめの半径を設定

	///////////半球表面
	int Nt;						//球表面の、θ方向の分割数
	double Lt;					//球表面の、θ方向の分割距離
	Model::Set_calc_N_and_L(PI/2*R,le*A,&Nt,&Lt);//半円の分割数は偶数・奇数どちらでもよい
	double d_theta=Lt/R;		//弧の長さがLtになる角度

	for(int k=0;k<Nt;k++)//loopはk<Ntで終わらせる。Ntに該当するところはすでに設置済み
	{
		double THETA=k*d_theta;	//θ
		double r=R*sin(THETA);	//その高さにおける円の半径
		double round=2*PI*r;//その高さにおける円周

		int Nf=Model :: Set_calc_division_N_circle(round,le);//球表面の、θ方向の分割数
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
			PART.push_back(PART0);
			PART[newN].Set(x,y,z,le);			
			newN++;
		}
	}
	if(Nt%2!=0)//Ntが奇数のときは、頂上に粒子が置かれなければならない。しかし上のloopはそれが不可能。よってここで追加
	{
		PART.push_back(PART0);
		PART[newN].Set(0,0,R,le);			
		newN++;
	}
	//////////////////////////////////////

	number=number+newN;

	if(flag!=HALFD_shell)
	{
	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=number;		//この関数呼び出し時における粒子数
	
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
					PART.push_back(PART0);
					PART[newN].Set(ii,jj,kk,le);			
					newN++;
				}
			}
		}
	}
	number=number+newN;

	newN=0;					//個の関数で新しく追加する粒子数
	beforeN=number;		//この関数呼び出し時における粒子数
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
					PART.push_back(PART0);
					PART[newN].Set(ii,jj,kk,le);			
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

	MD_3D(le,beforeN2,beforeN,newN,r,rigion);

	number=number+newN;
	}

	///上半球を下半球へｺﾋﾟｰし球を完成
	if(flag==FULL)
	{
		newN=0;					//新しく追加する粒子数
		beforeN=number;		//この時における粒子数

		for(int k=0;k<beforeN;k++)
		{
			if(PART[k].Get_Z()>0.4*le)
			{
				newN++;
				PART.push_back(PART0);
				PART[newN].Set(PART[k].Get_X(),PART[k].Get_Y(),-PART[k].Get_Z(),le);			
				newN++;
			}
		}
		number=number+newN;
	}///////////////////////////////
	if(flag==HALFD || flag==HALFD_shell)		//下半球がほしいときに、つくった上半球を上下反転させる
	{
		for(int k=beforeN2;k<number;k++) PART[k].Multiply(1,1,-1);
	}
	return (number);
}

void Model :: Set_point(int low_ID,int top_ID,int lefy_ID,int right_ID)
{
	
	point[0]=low_ID;
	point[1]=top_ID;
	point[2]=lefy_ID;
	point[3]=right_ID;
	
}

void Model :: Get_point(int &a,int &b,int &c,int &d)
{
	a=point[0];
	b=point[1];
	c=point[2];
	d=point[3];
}

void Model :: MD_2D(double le,int BstartID,int beforeN,int newN)
{
	//分子動力学によりnewN個の粒子の位置を最適化　IDがBstartIDからBendIDまでのは境界粒子なので動かさない

	double region[2][2];	//解析領域

	/////////////////////解析領域の決定
	region[A_X][0]=100; region[A_X][1]=-100;
	region[A_Y][0]=100; region[A_Y][1]=-100;
	for(int i=BstartID;i<beforeN;i++)
	{
		if(PART[i].Get_X()<region[A_X][0]) region[A_X][0]=PART[i].Get_X();
		else if(PART[i].Get_X()>region[A_X][1]) region[A_X][1]=PART[i].Get_X();

		if(PART[i].Get_Y()<region[A_Y][0]) region[A_Y][0]=PART[i].Get_Y();
		else if(PART[i].Get_Y()>region[A_Y][1]) region[A_Y][1]=PART[i].Get_Y();
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
//	cout<<"ここっぽい"<<endl;
//	vector<int> *MESH=new vector<int>[plane_SIZE];				//各メッシュに格納される粒子ID格納
	vector<vector<int> > MESH;
	MESH.resize(plane_SIZE);
	for(int i=BstartID;i<beforeN;i++)	//まずは境界粒子を格子に格納
	{
		int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);	//X方向に何個目の格子か 
		int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);	//Y方向に何個目の格子か
		int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
	{
		int i=beforeN+k;
		int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
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
			
			for(int i=BstartID;i<beforeN;i++)	//まずは境界粒子を格子に格納
			{
				int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);	//X方向に何個目の格子か 
				int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);	//Y方向に何個目の格子か
				int number=yn*grid_sizeX+xn;					//粒子iを含む格子の番号
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
			{
				int i=beforeN+k;
				int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
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
						double x=PART[j].Get_X()-PART[i].Get_X();
						double y=PART[j].Get_Y()-PART[i].Get_Y();
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
			PART[i].Add(dt*(U[k]+u)*0.5, dt*(V[k]+v)*0.5, 0);	
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
						double x=PART[j].Get_X()-PART[i].Get_X();
						double y=PART[j].Get_Y()-PART[i].Get_Y();
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
				double dX=PART[J].Get_X()-PART[i].Get_X();
				double dY=PART[J].Get_Y()-PART[i].Get_Y();
				PART[i].Add(-(dX/mindis*L), -(dY/mindis*L), 0);		
			}
			else if(J!=i && J>=beforeN)//leより近接している相手が内部粒子なら
			{
				double L=0.5*(le-mindis);//開くべき距離
				double dX=PART[J].Get_X()-PART[i].Get_X();
				double dY=PART[J].Get_Y()-PART[i].Get_Y();
				PART[i].Add(-(dX/mindis*L), -(dY/mindis*L), 0);	
				PART[J].Add(dX/mindis*L, dY/mindis*L, 0);	
			}
		}//////////*/
	}/////MD終了

	delete [] index;
//	delete [] MESH;
}

int Model :: Set_cylinder_face(int number,double le,double R,double height,int circle_start_id)
{
//半径R,高さheightの円柱の面を作成する。ただしこの関数呼び出し時において、すでに下面の円(Z=0)は作成済みとする
	//top_flag=ONなら円柱上面を作成する。OFFならしないが、側面だけは作成する。
	int beforeN=number;
	int circle_end_id=number;
	int newN=0;

	int Nv;				//水直の分割数
	double dL_V;		//水直の分割長さ
	double A=sqrt(3.0)/2;		//よく使う係数
	/////////直線分割///////////////////////////
	double temp_N=height/(le*A);			//仮の分割数。leで割り切れたら一番いいけど、そうもいかないときがある
	int Ns=(int) temp_N;				//真の分割数
	double difference=temp_N-Ns;		//仮と真の差
	if(difference>0.5) Ns++;
	dL_V=height/Ns;			//粒子の距離
	Nv=Ns;
	//////////////////////////////////////////////

	int Nr=(int)((2*PI*R)/le);//円周の分割数
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
					PART.push_back(PART0);
					PART[beforeN+newN].Set(R*cos(theta), R*sin(theta), jj, le);
					newN++;
				}
			}
		}
	}

	////////////////////////


	//上面作成（下面のコピー）ただしNvが奇数なら上面は下面と半格子ずれなければならない
	beforeN+=newN;
//	cout<<PART.size()<<" "<<beforeN<<endl;
	newN=0;
	if(Nv%2==0)	//偶数ならそのままｺﾋﾟｰ
	{
		
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				PART.push_back(PART0);
				PART[beforeN+newN].Set(PART[i].Get_X(), PART[i].Get_Y(), height, le);
				newN++;
			}
	}
	else
	{
		double d_theta=0.5*Lr/R;//この微小角度だけ回転させる。
			for(int i=circle_start_id;i<circle_end_id;i++)
			{
				double X2=PART[i].Get_X()*cos(d_theta)-PART[i].Get_Y()*sin(d_theta);//回転後の座標
				double Y2=PART[i].Get_X()*sin(d_theta)+PART[i].Get_Y()*cos(d_theta);
				PART.push_back(PART0);
				PART[beforeN+newN].Set(X2, Y2, height, le);
				newN++;
			}		
	}
	return beforeN+newN-number;
}

int Model :: Set_cylinder_in(int number,double le, double R, double height, int start_id)
{
	//半径R,高さheightの円柱内部を作成する。この関数呼び出し時に0<=i<numberの粒子で円柱表面が形成されているとする。
	int newN=0;					//個の関数で新しく追加する粒子数
	int beforeN=number;		//この関数呼び出し時における粒子数

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
						PART.push_back(PART0);
					//	cout<<PART.size()<<"   "<<beforeN+newN<<endl;
						PART[beforeN+newN].Set(ii, jj, kk, le);
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

	MD_3D(le,start_id,beforeN,newN,r,rigion);  //0をstert_idに変更
	return newN;
}

void Model :: MD_3D(double le,int BstartID,int beforeN,int newN,double r,double region[3][2])
{
//分子動力学によりnewN個の粒子の位置を最適化　IDがBstartIDからBendIDまでのは境界粒子なので動かさない
	double k0=1;
	double dt=0.001;
	int BendID=beforeN;
	
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
	double grid_width=le*((int)(r+1));								//格子の幅。rを含む整数*le
	int grid_sizeX=(int)((region[A_X][1]-region[A_X][0])/grid_width);	//X方向の格子の個数
	int grid_sizeY=(int)((region[A_Y][1]-region[A_Y][0])/grid_width);
	int grid_sizeZ=(int)((region[A_Z][1]-region[A_Z][0])/grid_width);
	int grid_SIZE=grid_sizeX*grid_sizeY*grid_sizeZ;
	int plane_SIZE=grid_sizeX*grid_sizeY;
	int *index=new int[newN];									//各内部粒子を含む格子番号
	vector<int> *MESH=new vector<int>[grid_SIZE];				//各メッシュに格納される粒子ID格納

	for(int i=BstartID;i<BendID;i++)	//まずは境界粒子を格子に格納
	{
		int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((PART[i].Get_Z()-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
		int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
		MESH[number].push_back(i);
	}
	for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
	{
		int i=beforeN+k;
		int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X方向に何個目の格子か 
		int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
		int zn=(int)((PART[i].Get_Z()-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
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
				int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
				int zn=(int)((PART[i].Get_Z()-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
				int number=zn*grid_sizeX*grid_sizeY+yn*grid_sizeX+xn;//粒子iを含む格子の番号
				MESH[number].push_back(i);
			}
			for(int k=0;k<newN;k++)	//つぎに内部粒子を格納
			{
				int i=beforeN+k;
				int xn=(int)((PART[i].Get_X()-region[A_X][0])/grid_width);//X方向に何個目の格子か 
				int yn=(int)((PART[i].Get_Y()-region[A_Y][0])/grid_width);//Y方向に何個目の格子か
				int zn=(int)((PART[i].Get_Z()-region[A_Z][0])/grid_width);//Z方向に何個目の格子か
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
								double x=PART[j].Get_X()-PART[i].Get_X();
								double y=PART[j].Get_Y()-PART[i].Get_Y();
								double z=PART[j].Get_Z()-PART[i].Get_Z();
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
								double x=PART[j].Get_X()-PART[i].Get_X();
								double y=PART[j].Get_Y()-PART[i].Get_Y();
								double z=PART[j].Get_Z()-PART[i].Get_Z();
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
			PART[i].Add(dt*(U[k]+u)*0.5, dt*(V[k]+v)*0.5, dt*(W[k]+w)*0.5);
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
							double x=PART[j].Get_X()-PART[i].Get_X();
							double y=PART[j].Get_Y()-PART[i].Get_Y();
							double z=PART[j].Get_Z()-PART[i].Get_Z();
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
				double dX=PART[J].Get_X()-PART[i].Get_X();
				double dY=PART[J].Get_Y()-PART[i].Get_Y();
				double dZ=PART[J].Get_Z()-PART[i].Get_Z();
				PART[i].Add(-dX/mindis*L, -dY/mindis*L, -dZ/mindis*L);
			}
			else if(J!=i && J>=beforeN)//leより近接している相手が内部粒子なら
			{
				double L=0.5*(le-mindis);//開くべき距離
				double dX=PART[J].Get_X()-PART[i].Get_X();
				double dY=PART[J].Get_Y()-PART[i].Get_Y();
				double dZ=PART[J].Get_Z()-PART[i].Get_Z();
				PART[i].Add(-dX/mindis*L, -dY/mindis*L, -dZ/mindis*L);
				PART[J].Add(dX/mindis*L, dY/mindis*L, dZ/mindis*L);
			}
		}//////////*/
	}/////MD終了

	delete [] index;
	delete [] MESH;
}

//直線を分割するさいの最適な分割数と分割距離の算出関数
void Model :: Set_calc_N_and_L(double dis,double le,int *N,double *L)
{
	double temp_N=dis/le;			//仮の分割数。leで割り切れたら一番いいけど、そうもいかないときがある
	int Ns=(int) temp_N;				//真の分割数
	double difference=temp_N-Ns;		//仮と真の差
	if(difference>0.5) Ns++;
	*L=dis/Ns;			//粒子の距離
	*N=Ns;
}

//円周分割数計算関数
int Model :: Set_calc_division_N_circle(double dis,double le)
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