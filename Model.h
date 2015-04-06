#pragma once
//
//モデルクラス
//

#include "Particle.h"	//継承元をinclude
class Model : Particle
{
	vector<Particle> PART;
	vector<point3D> node;
	vector<element3D> nelm;
	Particle PART0;					//初期化されたParticleクラスのデータ
	ofstream fq;
	int point[4];

public:
	Model();
	~Model();
	//////////分子動力学関数//////////
	virtual void MD_2D(double le,int BstartID,int beforeN,int newN);
	virtual void MD_3D(double le,int BstartID,int beforeN,int newN,double r,double region[3][2]);
		
	//////////円筒作成関数//////////
	virtual int Set_cylinder(int initial_ID, double le, double R, double height);	
	//////////球作成関数//////////
	virtual int Set_sphere_function(int initial_ID, double le, double R);
		//円作成関数
		virtual int Set_circle_edge(int initial_number, double le, double R);
		//円内部充填関数
		virtual int Set_circle_in(int number, double le, double R, int initial_number);
		//円内部充填6ピース対称関数
		virtual int Set_circle_in_using_6_pieces(int number,double le,double R,int  initial_ID);
		//球作成関数
		virtual int Set_sphere(int number,double le,double R,int flag);
		//円柱表面座標作成
		virtual int Set_cylinder_face(int number,double le,double R,double height, int circle_start_id);
		//円柱内部充填関数
		virtual int Set_cylinder_in(int number,double le, double R, double height, int circle_start_id);
		//直線を分割するさいの最適な分割数と分割距離の算出関数
		virtual void Set_calc_N_and_L(double dis,double le,int *N,double *L);
		//円周分割数計算関数
		virtual int Set_calc_division_N_circle(double dis,double le);

	//////////ファイル出力関数//////////
	void writedata(ofstream &fp, int id, double x, double y,double z, int type,int materialID,int surface,double val,double vx,double vy,double vz,double P,double h,int toBEM);
	
	//////////モデルセット関数//////////
	virtual int Model_set();
		//////////球モデル関数//////////
		virtual int Set_sphere_model();
		//////////円筒モデル関数//////////
		virtual int Set_cylinder_model();
		//////////立方体モデル関数//////////
		virtual int Set_cube_model();
		//////////引っ張り試験モデル//////////
		virtual int Set_tensiontest_model();
			virtual void Set_point(int low_ID,int top_ID,int lefy_ID,int right_ID);
			virtual void Get_point(int &a,int &b,int &c,int &d);
		//////////アクチュエータモデル関数//////////
		virtual int Set_actuator_model();
		//////////ベンチマーク問題モデル関数//////////
		virtual int Set_benchmark_model();
	
};