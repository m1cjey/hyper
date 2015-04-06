#pragma once

/*
Rigidbodyクラス

剛体用クラス
*/
#include "define.h"

class Rigidbody
{
	double Rg[DIMENSION];	//剛体の初期重心座標
	double R[DIMENSION];	//剛体の重心座標
	double V[DIMENSION];	//剛体の平均速度
	double theta[DIMENSION];//剛体の回転核
	double theta2[DIMENSION];//長いから分けた
	double I0[3][3];
	double I[3][3];	//剛体のモーメントテンソル	
	double J[3][3];	//Iの逆行列
//	double F[3];	//剛体に働く力
//	double L[3];	//角運動量
	int rigid_num;			//自分の固体番号

	vector<mpselastic> PARTr;//剛体を構成する粒子
	vector<double> ri[3];	//中心と粒子riの相対距離
	double Qj[4];		//クオータニオン
	double Wj[3];		//角速度
	
	mpsconfig CON;
public:
	Rigidbody();
	virtual ~Rigidbody();
	//重心取得関数
	virtual void Get_initial_particle(vector<mpselastic> PARTg);
		virtual void _Get_initial_center();
		virtual void _Get_initial_moment();
		virtual void _Get_new_center();
		virtual void _Get_velocity();
//	virtual void Get_center_update(vector<mpselastic> &PARTup);

	virtual void Renew_part_r_v(vector<mpselastic> &PARTg);//粒子の位置を移動
	virtual void Get_rigid_move(vector<mpselastic> &PART);
		virtual void  _Get_I(double (&E)[3][3]);
		virtual void _Rigid_force(vector<mpselastic> &PART);
		virtual void _Calc_quaternion(double q[4],double* ri);
	//クオータニオン積
	virtual void _Quaternion_crossproduct(double A[4],double B[4],double (&Ans)[4]);
	//テンソル積
	virtual void _Tensorial_product(double A[3][3],double B[3][3],double (&Ans)[3][3]);
	//逆テンソル
	virtual void _Inverse_matrix(double A[3][3],double (&Ans)[3][3]);
};

