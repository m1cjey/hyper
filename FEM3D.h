//#include "stdafx.h"
#pragma once
class point3D//節点クラス
{
public:
	double r[DIMENSION];
	int material; //材質
	int boundary_condition; //境界条件 0=普通　1,2=固定境界
	int particleID;			//対応する粒子番号 存在しないときは-1を格納
	int remesh;				//リメッシュ領域に属するか、しないか
};

class element3D//要素クラス
{
public:
    int node[5];//要素を構築するnodeの番号  反時計まわり→てっぺんの順に格納
	int elm[5];//要素と隣接する要素番号 //表面は0と格納
	int sides[6+1];//要素を構成する辺番号
	double volume;//体積の６倍
	double r[DIMENSION];//ボロノイ点座標(外接球中心)
	double RR;//外接球半径の二乗
	int map;//マッピング
	int material;//材質
};

class sides//辺クラス
{
public:
	int node[2+1];//辺を構成する節点番号格納
	int boundary_condition; //境界条件 0=普通　1,2=固定境界
};

void VOLT3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,int *jnb,double TIME,vector<mpselastic> &PART,int fluid_number,int **nei,double *RP);
void ELECTRO3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double *V,double **Ee);
void potential_calculation(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,double TIME,vector<mpselastic> &PART,int fluid_number,int **nei,double **F);
void conductive_approximation(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpselastic> &PART,double **F,int fluid_number);
void import_J0_density(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void check_J0(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,double **current);
void Avector3D_node_eddy2(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double **A,int *jnb,int *depth,int **nei2,int *branch_num,double **old_A,double II,double dt,int mps_num,double *V,int t);
void calc_transitional_EM_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F,int KTJ);
void Bflux3D_node(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **A,double **B);
int locate3D(vector<point3D> &NODE,vector<element3D> &ELEM,int nelm,double xp,double yp,double zp);

void plot_magnetic_flux_density(mpsconfig &CON, vector<mpselastic> &PART, vector<point3D> &NODE, vector<element3D> &ELEM, int nelm, double **B, int t); //磁束密度をプロット

void poly3D2(vector<point3D> &NODE,vector<element3D> &ELEM,int *iv,int *kv,int ip,int *nelm,mpsconfig *CON);
void smoothingFs3D(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,double *Fs);
void smoothingF3D(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,double *F[3],int t);

//辺要素生成関数
int make_edge_element(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,int *jnb,int **nei,sides *SIDE,int *branch_num,int **nei2,int KTE);
void calc_static_magnetic_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F,int KTJ);
void calc_variable_magnetic_field(mpsconfig &CON,int node,int nelm,vector<point3D> &NODE, vector<element3D> &ELEM,int *jnb,double dt,double TIME,int t,int **nei,int KTE,vector<mpselastic> &PART,int fluid_number,double **F);
void set_boundary_condition3D_edge(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,int *dn,int *NN,double *PHAT,double *A);
void Avector3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP);
void Avector3D2(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,int *jnb,int *branch_num,double **current,double *RP,double II,int *depth,int t);
void Bflux3D_side(mpsconfig &CON, vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *A,double **B,double *RP);

//力
void direct_divT3D(mpsconfig &CON,vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int *jnb,int **nei);
void kelvin_force3D(mpsconfig &CON,vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int*jnb,int **nei);
void coroid_pole(mpsconfig &CON,vector<mpselastic> &PART,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Be,int fluid_number,double **F,double *RP,int*jnb,int **nei);
void NODE_F3D(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,int node,int nelm,double **Ee,int *jnb,int **nei,double *RP,vector<mpselastic> &PART,double **F,int fluid_number);

//電流密度計算関数(辺要素用)
void calc_current(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,int *jnb,int *branch_num,double **current,int *depth,double II);
//辺要素電流密度計算関数
void denryu_side(mpsconfig &CON,vector<point3D> &NODE,vector<element3D> &ELEM,sides *SIDE,int node,int nelm,int side_num,double *T,double **current,double J1);
