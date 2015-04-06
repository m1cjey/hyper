#pragma once

#define DIMENSION 3
#define A_X 0
#define A_Y 1
#define A_Z 2
#define ON 1
#define OFF 0
#define PI 3.14159265358979323846264338327950288
#define POSITIVE 0
#define NEGATIVE 1
#define ERR 1.0e-12
#define EPS 1.0e-12

//////////////粒子素材///////////////
//ルール三桁の数字で置く
//特殊粒子 000_099
#define GHOST 000 //計算しない粒子
#define FACE_P 001//表面

//気体粒子 100_199
#define GAS 100
#define AIR 101   //空気

//流体粒子 200_299
#define FLUID 200
#define WATER 201

//固体粒子 300_399
#define SOLID 300 //固体
#define WALL 301 //固定壁
#define ELASTIC 302//弾性体
#define MAGELAST 303//磁性エラストマー
#define MAGNET 304
#define IRON 305
#define COIL 306
#define ELECTRODE 307	//電極
#define PLATE 308
#define ELECTRODE1 309	//電極1(nanoe円柱電極)
#define ELECTRODE2 310	//電極2(nanoe平板電極)
#define TERMINAL1 311	//引っ張り試験用の端子
#define TERMINAL2 312	//引っ張り試験用の端子
#define MAGELAST2 313 //磁性エラストマー２
#define HYPERELAST 314

//GPU CG法における係数行列の格納方法
#define CSR_scl 0
#define CSR_vec 1
#define ELL 2

#define Diric 0
#define Neumn 1
#define BOTH  2
#define CONSTANT 0
#define LINER 1

#define SIGNIFY 18		//ファイル出力の有効数字

#define UNDEFINED -1	//未定義