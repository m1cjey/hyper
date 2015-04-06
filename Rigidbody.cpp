#include "stdafx.h"
#include "Rigidbody.h"
#include "PART.h"
#include <fstream>
Rigidbody::Rigidbody(){}
Rigidbody::~Rigidbody(){
	PARTr.clear();
}

void Rigidbody::Get_initial_particle(vector<mpselastic> PARTg){
	//このPARTは1個の剛体を構成する粒子だけを含んでいる
	for(int i=0;i<PARTg.size();i++){	//剛体を構成する粒子情報の格納
		PARTr.push_back(PARTg[i]);
		ri[0].push_back(0);
		ri[1].push_back(0);
		ri[2].push_back(0);
	}
	_Get_initial_center();	//中心座標を取得
	_Get_initial_moment();	//初期形状のモーメントテンソル取得
	PARTg.clear();
	Qj[0]=1; Qj[1]=0; Qj[2]=0; Qj[3]=0;
	Wj[0]=0; Wj[1]=0; Wj[2]=0;
//	F[0]=0; F[1]=0; F[2]=0;
//	L[0]=0; L[1]=0; L[2]=0;
}

void Rigidbody::_Get_initial_center(){
	double X=0;	//剛体の中心座標
	double Y=0;
	double Z=0;
	
	for(int i=0;i<PARTr.size();i++){
		X+=PARTr[i].r[A_X];
		Y+=PARTr[i].r[A_Y];
		Z+=PARTr[i].r[A_Z];
	}
	Rg[A_X]=X/PARTr.size();	//中心座標格納
	Rg[A_Y]=Y/PARTr.size();
	Rg[A_Z]=Z/PARTr.size();
	R[A_X]=Rg[A_X];
	R[A_Y]=Rg[A_Y];
	R[A_Z]=Rg[A_Z];
}

void Rigidbody::_Get_initial_moment(){
	mpsconfig CON;
	double mass=CON.get_particle_mass();
	double X=0;	//原点から各粒子への距離
	double Y=0;
	double Z=0;
	for(int i=0; i<PARTr.size();i++){
		X+=PARTr[i].r[A_X]-Rg[A_X];
		Y+=PARTr[i].r[A_Y]-Rg[A_Y];
		Z+=PARTr[i].r[A_Z]-Rg[A_Z];		
	}
	I0[0][0]=mass*(pow(Y,2)+pow(Z,2)); I0[0][1]=-mass*Y*X;					I0[0][2]=-mass*Z*X;
	I0[1][0]=-mass*X*Z;				   I0[1][1]=mass*(pow(X,2)+pow(Z,2));	I0[1][2]=-mass*Z*Y;
	I0[2][0]=-mass*X*Z;				   I0[2][1]=-mass*Y*Z;					I0[2][2]=mass*(pow(X,2)*pow(Y,2));
}

void Rigidbody::_Get_new_center(){
	double X=0;	//剛体の中心座標
	double Y=0;
	double Z=0;
	
	for(int i=0;i<PARTr.size();i++){
		X+=PARTr[i].r[A_X];
		Y+=PARTr[i].r[A_Y];
		Z+=PARTr[i].r[A_Z];
	}
	R[A_X]=X/PARTr.size();	//中心座標格納
	R[A_Y]=Y/PARTr.size();
	R[A_Z]=Z/PARTr.size();

}
void Rigidbody::_Get_velocity(){
	double Vx=0;	//剛体の速度
	double Vy=0;
	double Vz=0;
	
	for(int i=0;i<PARTr.size();i++){
		Vx+=PARTr[i].u[A_X];
		Vy+=PARTr[i].u[A_Y];
		Vz+=PARTr[i].u[A_Z];
	}
	V[A_X]=Vx/PARTr.size();	//速度平均
	V[A_Y]=Vy/PARTr.size();
	V[A_Z]=Vz/PARTr.size();

}
/*void Rigidbody::Get_center_update(vector<mpselastic> &PARTup){
	double RA=0;
	double RI=0;
	double se=0.1e-9;
	PARTr.clear();	//粒子配列の初期化
	for(int i=0;i<PARTup.size();i++){	//剛体を構成する粒子情報の格納
		PARTr.push_back(PARTup[i]);
	}
	_Get_new_center();	//新しい中心座標を取得
	_Get_velocity();	//速度の平均
	////Aテンソル////
	A[0][0]=0.0;A[0][1]=0.0; A[0][2]=0.0;
	A[1][0]=0.0;A[1][1]=0.0; A[1][2]=0.0;
	A[2][0]=0.0;A[2][1]=0.0; A[2][2]=0.0;
	
	for(int i=0;i<PARTr.size();i++){
		A[0][0]+=(PARTup[i].r[A_X]-R[A_X])*(PARTup[i].initial_r[A_X]-Rg[A_X]); A[0][1]+=(PARTup[i].r[A_X]-R[A_X])*(PARTup[i].initial_r[A_Y]-Rg[A_Y]); A[0][2]+=(PARTup[i].r[A_X]-R[A_X])*(PARTup[i].initial_r[A_Z]-Rg[A_Z]);
		A[1][0]+=(PARTup[i].r[A_Y]-R[A_Y])*(PARTup[i].initial_r[A_X]-Rg[A_X]); A[1][1]+=(PARTup[i].r[A_Y]-R[A_Y])*(PARTup[i].initial_r[A_Y]-Rg[A_Y]); A[1][2]+=(PARTup[i].r[A_Y]-R[A_Y])*(PARTup[i].initial_r[A_Z]-Rg[A_Z]);
		A[2][0]+=(PARTup[i].r[A_Z]-R[A_Z])*(PARTup[i].initial_r[A_X]-Rg[A_X]); A[2][1]+=(PARTup[i].r[A_Z]-R[A_Z])*(PARTup[i].initial_r[A_Y]-Rg[A_Y]); A[2][2]+=(PARTup[i].r[A_Z]-R[A_Z])*(PARTup[i].initial_r[A_Z]-Rg[A_Z]);
//		cout<<"i="<<i<<"A[2][2]="<<A[2][2]<<endl;
	}
	if(A[0][0]<=se)A[0][0]=0; if(A[0][1]<=se)A[0][1]=0; if(A[0][2]<=se)A[0][2]=0;
	if(A[1][0]<=se)A[1][0]=0; if(A[1][1]<=se)A[1][1]=0; if(A[1][2]<=se)A[1][2]=0;
	if(A[2][0]<=se)A[2][0]=0; if(A[2][1]<=se)A[2][1]=0; if(A[2][2]<=se)A[2][2]=0;

//	cout<<"A[0][0]="<<A[0][0]<<" ,A[0][1]="<<A[0][1]<<" ,A[0][2]="<<A[0][2]<<endl;
//	cout<<"A[1][0]="<<A[1][0]<<" ,A[1][1]="<<A[1][1]<<" ,A[1][2]="<<A[1][2]<<endl;
//	cout<<"A[2][0]="<<A[2][0]<<" ,A[2][1]="<<A[2][1]<<" ,A[2][2]="<<A[2][2]<<endl;
	/////////////////
	//慣性モーメントの更新//
	I[0][0]=0.0;I[0][1]=0.0; I[0][2]=0.0;
	I[1][0]=0.0;I[1][1]=0.0; I[1][2]=0.0;
	I[2][0]=0.0;I[2][1]=0.0; I[2][2]=0.0;
	for(int i=0;i<PARTr.size();i++){
		I[0][0]+=pow(PARTup[i].r[A_Y]-R[A_Y],2)+pow(PARTup[i].r[A_Z]-R[A_Z],2); I[0][1]+=-(PARTup[i].r[A_X]-R[A_X])*(PARTup[i].r[A_Y]-R[A_Y]);			I[0][2]+=-(PARTup[i].r[A_X]-R[A_X])*(PARTup[i].r[A_Z]-R[A_Z]);
		I[1][0]+=-(PARTup[i].r[A_Y]-R[A_Y])*(PARTup[i].r[A_X]-R[A_X]);			I[1][1]+=pow(PARTup[i].r[A_X]-R[A_X],2)+pow(PARTup[i].r[A_Z]-R[A_Z],2); I[1][2]+=-(PARTup[i].r[A_Y]-R[A_Y])*(PARTup[i].r[A_Z]-R[A_Z]);
		I[2][0]+=-(PARTup[i].r[A_Z]-R[A_Z])*(PARTup[i].r[A_X]-R[A_X]);			I[2][1]+=-(PARTup[i].r[A_Z]-R[A_Z])*(PARTup[i].r[A_Y]-R[A_Y]);			I[2][2]+=pow(PARTup[i].r[A_X]-R[A_X],2)+pow(PARTup[i].r[A_Y]-R[A_Y],2);
	}
	////////////////////////
	//////Iの逆行列//////
	RI=I[0][0]*I[1][1]*I[2][2]+I[1][0]*I[2][1]*I[0][2]+I[2][0]*I[0][1]*I[1][2]-I[0][0]*I[2][1]*I[1][2]-I[2][0]*I[1][1]*I[0][2]-I[1][0]*I[0][1]*I[2][2];
		J[0][0]=(I[1][1]*I[2][2]-I[1][2]*I[2][1])/RI; J[0][1]=(I[0][2]*I[2][1]-I[0][1]*I[2][2])/RI; J[0][2]=(I[0][1]*I[1][2]-I[0][2]*I[1][1])/RI;
		J[1][0]=(I[1][2]*I[2][0]-I[1][0]*I[2][2])/RI; J[1][1]=(I[0][0]*I[2][2]-I[0][2]*I[2][0])/RI; J[1][2]=(I[0][2]*I[1][0]-I[0][0]*I[1][2])/RI;
		J[2][0]=(I[1][0]*I[2][1]-I[1][1]*I[2][0])/RI; J[2][1]=(I[0][1]*I[2][0]-I[0][0]*I[2][1])/RI; J[2][2]=(I[0][0]*I[1][1]-I[0][1]*I[1][0])/RI;
	//////////////////////

	///////回転角////////////
	for(int i=0;i<PARTr.size();i++){
		theta[A_X]+=(PARTup[i].u[A_Y]-V[A_Y])*(PARTup[i].r[A_Z]-R[A_Z])-(PARTup[i].u[A_Z]-V[A_Z])*(PARTup[i].r[A_Y]-R[A_Y]);
		theta[A_Y]+=(PARTup[i].u[A_Z]-V[A_Z])*(PARTup[i].r[A_X]-R[A_X])-(PARTup[i].u[A_X]-V[A_X])*(PARTup[i].r[A_Z]-R[A_Z]);
		theta[A_Z]+=(PARTup[i].u[A_X]-V[A_X])*(PARTup[i].r[A_Y]-R[A_Y])-(PARTup[i].u[A_Y]-V[A_Y])*(PARTup[i].r[A_X]-R[A_X]);
	}
	theta2[A_X]=J[0][0]*theta[A_X]+J[0][1]*theta[A_Y]+J[0][2]*theta[A_Z];
	theta2[A_Y]=J[1][0]*theta[A_X]+J[1][1]*theta[A_Y]+J[1][2]*theta[A_Z];
	theta2[A_Z]=J[2][0]*theta[A_X]+J[2][1]*theta[A_Y]+J[2][2]*theta[A_Z];
	/////////////////////////////

	//速度の更新//
	for(int i=0;i<PARTr.size();i++){
		PARTup[i].u[A_X]=((theta2[A_Y]*(PARTr[i].initial_r[A_Z]-Rg[A_Z]))-(theta2[A_Z]*(PARTr[i].initial_r[A_Y]-Rg[A_Y]))+V[A_X])*0.8;
		PARTup[i].u[A_Y]=((theta2[A_Z]*(PARTr[i].initial_r[A_X]-Rg[A_X]))-(theta2[A_X]*(PARTr[i].initial_r[A_Z]-Rg[A_Z]))+V[A_Y])*0.8;
		PARTup[i].u[A_Z]=((theta2[A_X]*(PARTr[i].initial_r[A_Y]-Rg[A_Y]))-(theta2[A_Y]*(PARTr[i].initial_r[A_X]-Rg[A_X]))+V[A_Z])*0.8;
	}//
	//////////////

	//位置の更新//
	B[0][0]=fabs(pow(A[0][0],2)+pow(A[1][0],2)+pow(A[2][0],2));    B[0][1]=fabs(A[0][0]*A[0][1]+A[1][0]*A[1][1]+A[2][0]*A[2][1]); B[0][2]=fabs(A[0][0]*A[0][2]+A[1][0]*A[1][2]+A[2][0]*A[2][2]);
	B[1][0]=fabs(A[0][1]*A[0][0]+A[1][1]*A[1][0]+A[2][1]*A[2][0]); B[1][1]=fabs(pow(A[0][1],2)+pow(A[1][1],2)+pow(A[2][1],2));	  B[1][2]=fabs(A[0][1]*A[0][2]+A[1][1]*A[1][2]+A[2][1]*A[2][2]);
	B[2][0]=fabs(A[0][2]*A[0][0]+A[1][2]*A[1][0]+A[2][2]*A[2][0]); B[2][1]=fabs(A[0][2]*A[0][1]+A[1][2]*A[1][1]+A[2][2]*A[2][1]); B[2][2]=fabs(pow(A[0][2],2)+pow(A[1][2],2)+pow(A[2][2],2));


	B[0][0]=pow(B[0][0],0.5); B[0][1]=pow(B[0][1],0.5); B[0][2]=pow(B[0][2],0.5);
	B[1][0]=pow(B[1][0],0.5); B[1][1]=pow(B[1][1],0.5); B[1][2]=pow(B[1][2],0.5);
	B[2][0]=pow(B[2][0],0.5); B[2][1]=pow(B[2][1],0.5); B[2][2]=pow(B[2][2],0.5);//
//	B[0][0]=A[0][0]; B[0][1]=A[0][1]; B[0][2]=A[0][2];
//	B[1][0]=A[1][0]; B[1][1]=A[1][1]; B[1][2]=A[1][2];
//	B[2][0]=A[2][0]; B[2][1]=A[2][1]; B[2][2]=A[2][2];//
//	cout<<"B[0][0]="<<B[0][0]<<" ,B[0][1]="<<B[0][1]<<" ,B[0][2]="<<B[0][2]<<endl;
//	cout<<"B[1][0]="<<B[1][0]<<" ,B[1][1]="<<B[1][1]<<" ,B[1][2]="<<B[1][2]<<endl;
//	cout<<"B[2][0]="<<B[2][0]<<" ,B[2][1]="<<B[2][1]<<" ,B[2][2]="<<B[2][2]<<endl;

	RA=B[0][0]*B[1][1]*B[2][2]+B[1][0]*B[2][1]*B[0][2]+B[2][0]*B[0][1]*B[1][2]-B[0][0]*B[2][1]*B[1][2]-B[2][0]*B[1][1]*B[0][2]-B[1][0]*B[0][1]*B[2][2];

	C[0][0]=(B[1][1]*B[2][2]-B[1][2]*B[2][1])/RA; C[0][1]=(B[0][2]*B[2][1]-B[0][1]*B[2][2])/RA; C[0][2]=(B[0][1]*B[1][2]-B[0][2]*B[1][1])/RA;
	C[1][0]=(B[1][2]*B[2][0]-B[1][0]*B[2][2])/RA; C[1][1]=(B[0][0]*B[2][2]-B[0][2]*B[2][0])/RA; C[1][2]=(B[0][2]*B[1][0]-B[0][0]*B[1][2])/RA;
	C[2][0]=(B[1][0]*B[2][1]-B[1][1]*B[2][0])/RA; C[2][1]=(B[0][1]*B[2][0]-B[0][0]*B[2][1])/RA; C[2][2]=(B[0][0]*B[1][1]-B[0][1]*B[1][0])/RA;
//	A[0][0]=fabs(A[0][0]); A[0][1]=fabs(A[0][1]); A[0][2]=fabs(A[0][2]);
//	A[1][0]=fabs(A[1][0]); A[1][1]=fabs(A[1][1]); A[1][2]=fabs(A[1][2]);
//	A[2][0]=fabs(A[2][0]); A[2][1]=fabs(A[2][1]); A[2][2]=fabs(A[2][2]);

	RA=A[0][0]*A[1][1]*A[2][2]+A[1][0]*A[2][1]*A[0][2]+A[2][0]*A[0][1]*A[1][2]-A[0][0]*A[2][1]*A[1][2]-A[2][0]*A[1][1]*A[0][2]-A[1][0]*A[0][1]*A[2][2];

	C[0][0]=(A[1][1]*A[2][2]-A[1][2]*A[2][1])/RA; C[0][1]=(A[0][2]*A[2][1]-A[0][1]*A[2][2])/RA; C[0][2]=(A[0][1]*A[1][2]-A[0][2]*A[1][1])/RA;
	C[1][0]=(A[1][2]*A[2][0]-A[1][0]*A[2][2])/RA; C[1][1]=(A[0][0]*A[2][2]-A[0][2]*A[2][0])/RA; C[1][2]=(A[0][2]*A[1][0]-A[0][0]*A[1][2])/RA;
	C[2][0]=(A[1][0]*A[2][1]-A[1][1]*A[2][0])/RA; C[2][1]=(A[0][1]*A[2][0]-A[0][0]*A[2][1])/RA; C[2][2]=(A[0][0]*A[1][1]-A[0][1]*A[1][0])/RA;
//	cout<<"C[0][0]="<<C[0][0]<<" ,C[0][1]="<<C[0][1]<<" ,C[0][2]="<<C[0][2]<<endl;
//	cout<<"C[1][0]="<<C[1][0]<<" ,C[1][1]="<<C[1][1]<<" ,C[1][2]="<<C[1][2]<<endl;
//	cout<<"C[2][0]="<<C[2][0]<<" ,C[2][1]="<<C[2][1]<<" ,C[2][2]="<<C[2][2]<<endl;//
//	cout<<"A[0][0]="<<A[0][0]<<" ,A[0][1]="<<A[0][1]<<" ,A[0][2]="<<A[0][2]<<endl;
//	cout<<"A[1][0]="<<A[1][0]<<" ,A[1][1]="<<A[1][1]<<" ,A[1][2]="<<A[1][2]<<endl;
//	cout<<"A[2][0]="<<A[2][0]<<" ,A[2][1]="<<A[2][1]<<" ,A[2][2]="<<A[2][2]<<endl;
//	cout<<"B[0][0]="<<B[0][0]<<" ,B[0][1]="<<B[0][1]<<" ,B[0][2]="<<B[0][2]<<endl;
//	cout<<"B[1][0]="<<B[1][0]<<" ,B[1][1]="<<B[1][1]<<" ,B[1][2]="<<B[1][2]<<endl;
//	cout<<"B[2][0]="<<B[2][0]<<" ,B[2][1]="<<B[2][1]<<" ,B[2][2]="<<B[2][2]<<endl;//
//	cout<<"C[0][0]="<<C[0][0]<<" ,C[0][1]="<<C[0][1]<<" ,C[0][2]="<<C[0][2]<<endl;//
	for(int i=0;i<PARTr.size();i++){
		PARTup[i].r[A_X]=(A[0][0]*C[0][0]+A[0][1]*C[1][0]+A[0][2]*C[2][0])*(PARTup[i].initial_r[A_X]-Rg[A_X])+(A[0][0]*C[0][1]+A[0][1]*C[1][1]+A[0][2]*C[2][1])*(PARTup[i].initial_r[A_Y]-Rg[A_Y])+(A[0][0]*C[0][2]+A[0][1]*C[1][2]+A[0][2]*C[2][2])*(PARTup[i].initial_r[A_Z]-Rg[A_Z])+R[A_X];
		PARTup[i].r[A_Y]=(A[1][0]*C[0][0]+A[1][1]*C[1][0]+A[1][2]*C[2][0])*(PARTup[i].initial_r[A_X]-Rg[A_X])+(A[1][0]*C[0][1]+A[1][1]*C[1][1]+A[1][2]*C[2][1])*(PARTup[i].initial_r[A_Y]-Rg[A_Y])+(A[1][0]*C[0][2]+A[1][1]*C[1][2]+A[1][2]*C[2][2])*(PARTup[i].initial_r[A_Z]-Rg[A_Z])+R[A_Y];
		PARTup[i].r[A_Z]=(A[2][0]*C[0][0]+A[2][1]*C[1][0]+A[2][2]*C[2][0])*(PARTup[i].initial_r[A_X]-Rg[A_X])+(A[2][0]*C[0][1]+A[2][1]*C[1][1]+A[2][2]*C[2][1])*(PARTup[i].initial_r[A_Y]-Rg[A_Y])+(A[2][0]*C[0][2]+A[2][1]*C[1][2]+A[2][2]*C[2][2])*(PARTup[i].initial_r[A_Z]-Rg[A_Z])+R[A_Z];//
//		PARTup[i].r[A_X]=(PARTup[i].initial_r[A_X]-Rg[A_X])+R[A_X];
//		PARTup[i].r[A_Y]=(PARTup[i].initial_r[A_Y]-Rg[A_Y])+R[A_Y];
//		PARTup[i].r[A_Z]=(PARTup[i].initial_r[A_Z]-Rg[A_Z])+R[A_Z];
	}
//	cout<<"D[0][0]="<<(A[0][0]*C[0][0]+A[0][1]*C[1][0]+A[0][2]*C[2][0])<<" ,D[0][1]="<<(A[0][0]*C[0][1]+A[0][1]*C[1][1]+A[0][2]*C[2][1])<<" ,D[0][2]="<<(A[0][0]*C[0][2]+A[0][1]*C[1][2]+A[0][2]*C[2][2])<<endl;
//	cout<<"D[1][0]="<<(A[1][0]*C[0][0]+A[1][1]*C[1][0]+A[1][2]*C[2][0])<<" ,D[1][1]="<<(A[1][0]*C[0][1]+A[1][1]*C[1][1]+A[1][2]*C[2][1])<<" ,D[1][2]="<<(A[1][0]*C[0][2]+A[1][1]*C[1][2]+A[1][2]*C[2][2])<<endl;
//	cout<<"D[2][0]="<<(A[2][0]*C[0][0]+A[2][1]*C[1][0]+A[2][2]*C[2][0])<<" ,D[2][1]="<<(A[2][0]*C[0][1]+A[2][1]*C[1][1]+A[2][2]*C[2][1])<<" ,D[2][2]="<<(A[2][0]*C[0][2]+A[2][1]*C[1][2]+A[2][2]*C[2][2])<<endl;

	///////////////////座標確認用
	for(int i=0;i<PARTup.size();i++){
		cout<<"i="<<i<<", rx="<<PARTup[i].r[A_X]<<", ry="<<PARTup[i].r[A_Y]<<", rz="<<PARTup[i].r[A_Z]<<endl;//
	}//
	
}//*/

void Rigidbody::Renew_part_r_v(vector<mpselastic> &PARTg){
	//このPARTは1個の剛体を構成する粒子だけを含んでいる
	double quate[4]={Qj[0],Qj[1],Qj[2],Qj[3]};	//前ステップのクオータニオン
	double w[3]={Wj[0],Wj[1],Wj[2]};	//前ステップの角速度
	double pR[3]={R[0],R[1],R[2]};		//前ステップの重心
	for(int i=0;i<PARTg.size();i++){
	PARTr[i]=PARTg[i];
	}
	_Get_new_center();	//新しい中心座標を取得
	_Get_velocity();	//速度の平均
	//相対位置ベクトル
	for(int i=0;i<PARTr.size();i++){
	ri[A_X][i]=PARTr[i].r[A_X]-R[A_X];
	ri[A_Y][i]=PARTr[i].r[A_Y]-R[A_Y];
	ri[A_Z][i]=PARTr[i].r[A_Z]-R[A_Z];
	}

	//回転後の座標
	for(int i=0;i<PARTr.size();i++){
		double a[3]={ri[A_X][i],ri[A_Y][i],ri[A_Z][i]};

		_Calc_quaternion(quate,a);

		ri[A_X][i]=a[0];
		ri[A_Y][i]=a[1];
		ri[A_Z][i]=a[2];
	}
	//各粒子の座標
	for(int i=0;i<PARTr.size();i++){
		PARTr[i].r[A_X]=pR[A_X]+ri[A_X][i];
		PARTr[i].r[A_Y]=pR[A_Y]+ri[A_Y][i];
		PARTr[i].r[A_Z]=pR[A_Z]+ri[A_Z][i];
	}
	//各粒子の速度
	for(int i=0;i<PARTr.size();i++){
		PARTr[i].u[A_X]=V[A_X]+w[A_Y]*PARTr[i].r[A_Z]-w[A_Z]*PARTr[i].r[A_Y];
		PARTr[i].u[A_Y]=V[A_Y]+w[A_Z]*PARTr[i].r[A_X]-w[A_X]*PARTr[i].r[A_Z];
		PARTr[i].u[A_Z]=V[A_Z]+w[A_X]*PARTr[i].r[A_Y]-w[A_Y]*PARTr[i].r[A_X];
	}
	PARTg=PARTr;
	R[0]=pR[0];
	R[1]=pR[1];
	R[2]=pR[2];
}

void Rigidbody::Get_rigid_move(vector<mpselastic> &PARTo){
	//このPARTはすべての粒子を含んでいる	
	//周辺粒子との力を計算,新しい重心の位置と回転を求める
	_Rigid_force(PARTo);
}



void Rigidbody::_Rigid_force(vector<mpselastic> &PART){
	double dis=CON.get_distancebp();
	double dt=CON.get_dt();
	double mass=CON.get_particle_mass();
	double F[3]={0,0,0},L[3]={0,0,0};	//力
	cout<<"剛体に働く力の計算"<<endl;
	double k=0.8,l=0.8,kt=0.8;	//ばね定数、減衰定数、摩擦定数
	double tr[3]={0,0,0};

	ofstream disp("disp.dat");
	ofstream disq("disq.dat");
	for(int i=0;i<PART.size();i++){
		double f[3]={0,0,0};
		if(PART[i].type==TERMINAL1){//自分の剛体番号を持つ粒子のみ計算
			int neighboursN0=static_cast<int>(PART[i].get_current_neighboursID().size());
			for(int k=0;k<neighboursN0;k++)
			{
				int j=PART[i].get_current_neighboursID()[k];
				if(PART[j].type!=TERMINAL1){/////////
					double rij[3]={PART[j].r[A_X]-PART[i].r[A_X] ,PART[j].r[A_Y]-PART[i].r[A_Y] ,PART[j].r[A_Z]-PART[i].r[A_Z]};
					double disij=pow(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2],0.5);
					disp<<"PART[j].r[A_Z]="<<PART[j].r[A_Z]<<", PART[i].r[A_Z]="<<PART[i].r[A_Z]<<endl;
					disp<<"disij="<<disij<<",dis="<<dis<<",i="<<i<<",j="<<j<<", rij[0]="<<rij[0]<<", rij[1]="<<rij[1]<<", rij[2]="<<rij[2]<<endl;
					if(disij<dis){
						//接触により働く力なので重力は除外
						f[0]+=(k*(dis-disij)*rij[0]/disij)+(l*(-PART[j].u[A_X]+PART[i].u[A_X]));
						f[1]+=(k*(dis-disij)*rij[1]/disij)+(l*(-PART[j].u[A_Y]+PART[i].u[A_Y]));
						f[2]+=(k*(dis-disij)*rij[2]/disij)+(l*(-PART[j].u[A_Z]+PART[i].u[A_Z]));
	//					L[0]+=rij[1]*f[2]-rij[2]*f[1];
	//					L[1]+=rij[2]*f[0]-rij[0]*f[2];
	//					L[2]+=rij[0]*f[1]-rij[1]*f[0];
					disq<<"接触中 f[0]="<<f[0]<<",f[1]="<<f[1]<<",f[2]="<<f[2]<<",disij="<<disij<<endl;
					}
				}
			}
	//		cout<<"f[0]="<<f[0]<<"f[1]="<<f[1]<<"f[2]="<<f[2]<<endl;
			//粒子iに働く力が求まった
			double r[3]={PART[i].r[A_X]-R[A_X],PART[i].r[A_Y]-R[A_Y],PART[i].r[A_Z]-R[A_Z]};
			tr[0]+=r[0];
			tr[1]+=r[1];
			tr[2]+=r[2];
			L[0]+=r[1]*f[2]-r[2]*f[1];
			L[1]+=r[2]*f[0]-r[0]*f[2];
			L[2]+=r[0]*f[1]-r[1]*f[0];
			f[2]+=(CON.get_g()*mass);
			F[0]+=f[0];
			F[1]+=f[1];
			F[2]+=f[2];
		}
	}
	disp.close();
	disq.close();
	cout<<"tr[0]="<<tr[0]<<",tr[1]="<<tr[1]<<",tr[2]="<<tr[2]<<endl;
	cout<<"F[0]="<<F[0]<<",F[1]="<<F[1]<<",F[2]="<<F[2]<<endl;
	cout<<"L[0]="<<L[0]<<",L[1]="<<L[1]<<",L[2]="<<L[2]<<endl;
	cout<<"剛体中心の並進と回転運動"<<endl;
	//剛体の並進と回転を計算
	R[A_X]+=F[A_X]*dt*dt/mass;
	R[A_Y]+=F[A_Y]*dt*dt/mass;
	R[A_Z]+=F[A_Z]*dt*dt/mass;
	cout<<"R[A_X]="<<R[A_X]<<", R[A_Y]="<<R[A_Y]<<", R[A_Z]="<<R[A_Z]<<endl;
	L[A_X]*=dt;
	L[A_Y]*=dt;
	L[A_Z]*=dt;
	double E[3][3];
	_Get_I(E);
	Wj[A_X]=E[0][0]*L[A_X]+E[0][1]*L[A_Y]+E[0][2]*L[A_Z];
	Wj[A_Y]=E[1][0]*L[A_X]+E[1][1]*L[A_Y]+E[1][2]*L[A_Z];
	Wj[A_Z]=E[2][0]*L[A_X]+E[2][1]*L[A_Y]+E[2][2]*L[A_Z];
	cout<<"Wj[A_X]="<<Wj[A_X]<<", Wj[A_Y]="<<Wj[A_Y]<<", Wj[A_Z]="<<Wj[A_Z]<<endl;
	double ax[3]={0,0,0}; //回転軸
	double axl=pow(pow(Wj[0],2)+pow(Wj[1],2)+pow(Wj[2],2),0.5);

	if(axl==0){ax[A_X]=0; ax[A_Y]=0; ax[A_Z]=0;} 
	else {ax[A_X]=Wj[A_X]/axl; ax[A_Y]=Wj[A_Y]/axl; ax[A_Z]=Wj[A_Z]/axl;}

	cout<<"ax[A_X]="<<ax[A_X]<<", ax[A_Y]="<<ax[A_Y]<<", ax[A_Z]="<<ax[A_Z]<<endl;
	double ang=pow(pow(Wj[A_X]*dt,2)+pow(Wj[A_Y]*dt,2)+pow(Wj[A_Z]*dt,2),0.5);
	cout<<"ang="<<ang<<endl;
	double p[4];
	p[0]=cos(ang/2);
	p[1]=ax[A_X]*sin(ang/2);
	p[2]=ax[A_Y]*sin(ang/2);
	p[3]=ax[A_Z]*sin(ang/2);
	cout<<"p[0]="<<p[0]<<", p[1]="<<p[1]<<", p[2]="<<p[2]<<", p[3]="<<p[3]<<endl;
	double qc[4]={Qj[0],Qj[1],Qj[2],Qj[3]};
	double qA[4]={0,0,0,0};
	_Quaternion_crossproduct(qc,p,qA);
	//剛体のクオータニオン更新
	Qj[0]=qA[0]; Qj[1]=qA[1]; Qj[2]=qA[2]; Qj[3]=qA[3];
	cout<<"Qj[0]="<<Qj[0]<<", Qj[1]="<<Qj[1]<<", Qj[2]="<<Qj[2]<<", Qj[3]="<<Qj[3]<<endl;
}

void Rigidbody::_Get_I(double (&E)[3][3]){
	double A[3][3];	//回転行列
	double B[3][3];	//回転行列の逆行列
	double C[3][3]; //初期慣性テンソルの逆行列
	double D[3][3];	//回転行列と慣性逆テンソルの積
	//回転行列
/*	A[0][0]=cos(theta[A_Z])*cos(theta[A_Y]); A[0][1]=-sin(theta[A_Z])*cos(theta[A_X])+cos(theta[A_Z])*sin(theta[A_Y])*sin(theta[A_X]); A[0][2]=sin(theta[A_Z])*sin(theta[A_X])+cos(theta[A_Z])*sin(theta[A_Y])*cos(theta[A_X]);
	A[1][0]=sin(theta[A_Z])*cos(theta[A_Y]); A[1][1]=cos(theta[A_Z])*cos(theta[A_X])+sin(theta[A_Z])*sin(theta[A_Y])*sin(theta[A_X]);  A[1][2]=-cos(theta[A_Z])*sin(theta[A_X])+sin(theta[A_Z])*sin(theta[A_Y])*cos(theta[A_X]);
	A[2][0]=-sin(theta[A_Y]);				 A[2][1]=cos(theta[A_Y])*sin(theta[A_X]);												   A[2][2]=cos(theta[A_Y])*cos(theta[A_X]);*/
	A[0][0]=1-2*(pow(Qj[2],2)+pow(Qj[3],2)); A[0][1]=2*(Qj[1]*Qj[2]-Qj[3]*Qj[0]);     A[0][2]=2*(Qj[3]*Qj[1]+Qj[2]*Qj[0]);
	A[1][0]=2*(Qj[1]*Qj[2]+Qj[3]*Qj[0]);	 A[1][1]=1-2*(pow(Qj[3],2)+pow(Qj[1],2)); A[1][2]=2*(Qj[2]*Qj[3]-Qj[1]*Qj[0]);
	A[2][0]=2*(Qj[3]*Qj[1]-Qj[2]*Qj[0]);	 A[2][1]=2*(Qj[2]*Qj[3]+Qj[1]*Qj[0]);	  A[2][2]=1-2*(pow(Qj[1],2)+pow(Qj[2],2));
	//回転行列の逆行列
	_Inverse_matrix(A,B);
//	double RA=A[0][0]*A[1][1]*A[2][2]+A[1][0]*A[2][1]*A[0][2]+A[2][0]*A[0][1]*A[1][2]-A[0][0]*A[2][1]*A[1][2]-A[2][0]*A[1][1]*A[0][2]-A[1][0]*A[0][1]*A[2][2];
//	B[0][0]=(A[1][1]*A[2][2]-A[1][2]*A[2][1])/RA; B[0][1]=(A[0][2]*A[2][1]-A[0][1]*A[2][2])/RA; B[0][2]=(A[0][1]*A[1][2]-A[0][2]*A[1][1])/RA;
//	B[1][0]=(A[1][2]*A[2][0]-A[1][0]*A[2][2])/RA; B[1][1]=(A[0][0]*A[2][2]-A[0][2]*A[2][0])/RA; B[1][2]=(A[0][2]*A[1][0]-A[0][0]*A[1][2])/RA;
//	B[2][0]=(A[1][0]*A[2][1]-A[1][1]*A[2][0])/RA; B[2][1]=(A[0][1]*A[2][0]-A[0][0]*A[2][1])/RA; B[2][2]=(A[0][0]*A[1][1]-A[0][1]*A[1][0])/RA;
	//慣性テンソルの逆行列
	_Inverse_matrix(I0,C);
//	double RC=I0[0][0]*I0[1][1]*I0[2][2]+I0[1][0]*I0[2][1]*I0[0][2]+I0[2][0]*I0[0][1]*I0[1][2]-I0[0][0]*I0[2][1]*I0[1][2]-I0[2][0]*I0[1][1]*I0[0][2]-I0[1][0]*I0[0][1]*I0[2][2];
//	C[0][0]=(I0[1][1]*I0[2][2]-I0[1][2]*I0[2][1])/RC; C[0][1]=(I0[0][2]*I0[2][1]-I0[0][1]*I0[2][2])/RC; C[0][2]=(I0[0][1]*I0[1][2]-I0[0][2]*I0[1][1])/RC;
//	C[1][0]=(I0[1][2]*I0[2][0]-I0[1][0]*I0[2][2])/RC; C[1][1]=(I0[0][0]*I0[2][2]-I0[0][2]*I0[2][0])/RC; C[1][2]=(I0[0][2]*I0[1][0]-I0[0][0]*I0[1][2])/RC;
//	C[2][0]=(I0[1][0]*I0[2][1]-I0[1][1]*I0[2][0])/RC; C[2][1]=(I0[0][1]*I0[2][0]-I0[0][0]*I0[2][1])/RC; C[2][2]=(I0[0][0]*I0[1][1]-I0[0][1]*I0[1][0])/RC;
	//回転行列と慣性逆テンソルの積
	_Tensorial_product(A,C,D);
//	D[0][0]=A[0][0]*C[0][0]+A[0][1]*C[1][0]+A[0][2]*C[2][0]; D[0][1]=A[0][0]*C[0][1]+A[0][1]*C[1][1]+A[0][2]*C[2][1]; D[0][2]=A[0][0]*C[0][2]+A[0][1]*C[1][2]+A[0][2]*C[2][2];
//	D[1][0]=A[1][0]*C[0][0]+A[1][1]*C[1][0]+A[1][2]*C[2][0]; D[1][1]=A[1][0]*C[0][1]+A[1][1]*C[1][1]+A[1][2]*C[2][1]; D[1][2]=A[1][0]*C[0][2]+A[1][1]*C[1][2]+A[1][2]*C[2][2];
//	D[2][0]=A[2][0]*C[0][0]+A[2][1]*C[1][0]+A[2][2]*C[2][0]; D[2][1]=A[2][0]*C[0][1]+A[2][1]*C[1][1]+A[2][2]*C[2][1]; D[2][2]=A[2][0]*C[0][2]+A[2][1]*C[1][2]+A[2][2]*C[2][2];
	//さらに逆行列との積
	_Tensorial_product(D,B,E);
//	E[0][0]=D[0][0]*B[0][0]+D[0][1]*B[1][0]+D[0][2]*B[2][0]; E[0][1]=D[0][0]*B[0][1]+D[0][1]*B[1][1]+D[0][2]*B[2][1]; E[0][2]=D[0][0]*B[0][2]+D[0][1]*B[1][2]+D[0][2]*B[2][2];
//	E[1][0]=D[1][0]*B[0][0]+D[1][1]*B[1][0]+D[1][2]*B[2][0]; E[1][1]=D[1][0]*B[0][1]+D[1][1]*B[1][1]+D[1][2]*B[2][1]; E[1][2]=D[1][0]*B[0][2]+D[1][1]*B[1][2]+D[1][2]*B[2][2];
//	E[2][0]=D[2][0]*B[0][0]+D[2][1]*B[1][0]+D[2][2]*B[2][0]; E[2][1]=D[2][0]*B[0][1]+D[2][1]*B[1][1]+D[2][2]*B[2][1]; E[2][2]=D[2][0]*B[0][2]+D[2][1]*B[1][2]+D[2][2]*B[2][2];
}
void Rigidbody::_Calc_quaternion(double q[4],double* ri){
	double QR[4]; //クオータニオンとベクトルの内積
	double RI[4]={0,ri[0],ri[1],ri[2]};
	double Vq[4]={q[0],-q[1],-q[2],-q[3]}; //クオータニオン共役
	double Ans[4]; //移動後の位置

	_Quaternion_crossproduct(q,RI,QR);

	//QR×クオータニオンの共役
	_Quaternion_crossproduct(QR,Vq,Ans);

	ri[0]=Ans[1];
	ri[1]=Ans[2];
	ri[2]=Ans[3];

}

void Rigidbody::_Quaternion_crossproduct(double A[4],double B[4],double (&Ans)[4]){
	Ans[0]=A[0]*B[0]-A[1]*B[1]-A[2]*B[2]-A[3]*B[3];
	Ans[1]=A[0]*B[1]+A[1]*B[0]+A[2]*B[3]-A[3]*B[2];
	Ans[2]=A[0]*B[2]-A[1]*B[3]+A[2]*B[0]+A[3]*B[1];
	Ans[3]=A[0]*B[3]+A[1]*B[2]-A[2]*B[1]+A[3]*B[0];
}

void Rigidbody::_Tensorial_product(double A[3][3],double B[3][3],double (&Ans)[3][3]){
	Ans[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0]; Ans[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1]+A[0][2]*B[2][1]; Ans[0][2]=A[0][0]*B[0][2]+A[0][1]*B[1][2]+A[0][2]*B[2][2];
	Ans[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0]; Ans[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1]; Ans[1][2]=A[1][0]*B[0][2]+A[1][1]*B[1][2]+A[1][2]*B[2][2];
	Ans[2][0]=A[2][0]*B[0][0]+A[2][1]*B[1][0]+A[2][2]*B[2][0]; Ans[2][1]=A[2][0]*B[0][1]+A[2][1]*B[1][1]+A[2][2]*B[2][1]; Ans[2][2]=A[2][0]*B[0][2]+A[2][1]*B[1][2]+A[2][2]*B[2][2];
}
void Rigidbody::_Inverse_matrix(double A[3][3],double (&Ans)[3][3]){
	double RA=A[0][0]*A[1][1]*A[2][2]+A[1][0]*A[2][1]*A[0][2]+A[2][0]*A[0][1]*A[1][2]-A[0][0]*A[2][1]*A[1][2]-A[2][0]*A[1][1]*A[0][2]-A[1][0]*A[0][1]*A[2][2];
	Ans[0][0]=(A[1][1]*A[2][2]-A[1][2]*A[2][1])/RA; Ans[0][1]=(A[0][2]*A[2][1]-A[0][1]*A[2][2])/RA; Ans[0][2]=(A[0][1]*A[1][2]-A[0][2]*A[1][1])/RA;
	Ans[1][0]=(A[1][2]*A[2][0]-A[1][0]*A[2][2])/RA; Ans[1][1]=(A[0][0]*A[2][2]-A[0][2]*A[2][0])/RA; Ans[1][2]=(A[0][2]*A[1][0]-A[0][0]*A[1][2])/RA;
	Ans[2][0]=(A[1][0]*A[2][1]-A[1][1]*A[2][0])/RA; Ans[2][1]=(A[0][1]*A[2][0]-A[0][0]*A[2][1])/RA; Ans[2][2]=(A[0][0]*A[1][1]-A[0][1]*A[1][0])/RA;
}