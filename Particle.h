#pragma once
//
//粒子クラス
//

class Particle 
{
	double X;	//粒子座標
	double Y;
	double Z;

	double R;	//粒子半径
	
public:
	

	virtual void Set(double x, double y, double z, double r);	//座標と半径をセット
	virtual void Add(double x_add, double y_add, double z_add);	//加減算
	virtual void Multiply(double x_mul, double y_mul, double z_mul);//積除算
	virtual double Get_X();
	virtual double Get_Y();
	virtual double Get_Z();
	
	
};