#include "stdafx.h"
#include "Particle.h"



void Particle :: Set(double x, double y, double z, double r)
{
	X=x;
	Y=y;
	Z=z;
	R=r;
}

void Particle :: Add(double x_add, double y_add, double z_add)
{
	X+=x_add;
	Y+=y_add;
	Z+=z_add;
}

void Particle :: Multiply(double x_mul, double y_mul, double z_mul)
{
	X*=x_mul;
	Y*=y_mul;
	Z*=z_mul;
}

double Particle ::  Get_X()
{
	return X;
}

double Particle ::  Get_Y()
{
	return Y;
}

double Particle ::  Get_Z()
{
	return Z;
}