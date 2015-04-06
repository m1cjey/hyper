#include "stdafx.h"
#include "Wall.h"

void Wall :: Set_boxmeshuWall(double size[3], double core[3])
{
	double **node=new double*[8];		//四角形の頂点
	for(int i=0;i<8;i++) node[i] = new double[3];
	for(int j=0;j<2;j++)
	{
		box_node[j*4][A_X]=core[A_X]-size[A_X]/2;
		box_node[j*4][A_Y]=core[A_Y]+size[A_Y]/2;
		box_node[j*4][A_Z]=core[A_Z]-((size[A_Z]/2)*(-1*j));
		box_node[j*4+1][A_X]=core[A_X]+size[A_X]/2;
		box_node[j*4+1][A_Y]=core[A_Y]+size[A_Y]/2;
		box_node[j*4+1][A_Z]=core[A_Z]-((size[A_Z]/2)*(-1*j));
		box_node[j*4+2][A_X]=core[A_X]-size[A_X]/2;
		box_node[j*4+2][A_Y]=core[A_Y]-size[A_Y]/2;
		box_node[j*4+2][A_Z]=core[A_Z]-((size[A_Z]/2)*(-1*j));
		box_node[j*4+3][A_X]=core[A_X]+size[A_X]/2;
		box_node[j*4+3][A_Y]=core[A_Y]-size[A_Y]/2;
		box_node[j*4+3][A_Z]=core[A_Z]-((size[A_Z]/2)*(-1*j));
	}
		
	for(int i=0;i<8;i++) delete[] node[i];
	delete[] node;
}

void Wall :: Set_flatWall()
{
	mpsconfig CON;					//初期化されたmpsconfigクラスのデータ
	double le=CON.get_distancebp();
	double x=30*le;	//x方向幅
	double y=30*le;	//y方向幅
	double z=-2*le;	//面の高さ

}