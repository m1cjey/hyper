
//#fasetは#facetの間違い？？？
//ファイル読み込みこれでできている？？？

#include "stdafx.h"

#define  A_R 0
#define  A_t 1 //θ


using namespace std;
//.nodeデータ取得関数
void tetgen_function::GetPointList(vector<tetgen_node> &NODE, tetgenio &in, tetgenio &out)
{
	NODE.clear();

	tetgen_node temp;
	for(int i=0;i<out.numberofpoints;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<3;n++)	temp.r[n]=out.pointlist[i*3+n];
		temp.attribute=(int)out.pointattributelist[i];
		temp.boundary=out.pointmarkerlist[i];

		NODE.push_back(temp);
	}
}


//.eleデータ取得関数(簡易版)
void tetgen_function::GetTetrahedronList(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<4;n++)	temp.node[n]=out.tetrahedronlist[i*4+n];
		temp.attribute=0; //恐らく、未定義で代入すると変な数値が入ってバグるのでここでは0にしとく

		ELEM.push_back(temp);
	}
}
//.eleデータ取得関数(MAGELAST用)
void tetgen_function::GetMTetrahedronList(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<4;n++)	temp.node[n]=out.tetrahedronlist[i*4+n];
		temp.attribute=MAGELAST; //恐らく、未定義で代入すると変な数値が入ってバグるのでここでは0にしとく

		ELEM.push_back(temp);
	}
}

void tetgen_function::Geteleattribute(vector<tetgen_node> &NODE,vector<tetgen_element> &ELEM,tetgenio &out)
{
 for(int i=0;i<out.numberoftetrahedra;i++)
    {
        int M1=(int)out.pointattributelist[out.tetrahedronlist[i*4+0]];
		int M2=(int)out.pointattributelist[out.tetrahedronlist[i*4+1]];
		int M3=(int)out.pointattributelist[out.tetrahedronlist[i*4+2]];
		int M4=(int)out.pointattributelist[out.tetrahedronlist[i*4+3]];
	
		///4頂点すべてが同じ材質なら要素もそれにならう。
		///ひとつでも異なっていたら空気と定義
		if(M1==M2 && M2==M3 && M3==M4) out.tetrahedronattributelist[i]=M1;
		else if(M1!=AIR && M2!=AIR && M3!=AIR && M4!=AIR) 
		{
			if(M1!=ELASTIC && M2!=ELASTIC && M3!=ELASTIC && M4!=ELASTIC) out.tetrahedronattributelist[i]=MAGELAST;
			//if(CON.get_model_number()==15) ELEM[i].material=AIR;
			else out.tetrahedronattributelist[i]=ELASTIC; //else ELEM[i].material=FLUID;
			if(M1==IRON || M2==IRON || M3==IRON || M4==IRON) out.tetrahedronattributelist[i]=IRON;//コイルの要素はコイル接点の内側になるようにする
		}
		else out.tetrahedronattributelist[i]=AIR;
    }
}


//.eleデータ取得関数(材質・要素要素関係含む)
void tetgen_function::GetTetrahedronList_full(vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	ELEM.clear();

	tetgen_element temp;
	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		temp.id=i+out.firstnumber;										//節点番号
		for(int n=0;n<4;n++) temp.node[n]=out.tetrahedronlist[i*4+n];	//構成節点
		for(int n=0;n<4;n++) temp.nei_elem[n]=out.neighborlist[i*4+n];	//要素-要素関係
		temp.attribute=(int)out.tetrahedronattributelist[i];				//材質
		//temp.volume=out.tetrahedronvolumelist[i];							//体積(取得できない)

		ELEM.push_back(temp);
	}
}


//.faceデータ取得関数
void tetgen_function::GetFacetList(vector<tetgen_facet> &FACE, tetgenio &in, tetgenio &out, int boundary)
{
	//※boundarymarkerは引数として与えている。PLCモードで作ったメッシュではないため境界が出力されない。

	FACE.clear();

	tetgen_facet temp;
	for(int i=0;i<out.numberoftrifaces;i++)
	{
		temp.id=i+out.firstnumber;
		for(int n=0;n<3;n++)	temp.node[n]=out.trifacelist[i*3+n];
	//	temp.boundary=out.trifacemarkerlist[i];
		temp.boundary=boundary;

		FACE.push_back(temp);
	}//*/

	//for(int i=0;i<3;i++)
	//for(int n=0;n<3;n++)	cout<<out.trifacelist[i*3+n]<<endl;
}


//.nodeファイル作成関数
void tetgen_function::MakeNodeFile(mpsconfig &CON, vector<tetgen_node> &NODE, char *filename)
{
	//cout<<filename<<" 出力"<<endl;
	ofstream fout(filename);
	
	fout<<"#node"<<endl;
	fout<<(int)NODE.size()<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	}
	fout.close();
}


//.nodeファイル作成関数
void tetgen_function::MakeNodeFile_NonAttributeAndBoundary(mpsconfig &CON, vector<tetgen_node> &NODE, char *filename)
{
	cout<<filename<<" 出力"<<endl;

	ofstream fout(filename);
	
	fout<<"#node"<<endl;
	fout<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<endl;
	}

	fout.close();
}


//.eleファイル作成関数
void tetgen_function::MakeElemFile(mpsconfig &CON, vector<tetgen_element> &ELEM, char *filename)
{
	ofstream fout(filename);

	fout<<(int)ELEM.size()<<"\t"<<"4"<<"\t"<<"1"<<endl;
	for(int i=0;i<(int)ELEM.size();i++)
	{
		fout<<ELEM[i].id<<"\t";
		for(int n=0;n<4;n++)	fout<<ELEM[i].node[n]<<"\t";
		fout<<ELEM[i].attribute<<endl;
	}

	fout.close();
}


//.faceファイル作成関数
void tetgen_function::MakeFaceFile(mpsconfig &CON, vector<tetgen_facet> &FACE, char *filename)
{
	ofstream fout(filename);

	fout<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fout<<FACE[i].id;
		for(int n=0;n<3;n++)	fout<<" "<<FACE[i].node[n];
		fout<<" "<<FACE[i].boundary;
		fout<<endl;
	}

	fout.clear();
}


//.polyファイル作成関数
void tetgen_function::MakePolyFile(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, char *filename)
{
	//cout<<filename<<" 出力"<<endl;

	ofstream fout(filename);

	//node list (ここでは出力しない)
	fout<<"#node"<<endl;
	fout<<"0"<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	//fout<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	//for(int i=0;i<(int)NODE.size();i++)
	//{
		//fout<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	//}

	//faset list
	fout<<"#faset"<<endl;
//	fout<<"#facet"<<endl;
	fout<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fout<<"1"<<" "<<"0"<<" "<<FACE[i].boundary<<endl;
		//fout<<"1"<<" "<<FACE[i].boundary<<endl;
		//fout<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<" "<<FACE[i].boundary<<endl;
		fout<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<endl;
	}

	//hole list
	fout<<"#hole"<<endl;
	fout<<"0"<<endl;
	fout<<endl;

	////////////////////////////region attributeの決定 (配列に格納してから出力する)/////////////////////////////
	//材質の指定は，境界内にある一点の座標を決め，そこの材質を指定することで，同じ境界内にある要素が全てその材質になる．
	//水は分裂を伴い，材質の指定が困難であるため，ここでは行わない．
	//後の材質の修正において，未定義となっている要素を水要素とする．

	vector<region_attribute_list> REGION;
	region_attribute_list temp;
	temp.id=0;
	temp.region_number=0;
	temp.region_attribute=0;
	
	if(CON.get_model_number()==2)
	{
		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//MAGNET
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//磁石の中心点
		temp.region_number=MAGNET;
		temp.region_attribute=MAGNET;
		REGION.push_back(temp);
	}

	//磁性エラストマー
	if(CON.get_model_number()==6)
	{
		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//空気(内部)
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=CON.get_magnet_r();
		temp.r[A_Z]=CON.get_magnet_Z()+CON.get_magnet_H()/2+0.0001;	//磁石の上部
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//コイル
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//磁石の中心点
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);
		int i_no;
		for(int i=0;i<NODE.size();i++){
			if(NODE[i].part_no==1)  i_no=i;//95
		}
		//MAGELAST
		temp.id+=1;
		temp.r[A_X]=NODE[i_no].r[A_X]-0.0001; //MAGELASTのノードが最初に追加される
		temp.r[A_Y]=NODE[i_no].r[A_Y];
		temp.r[A_Z]=NODE[i_no].r[A_Z]+0.0001;
		temp.region_number=MAGELAST;
		temp.region_attribute=MAGELAST;
		REGION.push_back(temp);
	}
	//磁性エラストマー
	if(CON.get_model_number()==7 && CON.get_model_number()==1 && CON.get_model_number()==11) 
	{
		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//空気2
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z()-CON.get_magnet_H()/2-0.001;	//磁石の中心点
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//コイル
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//磁石の中心点
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);

		int i_no;
		for(int i=0;i<NODE.size();i++)
		{
			if(NODE[i].part_no==1) i_no=i;
		}
		//MAGELAST

		temp.id+=1;
		temp.r[A_X]=NODE[i_no].r[A_X]-0.0001; //MAGELASTのノードが最初に追加される
		temp.r[A_Y]=NODE[i_no].r[A_Y];
		temp.r[A_Z]=NODE[i_no].r[A_Z]+0.0001;
		temp.region_number=MAGELAST;
		temp.region_attribute=MAGELAST;
		REGION.push_back(temp);
	}
	if(CON.get_model_number()==8)
	{
		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//コイル
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//磁石の中心点
		temp.region_number=COIL;
		temp.region_attribute=COIL;
		REGION.push_back(temp);

		int i_no;
		for(int i=0;i<NODE.size();i++){
			if(NODE[i].part_no==2)  i_no=i;//
		}
/*		//MAGELAST
		temp.id+=1;
		temp.r[A_X]=NODE[i_no].r[A_X]-0.00001; //MAGELASTのノードが最初に追加される
		temp.r[A_Y]=NODE[i_no].r[A_Y];
		temp.r[A_Z]=NODE[i_no].r[A_Z]+0.00001;
		temp.region_number=MAGELAST;
		temp.region_attribute=MAGELAST;
		REGION.push_back(temp);//*/
	}

	if(CON.get_model_number()==23	|| CON.get_model_number()==24)
	{
		double le=CON.get_distancebp();
		//空気
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_ZU()*0.9;	//解析領域上限の9割のところ
		temp.region_number=AIR;
		temp.region_attribute=AIR;
		REGION.push_back(temp);

		//MAGNET
		temp.id+=1;
		temp.r[A_X]=0;
		temp.r[A_Y]=0;
		temp.r[A_Z]=CON.get_magnet_Z();	//磁石の中心点
		temp.region_number=MAGNET;
		temp.region_attribute=MAGNET;
		REGION.push_back(temp);

		int i_no;
		for(int i=0;i<NODE.size();i++){
			if(NODE[i].part_no==2)  i_no=i;//
		}
		//MAGELAST
		temp.id+=1;
		temp.r[A_X]=NODE[i_no].r[A_X]+le; //MAGELASTのノードが最初に追加される
		temp.r[A_Y]=NODE[i_no].r[A_Y]+le;
		temp.r[A_Z]=NODE[i_no].r[A_Z]+le;
		temp.region_number=MAGELAST;
		temp.region_attribute=MAGELAST;
		REGION.push_back(temp);
	}


	////////////////////////////////////////////////////////////////////*/

	//region attribute list
	fout<<"#region attribute"<<endl;
	fout<<(int)REGION.size()<<endl;
	for(int i=0;i<(int)REGION.size();i++)
	{
		fout<<REGION[i].id<<" "<<REGION[i].r[A_X]<<" "<<REGION[i].r[A_Y]<<" "<<REGION[i].r[A_Z]<<" "<<REGION[i].region_number<<" "<<REGION[i].region_attribute<<endl;
	}

	fout.close();
}

//.smeshファイル作成関数
void tetgen_function::MakeSmeshFile(mpsconfig &CON, vector<tetgen_facet> &FACE, char *filename)
{
	double le=CON.get_distancebp();

	ofstream fsmesh(filename);

	//node list (ここでは出力しない)
	fsmesh<<"#node"<<endl;
	fsmesh<<"0"<<" "<<"3"<<" "<<"1"<<" "<<"1"<<endl;
	//fsmesh<<(int)NODE.size()<<" "<<"3"<<" "<<"0"<<" "<<"0"<<endl;
	//for(int i=0;i<(int)NODE.size();i++)
	//{
		//fsmesh<<NODE[i].id<<" "<<NODE[i].r[A_X]<<" "<<NODE[i].r[A_Y]<<" "<<NODE[i].r[A_Z]<<" "<<NODE[i].attribute<<" "<<NODE[i].boundary<<endl;
	//}

	//faset list
	fsmesh<<"#faset"<<endl;
//	fsmesh<<"#facet"<<endl;
	fsmesh<<(int)FACE.size()<<" "<<"1"<<endl;
	for(int i=0;i<(int)FACE.size();i++)
	{
		fsmesh<<"3"<<" "<<FACE[i].node[A_X]<<" "<<FACE[i].node[A_Y]<<" "<<FACE[i].node[A_Z]<<" "<<FACE[i].boundary<<endl;
	}

	//hole list
	fsmesh<<"#hole"<<endl;
	fsmesh<<"0"<<endl;
	fsmesh<<endl;

	//region attribute list
	fsmesh<<"#region attribute"<<endl;
	fsmesh<<"2"<<endl;
	fsmesh<<"1"<<" "<<"0"<<" "<<"0"<<" "<<le*2<<" "<<"1"<<" "<<"1"<<endl;
	fsmesh<<"2"<<" "<<"0"<<" "<<"0"<<" "<<-le*2<<" "<<"2"<<" "<<"2"<<endl;
	
	fsmesh.close();
}


//弾性体境界面作成
void tetgen_function::SetElastBoundary(mpsconfig &CON, vector<mpselastic> &PART, tetgen_config &TET, vector<tetgen_node> &NODEe, vector<tetgen_facet> &FACEe)
{
	cout<<"磁性エラストマー境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMe;
	vector<int> trans;

	tetgen_node temp;
	temp.id=0;

	double le=CON.get_distancebp();	//粒子間距離
	int type;
	int part_no;
	
	//磁性エラストマー
/*		for(int i=0;i<(int)PART.size();i++)
		{
			if(PART[i].type==MAGELAST)
			{
				part_no=i;
				for(int d=0;d<3;d++) temp.r[d]=PART[i].r[d];
				type=PART[i].type;
				temp.boundary=0;	
				if(PART[i].surface==0) temp.attribute=MAGELAST;
				else if(PART[i].surface==1) temp.attribute=FACE_P;
				trans.push_back(part_no);
				NODEe.push_back(temp);
				temp.id+=1;
			}
			
		}*/
	//磁性エラストマー
	int count=0;
		for(int i=0;i<(int)PART.size();i++)
		{
			if(PART[i].type==MAGELAST)
			{
				
				part_no=i;
				for(int d=0;d<3;d++) temp.r[d]=PART[i].r[d];
				type=PART[i].type;
				temp.boundary=0;	
				if(PART[i].surface==0)
				{
					temp.attribute=MAGELAST;
					count++;
				}
				else if(PART[i].surface==1) 
				{
					temp.attribute=FACE_P;
					//temp.boundary=1;
				}
				trans.push_back(part_no);
				NODEe.push_back(temp);
				temp.id+=1;
			}
			
		}
		cout<<"count"<<count<<endl;
	//nodeファイル作成
	cout<<"MREnode作成-----";
	MakeNodeFile(CON, NODEe, "MAGELAST.node");
	cout<<"OK"<<endl;

	//.nodeファイル読み取り
	in.load_node("MAGELAST");

	//まずは流体節点のみで分割
	cout<<"MREメッシュ分割-----";
	tetrahedralize("", &in, &out);
	cout<<"OK"<<endl;
		// i 要素内へ節点を追加(今の所無くても機能している)
		// f .faceファイルに境界ではない面も含める
		// e .edgeファイルの出力(ONにするとなぜか止まってしまう)
		// n .neighファイルの出力

	//出力
	out.save_nodes("MAGELAST_whole");
	out.save_elements("MAGELAST_whole");



	//////////////////ここまでで流体節点のみを使って、すべての要素が繋がった凸なメッシュができた(fluid_wholeで確認可能)*/


	///////////////不要な要素の削除

	//.nodeの取得
	GetPointList(NODEe, in, out);
	//.eleの取得
	GetMTetrahedronList(ELEMe, in, out);

	//長い要素の除去
	DelThinTetrahedron(CON, TET, NODEe, ELEMe, in, out);

	//節点-要素関係
	SetRelation_NodeElem(CON, NODEe, ELEMe);
	//要素-要素関係
	SetRelation_ElemElem(CON, NODEe, ELEMe);
	//要素-要素関係より弾性体表面取得
	if(CON.get_model_number()==23)
	{
		cout<<"要素除去後の流体表面定義";
		{
			//初期化
			FACEe.clear();	
			int id=0;	//id用(表面の数)
			for(int i=0;i<(int)ELEMe.size();i++)
			{
				for(int n=0;n<4;n++)
				{
					if(int j=0;j<(int)NODEe[ELEMe[i].node[n]].nei_elem.size();j++)
					{
						tetgen_facet temp;	
						int c=0;
						for(int f=0;f<4;f++)
						{
							if(n!=f)
							{
								temp.node[c]=ELEMe[i].node[f];
								c++;	
							}
						}
						temp.id=id;
						temp.boundary=MAGELAST;
						FACEe.push_back(temp);
						id++;
					}
				}	
			}
			cout<<"----------OK"<<endl;}
	}
	else{GetFacetList_from_neigh(CON, ELEMe, FACEe);}

	//粒子番号をNODEe側に代入
	for(int i=0;i<NODEe.size();i++)	NODEe[i].part_no=trans[i];

	//表面を構成する節点を選択し，配列番号を詰める 流体内部も粒子の節点を使う場合はコメントアウト
	//SelectFaceNode(CON, NODEe, FACEe);

	/////////////////要素確認用ファイル///////////////////////////////////
	out.save_nodes("boundary_MAGELAST");	//fluid.2.nodeと同じファイル
	MakeElemFile(CON, ELEMe, "boundary_MAGELAST.ele");
	MakeFaceFile(CON, FACEe, "boundary_MAGELAST.face");
	////////////////ここまででエラストマーのメッシュが切れた//////////////////////*/
//	out.save_elements("boundary_MAGELAST");
//	out.save_faces("boundary_MAGELAST");
}


//長い要素の除去
void tetgen_function::DelThinTetrahedron(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"不要な要素の削除  ";

	double le=CON.get_distancebp();
	double delL;
	int flag=0;

	//del_length.datがあれば読み込み
	ifstream del("del_length.dat");
	if(!del)
	{
		cout<<"tetgen_configより削除判定辺長さを決定  ";
		delL=TET.del_length;
	}
	else
	{
		cout<<"del_length.datより削除判定辺長さを決定  ";
		del>>delL;
		flag=1;		//ファイルから読み込んだらフラグON
	}
	del.close();

	if(flag==1)//ファイルの数字を戻しておく
	{
		ofstream del2("del_length.dat");
		del2<<TET.del_length<<endl;
		del2.close();
	}

	cout<<"del_length="<<delL<<endl;


	cout<<"除去前要素数: "<<(int)ELEM.size()<<endl;

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag1=UNDEFINED;
		int flag2=UNDEFINED;
		int flag3=UNDEFINED;	//1で削除
		int del=OFF;
		int count=0;

		//4点が表面節点で構成されていればフラグ1ON
		count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].attribute==FACE_P)	count+=1;
		}
		if(count==4)	flag1=ON;
		else			flag1=OFF;//*/

		/*//4点が内部節点で構成されていればフラグ2ON
		count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==FRFLUID)	count+=1;
		}
		if(count==4)	flag2=ON;
		else			flag2=OFF;//*/

		//一つでも長い辺があればフラグON
		for(int n1=0;n1<4;n1++)
		{
			for(int n2=n1+1;n2<4;n2++)
			{
				double dis=Distance(NODE[ELEM[i].node[n1]+out.firstnumber], NODE[ELEM[i].node[n2]+out.firstnumber]);
				if(dis>le*delL)	flag3=ON;
			}
		}//*/

		if(flag3==ON)	del=ON;
		if(flag1==ON)	del=ON;
		//if(flag3==2)	del=ON;

		//削除
		if(del==ON)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//イテレータ初期化
			it+=i;				//iを指定
			it=ELEM.erase(it);	//削除してイテレータを返す
			del_count++;
		}
		else i++;
	}

	//要素番号の振り直し
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"除去後要素数: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"の要素を削除 ----------OK"<<endl;
}


//不要な要素の除去(外側ダミー節点法用)
void tetgen_function::DelTetrahedron_OutsideDummy(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"不要な要素の削除"<<endl;
	cout<<"除去前要素数: "<<(int)ELEM.size()<<endl;

	double le=CON.get_distancebp();

	int del_count=0;
	int i=0;
	while(i<(int)ELEM.size())
	{
		int flag=0;	//1で削除
		
		//1つでもダミー(空気)節点があればフラグON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==AIR)
			{
				flag=1;
				break;
			}
		}//*/

		/*//4点が表面節点で構成されていればフラグON
		int count=0;
		for(int n=0;n<4;n++)
		{
			if(NODE[ELEM[i].node[n]].boundary==BOFLUID)	count+=1;
		}
		if(count==4)	flag=1;//*/


		if(flag==1)
		{
			vector<tetgen_element>::iterator it=ELEM.begin();	//イテレータ初期化
			it+=i;				//iを指定
			it=ELEM.erase(it);	//削除してイテレータを返す
			del_count++;
		}
		else i++;
	}

	//要素番号の振り直し
	for(int i=0;i<(int)ELEM.size();i++)	ELEM[i].id=i+out.firstnumber;

	cout<<"除去後要素数: "<<(int)ELEM.size()<<endl;
	cout<<del_count<<"の要素を削除 ----------OK"<<endl;
}


//節点-要素関係(tetgenioのedgeリストから取得) ※直前にtetgenioからデータを取得しておくこと
void tetgen_function::SetRelation_NodeNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	cout<<"節点-節点関係";

	//一応初期化
	for(int i=0;i<(int)NODE.size();i++)
	{
		NODE[i].nei_node.clear();
	}

	for(int i=0;i<out.numberofedges;i++)
	{
		int node1=out.edgelist[i*2+0];
		int node2=out.edgelist[i*2+1];
		
		NODE[node1].nei_node.push_back(node2);
		NODE[node2].nei_node.push_back(node1);
	}

	cout<<"----------OK"<<endl;
}


//節点-要素関係
void tetgen_function::SetRelation_NodeElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	cout<<"節点-要素関係";

	//一応初期化
	for(int i=0;i<(int)NODE.size();i++)
	{
		NODE[i].nei_node.clear();
		NODE[i].nei_elem.clear();
	}

	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)
		{
			NODE[ELEM[i].node[n]].nei_elem.push_back(i);	//要素を追加
		}
	}

	/*//出力 および最大数・最小数の出力
	//int max=0;
	//int min=(int)NODE[0].nei_elem.size();
	ofstream fout("neigh_node-elem.dat");
	for(int i=0;i<(int)NODE.size();i++)
	{
		fout<<NODE[i].id;
		for(int n=0;n<(int)NODE[i].nei_elem.size();n++)
		{
			fout<<" "<<NODE[i].nei_elem[n];
		}
		fout<<endl;
		//最大最小更新
		//if(max<(int)NODE[i].nei_elem.size())	max=(int)NODE[i].nei_elem.size();
		//if(min>(int)NODE[i].nei_elem.size())	min=(int)NODE[i].nei_elem.size();
	}
	fout.clear();
	//cout<<"最大数: "<<max<<endl;
	//cout<<"最小数: "<<min<<endl;
	//*/

	cout<<"----------OK"<<endl;
}


//要素-要素関係
void tetgen_function::SetRelation_ElemElem(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	cout<<"要素-要素関係";
	
	vector<int> nei_all;
	
	//初期化
	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)	ELEM[i].nei_elem[n]=-2;	//未定義を-2としておく
	}

	for(int i=0;i<(int)ELEM.size();i++)
	{
		//4つの節点の節点-要素関係にある要素を格納する
		nei_all.clear();
		
		for(int n=0;n<4;n++)
		{
			for(int j=0;j<(int)NODE[ELEM[i].node[n]].nei_elem.size();j++)
			{
				int elem=NODE[ELEM[i].node[n]].nei_elem[j];
				if(elem!=i)	nei_all.push_back(elem);	//要素iは除く
			}
		}//nei_allに格納完了(同じ要素番号が含まれている可能性あり)

		//確認用
		//if(i==10000)	for(int j=0;j<(int)nei_all.size();j++)	cout<<nei_all[j]<<endl;

		//面を探索
		for(int ni=0;ni<4;ni++)
		{
			//まずは要素iの面を指定
			int face[3];	//3つ番号で面を指定
			int c=0;//数え上げ変数

			for(int f=0;f<4;f++)
			{
				if(ni!=f)
				{
					face[c]=ELEM[i].node[f];
					c++;
				}
			}//face[3]にn番目の面が格納

			//面を探索
			int correct_nei=-1;
			for(int j=0;j<(int)nei_all.size();j++)
			{
				int count=0;//このカウントが3になれば確定

				for(int nj=0;nj<4;nj++)
				{
					int node_j=ELEM[nei_all[j]].node[nj];
					for(int f=0;f<3;f++)
					{
						if(node_j==face[f])	count++;
					}	
				}
				if(count==3)
				{
					correct_nei=nei_all[j];
					break;
				}
			}//もし見つからなかったらcorrecr_neiには-1が入っている

			ELEM[i].nei_elem[ni]=correct_nei;
		}
	}//*/

	/*//出力
	ofstream fout("neigh.dat");
	for(int i=0;i<(int)ELEM.size();i++)
	{
		fout<<ELEM[i].id;
		for(int n=0;n<4;n++)
		{
			fout<<" "<<ELEM[i].nei_elem[n];		
		}
		fout<<endl;
	}
	fout.clear();
	//*/

	cout<<"----------OK"<<endl;
}


//要素除去後の流体表面定義
void tetgen_function::GetFacetList_from_neigh(mpsconfig &CON, vector<tetgen_element> &ELEM, vector<tetgen_facet> &FACE)
{
	cout<<"要素除去後の流体表面定義";

	//初期化
	FACE.clear();
	
	int id=0;	//id用(表面の数)
	
	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int n=0;n<4;n++)
		{
			if(ELEM[i].nei_elem[n]==-1)==-1//対面が存在しない→表面
			{
				tetgen_facet temp;	
				int c=0;

				for(int f=0;f<4;f++)
				{
					if(n!=f)
					{
						temp.node[c]=ELEM[i].node[f];
						c++;
					}
				}

				temp.id=id;
				temp.boundary=MAGELAST;
				FACE.push_back(temp);
				id++;
			}
		}
	}

	cout<<"----------OK"<<endl;
}


//表面節点以外を削除，表面節点に番号を振りなおす
void tetgen_function::SelectFaceNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE)
{
	//boundaryをフラグに使わせてもらう。
	//boundary==-1は表面を構成していない節点→削除
	//boundary==-2は表面を構成している節点→新しい節点番号を与える
	
	//とりあえず初期化
	for(int i=0;i<(int)NODE.size();i++)	NODE[i].boundary=-1;
	
	//FACEにある節点はflagを-2に
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			NODE[FACE[i].node[n]].boundary=-2;
		}
	}

	//新しい節点番号の決定
	int id=0;
	for(int i=0;i<(int)NODE.size();i++)
	{
		if(NODE[i].boundary==-2)
		{
			NODE[i].boundary=id;
			id++;
		}
	}
	//表面節点にはboundaryに新しい節点番号が入る。内部節点には-1が入る

	//FACEの構成節点番号の変換
	for(int i=0;i<(int)FACE.size();i++)
	{
		for(int n=0;n<3;n++)
		{
			FACE[i].node[n]=NODE[FACE[i].node[n]].boundary;	//ここで-1や-2となるものはないはず。あれば上の処理が間違っている
		}
	}

	//内部節点の削除 NODEの節点番号の変換 boundaryを元に戻す
	int k=0;
	while(k<(int)NODE.size())
	{
		if(NODE[k].boundary==-1)
		{
			vector<tetgen_node>::iterator it=NODE.begin();	//イテレータ初期化
			it+=k;				//k番目を指定
			it=NODE.erase(it);	//削除してイテレータを返す
		}
		else
		{
			NODE[k].id=NODE[k].boundary;
			NODE[k].boundary=WATER;
			k++;
		}
	}

}


//ダミー節点を削除
void tetgen_function::DelDummyNode(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, int num_dummy)
{
	//NODEの中には前半に流体表面節点，後半にダミー節点が固まって格納されているので，後半のダミー節点の部分のみを消せばよい
	//ダミー要素を消してから表面データを取得しているので，表面を構成する節点番号の変更は不要

	//ダミー変数の数だけpopbackで末尾の要素から消す
	for(int i=0;i<num_dummy;i++)
	{
		NODE.pop_back();
	}
}


//空気境界面作成
void tetgen_function::SetAirBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	cout<<"空気境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

//	vector<tetgen_element> ELEMa;

	tetgen_node temp_n;
	temp_n.id=0;
	temp_n.attribute=AIR;
	temp_n.boundary=0;

	//分割数決定  1割だけオフセットしている
	double dL=TET.fine_air;
	int lx = int((CON.get_XL()-0.1*dL)/dL);
	int ux = int((CON.get_XR()+0.1*dL)/dL);
	int ly = int((CON.get_YD()-0.1*dL)/dL);
	int uy = int((CON.get_YU()+0.1*dL)/dL);
	int lz = int((CON.get_ZD()-0.1*dL)/dL);
	int uz = int((CON.get_ZU()+0.1*dL)/dL);

	//接点データ作成
	//静電無化
	if(CON.get_model_number()==14)
	{
		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//解析領域の端に来たとき節点を置く
					{
						temp_n.r[A_X]=x*dL;
						temp_n.r[A_Y]=y*dL;
						temp_n.r[A_Z]=z*dL;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}//*/

	//レイリー分裂
	if(CON.get_model_number()==20)
	{
		temp_n.boundary=ELECTRODE1;

		for(int z=lz;z<=uz;z++)
		{
			for(int y=ly;y<=uy;y++)
			{
				for(int x=lx;x<=ux;x++)
				{
					if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//解析領域の端に来たとき節点を置く
					{
						temp_n.r[A_X]=x*dL;
						temp_n.r[A_Y]=y*dL;
						temp_n.r[A_Z]=z*dL;
						NODEa.push_back(temp_n);
						temp_n.id+=1;
					}
				}
			}
		}
	}

	//磁性エラストマー
	//解析領域の端にきた時に節点を置いているがforループ使わなくても良いのでは？？？
	//あと円筒形に置いたほうが良い
	
	{
//		temp_n.boundary=ELECTRODE1;

		//for(int z=lz;z<=uz;z++)
		//{
		//	for(int y=ly;y<=uy;y++)
		//	{
		//		for(int x=lx;x<=ux;x++)
		//		{
		//			if(x==lx || x==ux || y==ly || y==uy || z==lz || z==uz)	//解析領域の端に来たとき節点を置く
		//			{
		//				temp_n.r[A_X]=x*dL;
		//				temp_n.r[A_Y]=y*dL;
		//				temp_n.r[A_Z]=z*dL;
		//				NODEa.push_back(temp_n);
		//				temp_n.id+=1;
		//			}
		//		}
		//	}
		//}
		/////////////////////////////////////////////////////////////////////////////
		int divN[3];
		divN[A_R]=20;
		divN[A_t]=20;
		divN[A_Z]=20;
		double regionR[2]={0.0, CON.get_RU()};
		double regionZ[2]={CON.get_ZD(), CON.get_ZU()};
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
					//中心点
	temp_n.r[A_X]=0;
	temp_n.r[A_Y]=0;
	temp_n.r[A_Z]=Zmin;					//解析領域の底面
	NODEa.push_back(temp_n);
	temp_n.id+=1;

	for(int n=1;n<=divN[A_R];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			double r=divL[A_R]*n;
			double theta=divL[A_t]*m;
			temp_n.r[A_X]=r*cos(theta);
			temp_n.r[A_Y]=r*sin(theta);
			temp_n.r[A_Z]=Zmin;					//解析領域の底面
			NODEa.push_back(temp_n);
			temp_n.id+=1;				//対応する粒子が存在しないから-1を格納
		}
	}
	//////////////////////////////////上面
					//中心点
	temp_n.r[A_X]=0;
	temp_n.r[A_Y]=0;
	temp_n.r[A_Z]=Zmax;					//解析領域の底面
	NODEa.push_back(temp_n);
			temp_n.id+=1;	
	for(int n=1;n<=divN[A_R];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			double r=divL[A_R]*n;
			double theta=divL[A_t]*m;
			temp_n.r[A_X]=r*cos(theta);
			temp_n.r[A_Y]=r*sin(theta);
			temp_n.r[A_Z]=Zmax;					//解析領域の底面
			NODEa.push_back(temp_n);
			temp_n.id+=1;					//対応する粒子が存在しないから-1を格納
		}
	}

	//側面
	double RR=divL[A_R]*divN[A_R];
	for(int n=1;n<divN[A_Z];n++)
	{
		for(int m=0;m<divN[A_t];m++)
		{
			double theta=divL[A_t]*m;
			temp_n.r[A_X]=RR*cos(theta);
			temp_n.r[A_Y]=RR*sin(theta);;
			temp_n.r[A_Z]=Zmin+n*divL[A_Z];		
			NODEa.push_back(temp_n);
			temp_n.id+=1;						//対応する粒子が存在しないから-1を格納
		}
	}

		/////////////////////////////////////////////////////////////////////////////
	}

	//nodeファイル作成
	MakeNodeFile(CON, NODEa, "NODEa1.node");
	//.eleの取得
//	GetTetrahedronList(ELEMa, in, out);

	in.load_node("NODEa1");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_air1");
	out.save_elements("boundary_air1");
	out.save_faces("boundary_air1");
	//*/

	//境界面データ取得
	GetFacetList(FACEa, in, out, AIR);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEa1, "FACEa1.face");
	//長い要素の除去
//	DelThinTetrahedron(CON, TET, NODEa, ELEMa, in, out);

	//節点-要素関係
//	SetRelation_NodeElem(CON, NODEa, ELEMa);
	//要素-要素関係
//	SetRelation_ElemElem(CON, NODEa, ELEMa);
	//要素-要素関係より弾性体表面取得
//	GetFacetList_from_neigh(CON, ELEMa, FACEa);
}


//水滴・エラストマー表面付近追加節点
void tetgen_function::SetAirFineBoundary(mpsconfig &CON, vector<mpselastic> &PART, tetgen_config &TET, vector<tetgen_node> &NODEa, vector<tetgen_facet> &FACEa)
{
	//水滴の表面の法線方向に何層かのメッシュ層を作成する

	cout<<"電磁石表面付近追加節点"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=AIR;
	temp.boundary=0;

	//分割数決定
	double le=CON.get_distancebp();
	double R=CON.get_magnet_r()+5*le;		//円筒半径
	double uz=CON.get_magnet_Z()+CON.get_magnet_H()/2+5*le;		//円筒最大高さ
	double lz=CON.get_magnet_Z()-CON.get_magnet_H()/2-5*le;		//円筒最小高さ -5*le
	double dx=le;					//円筒長さ方向メッシュ粗さ

		//上下側面でテーパーを付けたい時は分けるべきだがまとめてもいいかも・・・

	//中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=uz;
	temp.boundary=0;
	NODEa.push_back(temp);
	temp.id+=1;

	//上面
	for(double r=dx;r<R+0.1*dx;r+=dx)
	{
		int nr=static_cast<int>(2.0*PI*r/dx);//正何角形で近似するか（角度方向の分割数）
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=uz;
			temp.boundary=0;
			NODEa.push_back(temp);
			temp.id+=1;
		}
	}

	//側面	//R*=1.005;	//わずかに太くする(メッシュが繋がるのを防ぐため)・・・これはいらん
	for(double z=uz-dx;z>lz;z-=dx)
	{
			int nr=static_cast<int>(2.0*PI*R/dx);//r=dRの時はnr=7
			double d_theta=360.0/static_cast<double>(nr);

			for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
			{
					temp.r[A_X]=R*sin(theta*PI/180.0);
					temp.r[A_Y]=R*cos(theta*PI/180.0);
					temp.r[A_Z]=z;
					temp.boundary=0;
					NODEa.push_back(temp);
					temp.id+=1;
			}	
	}

	//下面
	

	for(double r=dx;r<R+0.1*dx;r+=dx)//	for(double r=R;r>le-0.1*le;r-=le)
	{
		int nr=static_cast<int>(2.0*PI*r/dx);
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
		{
				temp.r[A_X]=r*sin(theta*PI/180.0);
				temp.r[A_Y]=r*cos(theta*PI/180.0);
				temp.r[A_Z]=lz;
				temp.boundary=0;
				NODEa.push_back(temp);
				temp.id+=1;
		}
	}//*/

	//.stlファイル読み取り
//	in.load_stl("COIL");
	//nodeファイル作成
	MakeNodeFile(CON, NODEa, "air2.node");
	//.nodeファイル読み取り
	in.load_node("air2");

	//まずは流体節点のみで分割
	tetrahedralize("", &in, &out);
		// i 要素内へ節点を追加(今の所無くても機能している)
		// f .faceファイルに境界ではない面も含める
		// e .edgeファイルの出力(ONにするとなぜか止まってしまう)
		// n .neighファイルの出力
	//出力
	out.save_nodes("air2_whole");
	out.save_elements("air2_whole");
	out.save_faces("air2_whole");

	GetFacetList(FACEa, in, out, COIL);
}


//平板電極境界面作成
void tetgen_function::SetPlateElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEp, vector<tetgen_facet> &FACEp)
{
	cout<<"平板電極境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE2;
	temp.boundary=ELECTRODE2;

	//分割数決定
	double	h = TET.height_plate;						//平板高さ
	double	dh = TET.fine_plate_t;						//厚み方向メッシュ粗さ
	int		nh = int((TET.thickness_plate+0.1*dh)/dh);	//厚み方向分割数
	double	dL = TET.fine_plate_L;						//xy方向メッシュ粗さ
	int		nL = int((TET.length_plate/2+0.1*dL)/dL);	//xy方向分割数(片側)

	//接点データ作成
	for(int z=0;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				if(x==-nL || x==nL || y==-nL || y==nL || z==0 || z==nh)	//平板電極領域の端に来たとき節点を置く
				{
					temp.r[A_X]=x*dL;
					temp.r[A_Y]=y*dL;
					temp.r[A_Z]=z*dh+h;
					//temp.attribute=ELECTRODE;
					//temp.boundary=ELECTRODE;
					NODEp.push_back(temp);
					temp.id+=1;
				}
			}
		}
	}//*/

	//nodeファイル作成
	MakeNodeFile(CON, NODEp, "NODEp.node");

	in.load_node("NODEp");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_plate");
	out.save_elements("boundary_plate");
	out.save_faces("boundary_plate");
	//*/

	//境界面データ取得
	GetFacetList(FACEp, in, out, ELECTRODE2);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}


//円柱電極境界面作成
void tetgen_function::SetColumnElectrodeBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEc, vector<tetgen_facet> &FACEc)
{
	cout<<"円柱電極境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	//分割数決定
	double le=CON.get_distancebp();
	double rc=TET.radius_column;	//円柱半径
	double L=TET.length_column;		//円柱長さ
	double dL=TET.fine_column_L;	//円柱長さ方向メッシュ粗さ
	//double z0=0;					//上面の位置

	
	//上面
	//中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=0;
	NODEc.push_back(temp);
	temp.id+=1;

	for(double r=le;r<rc+0.1*le;r+=le)
	{
		int nr=int(2.0*PI*r/le);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			//double theta=nt*d_theta;
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=0;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}


	//側面
	int flag=0;
	double dz=le;
	double z=0;
	rc*=1.005;	//わずかに太くする(メッシュが繋がるのを防ぐため)
	int nr=int(2.0*PI*rc/(2*le));//leが2倍されていることに注意
	double d_theta=360.0/(double)nr;
	
	while(1)
	{
		//z方向への移動
		if(flag==0)			dz*=1.05;
		else if(flag==1)	dz=dL;
		z-=dz;

		if(dz>dL)	flag=1;
		if(z<-L+le)
		{
			z=-L+le;
			//flag=2;
			break;
		}

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=rc*sin(theta*PI/180.0);
			temp.r[A_Y]=rc*cos(theta*PI/180.0);
			temp.r[A_Z]=z;
			NODEc.push_back(temp);
			temp.id+=1;
		}

		//if(flag==2)	break;
	}

	//下面
	//中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=-L+le;
	NODEc.push_back(temp);
	temp.id+=1;

	//for(double r=le;r<rc+0.1*le;r+=le)
	for(double r=rc;r>le-0.1*le;r-=2*le)//leが2倍されていることに注意
	{
		int nr=int(2.0*PI*r/(2*le));//leが2倍されていることに注意
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L+le;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}//*/


	//nodeファイル作成
	MakeNodeFile(CON, NODEc, "NODEc.node");

	in.load_node("NODEc");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_column");
	out.save_elements("boundary_column");
	out.save_faces("boundary_column");
	//*/

	//境界面データ取得
	GetFacetList(FACEc, in, out, ELECTRODE1);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEc, "FACEc.face");
}


//土台境界面作成
void tetgen_function::SetBaseBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEb, vector<tetgen_facet> &FACEb)
{
	cout<<"土台境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=ELECTRODE1;
	temp.boundary=ELECTRODE1;

	//分割数決定
	double	le = CON.get_distancebp();
	double	rc = TET.radius_column;
	double	h = TET.length_column;						//平板高さ
	double	dh = TET.fine_base;						//厚み方向メッシュ粗さ
	int		nh = int((TET.thickness_base+0.1*dh)/dh);	//厚み方向分割数
	double	L = TET.length_base;
	double	dL = TET.fine_base;						//xy方向メッシュ粗さ
	int		nL = int((TET.length_base/2+0.1*dL)/dL);	//xy方向分割数(片側)

	//上面接続部分
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=-L;
	NODEb.push_back(temp);
	temp.id+=1;

	rc*=1.005;//円柱に合わせてわずかに太くする
	for(double r=rc;r>le-0.1*le;r-=2*le)//leが2倍されていることに注意
	{
		int nr=int(2.0*PI*r/(2*le));//leが2倍されていることに注意
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L;
			NODEb.push_back(temp);
			temp.id+=1;
		}
	}//*/

	//上面表面
	double r=rc+2*le;
	double s=2*le;
	while(r<sqrt(2.0)*L/2+dL)
	{
		int nr=int(2.0*PI*r/s);
		double d_theta=360.0/(double)nr;

		for(double theta=0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=-L;
			if(fabs(temp.r[A_X])<L/2 && fabs(temp.r[A_Y])<L/2)
			{
				NODEb.push_back(temp);
				temp.id+=1;
			}
		}
		r*=1.15;
		s*=1.15;
	}//*/

	//その他の面
	for(int z=0;z<=nh;z++)
	{
		for(int y=-nL;y<=nL;y++)
		{
			for(int x=-nL;x<=nL;x++)
			{
				if(x==-nL || x==nL || y==-nL || y==nL || z==nh)	//平板電極領域の端に来たとき節点を置く
				{
					//if(sqrt(temp.r[A_X]*temp.r[A_X]+temp.r[A_Y]*temp.r[A_Y])>rc+le || z>0)
					{
						temp.r[A_X]=x*dL;
						temp.r[A_Y]=y*dL;
						temp.r[A_Z]=-h-z*dL;
						NODEb.push_back(temp);
						temp.id+=1;
					}
				}
			}
		}
	}//*/

	//nodeファイル作成
	MakeNodeFile(CON, NODEb, "NODEb.node");

	in.load_node("NODEb");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_base");
	out.save_elements("boundary_base");
	out.save_faces("boundary_base");
	//*/

	//境界面データ取得
	GetFacetList(FACEb, in, out, ELECTRODE1);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEp, "FACEp.face");
}

void tetgen_function::SetMagnetBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEc, vector<tetgen_facet> &FACEc)
{
	cout<<"磁石境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	tetgen_node temp;
	temp.id=0;
	temp.attribute=MAGNET;
	temp.boundary=MAGNET;

	//分割数決定
	double le=CON.get_distancebp();
	double R=CON.get_magnet_r();	//磁石半径
	double L=CON.get_magnet_H();	//磁石高さ
	double Zmin=CON.get_magnet_Z()-0.5*CON.get_magnet_H();
	double Zmax=CON.get_magnet_Z()+0.5*CON.get_magnet_H();
	double dL=(Zmax-Zmin)/10;	//磁石長さ方向メッシュ粗さ
	double dR=le;//半径方向の分割単位

	//上下側面でテーパーを付けたい時は分けるべきだがまとめてもいいかも・・・

	//上面中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=Zmax;
	NODEc.push_back(temp);
	temp.id+=1;

	//上面
	for(double r=dR;r<R+0.1*dR;r+=dR)
	{
		int nr=static_cast<int>(2.0*PI*r/dR);//正何角形で近似するか（角度方向の分割数）
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=Zmax;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}

	//側面	//R*=1.005;	//わずかに太くする(メッシュが繋がるのを防ぐため)・・・これはいらん
	for(double z=Zmax-dL;z>Zmin;z-=dL)
	{
		//中心（忘れない）
		temp.r[A_X]=0;
		temp.r[A_Y]=0;	
		temp.r[A_Z]=z;
		NODEc.push_back(temp);
		temp.id+=1;

		for(double r=dR;r<R+0.1*dR;r+=dR)
		{
			int nr=static_cast<int>(2.0*PI*r/dR);//r=dRの時はnr=7
			double d_theta=360.0/static_cast<double>(nr);

			for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
			{
				temp.r[A_X]=r*sin(theta*PI/180.0);
				temp.r[A_Y]=r*cos(theta*PI/180.0);
				temp.r[A_Z]=z;
				NODEc.push_back(temp);
				temp.id+=1;
			}	
		}
	}

	//下面
	//中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=Zmin;
	NODEc.push_back(temp);
	temp.id+=1;

	for(double r=dR;r<R+0.1*dR;r+=dR)//	for(double r=R;r>le-0.1*le;r-=le)
	{
		int nr=static_cast<int>(2.0*PI*r/dR);
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
		{
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=Zmin;
			NODEc.push_back(temp);
			temp.id+=1;
		}
	}//*/

	//nodeファイル作成
	MakeNodeFile(CON, NODEc, "NODEc.node");

	in.load_node("NODEc");
	tetrahedralize("", &in, &out);
	out.save_nodes("boundary_column");
	out.save_elements("boundary_column");
	out.save_faces("boundary_column");
	//*/

	//境界面データ取得
	GetFacetList(FACEc, in, out, MAGNET);
	//.faceファイル作成
	//MakeFaceFile(CON, FACEc, "FACEc.face");
}

void tetgen_function::SetCOILBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEc, vector<tetgen_facet> &FACEc)
{
	cout<<"コイル境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

//	vector<tetgen_element> ELEMc;
	tetgen_node temp;
	temp.id=0;

	//分割数決定
	double le=CON.get_distancebp();
	double R=CON.get_magnet_r();	//磁石半径
	double L=CON.get_magnet_H();	//磁石高さ
	double Zmin=CON.get_magnet_Z()-0.5*CON.get_magnet_H();
	double Zmax=CON.get_magnet_Z()+0.5*CON.get_magnet_H();
	double dL=(Zmax-Zmin)/20;	//磁石長さ方向メッシュ粗さ
	double dR=0.001;//半径方向の分割単位

	//上下側面でテーパーを付けたい時は分けるべきだがまとめてもいいかも・・・

	//中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=Zmax;
	temp.attribute=IRON;
	temp.boundary=0;
	NODEc.push_back(temp);
	temp.id+=1;

	//上面
	for(double r=dR;r<R+0.1*dR;r+=dR)
	{
		int nr=static_cast<int>(2.0*PI*r/dR);//正何角形で近似するか（角度方向の分割数）
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360-0.1*d_theta;theta+=d_theta)
		{
			if(r<=dR*4){
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=Zmax;
			temp.attribute=IRON; //IRON
			temp.boundary=0;
			NODEc.push_back(temp);
			temp.id+=1;
			}
			else {
			temp.r[A_X]=r*sin(theta*PI/180.0);
			temp.r[A_Y]=r*cos(theta*PI/180.0);	
			temp.r[A_Z]=Zmax;
			temp.attribute=COIL; //COIL
			if(r>=R-0.1*dR){	//(r==R)だと誤差で通らない
				int a=int((Zmax-Zmin)/dL);		//コイルの表面にディレクレ形境界条件を挿入
				if(a%3==2) temp.boundary=23;
				else if(a%3==1) temp.boundary=22;
				else if(a%3==0) temp.boundary=21;
			}
			else temp.boundary=0;
			NODEc.push_back(temp);
			temp.id+=1;
			}
		}
	}

	//側面	//R*=1.005;	//わずかに太くする(メッシュが繋がるのを防ぐため)・・・これはいらん
	for(double z=Zmax-dL;z>Zmin;z-=dL)
	{
		//中心（忘れない）
		temp.r[A_X]=0;
		temp.r[A_Y]=0;	
		temp.r[A_Z]=z;
		temp.attribute=IRON; //IRON
		temp.boundary=0;
		NODEc.push_back(temp);
		temp.id+=1;

		for(double r=dR;r<R+0.1*dR;r+=dR)
		{
			int nr=static_cast<int>(2.0*PI*r/dR);//r=dRの時はnr=7
			double d_theta=360.0/static_cast<double>(nr);

			for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
			{
				if(r<=dR*4){
					temp.r[A_X]=r*sin(theta*PI/180.0);
					temp.r[A_Y]=r*cos(theta*PI/180.0);
					temp.r[A_Z]=z;
					temp.attribute=IRON; //IRON
					temp.boundary=0;
					NODEc.push_back(temp);
					temp.id+=1;
				}
				else {
					temp.r[A_X]=r*sin(theta*PI/180.0);
					temp.r[A_Y]=r*cos(theta*PI/180.0);	
					temp.r[A_Z]=z;
					temp.attribute=COIL; //COIL
					if(r>=R-0.1*dR){
						int a=int((z-Zmin)/dL);
						if(a%3==2) {temp.boundary=23;}
						else if(a%3==1) {temp.boundary=22;}
						else if(a%3==0) {temp.boundary=21;}
					}
					else {temp.boundary=0;}
					NODEc.push_back(temp);
					temp.id+=1;
				}
			}	
		}
	}

	//下面
	//中心
	temp.r[A_X]=0;
	temp.r[A_Y]=0;	
	temp.r[A_Z]=Zmin;
	temp.attribute=IRON; //IRON
	temp.boundary=0;
	NODEc.push_back(temp);
	temp.id+=1;

	for(double r=dR;r<R+0.1*dR;r+=dR)//	for(double r=R;r>le-0.1*le;r-=le)
	{
		int nr=static_cast<int>(2.0*PI*r/dR);
		double d_theta=360.0/static_cast<double>(nr);

		for(double theta=0.0;theta<360.0-0.1*d_theta;theta+=d_theta)
		{
			if(r<=dR*4){
				temp.r[A_X]=r*sin(theta*PI/180.0);
				temp.r[A_Y]=r*cos(theta*PI/180.0);
				temp.r[A_Z]=Zmin;
				temp.attribute=IRON; //IRON
				temp.boundary=0;
				NODEc.push_back(temp);
				temp.id+=1;
				}
			else {
				temp.r[A_X]=r*sin(theta*PI/180.0);
				temp.r[A_Y]=r*cos(theta*PI/180.0);	
				temp.r[A_Z]=Zmin;
				temp.attribute=COIL; //COIL
				if(r>=R-0.1*dR) {temp.boundary=21;}
				else {temp.boundary=0;}
				NODEc.push_back(temp);
				temp.id+=1;
			}
		}
	}//*/

	//.stlファイル読み取り
//	in.load_stl("COIL");
	//nodeファイル作成
	MakeNodeFile(CON, NODEc, "COIL.node");
	//.nodeファイル読み取り
	in.load_node("COIL");

	//まずは流体節点のみで分割
	tetrahedralize("", &in, &out);
		// i 要素内へ節点を追加(今の所無くても機能している)
		// f .faceファイルに境界ではない面も含める
		// e .edgeファイルの出力(ONにするとなぜか止まってしまう)
		// n .neighファイルの出力
	//出力
	out.save_nodes("COIL_whole");
	out.save_elements("COIL_whole");
	out.save_faces("COIL_whole");
	//////////////////ここまでで流体節点のみを使って、すべての要素が繋がった凸なメッシュができた(fluid_wholeで確認可能)*/
	//要素材質の決定
//	Geteleattribute(NODEc,ELEMc, out);
	///////////////不要な要素の削除

	//.nodeの取得
//	GetPointList(NODEc, in, out);
	//.eleの取得
//	GetTetrahedronList(ELEMc, in, out);
	//長い要素の除去
//	DelThinTetrahedron(CON, TET, NODEc, ELEMc, in, out);

	//節点-要素関係
//	SetRelation_NodeElem(CON, NODEc, ELEMc);
	//要素-要素関係
//	SetRelation_ElemElem(CON, NODEc, ELEMc);
	//境界面データ取得
	GetFacetList(FACEc, in, out, COIL);

/*	/////////////////要素確認用ファイル///////////////////////////////////
	out.save_nodes("boundary_COIL");	//fluid.2.nodeと同じファイル
	MakeElemFile(CON, ELEMc, "boundary_COIL.ele");
	MakeFaceFile(CON, FACEc, "boundary_COIL.face");
	////////////////ここまででエラストマーのメッシュが切れた//////////////////////*/

	//境界面データ取得
//	GetFacetList(FACEc, in, out, MAGNET);
}

void tetgen_function::SetIRONBoundary(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODEi, vector<tetgen_facet> &FACEi)
{
	cout<<"鉄心境界作成"<<endl;

	tetgenio in, out;
	in.initialize();
	out.initialize();

	vector<tetgen_element> ELEMi;
	tetgen_node temp;
	temp.id=0;

	//分割数決定
	int divN[3];
	divN[A_R]=20;//半径方向分割数
	divN[A_t]=20;//角度方向分割数（正何角形で近似するか）
	divN[A_Z]=20;//高さ方向分割数
	double Rmin=0.0;
	double Rmax=CON.get_magnet_r();

	double Zmin=CON.get_magnet_Z()-0.5*CON.get_magnet_H();
	double Zmax=CON.get_magnet_Z()+0.5*CON.get_magnet_H();
	double divL[3];
	divL[A_R]=(Rmax-Rmin)/divN[A_R];
	divL[A_t]=(2*PI)/divN[A_t];
	divL[A_Z]=(Zmax-Zmin)/divN[A_Z];
	double theta2=CON.get_magnet_angle()*PI*2/360;//磁石の回転角度

	//上下側面でテーパーを付けたい時は分けるべきだがまとめてもいいかも・・・

	for(int n=0;n<5;n++)
	{
		if(n==0)//中心
		{
			for(int k=0;k<=divN[A_Z];k++)
			{
				temp.r[A_X]=0;
				temp.r[A_Y]=0;
				temp.r[A_Z]=Zmin+divL[A_Z]*k;
				temp.attribute=IRON;				
				temp.boundary=0;			//境界条件
			
				if(fabs(theta2)>0.0175)//doubleを0と比較してはいけない（一致するわけがない）
				{
					double XX=temp.r[A_X];
					double ZZ=temp.r[A_Z]-CON.get_magnet_Z();
					double newX=XX*cos(theta2)-ZZ*sin(theta2);
					double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
					temp.r[A_X]=newX;
					temp.r[A_Z]=newZ;
				}
				NODEi.push_back(temp);
				temp.id+=1;	
			}
		}
		else 
		{
			for(int m=0;m<divN[A_t];m++)
			{
				for(int k=0;k<=divN[A_Z];k++)
				{
					double r=divL[A_R]*n;
					double theta=divL[A_t]*m;
					temp.r[A_X]=r*cos(theta);
					temp.r[A_Y]=r*sin(theta);
					temp.r[A_Z]=Zmin+divL[A_Z]*k;
					
					temp.attribute=IRON;	
					temp.boundary=0;			//境界条件

		//			temp.attribute=IRON;
		//			 if((n==divN[A_R]) && (k%3==0))temp.boundary=23;//外周側面
		//			else if((n==divN[A_R]) && (k%3==1))temp.boundary=22;
		//			else if((n==divN[A_R]) && (k%3==2))temp.boundary=21;
		/*			else if((n==6) && (k%3==0)) NODE[num].boundary_condition=21;	//鉄心との境界面
					else if((n==6) && (k%3==1))NODE[num].boundary_condition=22;
					else if((n==6) && (k%3==2))NODE[num].boundary_condition=23;  */
					/*temp.boundary=0;*/			//境界条件
	
					if(fabs(theta2)>0.0175)//磁石が傾いている場合(1degより大きい場合)・・・doubleを0と比較してはいけない
					//if(theta2!=0)これはダメ！！
					{
						double XX=temp.r[A_X];
						double ZZ=temp.r[A_Z]-CON.get_magnet_Z();
						double newX=XX*cos(theta2)-ZZ*sin(theta2);
						double newZ=XX*sin(theta2)+ZZ*cos(theta2)+CON.get_magnet_Z();
						temp.r[A_X]=newX;
						temp.r[A_Z]=newZ;
					}
					NODEi.push_back(temp);
					temp.id+=1;					//ID格納
				}
			}
		}
	}
	//nodeファイル作成
	MakeNodeFile(CON, NODEi, "IRON.node");

	//.nodeファイル読み取り
	in.load_node("IRON");

	//まずは流体節点のみで分割
	tetrahedralize("", &in, &out);
		// i 要素内へ節点を追加(今の所無くても機能している)
		// f .faceファイルに境界ではない面も含める
		// e .edgeファイルの出力(ONにするとなぜか止まってしまう)
		// n .neighファイルの出力

	//出力
	out.save_nodes("IRON_whole");
	out.save_elements("IRON_whole");

	//////////////////ここまでで流体節点のみを使って、すべての要素が繋がった凸なメッシュができた(fluid_wholeで確認可能)*/
	//要素材質の決定
//	Geteleattribute(NODEc,ELEMc, out);
	///////////////不要な要素の削除

	//.nodeの取得
	GetPointList(NODEi, in, out);
	//.eleの取得
	GetTetrahedronList(ELEMi, in, out);

	//長い要素の除去
//	DelThinTetrahedron(CON, TET, NODEc, ELEMc, in, out);

	//節点-要素関係
	SetRelation_NodeElem(CON, NODEi, ELEMi);
	//要素-要素関係
	SetRelation_ElemElem(CON, NODEi, ELEMi);
	//境界面データ取得
	GetFacetList(FACEi, in, out, IRON);

	/////////////////要素確認用ファイル///////////////////////////////////
	out.save_nodes("boundary_IRON");	//fluid.2.nodeと同じファイル
	MakeElemFile(CON, ELEMi, "boundary_IRON.ele");
	MakeFaceFile(CON, FACEi, "boundary_IRON.face");
	////////////////ここまででエラストマーのメッシュが切れた//////////////////////*/

	//境界面データ取得
//	GetFacetList(FACEc, in, out, MAGNET);
}

//境界節点・境界面データの結合
void tetgen_function::UniteBoundaryData(mpsconfig &CON, 
					   vector<tetgen_node> &NODE, vector<tetgen_node> &NODEa1, vector<tetgen_node> &NODEa2, vector<tetgen_node> &NODEp, vector<tetgen_node> &NODEc, vector<tetgen_node> &NODEb, vector<tetgen_node> &NODEw, 
					   vector<tetgen_facet> &FACE, vector<tetgen_facet> &FACEa1, vector<tetgen_facet> &FACEa2, vector<tetgen_facet> &FACEp, vector<tetgen_facet> &FACEc, vector<tetgen_facet> &FACEb, vector<tetgen_facet> &FACEw, 
					   vector<int> &TRANS)
{
	cout<<"境界節点・境界面データの結合"<<endl;

	NODE.clear();
	FACE.clear();
	tetgen_node temp_n;
	tetgen_facet temp_f;

	temp_n.boundary=0;
	temp_f.boundary=0;

	int offset_n=0;	//節点番号のオフセット量
	int offset_f=0;	//表面番号のオフセット量


	//水滴境界
	for(int i=0;i<(int)NODEw.size();i++)//節点
	{
		TRANS.push_back(NODEw[i].part_no);	//TRANSに粒子番号を格納(boundaryに粒子番号を格納している)

		temp_n=NODEw[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=WATER;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEw.size();i++)//表面
	{
		temp_f=FACEw[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=WATER;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEw.size();
	offset_f+=(int)FACEw.size();

	//空気境界
	for(int i=0;i<(int)NODEa1.size();i++)//節点
	{
		temp_n=NODEa1[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=AIR;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEa1.size();i++)//表面
	{
		temp_f=FACEa1[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=AIR;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEa1.size();
	offset_f+=(int)FACEa1.size();

	//空気境界 高解像度面
	for(int i=0;i<(int)NODEa2.size();i++)//節点
	{
		temp_n=NODEa2[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=AIR;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEa2.size();i++)//表面
	{
		temp_f=FACEa2[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=AIR;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEa2.size();
	offset_f+=(int)FACEa2.size();

	//平板電極境界
	for(int i=0;i<(int)NODEp.size();i++)//節点
	{
		temp_n=NODEp[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE2;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEp.size();i++)//表面
	{
		temp_f=FACEp[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE2;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEp.size();
	offset_f+=(int)FACEp.size();

	//円柱電極境界
	for(int i=0;i<(int)NODEc.size();i++)//節点
	{
		temp_n=NODEc[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE1;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEc.size();i++)//表面
	{
		temp_f=FACEc[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE1;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	offset_n+=(int)NODEc.size();
	offset_f+=(int)FACEc.size();

	//土台境界
	for(int i=0;i<(int)NODEb.size();i++)//節点
	{
		temp_n=NODEb[i];
		temp_n.id+=offset_n;
		//temp_n.boundary=ELECTRODE1;
		NODE.push_back(temp_n);
	}
	for(int i=0;i<(int)FACEb.size();i++)//表面
	{
		temp_f=FACEb[i];
		for(int n=0;n<3;n++)	temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=ELECTRODE1;
		FACE.push_back(temp_f);
	}

	//オフセット量更新
	//offset_n+=(int)NODEb.size();
	//offset_f+=(int)FACEb.size();

	//cout<<"----------OK"<<endl;
}

//境界節点・境界面データの追加
void tetgen_function::AddBoundaryData(mpsconfig &CON, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_node> &NODE, vector<tetgen_facet> &FACE, int attribute)
{
	//NODE,FACEに格納されている各部品のデータを，NODEall,FACEallに格納していく．

	tetgen_node temp_n;
	tetgen_facet temp_f;

	temp_n.boundary=0;
	temp_f.boundary=0;

	int offset_n=(int)NODEall.size();	//節点番号のオフセット量
	int offset_f=(int)FACEall.size();	//表面番号のオフセット量


	//節点の追加
	for(int i=0;i<(int)NODE.size();i++)//節点
	{
		temp_n=NODE[i];
		temp_n.id+=offset_n;
		temp_n.boundary=NODE[i].boundary;//なぜかコメントアウトされていた・・・
		NODEall.push_back(temp_n);
	}

	//面の追加
	for(int i=0;i<(int)FACE.size();i++)//表面
	{
		temp_f=FACE[i];
		for(int n=0;n<3;n++) temp_f.node[n]+=offset_n;
		temp_f.id+=offset_f;
		temp_f.boundary=FACE[i].boundary;
		FACEall.push_back(temp_f);
	}
}

//TRANS[]の格納
void tetgen_function::SetTRANS(vector<tetgen_node> &NODE, vector<int> &TRANS)
{
	//TRANS[i]には、節点番号iに対応する粒子番号を格納する。
	//FEM3D.cpp では、節点番号が1から始まるので、TRANS[0]には宣言後に-1を入れる。（ここでは既に入っている）

	for(int i=0;i<(int)NODE.size();i++)
	{
		TRANS.push_back(NODE[i].part_no);	//TRANSに粒子番号を格納
	}
}


//材質の決定（まだ未完成）
void tetgen_function::DecisionAttribute(vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	double M[4];	//節点の材質を格納
	
	for(int i=0;i<(int)ELEM.size();i++)
	{	
		if(ELEM[i].attribute==AIR)
		{
			for(int n=0;n<4;n++)	M[n]=NODE[ELEM[i].node[n]].attribute;

			if(M[0]==AIR || M[1]==AIR || M[2]==AIR || M[3]==AIR)	//1つでも空気節点があれば空気要素
			{
				ELEM[i].attribute=AIR;
				out.tetrahedronattributelist[i]=AIR;
			}
			else
			{
				ELEM[i].attribute=WATER;	//残りは水
				out.tetrahedronattributelist[i]=WATER;
			}
		}
	}
}


//材質の修正  この関数を使うときは直前にtetgenioからデータを取得しておくこと
void tetgen_function::ModifyAttribute(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	//tetgenioからデータを取得済みと仮定する
	//節点-節点関係が得られていると仮定する

	/*double le=CON.get_distancebp();
	double rc=TET.radius_column;
	double L=TET.length_column;*/

	//要素の材質の修正  この時点では水要素となる領域は材質番号が0(未定義)となっている
	for(int i=0;i<(int)ELEM.size();i++)
	{
		//未定義(0)の材質を水にする
		if(ELEM[i].attribute==0)
		{
			ELEM[i].attribute=WATER;
		}
	}

	//節点の材質を要素の材質と合わせる

	//まず全部空気にする
	for(int i=0;i<(int)ELEM.size();i++)
	{
		NODE[i].attribute=AIR;
	}
	//水節点の決定
	for(int i=0;i<(int)ELEM.size();i++)
	{
		if(ELEM[i].attribute==WATER)
		{
			for(int n=0;n<4;n++)	NODE[ELEM[i].node[n]].attribute=ELEM[i].attribute;
		}
	}
	//電極節点の決定
	for(int i=0;i<(int)ELEM.size();i++)
	{
		if(ELEM[i].attribute==WATER)
		{
			for(int n=0;n<4;n++)	NODE[ELEM[i].node[n]].attribute=ELEM[i].attribute;
		}
	}

	//節点-節点関係より節点の修正
	for(int i=0;i<(int)NODE.size();i++)
	{
		if(NODE[i].attribute==AIR)
		{
			int num_air=0;
			int num_ele=0;
		}
	}
}


//材質の修正  tetgenioを直接編集
void tetgen_function::ModifyAttribute_tetgenio(mpsconfig &CON, tetgen_config &TET, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	double le=CON.get_distancebp();
	double rc=TET.radius_column;
	double L=TET.length_column;

	//節点自動挿入された節点の材質をとりあえず空気とする
	for(int i=0;i<out.numberofpoints;i++)
	{
		if(out.pointattributelist[i]==0)
		{
			out.pointattributelist[i]=AIR;
		}
	}//*/

	//要素の材質の修正  この時点では水要素となる領域は材質番号がデフォルトで0(未定義)となっている
	if(CON.get_model_number()==14)
	{
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			//未定義の材質を水にする
			if(out.tetrahedronattributelist[i]==0)
			{
				out.tetrahedronattributelist[i]=WATER;//0だったものがWATERになる
			}
		}
	}

	if(CON.get_model_number()==2)
	{
			for(int i=0;i<out.numberofpoints;i++){
	if(out.pointattributelist[i]==FACE_P) out.pointattributelist[i]=MAGELAST;//要素を戻す
	}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			/////////////////要素を作る点の要素/////////////////////////
			int ai=(int)out.pointattributelist[out.tetrahedronlist[4*i]];
			int bi=(int)out.pointattributelist[out.tetrahedronlist[4*i+1]];
			int ci=(int)out.pointattributelist[out.tetrahedronlist[4*i+2]];
			int di=(int)out.pointattributelist[out.tetrahedronlist[4*i+3]];
			////////////////////////////////////////////////////////////
			
			if(ai==bi && bi==ci && ci==di) out.tetrahedronattributelist[i]=di;//すべて同じ素材なら要素もその素材
			
	//		if(ai==AIR || bi==AIR || ci==AIR || di==AIR) out.tetrahedronattributelist[i]=AIR;//ひとつでも空気接点を含んでいるなら空気	//これを入れるとMRE中に空気要素ができる
			
		}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			if(out.tetrahedronattributelist[i]==0) out.tetrahedronattributelist[i]=MAGELAST;
		}
	}
	else
	{
		for(int i=0;i<out.numberofpoints;i++){
	if(out.pointattributelist[i]==FACE_P) out.pointattributelist[i]=MAGELAST;//要素を戻す
	}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			/////////////////要素を作る点の要素/////////////////////////
			int ai=(int)out.pointattributelist[out.tetrahedronlist[4*i]];
			int bi=(int)out.pointattributelist[out.tetrahedronlist[4*i+1]];
			int ci=(int)out.pointattributelist[out.tetrahedronlist[4*i+2]];
			int di=(int)out.pointattributelist[out.tetrahedronlist[4*i+3]];
			////////////////////////////////////////////////////////////
			
			if(ai==bi && bi==ci && ci==di) out.tetrahedronattributelist[i]=di;//すべて同じ素材なら要素もその素材
			
			else if((ai==IRON || bi==IRON || ci==IRON || di==IRON) && !(ai==AIR || bi==AIR || ci==AIR || di==AIR)) {
				out.tetrahedronattributelist[i]=IRON;//ひとつでも鉄接点を含んでいるなら,かつ空気をひとつも含んでないなら鉄
			}
			//if(ai==AIR || bi==AIR || ci==AIR || di==AIR) out.tetrahedronattributelist[i]=AIR;//ひとつでも空気接点を含んでいるなら空気	//これを入れるとMRE中に空気要素ができる
			
		}
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			if(out.tetrahedronattributelist[i]==0) out.tetrahedronattributelist[i]=AIR;
		}
/*		//色の調整　小さい順に赤・緑・青・黄色 これすると解析が回らない　メッシュ材料にこの値が使用されるため
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			if(out.tetrahedronattributelist[i]==COIL) out.tetrahedronattributelist[i]=1;
			else if(out.tetrahedronattributelist[i]==MAGELAST) out.tetrahedronattributelist[i]=2;
			else if(out.tetrahedronattributelist[i]==AIR) out.tetrahedronattributelist[i]=3;
			else if(out.tetrahedronattributelist[i]==IRON) out.tetrahedronattributelist[i]=4;
			
		}*/
	}

	//自動追加された節点ののattributeとboundary_markerの修正
	//この時点では流体内部や電極内部や電極内部の節点は空気節点になっているのでそれぞれの材質に修正する
	//※電極の一番上の節点は電極の節点とするため、先に水の節点決めてから電極節点を決める
	
	if(CON.get_model_number()==2)
	{
		for(int i=0;i<out.numberoftetrahedra;i++)
		{
			if(out.tetrahedronattributelist[i]==MAGNET)
			{
				for(int n=0;n<4;n++)
				{
					out.pointattributelist[out.tetrahedronlist[i*4+n]]=MAGNET;
					out.pointmarkerlist[out.tetrahedronlist[i*4+n]]=MAGNET;
				}
			//	//材質はELECTRODEに戻しておく
			//	out.tetrahedronattributelist[i]=ELECTRODE;
			}
		}
		
	}
}


//要素の細分化
void tetgen_function::FineElement(mpsconfig &CON, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM, tetgenio &in, tetgenio &out)
{
	tetgenio add;
	add.initialize();

	vector<tetgen_node> NODEadd;
	tetgen_node temp;
	temp.id=0;
	temp.attribute=0;
	temp.boundary=0;

	for(int i=0;i<out.numberoftetrahedra;i++)
	{
		if(out.tetrahedronattributelist[i]==AIR)
		{
			double r[3]={0,0,0};
			//for(int d=0;d<3;d++)	r[d]=out.vpointlist[i*3+d];
			for(int d=0;d<3;d++)
			{
				for(int n=0;n<4;n++)	r[d]+=out.pointlist[3*(out.tetrahedronlist[4*i+n])+d];
				r[d]/=4;
			}

			if(fabs(r[A_X])<0.0005 && fabs(r[A_Y])<0.0005 && fabs(r[A_Z])<0.0005)
			{
				temp.r[A_X]=r[A_X];
				temp.r[A_Y]=r[A_Y];	
				temp.r[A_Z]=r[A_Z];
				NODEadd.push_back(temp);
				temp.id+=1;
			}
		}
	}

	MakeNodeFile(CON, NODEadd, "output-a.node");
}//*/

//節点間距離計算関数
double tetgen_function::Distance(tetgen_node &point1, tetgen_node &point2)
{
	double dis=0;

	for(int d=0;d<3;d++)	dis+=(point2.r[d]-point1.r[d])*(point2.r[d]-point1.r[d]);

	return sqrt(dis);
}


//要素重心座標計算関数
void tetgen_function::CalcBarycentricElement(mpsconfig&, vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM)
{
	double r[3]={0,0,0};

	for(int i=0;i<(int)ELEM.size();i++)
	{
		for(int d=0;d<3;d++)	for(int n=0;n<4;n++)	r[d]+=NODE[ELEM[i].node[n]].r[d];
		for(int d=0;d<3;d++)	ELEM[i].g[d]=r[d]/4;
	}
}
