///////////////////////注意///////////////////////
//TetViewでファイルを開いているとプロセスが終了せず
//outputファイルの読み書きが出来ないのでプログラムが進まない！
//必ずTetViewのプロセスを終了してから.exeを実行すること！！

//TetGen呼び出し関数  TetGenの入り口

#include "stdafx.h"


//TetGen入口  モデルにより分岐
void tetgen_function::call_TetGen(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	
	cout<<"TetGenによるメッシュ生成開始"<<endl;
	clock_t t1=clock();	//クロック数取得


	//////モデル毎に分岐/////////////////////////////
	TetGen_elastic(CON,PART,fluid_number,particle_number,NODEall,FACEall,ELEMall,TRANS); //磁性エラストマー用
	/////////////////////////////////////////////////


	clock_t t2=clock();	//クロック数取得　
	cout<<"メッシュ生成完了  CPU time="<<(t2-t1)/CLOCKS_PER_SEC<<endl;
}


//////これより下に、作りたいモデルのメッシュ生成プログラムを記述してください/////////////////////////////////////////////

//弾性体用　メッシュ生成
void tetgen_function::TetGen_elastic(mpsconfig &CON, vector<mpselastic> &PART, int fluid_number, int particle_number, vector<tetgen_node> &NODEall, vector<tetgen_facet> &FACEall, vector<tetgen_element> &ELEMall, vector<int> &TRANS)
{
	//tetgenioクラス宣言と初期化
	tetgenio in, out;
	in.initialize();
	out.initialize();

	//TetGen設定クラス
	tetgen_config TET;

	//全体配列初期化
	NODEall.clear();
	FACEall.clear();
	//各部品用配列（使い回し）
	vector<tetgen_node> NODE;
	vector<tetgen_facet> FACE;


	////////////エラストマー領域////////////
	NODE.clear();
	FACE.clear();
	SetElastBoundary(CON, PART, TET, NODE, FACE);
	SetTRANS(NODE, TRANS);	//粒子が節点に対応する時に実行？？？
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, MAGELAST);//*


	//////////////電磁石領域/////////////
	NODE.clear();
	FACE.clear();
	SetMagnetBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, MAGNET);//*/
	/////////////コイル領域//////////////
/*	NODE.clear();
	FACE.clear();
	SetCOILBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, COIL);
	////////////鉄心領域////////////////*/
/*	NODE.clear();
	FACE.clear();
	SetIRONBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, IRON);
	////////////空気領域////////////*/
	NODE.clear();
	FACE.clear();
	SetAirBoundary(CON, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);

/*	//////////////空気領域追加節点////////////
	NODE.clear();
	FACE.clear();
	SetAirFineBoundary(CON, PART, TET, NODE, FACE);
	AddBoundaryData(CON, NODEall, FACEall, NODE, FACE, AIR);
	//*/
	
	//.nodeファイル作成
	MakeNodeFile(CON, NODEall, "all_boundary.node");
	//.polyファイルの出力
	MakePolyFile(CON, TET, NODEall, FACEall, "all_boundary.poly");

	cout<<"polyファイル読み込み"<<endl;
	in.load_poly("all_boundary");	//.polyを読み込んだら自動的に.nodeも読み込むらしい
	cout<<"自動分割開始"<<endl;
	tetrahedralize("pq1.1a1.0e-4AYYn", &in, &out);	//("pq1.1a2.0e-7AYYn", &in, &out);//1.1未満では切れない デフォルトはrqa1.1AYYn  pq1.3a1.67e-7AYYn
	out.save_nodes("output");
	out.save_elements("output");
	//out.save_faces("output");
	//out.save_neighbors("output");
	//*/
	cout<<"材質の修正開始\n";
	//cout<<"材質修正"<<endl;
	//ModifyAttribute(CON, TET, NODEall, ELEMall, in, out);
	ModifyAttribute_tetgenio(CON, TET, NODEall, ELEMall, in, out);	//tetgenioを直接編集

	//出力
	cout<<"TetView用ファイル出力"<<endl;
	//MakeNodeFile(CON, NODEall, "output_final.node");
	//MakeElemFile(CON, ELEMall, "output_final.ele");
	out.save_nodes("output_final");
	out.save_elements("output_final");

	//節点・要素データ取得
	cout<<"TetGenより節点・要素データ取得"<<endl;
	GetPointList(NODEall, in, out);
	GetTetrahedronList_full(ELEMall, in, out);
}