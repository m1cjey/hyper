#ifndef CLASSH
#define CLASSH


//TetGenの設定
class tetgen_config
{
public:
	tetgen_config();
	tetgen_config(mpsconfig &CON);

	//磁石寸法
	double magnet_height;
	double magnet_radius;

	//静電霧化装置の寸法
	double radius_column;		//電極半径
	double length_column;		//電極長さ
	double height_plate;		//平板高さ
	double length_plate;		//平板の一辺の長さ
	double thickness_plate;		//平板の厚み
	double length_base;			//土台の一辺の長さ
	double thickness_base;		//土台の厚み

	//メッシュの粗さ(境界の1要素の辺の長さをどれぐらいにするか)
	double fine_air;
	double fine_plate_t;
	double fine_plate_L;
	double fine_column_L;
	double fine_base;

	//水滴表面のメッシュ層の設定
	int num_layer_out;			//流体外側メッシュ層数
	int num_layer_in;			//流体内側メッシュ層数
	double thick_layer;			//境界メッシュ1層の厚さ
	
	//長い要素を削除する閾値
	double del_length;			//leの何倍以上の辺を持つ要素を消すか
};

//ノード
class tetgen_node
{
public:
	int id;
	double r[3];
	int attribute;//性質
	int boundary;//境界
	vector<int> nei_node;
	vector<int> nei_elem;
	int part_no;
};

//ファセット
class tetgen_facet
{
public:
	int id;
	int node[3];
	int boundary;
};

//要素
class tetgen_element
{
public:
	int id;
	int node[4];
	double g[3];	//重心座標
	int attribute;
	double volume;
	int nei_elem[4];
};

//
class region_attribute_list
{
public:
	int id;
	double r[3];
	int region_number;
	int region_attribute;
};


#endif