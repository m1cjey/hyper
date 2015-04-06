/*------------------------------------------------------------------------------------------------------------------------------

【自作TetGen用関数】
TetGenでメッシュ生成を行うために用いる関数は、全てtetgen_functionクラスに属しています。
まず、tetgen_functionクラスの適当なインスタンス(オブジェクト)を生成してcall_TetGenに入ってください。
以後、このクラスの関数内であれば、インスタンスの生成は不要です。

------------------------------------------------------------------------------------------------------------------------------*/

#ifndef TETFUNCH
#define TETFUNCH


class tetgen_function
{
public:

	///////////tetcall.cpp/////////////////////////////////////////////////////////////////////////////////

	void TetGen_elastic(mpsconfig &CON, vector<mpselastic>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//弾性体用　メッシュ生成
	void call_TetGen(mpsconfig &CON, vector<mpselastic>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//TetGenメッシュ生成の入口
	void TetGen_nanoe(mpsconfig &CON, vector<mpselastic>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//静電霧化用  メッシュ生成
	void TetGen_rayleigh(mpsconfig &CON, vector<mpselastic>&, int, int, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_element>&, vector<int>&);		//静電霧化用  メッシュ生成

	///////////tetfunc.cpp/////////////////////////////////////////////////////////////////////////////////

	void GetPointList(vector<tetgen_node>&, tetgenio&, tetgenio&);
	void Geteleattribute(vector<tetgen_node> &NODE, vector<tetgen_element> &ELEM,tetgenio &out);
	//void GetPointList_Fluid(vector<tetgen_node>&, tetgenio&, tetgenio&);
	void GetTetrahedronList(vector<tetgen_element>&, tetgenio&, tetgenio&);
	void GetMTetrahedronList(vector<tetgen_element>&, tetgenio&, tetgenio&);
	void GetTetrahedronList_full(vector<tetgen_element>&, tetgenio&, tetgenio&);
	void GetFacetList(vector<tetgen_facet>&, tetgenio&, tetgenio&, int);
	void GetFacetList_from_neigh(mpsconfig &CON, vector<tetgen_element>&, vector<tetgen_facet>&);
	void SelectFaceNode(mpsconfig &CON, vector<tetgen_node>&, vector<tetgen_facet>&);
	void DelDummyNode(mpsconfig &CON, vector<tetgen_node>&, vector<tetgen_facet>&, int);

	void DelThinTetrahedron(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void DelTetrahedron_OutsideDummy(mpsconfig &CON, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void SetRelation_NodeNode(mpsconfig &CON, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void SetRelation_NodeElem(mpsconfig &CON, vector<tetgen_node>&, vector<tetgen_element>&);
	void SetRelation_ElemElem(mpsconfig &CON, vector<tetgen_node>&, vector<tetgen_element>&);

	void MakeNodeFile(mpsconfig &CON, vector<tetgen_node>&, char*);
	void MakeNodeFile_NonAttributeAndBoundary(mpsconfig &CON, vector<tetgen_node>&, char*);
	void MakeElemFile(mpsconfig &CON, vector<tetgen_element>&, char*);
	void MakeFaceFile(mpsconfig &CON, vector<tetgen_facet>&, char*);
	void MakePolyFile(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&, char*);
	void MakeSmeshFile(mpsconfig &CON, vector<tetgen_facet>&, char*);

	void SetTRANS(vector<tetgen_node>&, vector<int>&);

	
	void SetElastBoundary					(mpsconfig &CON, vector<mpselastic>&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetAirBoundary(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetAirFineBoundary(mpsconfig &CON, vector<mpselastic>&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetPlateElectrodeBoundary(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetColumnElectrodeBoundary(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetBaseBoundary(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetMagnetBoundary(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetCOILBoundary(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);
	void SetIRONBoundary(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_facet>&);

	void UniteBoundaryData(mpsconfig&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<tetgen_facet>&, vector<int>&);
	void AddBoundaryData(mpsconfig &CON, vector<tetgen_node>&, vector<tetgen_facet>&, vector<tetgen_node>&, vector<tetgen_facet>&, int);

	void DecisionAttribute(vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void ModifyAttribute(mpsconfig&, tetgen_config&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void ModifyAttribute_tetgenio(mpsconfig &CON, tetgen_config&, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);
	void FineElement(mpsconfig &CON, vector<tetgen_node>&, vector<tetgen_element>&, tetgenio&, tetgenio&);

	double Distance(tetgen_node&, tetgen_node&);
	void CalcBarycentricElement(mpsconfig&, vector<tetgen_node>&, vector<tetgen_element>&);

};

#endif