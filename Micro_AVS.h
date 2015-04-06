#pragma once
/*Micro_AVSファイル出力クラス*/
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

class Micro_AVS{
	std::string Filename;
	int Step;
	int Dimension; 

	std::vector<double> Xr;
	std::vector<double> Yr;
	std::vector<double> Zr;

	std::vector<double> Xdata;
	std::vector<double> Ydata;
	std::vector<double> Zdata;

	static std::vector<std::string> Filelist;
	static bool fast_flag;
public:
	Micro_AVS::Micro_AVS();

	//Micro＿AVSで出力する値のリスト作成
	virtual void make_list(double xr,double yr,double xdata,double ydata);
	virtual void make_list(double xr,double yr,double zr,double xdata,double ydata,double zdata);
	//Micro_AVS出力
	virtual void Output_vector_MicroAVS(std::string fname,int step);
	virtual void Output_mgf_MicroAVS(std::string fname,int max_step);
	//今までに作成したフォルダ名を取得
	virtual void Get_folder_name();

protected:
	//ベクトルファイル
	//2次元
	virtual void _Wright_vector_2d();
	//3次元
	virtual void _Wright_vector_3d();
	
	//mgfファイル
	//2次元
//	virtual void _Wright_mgf_2d(); //実体はない
	//3次元
	virtual void _Wright_mgf_3d();

	//新しく作成したフォルダ名を記憶する
	virtual void _Foldername_storing();




};