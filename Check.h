#pragma once
/*
Checkクラス

後から見直したときにわかるように解析条件・計算情報を記録するクラス

*/
class Check
{
	string checkf;

public:
	Check();
	virtual ~Check();
	////////////////年月日取得//////////////////
	virtual string Set_y_m_d();

	///////////解析条件ファイル出力/////////////
	virtual void Out_put_config();

	/////////////クーラン数の確認///////////////
	virtual void Courant_condition(vector<mpselastic> &PART);
};