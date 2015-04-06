#ifndef MPSPARTICLE //ifndefでコンパイルしない
#define MPSPARTICLE

class vector3D
{
	double component[3];
public:
	vector3D(double *r0);
	void operator=(double vc[3]);
	double *get_comp(){return component;}
	double get_comp(int D){return component[D];}
};

class mpsparticle
{
public:
	unsigned ID;
	double r[DIMENSION];
	double initial_r[DIMENSION];
	double q0[DIMENSION];

	double u[DIMENSION];

	double r_temp[DIMENSION];
	double u_temp[DIMENSION];

	double n[DIMENSION];	//単位法線ベクトル(TetGen用)　どこで決定？

	double ang[DIMENSION+1];//角度・クォータニオン
	double ang_u[DIMENSION];
	double ang_temp[DIMENSION+1];//角度・クォータニオン
	double ang_u_temp[DIMENSION];
	double P;				//圧力
	double volumetric_strain;//体積ひずみ

	double p[DIMENSION];	//物理量

	bool contact;

	double dirP_st;	//表面張力による圧力ディリクレ値
	double dirP_em;	//電磁力による圧力ディリクレ値

	double h;				//エンタルピー
	double PND;				//reを用いた粒子数密度
	double PND2;
	double PND0;		
	double val;					//その都度適当な値の格納に利用
	double eforce[DIMENSION];	//電磁力

	int fly;				//飛散タイプ(TOUCH, GROUP, ISOLATION)
	int type;				//FLUID INWALL OUTWALL
	int materialID;			//材質番号 1ならdensity,2ならdensity2などを使用
	int surface;			//0:内部 1:表面
	int index;				//格納されている格子の番号
	int N;					//re内に存在する周辺粒子数・・・これがpublicになっているので添字にも代入されてしまいうる！！！
	int N2;					//re2内に存在する周辺粒子数
	int N3;					//re3内に存在する周辺粒子数
	int NEI[300];			//re内に存在する周辺粒子番号
	int NEI2[300];
	int NEI3[450];
	double PAcc[3];			//[murao] Permanent acceleration

	unsigned int toFEM;		//FEMに転送するかどうか。ON or OFF
	double dir_Pst;			//表面張力による圧力
	double dir_Pem;			//電磁力による圧力

	double ave_lambda;		//粒子の平均λ
	//PND=particle_number_density
	void set_PND(const double pnd){PND=pnd;}
	double get_PND() const{return PND;}

	//超弾性
	int N0;
	int NEI0[300];



};

//弾性体計算用派生クラス
class mpselastic: public mpsparticle
{
	double youngs_modulus;	//ヤング率
	double poisson_ratio;	//ポアソン比

	//要は周辺にある粒子のIDだけ持っておけば良いのだが、その都度位置関係を計算していたら面倒
	vector<int> initial_neighboursID;	//初期の周辺粒子ID
	vector<int> current_neighboursID;	//現在の周辺粒子ID

	//位置座標
	vector<double> initial_neighbourX;	//freeonで初期化する/ size==initial_neighboursN
	vector<double> initial_neighbourY;
	vector<double> initial_neighbourZ;
	vector<double> current_neighbourX;	//freeonで初期化する
	vector<double> current_neighbourY;
	vector<double> current_neighbourZ;

	vector<double> initial_distancebps;
	vector<double> current_distancebps;

	vector<vector3D> r0_ij;
	vector<vector3D> r0_ji;

	//力学パラメータ(iに属するものはelasticクラスからこちらへ移動する)
	double elastic_density;					//弾性体の密度（PNDから求める)

	double stress_accel[DIMENSION];			//応力加速度
	double stress_visco_accel[DIMENSION];	//せん断応力加速度
	double pressure_accel[DIMENSION];		//圧力加速度
	double total_accel[DIMENSION];			//加速度合計
	
	bool stop_on_floor; //床との接触判定
	bool acceleration_upward;//加速度の判定（上向きか下向きか）

public:

	mpselastic();

	//イテレータを戻り値に使うことは可能？

	vector<int>& get_initial_neighboursID() {return initial_neighboursID;}	//参照を返すほうが早いはず（直接initial_neighboursIDにアクセスしているのと同じになるので）
	vector<int>& get_current_neighboursID() {return current_neighboursID;}	//参照を返す方が早いはず（直接initial_neighboursIDにアクセスしているのと同じになるので）

	void set_current_neighboursID(const int j){current_neighboursID.push_back(j);}//current_neighboursID.size()==N
	void reset_current_neighboursID(){current_neighboursID.clear();}

	void get_initial_neighbours_position(int k0, double &X0, double &Y0, double &Z0){X0=initial_neighbourX[k0]; Y0=initial_neighbourY[k0]; Z0=initial_neighbourZ[k0];}
	void get_initial_neighbours_position(int k0, double *r0){r0[0]=initial_neighbourX[k0]; r0[1]=initial_neighbourY[k0]; r0[2]=initial_neighbourZ[k0];}
	void get_current_neighbours_position(int k, double &X, double &Y, double &Z){X=current_neighbourX[k]; Y=current_neighbourY[k]; Z=current_neighbourZ[k];}
	void get_current_neighbours_position(int k, double *r){r[0]=current_neighbourX[k]; r[1]=current_neighbourY[k]; r[2]=current_neighbourZ[k];}
	void set_current_neighbours_position(double X,double Y, double Z){current_neighbourX.push_back(X); current_neighbourY.push_back(Y); current_neighbourZ.push_back(Z);}
	void reset_current_neighbours_position(){current_neighbourX.clear(); current_neighbourY.clear(); current_neighbourZ.clear();}

	void set_current_distancebps(const double dis){current_distancebps.push_back(dis);}//これいらない。initial_distancebpsだけで十分
	double get_initial_distancebps(const int k){return initial_distancebps[k];}
	double get_current_distancebps(const int k){return current_distancebps[k];}
	vector<double>& get_initial_distancebps(){return initial_distancebps;}
	vector<double>& get_current_distancebps(){return current_distancebps;}

	void set_r0_ij(double *r0);
	void set_r0_ji(double *r0);
	vector<vector3D>& get_r0_ij(){return r0_ij;}
	vector<vector3D>& get_r0_ji(){return r0_ji;}

	double get_density() const{return elastic_density;}
	void set_density(const double d){elastic_density=d;}
	void reset_current_distancebps(){current_distancebps.clear();}//claer()しないと配列が増え続ける

	double get_pressure_accel(int D){return pressure_accel[D];}
	void set_pressure(const double *pressure){for(int D=0;D<3;D++) pressure_accel[D]=pressure[D];}
	
	double get_stress_accel(int D){return stress_accel[D];}
	double get_stress_visco_accel(int D){return stress_visco_accel[D];}
/*	void reset_stress(){for(int D=0;D<3;D++){ 
		stress_accel[D]=0;
		stress_visco_accel[D]=0;
	}
	}//*/

	void add_pressure(const double* press){for(int D=0;D<3;D++) pressure_accel[D]+=press[D];}

	void add_stress_accel(const double* st){for(int D=0;D<3;D++) stress_accel[D]+=st[D];}
	void add_stress_visco_accel(const double* sv){for(int D=0;D<3;D++) stress_visco_accel[D]+=sv[D];}

	void mul_pressure(const double coef){for(int D=0;D<3;D++) pressure_accel[D]*=coef;}

	void mul_stress_accel(const double coef){for(int D=0;D<3;D++) stress_accel[D]*=coef;}
	void mul_stress_visco_accel(const double coef){for(int D=0;D<3;D++) stress_visco_accel[D]*=coef;}

	void initialize_particles(elastic &ELAST, const int t);//これはelasticのメンバで良いのでは？
	void reset_particle_acceleration();

	double get_total_accel(int D){return total_accel[D];}
	void set_total_accel(const double *acceleration){for(int D=0;D<3;D++) total_accel[D]=acceleration[D];}


	//弾性定数のアクセサ
	double get_youngs_modulus(){return youngs_modulus;}
	double get_poisson_ratio(){return poisson_ratio;}
	void set_youngs_modulus(double E){youngs_modulus=E;}
	void set_poisson_ratio(double nu){poisson_ratio=nu;}

	//床反力
	bool get_stop_on_floor(){return stop_on_floor;}
	void set_stop_on_floor(bool flag){stop_on_floor=flag;}

	//加速度の判定
	bool get_acceleration_upward(){return acceleration_upward;}
	void set_acceleration_upward(bool flag){acceleration_upward=flag;}

	//加速度チェック
	void check_acceleration();

	//弾性係数チェック
	void check_elastic_constant();

	//mapに入れるとソートされてしまうので注意
	//接触判定で使う
	//freeonでneighbour_distanceを求めてinitialize_particlesでmapにコピーする
	map<int, double> distance0;		//初期の周辺粒子との距離(sizeはneighbours0), IDをキーにして粒子間距離を参照する //NEIの置き換えを図る
	map<int, double> distance;		//現在の周辺粒子との距離(sizeはneighbours), IDをキーにして粒子間距離を参照する
};

#endif