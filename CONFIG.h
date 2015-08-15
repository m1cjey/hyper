////新しい粒子を追加したときいじらないといけないのは①u_laplacian_fのノンスリップ②P_gradient③P_gradient2④calc_Temperatureのk[i]
////⑤ freeon
#ifndef config
#define config

using namespace std;

class mpsconfig
{    
//解析条件

	double dt;//時間刻み幅
	double dt_for_FEM; //FEM用時間刻み
	int step;
	int current_step;
	double current_time;

	int avoid_step;
	int avoid_step2;
	int avoid_step3;
	int avoid_step4;
	int avoid_step5;
	int avoid_step6;
	int avoid_step7;

	bool flag_cut_x_movie;
	bool flag_cut_y_movie;

	int dimension;//次元
	double maxX;//解析領域
	double minX;
	double maxY;
	double minY;
	double maxZ;
	double minZ;
	
//弾性体（流体）の物性値
	
	double MRE_density;		//MRE密度
	double Silicone_density;	//Silicone密度
	double nensei;		//粘性係数
	double sigma;		//表面張力係数
	double Cp;			//定圧比熱
	double k;			//熱伝導率
	double latent_H;	//潜熱(latent_heat)
	double MP;			//融点(melting point)
	double CTE;			//熱膨張係数[1/K] 
	double E_m;			//ヤング率
	double v_m;			//ポアソン比
	double E_e;
	double v_e;
	bool nonlinear_elastic;
	
//粒子配置用
	
	int fluidwidth;		//流体の代表長さ
	double distancebp;	//初期粒子間距離
	double wlength;		//左右の壁の距離。流体の何倍か
	double height;		//流体の重心高さ
	double ground;		//床上面の位置
	
//粒子法用パラメータ
	
	double re;		//一般的な粒子半径re
	double re2;		//ラプラシアン用のre
	double re3;		//表面張力用のre
	double re4;		//freeon
	double re_elastic;
	double beta;	//β
	int dx;			//index用格子幅。最大粒子半径を含む整数にすること
	double times;	//ベクトルプロット用の倍率

	int motion_interval;
	
	//表面張力関係
	
	int surface_tension; //表面張力スイッチ　0=OFF 1=教科書 2=粒子間引力
	int smooth;			//スムージングスイッチ　1=ON 0=OFF
	int SFT;				//1=ON 0=OFF 表面張力(surface tension)の温度(T)依存性スイッチ
	int smn;				//スムージングの回数
	int smth_sumP;    //for ver3. curvのスムージング 0=OFF 1=ON
	int wall_direct;		//for ver2. BDWALLの法線ベクトル計算方法　0=無視 1=通常 2=水平方向 3=垂直方向
	int non_slip;		//for ver.5 粘着条件 ONならBDWALL上に仮想的に液体が存在する設定。その際はBDWALLの配置に留意
	int suf_modify;		//曲率計算で得られた曲面に一致するように粒子位置を修正するか、しないか
	int interpolate_curv;//曲率が計算されなかった粒子の曲率を、周辺粒子の値から補間
	
	//圧力関係

	int iteration_count;//圧力計算を行う回数 通常は1に設定
	int solution;		//連立方程式の解放　0=CG method 1=ICCG method
	int B_of_P;			//圧力計算の解行列 0=教科書 1=速度発散 2=0+1
	double w_div;		//発散の重み 
	int HL_sw;			//ラプラシアンに対し高次の離散化を行う 1=ON 0=OFF
	int	dir_for_P;		//圧力計算の際、表面に非ゼロのディリクレ値を代入するか、しないか
	int initialP;		//圧力計算時に初期値として現圧力を代入するかどうか　1=ON 0=OFF
	double CGep;		//圧力計算時における収束判定
	int	omp_P;			//圧力計算時のCG法でopenMPを使用するか 1=ON 0=OFF
	int negativeP;		//負圧　1=許可 0=不許可
	int set_P;			//圧力強制決定スイッチ 1=ON 0=OFF
	int Pgrad;			//圧力勾配計算法 0=教科書 1=表面のみ法線
	int minP;			//最少圧力スイッチ 1=ON 0=OFF
	int ave_Pdirect;	//圧力勾配用法線ベクトルの平均化 1=ON 0=OFF
	int Pgrad_order;	//圧力勾配ver.4における精度 1=線形 2=2次
	int artP_sw;		//人工圧力計算フラグ 0=OFF 1=artP等方 2=Monaghan
	double artP;		//人工圧力値[Pa]
	double Pgrad_times;	//圧力勾配をプロットする際の、表面張力に対する倍率 通常は1に設定
	int P_AVS;		//microAVS用の圧力ファイルを出力するstep間隔。0ならOFF

	//温度場関係
	
	int T_field;			//温度場解析　　1=ON 0=OFF
	int insulate;			//壁との断熱状態 0=断熱　1=非断熱
	int T_laplacian;		//温度のラプラシアン。　0=教科書　1=λ[i] 2=発散・勾配
	double wall_density;	//壁の密度[kg/m^3]
	double wall_Cp;			//壁の比熱[J/(kg・K)]
	double wall_k;			//壁の熱伝導率[W/(m・K)]
	double roomT;			//室温 int型でよい
	int air_cool;			//空気との熱伝達を考慮するかしないか 1=ON 0=OFF
	int T_expansion;		//熱膨張計算　1=ON 0=OFF
	int buoyant;			//浮力(密度のブシネスク近似)  1=ON 0=OFF
	double TplotZ;			//3D解析において、XY平面の温度を出力するときのＺ座標
	int T_AVS;				//microAVS用の温度ファイルを出力するstep間隔。0ならOFF

	//電磁場の解法

	int	EM_method;		//電磁場の解法 0=OFF 1=FEM 2=BEM 3=磁気モーメント法
	int	EM_calc_type;	//0=OFF 1=電場 2=磁場 3=磁場(渦電流)
	int	EM_interval; //電磁場計算を何ステップに一回行うか。通常は1に設定
	int region_shape;	//解析領域形状　0=立方体 1=円筒
	double XL;			//解析領域左端
	double XR;    //解析領域右端
	double YU;    //解析領域上端
	double YD;    //解析領域下端
	double ZU;			//解析領域Z方向上端
	double ZD;			//解析領域Z方向下端
	double RU;			//解析領域が円筒形となるときのその半径

	//メッシュ
	int mesh_input;		//meshの作成方法 0:自分 1:Magnet 2:TetGen
	int remesh_sw;		//ON:remesh領域を考慮し、そこだけremesh OFF:FULL領域を毎回分割
	double boxalpha;		//スーパーボックスの大きさ 1～10の実数にすること。通常は2に設定
	int fine;			//節点の自動追加を行うかどうか 0:OFF 1:ON
	double co_fine;			//節点自動追加における係数.大きいほど物体から離れたら粗くなる
	int add_points;		//自動節点追加数
	int poly_memory;		//poly3D()で確保するメモリの数
	int air_layer;		//物体の周囲に空気層を生成する層数(0なら生成しない)
	double layer_depth;		//空気層の幅。初期粒子間距離の何倍か
	int mesh_output_interval;//メッシュ情報を、有限要素法のステップに対して、何ステップに一度出力するか。
	double FEMCGep;
	double MRTRep;

	//BEM関係
	int BEM_elm_type;	//要素タイプ 0:節点要素 1:辺要素
	int BEM_M_solver;	//境界要素法における行列解法 0:ガウスの消去法 1:BiCGStab法
	int gauss_N;			//Gauss積分の評価点数 3,4,7のどれか
	int tree_sw;			//BEMでtree法を使うか、使わないか
	double tree_mesh_size;	//tree法で使用する最大レベルのセルのサイズ。leの何倍か
	int p_max;			//tree法の無限級数の最大項数。0から6の数字に対応

	//FEM関係
	int FEM_elm_type;	//要素タイプ 0:節点要素 1:辺要素
	int FEM_smn;			//電磁力スムージング回数　0ならOFF
	int max_DN;			//ディリクレ型境界条件をとる最大節点数
	int FEMCG;			//FEMにおける行列解法 0:CG 1:ICCG
	double CGaccl;			//CG,ICCG法における加速ファクタ　CGaccelerator
	double EMCGep;			//電磁場のICCGの収束判定
	double FEMtimes;		//電磁力をプロットする際の、表面張力に対する倍率 通常は1に設定
	double legend_F;		//F.datの凡例に出力する力[N]
	int plot_F_type;	//F.datの出力タイプ  0=[N]表示 1=[Pa]表示
	double legend_F_Pa;	//F.datの凡例に出力する力[Pa]

	///////FEM ver2関係
	double surface_depth;			//表面の厚みの半分(これにleをかけたものが本当の厚み(の半分))
	
	///////電界計算
	double V;			//電圧
	int V_step;			//電圧印加開始ステップ
	double r_perm;		//比誘電率 (relative permittivity)
	int V_con;			//電圧条件　0:パルス　1:リニア　2:時定数
	double initial_V;	//電圧条件リニアのときの初期電圧
	double E_times;		//ファイル出力する際の、電界の倍率
	int eleforce;		//静電力計算方法 1:節点力法 2:積分面
	int charge;			//電荷を考慮するかしないか。0=OFF 1=ON
	int plot_E_type;	//電界出力タイプ 1=ベクトル 2=スカラー

	//磁場計算
	
	int J_input_way;		//電流密度入手方法 0:自分 1:ソフト
	double J0;				//強制電流密度[A/m2]
	double I0;				//強制電流値[A]
	double RP;				//比透磁率(Relative Permeability)
	double ele_conduc;		//電気伝導率
	int Hz;					//交流の周波数
	int div_Hz;				//１周期の分割数(解析精度)
	int jheat;				//渦電流による発熱を考慮するか　0=OFF 1=ON
	int	m_force;			//電磁力計算方式 0=マクスウェルの応力 1=体積力
	int NLBHI;				//体積力において、要素Ｂから要素Ｈを求める際に非線形性を考慮するか、しないか(non linier B H inverter)
	int NLMH;				//Ｍの算出に非線形性を考慮するか、しないか
	double magnet_H;		//永久磁石の大きさ
	double magnet_B;		//永久磁石の強さ[T]
	double magnet_angle;	//永久磁石の着磁方向 0なら+Z方向となる。そこから角度をつけたいなら、その角度[deg]を入力する
	double magnet_r;		//永久磁石の半径
	double magnet_Z;		//永久磁石の中心のZ座標
	int uniform_B_sw;		//解析領域中に一様磁場を発生させるか否か 0=OFF 1=ON
	int magnetic_layers;	//永久磁石周辺の空気層の数
	double uniform_B;		//一様磁場の大きさ[T]
	double B_times;			//ファイル出力する際の、磁束密度の倍率
	int plot_B_type;		//磁束密度出力タイプ 1=ベクトル 2=スカラー

	int FEM_calc_type;		//0=OFF 1=電場 2=磁場 3=磁場(渦電流)
	int ele_type;
	//各種スイッチ　どのような計算を考慮するか、しないか
	
	double g;					//重力加速度　通常は-9.8とすること
	int restart;				//restart用スイッチ
	int autosave;			//オートセーブ間隔。無効にしたいときは大きな数字を代入しておく
	double courant;			//クーラン数条件 0ならOFF
	int modify_position;	//粒子間距離がleより小さい場合これを修正するか、しないか
	int vis_calc_type;		//粘性項計算手法　0=陽解法 1=陰解法
	int wall_adheision;		//壁の粘性状態　1=ノンスリップ 0=ﾌﾘｰｽﾘｯﾌﾟ
	int laplacian;			//ラプラシアン。　0=教科書　1=λ[i] 2=発散・勾配
	int	vis_solver;			//粘性項を陰解析で解く際の行列ｿﾙﾊﾞｰ 0:CG 1:ICCG
	int initial_u;			//粘性項を陰的に解く際に、初期値として現在速度を入力するか、しないか
	int temporary_r;		//陽解析後に仮の位置を計算するかしないか。 1=ON 0=OFF 通常はONに設定
	
	int fix_center;			//1=ON 0=OFF
	int freeon;				//粒子依存関係関数　1:並列化可能 2:並列不可
	int freeon3sw;			//freeon3を計算するかしないか 1=ON 0=OFF
	int surface_judge2;		//surface_judge2()を使用するかしないか  1=ON 0=OFF
	int move_prtcl;			//移動粒子を考慮するかしないか 1=ON 0=OFF
	int move_u_dirct;		//移動粒子を移動する方向　現在は±X方向=±1,±Y方向=±2,±Z方向=±3
	double move_speed;		//移動粒子の移動速度[m/s]
	int check_something;	//check_something()を実行するかしないか 1=ON 0=OFF
	
	int model_number;
	int model_set_way;		//modelをセットする方法　0=正方格子 1=MD

	double R1;

	bool switch_FEM;			//FEMを実行するかしないか OFF: 0, ON: 1
	bool switch_vis;			//粘弾性計算するかしないか OFF: 0, ON: 1

	///////速度プロット変数
	int speed_plot_particle;	//速度をプロットする粒子の種類 1=すべて 2=fluid
	double speedtimes;			//速度プロット時の、座標に対する速度の倍率
	int speed_face;				//speed.datの出力面 0=YZ平面 1=XZ 2=XY
	double speed_face_p;		//3D解析時のspeed.datの出力面の座標
	int ax_sym_modify;			//3D時のspeed.datに関して、軸対称による出力修正を行うか否か　1=ON 0=OF
	int flat_speed_plot;		//水平方向の速度をプロットするかしないか
	double flat_speed_p;		//flat_speed.datの出力面の座標
	int relative_speed;			//重心に対する相対速度を出力するかしないか 1=ON 0=OFF
	int speed_AVS;				//microAVSによる3D速度分布出力するかしないか 1=ON 0=OFF
	double legend_speed;		//speed.datの凡例に出力する速度[m/s]
	int set_zero_speed;			//restart時に速度をゼロセットするかしないか  1=ON 0=OFF

	//圧力チェック用
	int pressure_face;
	double pressure_face_p;
	double pressure_times;
	double legend_pressure;		//speed.datの凡例に出力する速度[m/s]

	///////GPU関係
	int M_form;				//CG法における係数行列の格納方法 CSR_scl,CSR_vec,ELLの3つをｻﾎﾟｰﾄ
	int MAX_thread;			//ひとつのSMあたりの最大スレッド数　ふつうは512
	
	//ファイル出力変数
	
	int interval;			//AVS用ステップ間隔
	int AVS;				//0:普通　1:圧力　2:温度

	//////AVSファイル出力
	int avs_eforce_interval;//AVS電磁力ファイルを何回に1回出力するか 0=出力しない
	int F_interval;//電磁力F.datを何ステップ毎に出力するか
	int avs_mesh1_interval;	//AVS電位ファイル(断面)を何回に1回出力するか 0=出力しない
	int avs_mesh2_interval;	//AVSメッシュファイル(断面)を何回に1回出力するか 0=出力しない
	int avs_mesh3_interval;	//AVSメッシュファイル(材質)を何回に1回出力するか 0=出力しない

	double P_size_AVS;		//出力する粒子のサイズ(leの何倍か)
	double maxT;			//AVS(2)における最大温度
	double minT;			//AVS(2)における最小温度

	double max_pressure;
	double min_pressure;

	bool flag_modify_density;

	double times_Pa;

	double ave_P_for_FEM_flag;
	bool Ferror;
	bool poise_flag;
	double length;


	//超弾性計算
	int flag_ELAST;
	int flag_HYPER;
	double hyper_density;
	double c01;
	double c10;
	int flag_wall;
	double r_z_wall;
	double h_dis;
	double h_vis;
	int flag_vis;
	int tension_test;	//引っ張り試験解析用15/2/8
	int nr_time;

public:
	mpsconfig();
	void Set_length(double len){length=len;}
	double Get_length(){return length;}
	double get_dt(){return dt;}
	int get_step(){return step;}
	int get_current_step(){return current_step;}
	void set_current_step(int cs){current_step=cs;}
	double get_current_time(){return current_time;}
	void set_current_time(const double ct){current_time=ct;}
	int get_avoid_step(){return avoid_step;}
	int get_avoid_step2(){return avoid_step2;}
	int get_avoid_step3(){return avoid_step3;}
	int get_avoid_step4(){return avoid_step4;}
	int get_avoid_step5(){return avoid_step5;}
	int get_avoid_step6(){return avoid_step6;}
	int get_avoid_step7(){return avoid_step7;}

	bool get_cut_x(){return flag_cut_x_movie;}
	bool get_cut_y(){return flag_cut_y_movie;}

	int get_dimension(){return dimension;}
	double get_maxX(){return maxX;}
	double get_minX(){return minX;}
	double get_maxY(){return maxY;}
	double get_minY(){return minY;}
	double get_maxZ(){return maxZ;}
	double get_minZ(){return minZ;}
	
	double get_density(){return MRE_density;}
	double Get_MRE_density(){return MRE_density;}
	double Get_Silicone_density(){return Silicone_density;}
	bool get_modify_density(){return flag_modify_density;}
	double get_nensei(){return nensei;}
	double get_sigma(){return sigma;}
	double get_vis(){return nensei/MRE_density;}	//動粘性係数
	double Get_MRE_vis(){return nensei/MRE_density;}
	double Get_Silicone_vis(){return nensei/Silicone_density;}
	double get_Cp(){return Cp;}
	double get_k(){return k;}
	double get_latent_H(){return latent_H;}
	double get_MP(){return MP;}
	double get_CTE(){return CTE;}
	double get_E_m(){return E_m;}
	double get_v_m(){return v_m;}
	double get_E_e(){return E_e;}
	double get_v_e(){return v_e;}

	bool get_Ferror(){return Ferror;}//Feエラー
	void set_Ferror(bool ero){Ferror=ero;}
	int get_fluidwidth(){return fluidwidth;}
	double get_distancebp(){return distancebp;}
	double get_wlength(){return wlength;}
	double get_height(){return height;}
	
	double get_re(){return re;}
	double get_re2(){return re2;}
	double get_re3(){return re3;}
	double get_re4(){return re4;}
	double get_re_elastic(){return re_elastic;}
	double get_beta(){return beta;}
	int get_dx(){return dx;}
	double get_times(){return times;}

	int get_motion_interval(){return motion_interval; }

	int get_X_mesh(){return (int)((maxX-minX)/(distancebp*dx)+0.00000000000001);}		//X軸方向の格子数 丸め誤差を防ぐために0.001を足している（？）
	int get_Y_mesh(){return (int)((maxY-minY)/(distancebp*dx)+0.00000000000001);}		//Y軸方向の格子数 丸め誤差を防ぐために0.001を足している
	int get_Z_mesh(){return (int)((maxZ-minZ)/(distancebp*dx)+0.00000000000001);}		//Z軸方向の格子数 丸め誤差を防ぐために0.001を足している
	int get_number_of_mesh(){return (int)((maxX-minX)/(distancebp*dx)*(maxY-minY)/(distancebp*dx)*(maxZ-minZ)/(distancebp*dx)+0.001);}//格子数：X_mesh*Y_mesh*Z_mesh
	
	int get_surface_tension(){return surface_tension;}
	int get_smooth(){return smooth;}
	int get_SFT(){return SFT;}
	int get_smn(){return smn;}
	int get_smth_sumP(){return smth_sumP;}
	int get_wall_direct(){return wall_direct;}
	int get_non_slip(){return non_slip;}
	int get_suf_modify(){return suf_modify;}
	int get_interpolate_curv(){return interpolate_curv;}

	int get_iteration_count(){return iteration_count;}
	int get_
		(){return Pgrad;}
	int get_ave_Pdirect(){return ave_Pdirect;}
	int get_solution(){return solution;}
	int get_B_of_P(){return B_of_P;}
	double get_w_div(){return w_div;}
	int get_HL_sw(){return HL_sw;}
	int get_dir_for_P(){return dir_for_P;}
	int get_initialP(){return initialP;}
	double get_CGep(){return CGep;}
	int get_omp_P(){return omp_P;}
	int get_negativeP(){return negativeP;}
	int get_minP(){return minP;}
	int get_set_P(){return set_P;}
	int get_artP_sw(){return artP_sw;}
	double get_artP(){return artP;}
	double get_Pgrad_times(){return Pgrad_times;}
	int get_Pgrad_order(){return Pgrad_order;}
	int get_P_AVS(){return P_AVS;}
	
	int get_T_field(){return T_field;}
	int get_insulate(){return insulate;}
	int get_T_laplacian(){return T_laplacian;}
	double get_wall_density(){return wall_density;}
	double get_wall_Cp(){return wall_Cp;}
	double get_wall_k(){return wall_k;}
	double get_roomT(){return roomT;}
	int get_air_cool(){return air_cool;}
	int get_T_expansion(){return T_expansion;}
	int get_buoyant(){return buoyant;}
	double get_TplotZ(){return TplotZ;}
	int get_T_AVS(){return T_AVS;}

	int get_EM_method(){return EM_method;}
	int get_EM_calc_type(){return EM_calc_type;}
	int get_EM_interval(){return EM_interval;}
	int get_region_shape(){return region_shape;}
	double get_XL(){return XL;}
	double get_XR(){return XR;}
	double get_YU(){return YU;}
	double get_YD(){return YD;}
	double get_ZU(){return ZU;}
	double get_ZD(){return ZD;}
	double get_RU(){return RU;}

	int get_mesh_input(){return mesh_input;}
	int get_remesh_sw(){return remesh_sw;}
	double get_boxalpha(){return boxalpha;}
	int get_fine(){return fine;}
	double get_co_fine(){return co_fine;}
	int get_add_points(){return add_points;}
	int get_poly_memory(){return poly_memory;}
	int get_air_layer(){return air_layer;}
	double get_layer_depth(){return layer_depth;}
	int get_mesh_output_interval(){return mesh_output_interval;}
	
	int get_BEM_elm_type(){return BEM_elm_type;}
	int get_BEM_M_solver(){return BEM_M_solver;}
	int get_gauss_N(){return gauss_N;}
	int get_FEM_elm_type(){return FEM_elm_type;}
	int get_FEM_smn(){return FEM_smn;}
	int get_max_DN(){return max_DN;}
	int get_FEMCG(){return FEMCG;}
	double get_CGaccl(){return CGaccl;}
	double get_EMCGep(){return EMCGep;}
	double get_FEMtimes(){return FEMtimes;}
	double get_legend_F(){return legend_F;}
	int get_tree_sw(){return tree_sw;}
	double get_tree_mesh_size(){return tree_mesh_size;}
	int get_p_max(){return p_max;}
	int    get_plot_F_type(){return plot_F_type;}
	double get_legend_F_Pa(){return legend_F_Pa;}

	double get_surface_depth(){return surface_depth;}
	
	double get_V(){return V;}
	int get_V_step(){return V_step;}
	double get_r_perm(){return r_perm;}
	int get_V_con(){return V_con;}
	double get_initial_V(){return initial_V;}
	double get_E_times(){return E_times;}
	int	get_eleforce(){return eleforce;}
	int get_charge(){return charge;}
	int get_plot_E_type()	{return plot_E_type;}
	
	int get_J_input_way(){return J_input_way;}
	double get_J0(){return J0;}
	double get_I0(){return I0;}
	double get_RP(){return RP;}
	double get_ele_conduc(){return ele_conduc;}
	int get_Hz(){return Hz;}
	int get_div_Hz(){return div_Hz;}
	int get_jheat(){return jheat;}
	int get_m_force(){return m_force;}
	int get_NLBHI(){return NLBHI;}
	int get_NLMH(){return NLMH;}
	double get_magnet_H(){return magnet_H;}
	double get_magnet_r(){return magnet_r;}
	double get_magnet_Z(){return magnet_Z;}
	double get_magnet_angle(){return magnet_angle;}
	double get_magnet_B(){return magnet_B;}
	int get_uniform_B_sw(){return uniform_B_sw;}
	int get_magnetic_layers(){return magnetic_layers;}
	double get_uniform_B(){return uniform_B;}
	double get_B_times(){return B_times;}
	int get_plot_B_type(){return plot_B_type;}

    int get_FEM_calc_type()	{return FEM_calc_type;}
	int get_ele_type(){return ele_type;}
	double get_FEMCGep() {return FEMCGep;}
	double get_MRTRep() {return MRTRep;}

	double get_g(){return g;}
	int get_restart(){return restart;}
	int get_autosave(){return autosave;}
	double get_courant(){return courant;}
	int get_modify_position(){return modify_position;}
	int get_vis_calc_type(){return vis_calc_type;}
	int get_wall_adheision(){return wall_adheision;}
	int get_laplacian(){return laplacian;}
	int get_vis_solver(){return vis_solver;}
	int get_initial_u(){return initial_u;}
	int get_temporary_r(){return temporary_r;}
	
	int get_fix_center(){return fix_center;}
	int get_freeon(){return freeon;}
	int get_freeon3sw(){return freeon3sw;}
	int get_surface_judge2(){return surface_judge2;}
	int get_move_prtcl(){return move_prtcl;}
	int get_move_u_dirct(){return move_u_dirct;}
	double get_move_speed(){return move_speed;}
	int get_check_something(){return check_something;}
	int get_set_zero_speed(){return set_zero_speed;}

	int get_model_number(){return model_number;}
	int get_model_set_way(){return model_set_way;}

	double get_R1(){return R1;}

	bool get_FEM_flag(){return switch_FEM;}
	void set_FEM_flag(bool sw){switch_FEM=sw;}
	bool get_vis_flag(){return switch_vis;}

	int get_speed_plot_particle(){return speed_plot_particle;}
	double get_speedtimes(){return speedtimes;}
	int get_speed_face(){return speed_face;}
	double get_speed_face_p(){return speed_face_p;}
	int get_ax_sym_modify(){return ax_sym_modify;}
	int get_flat_speed_plot(){return flat_speed_plot;}
	double get_flat_speed_p(){return flat_speed_p;}
	int get_relative_speed(){return relative_speed;}
	int get_speed_AVS(){return speed_AVS;}
	double get_legend_speed(){return legend_speed;}

	int get_pressure_face(){return pressure_face;}			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ
	double get_pressure_face_p(){return pressure_face_p;}		//3D解析時のspeed.datの出力面の座標
	double get_pressure_times(){return pressure_times;}
	double get_legend_pressure(){return legend_pressure;}


	int get_M_form(){return M_form;}
	int get_MAX_thread(){return MAX_thread;}
	
	int get_interval(){return interval;}
	int get_AVS(){return AVS;}
	double get_P_size_AVS()	{return P_size_AVS;}

	double get_maxT(){return maxT;}
	double get_minT(){return minT;}
	double get_max_pressure(){return max_pressure;}
	double get_min_pressure(){return min_pressure;}

	double get_particle_mass();
	void modify_density(double coefficient){MRE_density*=coefficient;}

	double get_particle_volume(){ return (get_particle_mass()/get_density());}

	int get_F_interval(){return F_interval;}

	int get_avs_eforce_interval(){return avs_eforce_interval;}
	int get_avs_mesh1_interval(){return avs_mesh1_interval;}
	int	get_avs_mesh2_interval(){return avs_mesh2_interval;}
	int	get_avs_mesh3_interval(){return avs_mesh3_interval;}

	double get_times_Pa(){return times_Pa;}

	double get_ave_P_for_FEM_flag(){return ave_P_for_FEM_flag;}

	bool get_nonlinear_elastic_flag(){return nonlinear_elastic;}
	void set_nonlinear_elastic_flag(bool flag){nonlinear_elastic=flag;}

	void change_step_size(){dt=dt_for_FEM;} //時間刻み幅をFEM計算用に変更する

	double get_ground_position(){return ground;}
	void set_ground_position(double g){ground=g;}

	bool get_poise_flag(){return poise_flag;}
	void set_poise_flag(bool flg){poise_flag=flg;}



	//超弾性計算
	int get_flag_ELAST(){return flag_ELAST;}
	int get_flag_HYPER(){return flag_HYPER;}
	double get_hyper_density(){return hyper_density;}
	double get_c10(){return c10;}
	double get_c01(){return c01;}
	int get_flag_wall(){return flag_wall;}
	double get_h_dis(){return h_dis;}
	double get_h_viscousity(){return h_vis;}
	int get_flag_vis(){return flag_vis;}
	int get_nr(){return nr_time;}
};

class elastic; //前方宣言 これがないとエラーが出る?

#endif