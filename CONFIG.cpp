#include "stdafx.h"	
//#include <fstream>

//コンストラクタ
//解析条件
mpsconfig::mpsconfig()
{
/*	ifstream fin("config.txt");
	string buffer;
	stringstream ss;
	if(fin.fail()){
		cout<<"コンフィグファイルを開けませんでした";
		exit(1);
	}
	//boost::formatで実装するのがラクでスマート・・・
	//get_char()で実装して#(コメント文)を読み飛ばしても良い
	getline(fin, buffer);
	ss<<buffer;
	ss>>step;
	fin.close();*/
	
	step=40000;				//全step数	step=20000;//40000;	//30000;//10000;;	//79*20+1;
	switch_FEM=false;		//FEMを実行するかしないか false
	nonlinear_elastic=false;	//弾性体非線形計算するかtrue
	switch_vis=OFF;			//粘性項計算するかしないか・・・これはあとで消す
	FEMCG=2;				//FEMにおける行列解法 0:CG 1:ICCG 2:並列ICCG 3:MRTR 4:ICMRTR

//	dt= (switch_FEM==OFF) ? 1.0e-5: 5.0e-6; //0.0001;不安定要因！ 0.00001:推奨(Courant数考えて) //Cf. dt_for_FEM=0.000001/2;
	dt=2.5e-5;
	dt_for_FEM=1.0e-3;
	//FEMだと0.000001で止まる・・・
	current_step=1;
	current_time=0.0;
	dimension=3;

	interval=4; //10	//particle_movie.mgfの出力間隔。2以上の整数にすること
	EM_interval=1;//1	//電磁場計算を何ステップに一回行うか。通常は1に設定
	motion_interval=1;	//運動方程式を何回に一回解くか
	
	//この圧力以上になったらFEMスタート cf. PostProcessing.cpp
	//物性値を少し変えてギリギリを狙うと到達しない可能性がある・・・
	ave_P_for_FEM_flag=8000000000;//80.0;//75.0;//70.0;

//モデル
	model_number=21;			//4:引っ張り試験片 7:MREアクチュエータ 12:剛体
	model_set_way=1;		//modelをセットする方法　0=正方格子 1=MD

//モデル１,11専用
/*	R1=9*0.001;
	avoid_step=0;
	avoid_step2=2425;
	avoid_step3=2435;
	avoid_step4=2445;
	avoid_step5=5935;
	avoid_step6=8410;
	avoid_step7=10710;*/

//出力動画
	flag_cut_x_movie=OFF;
	flag_cut_y_movie=OFF;

//電磁力計算
	region_shape=1;//1;		//解析領域形状　0=立方体 1=円筒
	EM_method=1;			//電磁場の解法 0=OFF 1=FEM 2=BEM 3=磁気ﾓｰﾒﾝﾄ法 4=FEM2
	EM_calc_type=2;			//0=デローニのみ 1=電場 2=静磁場 3=動磁場 4=磁位
//	EM_interval=1;//1		//電磁場計算を何ステップに一回行うか。通常は1に設定
	//解析領域
	XR=2.0;
	XL=-2.0;
	YU=2.0;
	YD=-2.0;
	/*
	XR=0.1;//0.01;		
	XL=-0.1;//-0.01;
	YU=0.1;//0.01;
	YD=-0.1;//-0.01;*/
	//円筒領域はこれだけ決める
/*	ZU=distancebp*20;
	ZD=-distancebp*20;
	RU=distancebp*10;*/
	
	//FRMcheck用	15/2/10
	ZU=2.5;//0.10; //0.2
	ZD=-2.5;//0.10; //0.2 				//液滴 -0.01 コイル:-0.15 るつぼ:-0.0002
	RU=2.0;//0.10;//0.1;				//解析領域が円筒形となるときのその半径

//流体の物性値
	MRE_density=1826;          //water:997.04  エタノール:798[kg/m3]
	Silicone_density=980;
	flag_modify_density=OFF;	//密度補償するかどうか
	nensei=1.00; //[Pa・s]//water:0.001 エタノール:0.001084 nensei 8.0
	sigma=0.07196;			//water:0.07196 エタノール:0.02361 SUS404:1.85 表面張力係数
	Cp=640/10;     			//water:4.2[kJ/(kg・K)] 鋼:800J SUS404:645J/(kgK)　
	k=28;       			//water:0.6[W/(m・K)] //熱伝導率
	latent_H=0;        		//water:334000J/kg  鋼209.3J　潜熱（エンタルピー）
	MP=273;		//融点[K]
	CTE=0.002;//2.1e-4;//2.1e-4;	//線膨張係数 SUS410:10.4e-6  水:2.1e-4

//弾性体の物性値
	/////MAGELAST/////
	E_m=20000.0;//1000; //200000		//ヤング率；500kPaだと圧力が異常な値になって止まる(dtのせいだと思われる)
	v_m=0.35;//0.49	//ポアソン比
	///////////////////
	/////ELASTIC//////
	E_e=18500.0;
	v_e=0.28;
	///////////////////
//粒子配置用
	fluidwidth=20; //30;//40//15[個]	//fluidwidth=20*2;
	distancebp=0.025;///0.001/2;//0.005; //distancebp=0.0125;[mm]
	wlength=2;
	height=0.0;//0.005;    

//解析領域
	maxX=2.0;//0.1;	//0.1/2;	//1
	minX=-2.0;//-0.1;	//-0.1/2;
	maxY=2.0;//0.1;	//0.1/2;	//0.4;
	minY=-2.0;//-0.1;	//-0.1/2;	//-0.6; //-1.0
	maxZ=2.5;//0.1;	//0.1/2;	//0.3;
	minZ=-2.5;//-0.1;	//-0.1/2;	//-0.6;  //indexの関係上、Z方向には余裕をもつこと。

	//FEMcheck用15/2/10
/*	maxX=0.2;	//0.1/2;	//
	minX=-0.2;	//-0.1/2;
	maxY=0.2;	//0.1/2;	//0.4;
	minY=-0.2;	//-0.1/2;	//-0.6; //-1.0
	maxZ=0.2;	//0.1/2;	//0.3;
	minZ=-0.2;	//-0.1/2;	//-0.6;  //indexの関係上、Z方向には余裕をもつこと。*/


//粒子法用パラメータ
	re=2.1;//2.1; //勾配・発散用の計算に使う
	re2=2.1; //laplacianに使う
	re3=3.0; //表面判定
	re4=3.1; //
	re_elastic=re;//弾性体計算用
	beta=0.8;
	dx=4;					//index用格子幅。最大粒子半径を含む整数にすること
    times=1;
	
//表面張力関係                                      
	surface_tension=0;      //0=OFF 1=形状
	smooth=OFF;				//スムージング　1=ON 0=OFF
	SFT=0;					//1=ON 0=OFF 表面張力(surface tension)の温度(T)依存性スイッチ
	smn=1;					//スムージングの回数
	smth_sumP=0;			//curvのスムージング回数 0ならOFF
    wall_direct=1;			//for ver2~5. BDWALLの法線ベクトル計算方法　0=無視 1=通常 2=水平方向 3=垂直方向
	non_slip=OFF;			//for ver.1 粘着条件 ONならBDWALL上に仮想的に液体が存在する設定。その際はBDWALLの配置に留意
	suf_modify=ON;			//曲率計算で得られた曲面に一致するように粒子位置を修正するか、しないか
	interpolate_curv=ON;	//曲率が計算されなかった粒子の曲率を、周辺粒子の値から補間するか、しないか

//圧力計算
	iteration_count=0;//1	//圧力計算を行う反復回数 通常は1に設定	0でFEM	
	solution=1;             //0=CG method 1=ICCG 2=ICCG2はメモリ節約
	B_of_P=1;				//圧力計算の解行列 0=教科書 1=速度発散 2=0+1 3=速度発散2 4=3+PND 5=(ni-nk)+(nk-n0)
	w_div=100;				//発散の重み
	HL_sw=OFF;				//ラプラシアンに対し高次の離散化を行う 1=ON 0=OFF
	dir_for_P=0;			//圧力計算の際、表面に非ゼロのDirichlet値を代入するか、しないか 0=OFF 1=表面張力 2=電磁力 3=両方
	initialP=ON;            //圧力計算時に初期値として現圧力を代入するかどうか　1=ON 0=OFF
	CGep=1e-3;				//圧力計算時における収束判定
	omp_P=OFF;				//圧力計算時のCG法でopenMPを使用するか 1=ON 0=OFF
	negativeP=ON;			//負圧　1=許可 0=不許可
	set_P=OFF;              //0=OFF 1=静水圧 2=関数
	Pgrad=3;				//圧力勾配計算法 1=教科書 2=表面のみ法線 3=2+minP=0のとき表面反発 4:WLSM 5=CMPS 6=Pj+Pi 7=MPS-AS 
	minP=ON;				//最少圧力スイッチ 1=ON 0=OFF
	ave_Pdirect=OFF;		//圧力勾配用法線ベクトルの平均化 1=ON 0=OFF
	Pgrad_order=2;			//圧力勾配ver.4における精度 1=線形 2=2次
	artP_sw=OFF;			//人工圧力計算フラグ 0=OFF 1=artP等方 2=Monaghan
	artP=100;				//人工圧力値 通常は0に設定
	Pgrad_times=20;			//圧力勾配をプロットする際の、表面張力に対する倍率 通常は1に設定
	P_AVS=10000;				//microAVS用の圧力ファイルを出力するstep間隔。0ならOFF
	
//温度場関係
	T_field=0;              //温度場解析　　1=ON 0=OFF
	insulate=1;             //壁との断熱状態 0=断熱　1=非断熱
	T_laplacian=0;          //温度のラプラシアン　0=教科書　1=λ[i] 2=発散・勾配
	wall_density=7850;		//壁の密度[kg/m^3]
	wall_Cp=460;			//壁の比熱[J/(kg・K)]
	wall_k=24;				//壁の熱伝導率[W/(m・K)]
	roomT=313;				//室温(20℃) int型でよい
	air_cool=ON;			//空気との熱伝達を考慮するかしないか 1=ON 0=OFF
	T_expansion=OFF;		//熱膨張計算 1=ON 0=OFF
	buoyant=OFF;			//浮力(密度のBoussinesq近似)  1=ON 0=OFF
	TplotZ=0.004;			//3D解析において、XY平面の温度を出力するときのＺ座標
	T_AVS=20;				//microAVS用の温度ファイルを出力するstep間隔。0ならOFF

//メッシュ
	mesh_input=2;//2;				//meshの作成方法 0:自分 1:Magnet 2:TetGen
	remesh_sw=OFF; //OFF;2012/02/23	//ON:remesh領域を考慮し、そこだけremesh OFF:FULL領域を毎回分割
	boxalpha=2.0;				//スーパーボックスの大きさ 1～10の実数にすること。通常は2に設定
	fine=ON;//ON2012/07/08		//節点の自動追加を行うかどうか 0:OFF 1:辺ベース 2:重心ベース(現在は消去)
	co_fine=3.0;				//節点自動追加における係数.辺ベースの場合は最少辺と最大辺の比率のしきい値
	add_points=50000;//50000	//自動節点追加数
	poly_memory=100000;//50000	//poly3D()で確保するメモリの数
	air_layer=2;				//物体の周囲に空気層を生成する層数(0なら生成しない)
	layer_depth=0.5;//0.5(2012/03/03);		//空気層の幅。初期粒子間距離の何倍か
	mesh_output_interval=1;//1;	//メッシュ情報を、有限要素法のステップに対して、何ステップに一度出力するか。
	FEMCGep=5.0e-5;			//1.3e-4PICCG 元々5.0e-6	//本プログラムでは現在不使用15/5/24
	MRTRep=6.8e-4;		//MRTR(&ICMRTR)法の収束判定	//本プログラムでは現在不使用15/5/24
	FEM_calc_type=2;	//15/5/24	//3;		//0=OFF 1=電場 2=磁場 3=磁場(渦電流) 4=磁位 5=非線形静磁場
	ele_type=1;				//(mesher=0の場合) 要素ﾀｲﾌﾟ 0:節点要素 1:辺要素

//磁場計算
	J_input_way=0;			//電流密度入手方法 0:OFF 1:自分 2:ソフト
	J0=0.0;//100;//180000000;		//強制電流密度[A/m2]
	I0=450.0;//200; mA?
	RP=2.0;//2//1.5;//1.28;	//比透磁率
	ele_conduc=1e7;			//電気伝導率
	Hz=10;			//周波数
	div_Hz=4;				//１周期の分割数(解析精度) 4の倍数がいい 40?
	jheat=0;				//渦電流による発熱を考慮するか　0=OFF 1=ON
	m_force=1;//1			//電磁力計算方式 1=節点力法 2=kelvin 3=積分面 4=divT(マクスウェルの応力テンソル) 5=VAG 6=積分つき節点力法 7=MC
	NLBHI=0;				//体積力において、要素Ｂから要素Ｈを求める際に非線形性を考慮するか、しないか(non linier B H inverter)
	NLMH=OFF;				//Ｍの算出に非線形性を考慮するか、しないか
	magnet_H=5.0*distancebp;			//永久磁石の高さ0.005
	magnet_r=3.0*distancebp;//0.01	//永久磁石の半径0.005 　　　　　　　　　　　//J_input_way=2:半径ではなく直径、J_put_way=0:半径　と思われる。
	magnet_Z=-(6*sqrt(2.0)+3.0)*distancebp;//-8*distancebp; //-45*0.0005-0.005; //-(fluidwidth)*distancebp-0.01; //-0.0125*0.8		//永久磁石の中心のZ座標 42*0.0005 //-magnet_H/2-(15*distancebp+0.002) モデル5:-0.035
	magnet_angle=0.0;			//永久磁石の着磁方向 0なら+Z方向となる。そこから角度をつけたいなら、その角度[deg]を入力する
	magnet_B=1.25;//0.145;//1.20;			//永久磁石の強さ[T] Avector3D()で指定
	magnetic_layers=1;	//永久磁石周辺の空気層の数1層はすでにある 1+
	uniform_B=0.00; //0.01;	//一様磁場の大きさ[T]?
	B_times=1;//0.1			//ファイル出力する際の、磁束密度の倍率
	plot_B_type=1;			//磁束密度出力タイプ 1=ベクトル 2=スカラー 0=OFF

	//2012-11-02 0:38 辺要素を分岐させる！！！これがないとディリクレ値が全てゼロ！！
	//コンストラクタで値が未定義で実行されていない関数未だあるのでは？？
	uniform_B_sw=OFF;		//解析領域中に一様磁場を発生させるか否か 0=OFF 1=ON　

//BEM関係
	BEM_elm_type=CONSTANT;	//要素タイプ 0:一定要素 1:線形要素
	BEM_M_solver=2;			//境界要素法における行列解法 0:ガウスの消去法 1:BiCGStab法 2:BiCGStab2法
	gauss_N=7;				//Gauss積分の評価点数 3,4,7のどれか
	FEM_elm_type=1;			//要素タイプ 0:節点要素 1:辺要素
	FEM_smn=1;				//電磁力スムージング回数　0ならOFF マイナスなら表面のみ プラスなら内部も。
	max_DN=25000;           //Dirichlet型境界条件をとる最大節点(辺)数
//	FEMCG=1;				//FEMにおける行列解法 0:CG 1:ICCG 2:並列ICCG	//なぜコメントアウトされているのか15/5/24
	CGaccl=1.3;				//ICCG法における加速ファクタ　1のときファクタOFF
	EMCGep=1.0e-12;			//電磁場のICCGの収束判定1e-5  アクチュエータ2.5e-5　ICCG
	FEMtimes=5;				//電磁力をプロットする際の、表面張力に対する倍率 通常は1に設定
	legend_F=2e-7;			//F.datの凡例に出力する力[N]
	tree_sw=ON;				//BEMでtree法を使うか、使わないか
	tree_mesh_size=4;		//tree法で使用する最大レベルのセルのサイズ。leの何倍か
	p_max=3;				//tree法の無限級数の最大項数。0から6の数字に対応
	plot_F_type=1;			//F.datの出力タイプ  0=[N]表示 1=[Pa]表示
	legend_F_Pa=300;		//F.datの凡例に出力する単位[Pa]

//FEMver2関係
	surface_depth=1;		//表面の厚みの半分(これにleをかけたものが本当の厚み(の半分))
	
//電界計算
	V=40000;				//電界計算用電圧　もうひとつは０とおいている。
	V_step=100;//1500;
	r_perm=80;				//比誘電率
	V_con=0;				//電圧条件　0:パルス　1:リニア　2:時定数
	initial_V=2000;			//電圧条件リニアのときの初期電圧
	E_times=1e-10;			//電界出力倍率 0なら出力しない
	eleforce=4;				//静電力計算方法 1:節点力法 2:積分面 3:表面力 4:divT
	charge=0;				//電荷考慮 0=OFF 1=電荷密度 2=クーロン力
	plot_E_type=1;			//電界出力タイプ 1=ベクトル 2=スカラー 0=OFF
	
//各種スイッチ　どのような計算を考慮するか、しないか
	g=-9.8;					//-9.8
	restart=OFF;				//1=ON 0=OFF
	autosave=100;			//オートセーブ間隔。無効にしたいときは大きな数字を代入しておく
	courant=1.0;			//Courant数条件　0ならＯＦＦ
	modify_position=OFF;		//粒子間距離がleより小さい場合これを修正するか、しないか
	vis_calc_type=1;		//粘性項計算手法　0=POSITIVE=陽解法 1=NEGATIVE=陰解法
	wall_adheision=2;		//0=フリースリップ 1=ノンスリップ  2=2*ノンスリップ //初期では：2
	laplacian=0;//0;        //0=教科書　1=λ[i] 2=発散・勾配 //初期では2
	vis_solver=0;			//粘性項を陰解析で解く際の行列ソルバー 0:CG 1:ICCG
	initial_u=OFF;			//粘性項を陰的に解く際に、初期値として現在速度を入力するか、しないか
	temporary_r=OFF;		//陽解析後に仮の位置を計算するかしないか。 1=ON 0=OFF 通常はONに設定
	fix_center=0;			//1=ON 0=OFF
	freeon=1;				//粒子依存関係関数　1:並列化可能 2:並列不可 4:GPU
	freeon3sw=1;			//freeon3を計算するかしないか 1=ON 0=OFF
	surface_judge2=OFF;		//surface_judge2()を使用するかしないか  1=ON 0=OFF
	move_prtcl=OFF;			//移動粒子を考慮するかしないか 1=ON 0=OFF
	move_u_dirct=-3;//-3;//2;			//移動粒子を移動する方向　現在は±X方向=±1,±Y方向=±2,±Z方向=±3
	move_speed=1.5e-3;//1.5e-3;//1e-3;//8.333*1e-3;//12.5*1e-3;	//移動粒子の移動速度[m/s]
	check_something=OFF;		//check_something()を実行するかしないか 1=ON 0=OFF

	
//速度プロット変数
	speed_plot_particle=2;	//速度をプロットする粒子の種類 1=すべて 2=fluid 3=壁
	speedtimes=5.0e-3;		//速度プロット時の、座標に対する速度の倍率
	speed_face=0;			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ
	speed_face_p=0.0;		//3D解析時のspeed.datの出力面の座標
	ax_sym_modify=OFF;		//3D時のspeed.datに関して、軸対称による出力修正を行うか否か　1=ON 0=OF
	flat_speed_plot=OFF;	//水平方向の速度(XY面)をプロットするかしないか1=ON 0=OFF
	flat_speed_p=0.004;		//flat_speed.datの出力面の座標
	relative_speed=OFF;		//重心に対する相対速度を出力するかしないか 1=ON 0=OFF
	speed_AVS=OFF;			//microAVSによる3D速度分布出力するかしないか 1=ON 0=OFF
	legend_speed=0.1;		//speed.datの凡例に出力する速度[m/s]
	set_zero_speed=OFF;		//restart時に速度をゼロセットするかしないか  1=ON 0=OFF
	
//圧力による加速度のプロット変数
	pressure_face=1;			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ
	pressure_face_p=0.0;		//3D解析時のspeed.datの出力面の座標
	pressure_times=1.0e-6;
	legend_pressure=0.1;

//GPU関係
	M_form=ELL;				//CG法における係数行列の格納方法 CSR_scl,CSR_vec,ELLの3つをサポート
	MAX_thread=512;			//ひとつのSMあたりの最大スレッド数　ふつうは512

//ファイル出力変数 
//	interval=50; //10		//2以上の整数にすること
	AVS=0;                  //0:普通　1:圧力　2:温度 3:壁非表示 4：表面のみ 5:壁 6:特定 
	P_size_AVS=1;			//出力する粒子のサイズ(leの何倍か)
	maxT=1000;
	minT=293;

//AVS関連出力ファイル
	F_interval=2;			//F.datのログを何ステップ毎に出力するか
	avs_eforce_interval=0;			//AVS電磁力ファイルを何ステップ毎に出力するか
	avs_mesh1_interval=200;			//AVS電位ファイル(断面)を何ステップ毎に出力するか
	avs_mesh2_interval=200;			//AVSメッシュファイ有ル(断面)を何ステップ毎に出力するか
	avs_mesh3_interval=200;			//AVSソリッドモデルを何ステップ毎に出力するか

	times_Pa=3E-7;	//gnuplot用ファイルのベクトル長さの倍率 面積力[Pa]

	max_pressure=10000;
	min_pressure=0;


//超弾性計算 
	flag_ELAST=OFF;
	flag_HYPER=ON;
	hyper_density=1000;          //water:997.04  エタノール:798[kg/m3]
	c10=30000;//30000;
	c01=20000;//20000;
	flag_wall=ON;
	h_dis=1.9*distancebp;
	h_vis=1.0;
	flag_vis=OFF;
	nr_time=1000;	//15/2/8
}


//質量: 最初に分かるのは粒子一個あたりの体積だけ。あとは初期配置に拠って密度が補正されるのでそれをかける
//充填率: 正方 68%　最密 74%
double mpsconfig::get_particle_mass()
{
	double mass;
	if(model_set_way==0)//正方格子のとき充填率＝PI/6=52.3%
	{
		if(dimension==2){
			mass=MRE_density*distancebp*distancebp;
		}else{
			mass=(0.523598775)*MRE_density*distancebp*distancebp*distancebp;
		}
	}
	else if(model_set_way==1)//最密格子のとき(充填率はsqrt(2)*PI/6=0.74048049)
	{
		if(dimension==2){
			mass=MRE_density*sqrt(3.0)/4*distancebp*distancebp;
		}else{
			//mass=density*sqrt(2.0)/12*distancebp*distancebp*distancebp;
			//mass=density*distancebp*distancebp*distancebp/sqrt(2.0);
			//(4*PI/3)/2^3=4.1887902/8=0.523598776
			//(sqrt(2)*PI/6)*(4*PI/3)/2^3=0.387714678
		//	mass=(0.74048049)*distancebp*distancebp*distancebp*density;
			mass=distancebp*distancebp*distancebp*MRE_density;
		}
	}
	else{
		cout<<"モデルの積み方が不定です 質量を計算できません"<<endl;
		exit(1);
	}

	return mass;
}