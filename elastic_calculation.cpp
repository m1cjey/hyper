#include "stdafx.h"
#include "Rigidbody.h"
#include "Micro_AVS.h"
//#include<boost\array.hpp>
//#include <boost\shared_ptr.hpp>
using namespace std;

//粘弾性体計算
//エネルギーは関数を作って計算させる
void calc_elastic(vector<mpselastic> &PART, elastic &ELAST, int t, double **F)
{
	cout<<"粒子移動計算開始"<<endl;

	int dimension=ELAST.get_dimension();

/*	if(t==1){
		rigid_calc(PART, ELAST);
	}*/

	if(ELAST.get_modify_density()==ON) calc_modified_density(PART, ELAST); //密度の修正

	if(ELAST.get_nonlinear_elastic_flag()==false)
	{
		calc_pre_velocity_and_position(PART, ELAST, F);//前回のステップの内力（加速度）を用いて仮の位置・速度を計算

		#pragma omp parallel for
		//剛体の計算
//		calc_rigid(PART,rigids);

		for(int i=0;i<PART.size();i++) PART[i].reset_particle_acceleration(); //加速度のリセット

		calc_quaternion_using_vectorSTL(ELAST, PART);
		calc_angular_velocity_using_vectorSTL(ELAST, PART);

		//初期位置ベクトルの回転（現在位置ベクトルは回転させなくてよい）
		calc_r0_ij(PART); 

		//
		calc_accel_for_3D(PART, ELAST); 

		//圧力と接触反力による加速度を計算
		contact_judge(PART, ELAST);

		//各加速度の合計
		modify_acceleration(PART, ELAST, F);

		//仮の位置・速度から求めた内力（加速度）でu, rを修正
		calc_post_velocity_and_position(PART, ELAST,t);




/*		for(int i=0;i<PART.size();i++){
		cout<<"i="<<i<<", rx="<<PART[i].r[A_X]<<", ry="<<PART[i].r[A_Y]<<", rz="<<PART[i].r[A_Z]<<endl;//
	}//*/
		//ハミルトニアンのファイル出力
		calc_hamiltonian(ELAST, t);
	}
	else
	{
		calc_pre_velocity_and_position(PART, ELAST, F);//前回のステップの内力（加速度）を用いて仮の位置・速度を計算

		#pragma omp parallel for
		for(int i=0;i<PART.size();i++) PART[i].reset_particle_acceleration(); //加速度のリセット

		calc_quaternion_using_vectorSTL(ELAST, PART);
		calc_angular_velocity_using_vectorSTL(ELAST, PART);

		//初期位置ベクトルの回転
		calc_r0_ij(PART); 


		calc_nonlinear_accel_for_3D_test(PART, ELAST);

		//圧力と接触反力による加速度を計算
		contact_judge(PART, ELAST);	

		//各加速度の合計
		modify_acceleration(PART, ELAST, F);

		//仮の位置・速度から求めた内力（加速度）でu, rを修正
		calc_post_velocity_and_position(PART, ELAST,t);

		//剛体の位置修正
//		calc_rigid(PART,rigids);

		//ハミルトニアンのファイル出力
		calc_hamiltonian(ELAST, t);
	}
}

void calc_accel_for_3D(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double mag_lambda=ELAST.get_mag_lambda();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double elas_lambda=ELAST.get_elas_lambda();
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //これは既にre*le
	double vis=ELAST.get_nensei(); //粘度（動粘度ではない）
	double g=ELAST.get_g();

	double KE=0.0;	//運動エネルギー
	double EE1=0.0;	//ひずみエネルギー
	double EE2=0.0;	//体積ひずみによる
	double PE=0.0;	//ポテンシャル
	double ground=ELAST.get_ground_position(); //床のｚ座標の取得

	cout<<"ground: "<<ground<<endl;
	double mindis=100;
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
		//WALLは除外してはいけない

		double EE1_temp=0.0;

		//運動エネルギーの更新
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST){
			for(int D=0;D<3;D++) KE+=PART[i].u_temp[D]*PART[i].u_temp[D];
		}

		//位置ポテンシャルの更新
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST){
			PE+=(PART[i].r[2]-ground);
		}

		int neighboursN0=static_cast<int>(PART[i].get_initial_neighboursID().size());

		//////////////////初期配置で周辺にあるものだけを追跡する・・・/////////////////////////
		for(int k=0;k<neighboursN0;k++)
		{
			int j=PART[i].get_initial_neighboursID()[k];

			//iからみた初期配置での相対座標
			double r_ij_Init[3], r_ji_Init[3]; 
			for(int D=0;D<3;D++){
				r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
				r_ji_Init[D]=-r_ij_Init[D];
			}

			//現在配置での相対座標
			double r_ij[3], r_ji[3];
			for(int D=0;D<3;D++)
			{
				r_ij[D]=PART[j].r_temp[D]-PART[i].r_temp[D];
				r_ji[D]=-r_ij[D];
			}

			double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離
			double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離
			if(dis<=mindis) mindis=dis;
			double r_ij_zero[3], r_ji_zero[3];

			//calc_r0_ij(PART)で"現在位置のクォータニオンを使って"回転させた初期位置ベクトルを取得
			for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
			for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

			double w=kernel(re, dis0);

			//回転後の初期配置相対位置の単位ベクトル
			double n_ij_zero[3], n_ji_zero[3];
			for(int D=0;D<3;D++) 
			{
				n_ij_zero[D]=r_ij_zero[D]/dis0;
				n_ji_zero[D]=r_ji_zero[D]/dis0;
			}

			//現在配置相対位置単位ベクトル
			double n_ij_current[3], n_ji_current[3];
			for(int D=0;D<3;D++) 
			{
				n_ij_current[D]=r_ij[D]/dis;
				n_ji_current[D]=r_ji[D]/dis;
			}

			//変位ベクトル
			double U_ij[3], U_ji[3]; 
			for(int D=0;D<3;D++)
			{
				U_ij[D]=r_ij[D]-r_ij_zero[D];
				U_ji[D]=r_ji[D]-r_ji_zero[D];
			}

			//ひずみベクトル
			double E_ij[3], E_ji[3];
			for(int D=0;D<3;D++)
			{
				E_ij[D]=U_ij[D]/dis0;
				E_ji[D]=U_ji[D]/dis0;
			}
				
			//圧力計算の準備・・・e_vol_iのΣ[]（変位の発散）を計算
			//垂直ひずみの和を取る（応力テンソルの対角成分の和に重み付け）
			double volumetric_strain=0.0;
			for(int D=0;D<3;D++) volumetric_strain+=E_ij[D]*n_ij_current[D];
//			for(int D=0;D<3;D++) PART[i].P+=E_ij[D]*n_ij_zero[D]*w; 
			volumetric_strain*=w;
			PART[i].P+=volumetric_strain; //＝ではないので注意！

			//ひずみエネルギー
			for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;

			//応力ベクトル（取り扱いを変更・・・重み関数は粒子iで常に加重平均する。そうしないと圧力の重み関数との整合性が取れない）
			double sigma_ij[3], sigma_ji[3];
			for(int D=0;D<3;D++)
			{
				if(PART[i].type==MAGELAST){
				sigma_ij[D]=mag_shear_modulus*w*(E_ij[D]-E_ji[D])/dis0;
				}
				else if(PART[i].type==ELASTIC){
					sigma_ij[D]=elas_shear_modulus*w*(E_ij[D]-E_ji[D])/dis0;
				}
			}
			PART[i].add_stress_accel(sigma_ij);

//			if(fabs((dis-dis0)/dis0)>1.04) w=0;//破壊条件

			//ひずみ速度の計算
			double strain_vi[3], strain_vj[3];

			//iからみたひずみ速度
			strain_vi[0]=((PART[j].u_temp[0]-PART[i].u_temp[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
			strain_vi[1]=((PART[j].u_temp[1]-PART[i].u_temp[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y軸周りのひずみ速度
			strain_vi[2]=((PART[j].u_temp[2]-PART[i].u_temp[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z軸周りのひずみ速度
			//anglar_u1..3は配列に書き換え

			//jからみたひずみ速度
			strain_vj[0]=((PART[i].u_temp[0]-PART[j].u_temp[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
			strain_vj[1]=((PART[i].u_temp[1]-PART[j].u_temp[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y軸周りのひずみ速度
			strain_vj[2]=((PART[i].u_temp[2]-PART[j].u_temp[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z軸周りのひずみ速度
			
			double sigma_v_ij[3], sigma_v_ji[3];
			for(int D=0;D<3;D++)
			{
				sigma_v_ij[D]=vis*w*(strain_vi[D]-strain_vj[D])/dis;
//				sigma_v_ji[D]=2*vis*w*strain_vj[D]/PART[j].PND/dis;
//				sigma_v_ij[D]-=sigma_v_ji[D];
			}

			PART[i].add_stress_visco_accel(sigma_v_ij);
		}//for(int k=0;k<neighboursN0;k++)ループ終了 初期粒子配置での周囲にある粒子
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//		PART[i].P*=dimension/PART[i].PND;//体積ひずみが求められた（これはエネルギー計算する前に求めておくこと）
		PART[i].P*=dimension/PART[i].PND0;//体積ひずみが求められた（これはエネルギー計算する前に求めておくこと）

//		double coef=dimension/PART[i].get_density()/PART[i].PND;
		double coef=dimension/PART[i].get_density()/PART[i].PND0;

		PART[i].mul_stress_accel(coef);
		PART[i].mul_stress_visco_accel(coef);

		//ひずみエネルギー
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
//			EE1+=EE1_temp/PART[i].PND/PART[i].get_density();//第二項のΣ 密度はそれぞれ違う
	//		if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
			if(PART[i].type==MAGELAST){
			EE1+=(EE1_temp/PART[i].PND0/PART[i].get_density())*mag_shear_modulus*mass*dimension;//第二項のΣ 密度はそれぞれ違う
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*mag_lambda;//第三項のΣ 密度はそれぞれ違う
			}
			else if(PART[i].type==ELASTIC){
			EE1+=(EE1_temp/PART[i].PND0/PART[i].get_density())*elas_shear_modulus*mass*dimension;//第二項のΣ 密度はそれぞれ違う
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*elas_lambda;//第三項のΣ 密度はそれぞれ違う
			}
		}

//		if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
		if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
		PART[i].P*=-mag_lambda;//体積ひずみによる圧力が求められた
		}
		else if(PART[i].type==ELASTIC){
		PART[i].P*=-elas_lambda;
		}
		else if(PART[i].type==WALL){
		PART[i].P*=-elas_lambda;
		}
//		これで置き換えずにWALL込で計算して上のブロックで関数で有無を判定して符号を逆転させる方が確実では・・・
		
		//粒子数密度が増加した場合はPを置き換える
		if(PART[i].PND>PART[i].PND0)
		{
			if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2)
			{
				double contact_pressure=mag_lambda*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
				if(contact_pressure>PART[i].P)
				{
					PART[i].P=contact_pressure;//密度が非線形なのにここで粒子数密度が線形として良い？→OK（n0とρは比例）
				}
				PART[i].contact=true; //フラグを立てておく
			}
			else if(PART[i].type==ELASTIC){
			double contact_pressure=elas_lambda*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
			if(contact_pressure>PART[i].P)
			{
				PART[i].P=contact_pressure;//密度が非線形なのにここで粒子数密度が線形として良い？→OK（n0とρは比例）
			}
			PART[i].contact=true; //フラグを立てておく
			}
			else if(PART[i].type==WALL){
			double contact_pressure=elas_lambda*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
			if(contact_pressure>PART[i].P)
			{
				PART[i].P=contact_pressure;//密度が非線形なのにここで粒子数密度が線形として良い？→OK（n0とρは比例）
			}
			PART[i].contact=true; //フラグを立てておく
			}
		}


	}//for(i=0;i<PART.size();i++)ループ終了
	cout<<"最小粒子間距離"<<mindis<<endl;
	KE*=0.5*mass;
	PE*=-g*mass;
	
//	EE1*=shear_modulus*mass*dimension;
//	EE2*=0.5*mass*lambda;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//全系の弾性エネルギー
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);
}

void calc_nonlinear_accel_for_3D_test(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //これは既にre*le
	double vis=ELAST.get_nensei(); //粘度（動粘度ではない）
	double g=ELAST.get_g();

//	vector<double> lambda;

	double KE=0.0;	//運動エネルギー
	double EE1=0.0;	//ひずみエネルギー
	double EE2=0.0;	//体積ひずみによる
	double PE=0.0;	//ポテンシャル
	double ground=ELAST.get_ground_position(); //床のｚ座標の取得


//	cout<<"ground: "<<ground<<endl;
	
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2)//WALLは弾性変形しないので内力計算に含めない
		{

			double EE1_temp=0.0;
			double EE2_temp=0.0;
			double volumetric_strain=0.0;//体積ひずみ（iごとに異なる）

			//運動エネルギーの更新
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST){
				for(int D=0;D<3;D++) KE+=PART[i].u[D]*PART[i].u[D];
			}

			//位置ポテンシャルの更新
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST){
				PE+=(PART[i].r[2]-ground);
			}

			int neighboursN0=static_cast<int>(PART[i].get_initial_neighboursID().size());

			//着目した粒子iと周辺粒子jとの「関係」を計算
			//初期配置で周辺にあるものだけを追跡する・・・
			for(int k=0;k<neighboursN0;k++) //pressureと合わせないと狂う？
			{
				int j=PART[i].get_initial_neighboursID()[k];

				//iからみた初期配置での相対座標
				double r_ij_Init[3], r_ji_Init[3]; 
				for(int D=0;D<3;D++){
					r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
					r_ji_Init[D]=-r_ij_Init[D];
				}

				//現在配置での相対座標
				double r_ij[3], r_ji[3];
				for(int D=0;D<3;D++)
				{
					r_ij[D]=PART[j].r_temp[D]-PART[i].r_temp[D];
					r_ji[D]=-r_ij[D];
				}

				double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離
				double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離
			
				double r_ij_zero[3], r_ji_zero[3];

				//calc_r0_ij(PART)で"現在位置のクォータニオンを使って"回転させた初期位置ベクトルを取得
				for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
				for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

				double w=kernel(re, dis0);

				//回転後の初期配置相対位置の単位ベクトル
				double n_ij_zero[3], n_ji_zero[3];
				for(int D=0;D<3;D++) 
				{
					n_ij_zero[D]=r_ij_zero[D]/dis0;
					n_ji_zero[D]=r_ji_zero[D]/dis0;
				}

				//現在配置相対位置単位ベクトル
				double n_ij_current[3], n_ji_current[3];
				for(int D=0;D<3;D++) 
				{
					n_ij_current[D]=r_ij[D]/dis;
					n_ji_current[D]=r_ji[D]/dis;
				}

				//変位ベクトル
				double U_ij[3], U_ji[3]; 
				for(int D=0;D<3;D++)
				{
					U_ij[D]=r_ij[D]-r_ij_zero[D];
					U_ji[D]=r_ji[D]-r_ji_zero[D];
				}

				//ひずみベクトル
				double E_ij[3], E_ji[3];
				for(int D=0;D<3;D++)
				{
					E_ij[D]=U_ij[D]/dis0;
					E_ji[D]=U_ji[D]/dis0;
				}

				//ヤング率とポアソン比の計算
				double En_ij[3];
				double volumetric_temp=0.0;
				for(int D=0;D<3;D++) volumetric_temp+=E_ij[D]*n_ij_current[D];	//E_ijとn_ij_currentの内積
				for(int D=0;D<3;D++) En_ij[D]=volumetric_temp*n_ij_current[D];	//ひずみベクトルをr_ij方向へ正射影
				double ns=sqrt(En_ij[0]*En_ij[0]+En_ij[1]*En_ij[1]+En_ij[2]*En_ij[2]); //縦ひずみ
				ns*=100;	//[%]
				//シリコーンの場合
				//youngs_modulus=-6.11619E-11*x^9+1.66522E-8*x^8-1.93253E-6*x^7+1.2459E-4*x^6-4.87141E-3*x^5+0.118275*x^4-1.75252*x^3+14.9237*x^2-64.3273*x+123.411;
				double elas_youngs_modulus=(((((((-0.000004092500180)*ns+0.000860882692091)*ns-0.073155040824811)*ns+3.215486271963979)*ns-77.712461766914529)*ns+1012.872595336284500)*ns-6495.847221601891800)*ns+18812.474059104068000;
				//-1.09176E-12*x^10+1.56272E-10*x^9-9.70758E-9*x^8+3.42745E-7*x^7-7.56378E-6*x^6+0.00010807*x^5+-0.00100264*x^4+5.93234E-3*x^3-0.0219534*x^2+0.0537726*x^1+0.061154
				double elas_poisson_ratio=(((((((+0.000000000918578)*ns-0.000000090882723)*ns+0.000003665373271)*ns-0.000078101785064)*ns+0.000960394839085)*ns-0.007112271340844)*ns+0.033194941391016)*ns+0.069899017991200
;

				//MREの場合
				//youngs_modulus=-4.91157E-8*x^7+9.42896E-6*x^6-7.31665E-4*x^5+2.95145E-2*x^4-0.663258*x^3+8.25586*x^2-53.0674*x+190.48
				double mag_youngs_modulus=(((((((-0.000002526899825)*ns+0.000526467419950)*ns-0.044510412964222)*ns+1.962528096212004)*ns-48.344643329754824)*ns+664.648197882870590)*ns-4894.930152042784400)*ns+23399.990870638816000;
				//-9.6723E-12*x^10+1.22758E-9*x^9-6.73685E-8*x^8+2.09375E-7*x^7-4.05552E-5*x^6+5.08152E-4*x^5-0.00414923*x^4+0.0218461*x^3-0.072922*x^2+0.151892*x^1+0.17713
				double mag_poisson_ratio=(((((((+0.000000006114133)*ns-0.000000557177269)*ns+0.000020762913954)*ns-0.000408515956157)*ns+0.004578221540415)*ns-0.029448326315701)*ns+0.103023742143433)*ns+0.192042424941781;

				if(i==1511 && j==1513){
					ofstream pr("poason.dat", ios::app);
					cout<<"youngs_modulus: "<<mag_youngs_modulus<<", poisson_ratio: "<<mag_poisson_ratio<<endl;
					pr<<"youngs_modulus="<<mag_youngs_modulus<<", poisson_ratio="<<mag_poisson_ratio<<", strain="<<ns<<endl;
					pr.close();
				}
				//ラメ定数mu
				double elas_shear_modulus=elas_youngs_modulus/(2.0*(1.0+elas_poisson_ratio));
				double mag_shear_modulus=mag_youngs_modulus/(2.0*(1.0+mag_poisson_ratio));
				//ラメ定数lambda
				double elas_lambda_temp=(elas_poisson_ratio*elas_youngs_modulus)/((1.0+elas_poisson_ratio)*(1.0-2.0*elas_poisson_ratio));
				double mag_lambda_temp=(mag_poisson_ratio*mag_youngs_modulus)/((1.0+mag_poisson_ratio)*(1.0-2.0*mag_poisson_ratio));
				
	//			cout<<"shear_modulus: "<<shear_modulus<<", lambda: "<<lambda<<endl;

				//圧力計算の準備・・・e_vol_iのΣ[]（変位の発散）を計算
				//垂直ひずみの和を取る（応力テンソルの対角成分の和に重み付け）
				if(PART[i].type==ELASTIC){
				volumetric_temp*=-elas_lambda_temp*w;
				PART[i].P+=volumetric_temp;
				PART[i].ave_lambda+=elas_lambda_temp;
				//ひずみエネルギー(重みはいらない)
				for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;
				EE2_temp+=elas_lambda_temp*volumetric_temp*volumetric_temp;

				//応力ベクトル（取り扱いを変更・・・重み関数は粒子iで常に加重平均する。else, 圧力の重み関数との整合性が取れない）
				double sigma_ij[3], sigma_ji[3];
				for(int D=0;D<3;D++) sigma_ij[D]=elas_shear_modulus*w*(E_ij[D]-E_ji[D])/dis0;
				PART[i].add_stress_accel(sigma_ij);
				}
				///////////////////////////////////////////////////

				else if(PART[i].type==MAGELAST|| PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
				volumetric_temp*=-mag_lambda_temp*w;
				PART[i].P+=volumetric_temp;
				PART[i].ave_lambda+=mag_lambda_temp;
				//ひずみエネルギー(重みはいらない)
				for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;
		//		EE2_temp+=mag_lambda_temp*volumetric_temp*volumetric_temp;
				EE2_temp+=mag_lambda_temp*volumetric_temp*volumetric_temp/w*w; //重みをつけなければ

				//応力ベクトル（取り扱いを変更・・・重み関数は粒子iで常に加重平均する。else, 圧力の重み関数との整合性が取れない）
				double sigma_ij[3], sigma_ji[3];
				for(int D=0;D<3;D++) sigma_ij[D]=mag_shear_modulus*w*(E_ij[D]-E_ji[D])/dis0;
				PART[i].add_stress_accel(sigma_ij);
				}
				

	//			if(fabs((dis-dis0)/dis0)>1.04) w=0;//破壊条件

				//ひずみ速度の計算
				double strain_vi[3], strain_vj[3];

				//iからみたひずみ速度
				strain_vi[0]=((PART[j].u[0]-PART[i].u[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
				strain_vi[1]=((PART[j].u[1]-PART[i].u[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y軸周りのひずみ速度
				strain_vi[2]=((PART[j].u[2]-PART[i].u[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z軸周りのひずみ速度
				//anglar_u1..3は配列に書き換え

				//jからみたひずみ速度
				strain_vj[0]=((PART[i].u[0]-PART[j].u[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
				strain_vj[1]=((PART[i].u[1]-PART[j].u[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y軸周りのひずみ速度
				strain_vj[2]=((PART[i].u[2]-PART[j].u[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z軸周りのひずみ速度
			
				//粘性応力ベクトル
				double sigma_v_ij[3], sigma_v_ji[3];
				for(int D=0;D<3;D++)
				{
					sigma_v_ij[D]=vis*w*(strain_vi[D]-strain_vj[D])/dis;
	//				sigma_v_ji[D]=vis*w*strain_vj[D]/PART[j].PND/dis;
	//				sigma_v_ij[D]-=sigma_v_ji[D];
				}

				PART[i].add_stress_visco_accel(sigma_v_ij);
			}//for(int k=0;k<neighboursN0;k++)ループ終了

			//λの平均
			PART[i].ave_lambda/=neighboursN0;	//あくまで初期粒子間のλ
			//圧力
			PART[i].P*=dimension/PART[i].PND0;

	//		double coef=dimension/PART[i].get_density()/PART[i].PND;
			double coef=dimension/PART[i].get_density()/PART[i].PND0;

			PART[i].mul_stress_accel(coef);
			PART[i].mul_stress_visco_accel(coef);

			//ひずみエネルギー
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST || PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
				if(PART[i].type==MAGELAST){
				EE1+=(EE1_temp/PART[i].PND0/PART[i].get_density())*mag_shear_modulus*mass*dimension;//第二項のΣ 密度はそれぞれ違う
				}
				else if(PART[i].type==ELASTIC){
				EE1+=(EE1_temp/PART[i].PND0/PART[i].get_density())*elas_shear_modulus*mass*dimension;//第二項のΣ 密度はそれぞれ違う
				}
				EE2+=EE2_temp/PART[i].get_density();//第三項のΣ 密度はそれぞれ違う
			}

			//粒子数密度が増加した場合はPを置き換える
	

			if(PART[i].PND>PART[i].PND0)
		{
			if(PART[i].type==MAGELAST || PART[i].type==TERMINAL1||PART[i].type==TERMINAL2){
			double contact_pressure=ELAST.get_mag_lambda()*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
			if(contact_pressure>PART[i].P)
			{
				PART[i].P=contact_pressure;//密度が非線形なのにここで粒子数密度が線形として良い？→OK（n0とρは比例）
			}
			PART[i].contact=true; //フラグを立てておく
			}
			else if(PART[i].type==ELASTIC){
			double contact_pressure=ELAST.get_elas_lambda()*(PART[i].PND-PART[i].PND0)/PART[i].PND0;
			if(contact_pressure>PART[i].P)
			{
				PART[i].P=contact_pressure;//密度が非線形なのにここで粒子数密度が線形として良い？→OK（n0とρは比例）
			}
			PART[i].contact=true; //フラグを立てておく
			}
			
		}
		}
	}//for(i=0;i<PART.size();i++)ループ終了

	//壁からの反力
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++){
		if(PART[i].type==WALL){
			int neighbours=static_cast<int>(PART[i].get_current_neighboursID().size());
			for(int k=0;k<neighbours;k++) //pressureと合わせないと狂う？
			{
				int j=PART[i].get_current_neighboursID()[k];
				if(PART[i].type!=WALL){ //壁の影響半径内に入った素材の圧力の合計を壁からの反力とする。
				PART[i].P+=PART[j].P;
				}
			}
		}
	}
	

	KE*=0.5*mass;
	PE*=-g*mass;
//	EE1*=shear_modulus*mass*dimension;
	EE2*=0.5*mass;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//全系の弾性エネルギー
	ELAST.set_elastic_energy1(EE1);
	ELAST.set_elastic_energy2(EE2);
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);
}

void contact_judge(vector<mpselastic> &PART, elastic &ELAST)
{
	//アルゴリズム
	// 0. i周辺の粒子数密度が増加した場合，影響半径内にある粒子を探索し，以下を行う
	// 1. 「接触の可能性がある粒子」（(PART[j].PND>PART[j].PND0)が真？）を調べる
	// 2. 圧力を置換
	// 3. 初期配置の粒子と重複しないように接触の可能性がある粒子との間で力を計算する
	double dimension=static_cast<double>(ELAST.get_dimension());
	double re=ELAST.get_r();
	double le=ELAST.get_distancebp();
	double ground=ELAST.get_ground_position();

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
//		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			//初期配置の粒子とは普通に圧力勾配を計算（内力計算）
			if(PART[i].contact==false)
			{
				size_t neighbourN0=PART[i].get_initial_neighboursID().size();//初期配置での周辺粒子数を取得
				double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};

				for(int k=0;k<neighbourN0;k++)
				{
					int j=PART[i].get_initial_neighboursID()[k];

					//現在配置での相対座標
					double r_ij[3], r_ji[3];
					for(int D=0;D<3;D++)
					{
						r_ij[D]=PART[j].r_temp[D]-PART[i].r_temp[D];
						r_ji[D]=-r_ij[D];
					}

					//現在粒子間距離
					double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis);

					//弾性体の内力では重み付けにdis0を使うべき（初期配置からの変形が重要なので）
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++)
					{
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					double r_ij_zero[3], r_ji_zero[3];
					//calc_r0_ij(PART)で"現在位置のクォータニオンを使って"回転させた初期位置ベクトルを取得
					for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
					for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離

					double w=kernel(re, dis0);

					double press_accel[3];
					for(int D=0;D<3;D++){
						//重み関数の加重平均のとり方を変更。iの粒子数密度で計算する
						//重みは初期の距離を使うが勾配は現在配置の物を使う
						//press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis0/dis0;
						//勾配にdis0を使ってはいけない？
						press_accel[D]=(PART[i].P*r_ij_zero[D]-PART[j].P*r_ji_zero[D])*w/dis/dis;
					}

					PART[i].add_pressure(press_accel);
				}
	//			double coef=-dimension/PART[i].get_density()/PART[i].PND;
				double coef=-dimension/PART[i].get_density()/PART[i].PND0;
				PART[i].mul_pressure(coef);
			}
			else if(PART[i].contact==true)
			{
				//現在位置での周辺粒子数を取得
				size_t neighbourN=PART[i].get_current_neighboursID().size();
				double press_accel_temp[3]={0.0, 0.0, 0.0};

				for(int k=0;k<neighbourN;k++)
				{
					int j=PART[i].get_current_neighboursID()[k];

					//初期配置では近くに存在しない粒子との圧力勾配を計算する
					//現在配置での相対座標
					double r_ij[3], r_ji[3];
					for(int D=0;D<3;D++)
					{
						r_ij[D]=PART[j].r_temp[D]-PART[i].r_temp[D];
						r_ji[D]=-r_ij[D];
					}

					double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離

					double w=kernel(re, dis);
	//				if(PART[j].type!=WALL){
					for(int D=0;D<3;D++){
						press_accel_temp[D]+=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;
					}
	//				}
				}//for(int k=0;k<neighbourN;k++)

				double coef=-dimension/PART[i].get_density()/PART[i].PND;//これは弾性変形とは限らないのでPNDで割る
				for(int D=0;D<3;D++) press_accel_temp[D]*=coef;

				PART[i].add_pressure(press_accel_temp);
			}
//			if(PART[i].r_temp[A_Z]<ground+0.5*le) PART[i].set_stop_on_floor(true);  //接触確認フラグ
//			else PART[i].set_stop_on_floor(false);
		}
	}
}

void modify_acceleration(vector<mpselastic> &PART, elastic &ELAST, double **F)
{
	//加速度をゼロにしても速度成分はゼロにならない（慣性で動き続けるため）
	//加速度ゼロ＝力の釣り合い
	//いずれ速度の微分値がゼロになる・・・？
	//ならなければゼロになるような条件を付加

	double dt=ELAST.get_dt();//OK?確認する・・・
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	double total_accel[3];
	double dimension=static_cast<double>(ELAST.get_dimension());
	double re=ELAST.get_r();
	double le=ELAST.get_distancebp();

	g[A_Z]=ELAST.get_g();

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		for(int D=0;D<3;D++) total_accel[D]=(PART[i].get_stress_accel(D)+PART[i].get_stress_visco_accel(D)+PART[i].get_pressure_accel(D)+PART[i].PAcc[D]+g[D]+F[D][i]/mass);

		//下向き加速度が増加した場合←これは過拘束では？
		if((PART[i].get_stop_on_floor()==true) && (total_accel[A_Z]<0.0))
		{
			total_accel[A_Z]+=-2.0*PART[i].u_temp[A_Z]/dt;//完全弾性衝突を仮定
			PART[i].set_total_accel(total_accel);
			PART[i].set_acceleration_upward(false);
		}
		else
		{
			PART[i].set_total_accel(total_accel);
			PART[i].set_acceleration_upward(true);
		}
	}
}

void calc_pre_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, double **F)
{
	//F[D][i]これは粒子一つ一つに対してx, y, z方向それぞれの成分を持つ(節点力法を利用すべき) //Fの単位は[N]!!!
	int dimension=ELAST.get_dimension();
	double dt=ELAST.get_dt();
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	int symplectic=ELAST.get_symp_flag();
	int symplectic_order=ELAST.get_symp_order();
	
	if(dimension==2) g[A_Y]=ELAST.get_g(); 
	if(dimension==3) g[A_Z]=ELAST.get_g();

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		//仮の位置・仮の速度の更新
		for(int D=0;D<3;D++)
		{
			PART[i].r_temp[D]=PART[i].r[D];
			PART[i].u_temp[D]=PART[i].u[D];
		}
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST|| PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2)
		{
			if(dimension==3)
			{
				//3Dではシンプレクティック条件を満たさないのに適用していた！オイラー法で計算すると・・・
				if(symplectic==OFF)
				{
					for(int D=0;D<dimension;D++){
//							double u_temp=PART[i].u[D];
//							PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+ELAST.get_P_visco_stress(D, i)+g[D]+F[D][i]/mass);
//							PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+g[D]+F[D][i]/mass);
					}
				}
				else
				{
//					check_velocity_and_position(PART, i, mass, g, F);

					for(int D=0;D<dimension;D++){
				//		PART[i].u_temp[D]+=dt*(PART[i].get_stress_accel(D)+PART[i].get_pressure_accel(D)+PART[i].get_stress_visco_accel(D)+g[D]+F[D][i]/mass);
				//		PART[i].r_temp[D]+=dt*PART[i].u[D];
						PART[i].u_temp[D]=PART[i].u[D];
						PART[i].r_temp[D]=PART[i].r[D];
					}
				}
			}
		}
	}
}

void calc_post_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST,int t)
{
	cout<<"calc_post_velocity_and_position"<<endl;
	//F[D][i]これは粒子一つ一つに対してx, y, z方向それぞれの成分を持つ(節点力法を利用すべき) //Fの単位は[N]!!!
	int dimension=ELAST.get_dimension();
	double dt=ELAST.get_dt();
	double density=ELAST.get_density();
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	int symplectic=ELAST.get_symp_flag();
	int symplectic_order=ELAST.get_symp_order();

	double accleration_epsilon=1.0e-12;
	double ground=ELAST.get_ground_position();
	double le=ELAST.get_distancebp();

	if(dimension==2) g[A_Y]=ELAST.get_g(); 
	if(dimension==3) g[A_Z]=ELAST.get_g();

	int tei=0;
	int tei2=120;
	int ue=2904;

if(ELAST.get_model_number()==4){
	if(ELAST.get_poise_flag()==ON){	//初めの方はできるだけゆっくり引っ張りたい
		for(int i=tei;i<tei2;i++){
			PART[i].u[A_Z]=-0.05;//-0.05
			PART[i+ue].u[A_Z]=0.05;//0.05 
		}
	}
	else {
		for(int i=tei;i<tei2;i++){
			PART[i].u[A_Z]=0.0;
			PART[i+ue].u[A_Z]=0.0;
		}
	}
}
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type!=WALL && PART[i].type!=TERMINAL1 && PART[i].type!=TERMINAL2){	//動かない粒子
			if(symplectic==OFF)
			{
				for(int D=0;D<dimension;D++){
	//						double u_temp=PART[i].u[D];
	//						PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+ELAST.get_P_visco_stress(D, i)+g[D]+F[D][i]/mass);
	//						PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+g[D]+F[D][i]/mass);
				}
			}
			else
			{
	//			check_velocity_and_position(PART, i, mass, g, F);

				//速度の更新
				if(ELAST.get_model_number()==4){
				if(!((i>=tei && i<tei2) || (i>=ue && i<ue+tei2))){
				for(int D=0;D<dimension;D++) PART[i].u[D]+=dt*PART[i].get_total_accel(D);
				}
				}
				else {
					for(int D=0;D<dimension;D++) PART[i].u[D]+=dt*PART[i].get_total_accel(D);
				}
				//速度の修正
				//加速度をゼロにしても慣性で動き続けるのでその場合は速度をゼロにする
				if(PART[i].get_stop_on_floor()==true)//床の近くにないなら特に何もしない
				{
						//運動量保存則を満たすようにF=Δmv/Δtを決める
						//完全弾性衝突を仮定する
						//cout<<"PART["<<i<<"].u[A_Z]="<<PART[i].u[A_Z]<<endl;
						//PART[i].u[A_Z]*=0.0;これでもOK
						PART[i].u[A_Z]*=-1.0;
				}
				
				//位置の更新	
				for(int D=0;D<dimension;D++) {
					PART[i].r[D]+=dt*PART[i].u[D];//次のステップの位置を格納
				}
			}
		}
	}
}

void calc_rigid(vector<mpselastic> &PART,vector<Rigidbody> &rigids){
	
	vector<mpselastic> PARTa;
	vector<mpselastic> PARTb;
	int j=0, k=0;
	bool calcf=false;
	//現在の粒子位置を移動
	for(int i=0;i<PART.size();i++){
		if(PART[i].type==TERMINAL1){
			PARTa.push_back(PART[i]);
			calcf=true;
			}		
		else if(PART[i].type==TERMINAL2){
			PARTb.push_back(PART[i]);
			calcf=true;
			}
	}
	if(calcf==true){
		cout<<"剛体計算開始"<<endl;
		rigids[0].Renew_part_r_v(PARTa);
//		rigids[1].Renew_part_r_v(PARTb);
	//全体粒子リストに剛体粒子を戻す
	//粒子を格納する順番は見つけてきた順番と同じはず
		cout<<"剛体粒子を全粒子リストに上書き"<<endl;
	for(int i=0;i<PART.size();i++){
		if(PART[i].type==TERMINAL1){
			PART[i].r[A_X]=PARTa[j].r[A_X];
			PART[i].r[A_Y]=PARTa[j].r[A_Y];
			PART[i].r[A_Z]=PARTa[j].r[A_Z];
			j++;
			}		
		else if(PART[i].type==TERMINAL2){
			PART[i].r[A_X]=PARTb[k].r[A_X];
			PART[i].r[A_Y]=PARTb[k].r[A_Y];
			PART[i].r[A_Z]=PARTb[k].r[A_Z];
			k++;
			}
	}
	//剛体の重心を移動
	cout<<"重心移動計算開始"<<endl;
	rigids[0].Get_rigid_move(PART);
//	rigids[1].Get_rigid_move(PART);

	Micro_AVS avs;
	for(int i=0;i<PART.size();i++){
		avs.make_list(PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z],0,0,0);
	}
	avs.Output_mgf_MicroAVS("move_particle",1);
	cout<<"剛体計算終了"<<endl;
	}
//	getchar();
	PARTa.clear();
	PARTb.clear();
	
}

//クォータニオンによるset_r0_ij()の計算
//初期周辺粒子をすべてRiで回転させて保持する・・・使いまわしたいので
void calc_r0_ij(vector<mpselastic> &PART)
{
	//PART[i].ang[D]はすでに計算されているという前提で←calc_quaternion(), calc_angular_velocity();
	for(int i=0;i<PART.size();i++)
	{
//		if(PART[i].type==ELASTIC)
		{
			size_t neighboursN0=PART[i].get_initial_neighboursID().size();
			double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};

			for(int k=0;k<neighboursN0;k++)
			{
				int j=PART[i].get_initial_neighboursID()[k];
				double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

				double r_ij_Init[3]; for(int D=0;D<3;D++) r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
				double r_ji_Init[3]; for(int D=0;D<3;D++) r_ji_Init[D]=-r_ij_Init[D];
				double r0_ij[3], r0_ji[3];

				rotate_r0(r_ij_Init, qi, r0_ij); //qiを用いて回転
				rotate_r0(r_ji_Init, qj, r0_ji); //qjを用いて回転

				PART[i].set_r0_ij(r0_ij);
				PART[i].set_r0_ji(r0_ji);
			}
		}
	}
//	cout<<"回転成分の計算完了"<<endl;
}

//クォータニオンを用いた回転行列の計算
void rotate_r0(double const *rInit, double const *q, double *result)
{
	result[0]=rInit[0]*(1-2*q[1]*q[1]-2*q[2]*q[2])+rInit[1]*(2*q[0]*q[1]-2*q[3]*q[2])+rInit[2]*(2*q[0]*q[2]+2*q[3]*q[1]);
	result[1]=rInit[0]*(2*q[0]*q[1]+2*q[3]*q[2])+rInit[1]*(1-2*q[0]*q[0]-2*q[2]*q[2])+rInit[2]*(2*q[1]*q[2]-2*q[3]*q[0]);
	result[2]=rInit[0]*(2*q[0]*q[2]-2*q[3]*q[1])+rInit[1]*(2*q[1]*q[2]+2*q[3]*q[0])+rInit[2]*(1-2*q[0]*q[0]-2*q[1]*q[1]);
}

//角速度の更新・・・LU分解にポインタでなくvectorの参照を渡す
//shared_ptrでも良い
void calc_angular_velocity_using_vectorSTL(elastic &ELAST, vector<mpselastic> &PART)
{
	double re=ELAST.get_r();
	
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		//行列をゼロに初期化
		vector<vector<double>> A(3, vector<double>(3, 0.0));
		vector<double> b(3, 0.0);
		
		size_t neighboursN0=PART[i].get_initial_neighboursID().size(); //この方が意味が通りやすい

//			if(PART[i].type==static_cast<int>(ELASTIC)) //WALLの角速度は計算しない
		{
			for(int k=0;k<neighboursN0;k++)
			{
				int j=PART[i].get_initial_neighboursID()[k];//参照を取ってきたほうが早い

				double X=PART[j].r_temp[A_X]-PART[i].r_temp[A_X];
				double Y=PART[j].r_temp[A_Y]-PART[i].r_temp[A_Y];
				double Z=PART[j].r_temp[A_Z]-PART[i].r_temp[A_Z];
			
				double vX=PART[j].u_temp[A_X]-PART[i].u_temp[A_X];
				double vY=PART[j].u_temp[A_Y]-PART[i].u_temp[A_Y];
				double vZ=PART[j].u_temp[A_Z]-PART[i].u_temp[A_Z];

				double dis0=PART[i].get_initial_distancebps()[k];	//初期粒子間距離
//				double dis=PART[i].get_current_distancebps()[k]; //NG！影響半径内にある粒子がjとは限らない！・・・この関数必要？
				double dis=sqrt(X*X+Y*Y+Z*Z);//これはOK。iとjとの関係を考えているので

				double w=kernel(re, dis0);

				//係数行列
				A[0][0]+=(Y*Y+Z*Z)*w/dis/dis;
				A[1][1]+=(Z*Z+X*X)*w/dis/dis;
				A[2][2]+=(X*X+Y*Y)*w/dis/dis;			//A[2][2]=(X*X+Y*Y)*w/dis/dis;＋＝が正解・・・2012-09-13
				A[0][1]-=(X*Y)*w/dis/dis;
				A[0][2]-=(Z*Y)*w/dis/dis;
				A[1][2]-=(Y*Z)*w/dis/dis;

				//右辺ベクトル
				b[0]+=(Y*vZ-Z*vY)*w/dis/dis;
				b[1]+=(Z*vX-X*vZ)*w/dis/dis;
				b[2]+=(X*vY-Y*vX)*w/dis/dis;
//				vector_product(r, v, b);				//この方式だと+=が合計されない（bがその都度リセットされる）
//				for(int D=0;D<3;D++) b[D]*=w/dis/dis;

			}

			//要素が対称なので数え上げが終わってから代入
			A[1][0]=A[0][1];
			A[2][0]=A[0][2];
			A[2][1]=A[1][2];

			//LU分解法でAω=bを解く。解はbに格納される
			lu_decomposition(A, b);

			for(int D=0;D<3;D++)
			{
				PART[i].ang_u_temp[D]=b[D];
//				cout<<"PART["<<i<<"].ang_u["<<D<<"]="<<PART[i].ang_u[D]<<" ";
			}
//			cout<<endl;
		}
	}

//	cout<<"角速度計算完了"<<endl;
}

//並列計算用・・・LU分解にポインタでなくvectorを渡す
void calc_quaternion_using_vectorSTL(elastic &ELAST, vector<mpselastic> &PART)
{
	double re=ELAST.get_r();
	double r0[3];

	const int MAX_ITERATION=1000;		//最大繰り返し回数
	const double EPSILON=1.0e-6;//pow(10.0, -12); //トレランス。マシンイプシロンはDBL_EPSILONで得る 約1.0e-16
	int particles_over_max_iteration=0; //MAX_ITERATIONを超えた粒子数をカウントする

	//クォータニオン用の係数行列用関数ポインタテーブル
	double (* const rotation[])(const double *r, const double *r0, const vector<double> &q)={rotx, roty, rotz};
	double (* const calc_jacobi_matrix[3][4])(const double *r, const double *r0, const vector<double> &q)={
		{rotx_x, rotx_y, rotx_z, rotx_s},
		{roty_x, roty_y, roty_z, roty_s},
		{rotz_x, rotz_y, rotz_z, rotz_s},
	};
	double (* const jacobi_norm[4])(const vector<double> &q)={rotn_x, rotn_y, rotn_z, rotn_s};

	vector<vector<double>> jacobi_matrix(4, vector<double>(4, 0.0));
	vector<double> d(4, 0.0);
	vector<double> qi(4);

	for(int i=0;i<PART.size();i++)
	{
		{
			int itr=0; //反復回数のリセット
			size_t neighboursN0=PART[i].get_initial_neighboursID().size();//re内に含まれる初期粒子数・・・PART[i].Nは初期ではない！
			
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					jacobi_matrix[j][k]=0.0;
				}
				d[j]=0.0;
				qi[j]=PART[i].ang[j];
			}

			do{
				for(int k=0;k<neighboursN0;k++)//re内で和を取る
				{
					//j: ID, k: 配列の添字！！ IDと粒子の添字の対照には気をつけること！！
					int j=PART[i].get_initial_neighboursID()[k];

					if(j!=i)//このifはいらない
					{
						//iからみた相対座標・・・現在の影響半径内の粒子ではなく、初期に近傍にあった粒子のIDを用いる。従ってPART[i].get_current_position();は無意味！！！
						double X=PART[j].r_temp[A_X]-PART[i].r_temp[A_X];
						double Y=PART[j].r_temp[A_Y]-PART[i].r_temp[A_Y];
						double Z=PART[j].r_temp[A_Z]-PART[i].r_temp[A_Z];

						double dis0=PART[i].get_initial_distancebps()[k];//r_init//これは[i][j]ではなく[i][k]！！
						PART[i].get_initial_neighbours_position(k, r0);

						double w=kernel(re, dis0);

						//右辺ベクトルの作成
						double r[3]={X, Y, Z};//r_ij
						for(int D=0;D<3;D++) d[D]-=rotation[D](r, r0, qi)*w/dis0/dis0;//Σに注意！(-1)を忘れないように.和を取る前に初期化！

						//ヤコビ行列の作成(改善→係数の和を渡すようにする！)
						for(int row=0;row<3;row++)
							for(int col=0;col<4;col++)
								jacobi_matrix[row][col]+=calc_jacobi_matrix[row][col](r, r0, qi)*w/dis0/dis0;
					}
				}

				//規格化条件
				d[3]=-rotn(qi);
				for(int col=0;col<4;col++) jacobi_matrix[3][col]=jacobi_norm[col](qi); //規格化条件の偏微分
				
				lu_decomposition(jacobi_matrix, d);
				for(int D=0;D<4;D++) qi[D]+=d[D]; //qi[]の更新
		
				//反復回数のインクリメント
				itr++;
			}while(check_q_norm(qi, EPSILON) && itr<MAX_ITERATION);

			if(itr==MAX_ITERATION){
				particles_over_max_iteration++;
			}else{
				for(int D=0;D<4;D++) PART[i].ang[D]=qi[D];
//				std::cout<<"PART["<<i<<"], finished, iteration: "<<itr<<", matrix converged"<<std::endl;
//				cout<<"q convergence: "<<boolalpha<<check_q_residue(qi, EPSILON)<<"; ";
//				for(int m=0;m<4;m++) cout<<"q["<<m<<"]="<<qi[m]<<" ";
//				cout<<"iteration: "<<itr<<endl;
			}
		}
	}
	cout<<"particles_over_max_iteration: "<<particles_over_max_iteration<<endl;
}

//回転行列と規格化条件の偏導関数
inline double rotx(double *r, double *r0, double *q)
{
	return (r[1]*(2*(r0[0]*(q[2]*q[0]-q[1]*q[3])+r0[1]*(q[1]*q[2]+q[0]*q[3]))+r0[2]*(q[2]*q[2]+q[3]*q[3]-q[0]*q[0]-q[1]*q[1]))-r[2]*(2*(r0[0]*(q[0]*q[1]+q[2]*q[3])+r0[2]*(q[1]*q[2]-q[0]*q[3]))+r0[1]*(q[1]*q[1]+q[3]*q[3]-q[2]*q[2]-q[0]*q[0])));
}

inline double roty(double *r, double *r0, double *q)
{
	return (r[2]*(2*(r0[1]*(q[0]*q[1]-q[2]*q[3])+r0[2]*(q[2]*q[0]+q[1]*q[3]))+r0[0]*(q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2]))-r[0]*(2*(r0[0]*(q[2]*q[0]-q[1]*q[3])+r0[1]*(q[1]*q[2]+q[0]*q[3]))+r0[2]*(q[2]*q[2]+q[3]*q[3]-q[0]*q[0]-q[1]*q[1])));
}

inline double rotz(double *r, double *r0, double *q)
{
	return (r[0]*(2*(r0[0]*(q[0]*q[1]+q[2]*q[3])+r0[2]*(q[1]*q[2]-q[0]*q[3]))+r0[1]*(q[1]*q[1]+q[3]*q[3]-q[2]*q[2]-q[0]*q[0]))-r[1]*(2*(r0[1]*(q[0]*q[1]-q[2]*q[3])+r0[2]*(q[2]*q[0]+q[1]*q[3]))+r0[0]*(q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2])));
}

inline double rotn(double *q)
{
	double norm=0.0;
	for(int i=0;i<4;i++) norm+=q[i]*q[i];
	return norm-1.0;
}

//回転行列と規格化条件の偏導関数(オーバーロード)
inline double rotx(const double *r, const double *r0, const vector<double> &q)
{
	return (r[1]*(2*(r0[0]*(q[2]*q[0]-q[1]*q[3])+r0[1]*(q[1]*q[2]+q[0]*q[3]))+r0[2]*(q[2]*q[2]+q[3]*q[3]-q[0]*q[0]-q[1]*q[1]))-r[2]*(2*(r0[0]*(q[0]*q[1]+q[2]*q[3])+r0[2]*(q[1]*q[2]-q[0]*q[3]))+r0[1]*(q[1]*q[1]+q[3]*q[3]-q[2]*q[2]-q[0]*q[0])));
}

inline double roty(const double *r, const double *r0, const vector<double> &q)
{
	return (r[2]*(2*(r0[1]*(q[0]*q[1]-q[2]*q[3])+r0[2]*(q[2]*q[0]+q[1]*q[3]))+r0[0]*(q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2]))-r[0]*(2*(r0[0]*(q[2]*q[0]-q[1]*q[3])+r0[1]*(q[1]*q[2]+q[0]*q[3]))+r0[2]*(q[2]*q[2]+q[3]*q[3]-q[0]*q[0]-q[1]*q[1])));
}

inline double rotz(const double *r, const double *r0, const vector<double> &q)
{
	return (r[0]*(2*(r0[0]*(q[0]*q[1]+q[2]*q[3])+r0[2]*(q[1]*q[2]-q[0]*q[3]))+r0[1]*(q[1]*q[1]+q[3]*q[3]-q[2]*q[2]-q[0]*q[0]))-r[1]*(2*(r0[1]*(q[0]*q[1]-q[2]*q[3])+r0[2]*(q[2]*q[0]+q[1]*q[3]))+r0[0]*(q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2])));
}

inline double rotn(const vector<double> &q)
{
	double norm=0.0;
	for(int i=0;i<4;i++) norm+=q[i]*q[i];
	return norm-1.0;
}

//クォータニオンの残差をチェックする
bool check_q_norm(const double *q, const double EPSILON)
{
	double norm=0.0;
	for(unsigned i=0;i<4;i++) norm+=q[i]*q[i];
	return (fabs(1.0-norm)>EPSILON) ? true: false;
}

bool check_q_norm(const vector<double> &q, const double EPSILON)
{
	double norm=0.0;
	for(unsigned i=0;i<4;i++) norm+=q[i]*q[i];
	return (fabs(1.0-norm)>EPSILON) ? true: false;
}

//エネルギーの確認
void calc_hamiltonian(elastic &ELAST, int t)
{
//	double hamiltonian=ELAST.get_elastic_energy()+ELAST.get_kinetic_energy()+ELAST.get_potential_energy();
	ofstream fout1("./Elastic/hamiltonian.dat", ios::app);
	ofstream fout2("./Elastic/kinetic_energy.dat", ios::app);
	ofstream fout3("./Elastic/elastic_energy.dat", ios::app);
	ofstream fout4("./Elastic/elastic_energy1.dat", ios::app);
	ofstream fout5("./Elastic/elastic_energy2.dat", ios::app);
	ofstream fout6("./Elastic/potential_energy.dat", ios::app);
	
	if(fout1.fail() || fout2.fail() || fout3.fail() || fout4.fail() || fout5.fail() || fout6.fail()){
		system("mkdir Elastic");
		ofstream fout1("./Elastic/hamiltonian.dat", ios::app);
		ofstream fout2("./Elastic/kinetic_energy.dat", ios::app);
		ofstream fout3("./Elastic/elastic_energy.dat", ios::app);
		ofstream fout4("./Elastic/elastic_energy1.dat", ios::app);
		ofstream fout5("./Elastic/elastic_energy2.dat", ios::app);
		ofstream fout6("./Elastic/potential_energy.dat", ios::app);

		if(fout1.fail() || fout2.fail() || fout3.fail() || fout4.fail())//再試行
		{
			cout<<"エネルギー確認用のディレクトリが開けません"<<endl;
			cout<<"プログラムを終了します"<<endl;
			exit(1);
		}
	}

	fout1<<t<<"\t"<<ELAST.get_hamiltonian()<<endl;
	fout2<<t<<"\t"<<ELAST.get_kinetic_energy()<<endl;
	fout3<<t<<"\t"<<ELAST.get_elastic_energy()<<endl;
	fout4<<t<<"\t"<<ELAST.get_elastic_energy1()<<endl;
	fout5<<t<<"\t"<<ELAST.get_elastic_energy2()<<endl;
	fout6<<t<<"\t"<<ELAST.get_potential_energy()<<endl;

	fout1.close();
	fout2.close();
	fout3.close();
	fout4.close();
	fout5.close();
	fout6.close();

	//ofstream fout5("difference_of_EE.dat", ios::app);
	//fout5<<t<<"\t"<<(ELAST.get_elastic_energy()-ELAST.get_last_elastic_energy())<<endl;
	//fout5.close()

/*
	if(ELAST.get_FEM_flag()==ON)
	{
		//本来はELASTIC以外にもエネルギーを計算すべき
		//<0だとちょっとした誤差でスイッチが入るので<
//		if((ELAST.get_elastic_energy()-ELAST.get_last_elastic_energy())<0)
		if((ELAST.get_elastic_energy()-ELAST.get_last_elastic_energy())<-1e-10)
		{
			ELAST.set_FEM_switch(ON);
			cout<<"step: "<<t<<"FEM switch: "<<boolalpha<<ELAST.get_FEM_switch()<<endl;
		}else{
			cout<<"FEM switch: "<<boolalpha<<ELAST.get_FEM_switch()<<endl;
		}
	}
*/
	ELAST.set_last_elastic_energy(ELAST.get_elastic_energy());
	ELAST.set_last_kinetic_energy(ELAST.get_kinetic_energy());
	ELAST.set_last_potential(ELAST.get_potential_energy());

}

//ベクトル積の計算
void vector_product(double *a, double *b, double *result)
{
	result[0]=a[A_Y]*b[A_Z]-a[A_Z]*b[A_Y];
	result[1]=a[A_Z]*b[A_X]-a[A_X]*b[A_Z];
	result[2]=a[A_X]*b[A_Y]-a[A_Y]*b[A_X];
}

//加速度の確認
void check_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, const int i, double **F)
{
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	if(2==ELAST.get_dimension()) g[1]=ELAST.get_g();
	else g[2]=g[1]=ELAST.get_g();

	cout<<"checking 3d velocity and position"<<endl;
	for(int D=0;D<3;D++) cout<<"PART["<<i<<"].u["<<D<<"]="<<PART[i].u[D]<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"PART["<<i<<"].r["<<D<<"]="<<PART[i].r[D]<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"PART["<<i<<"].pr["<<D<<"]="<<PART[i].get_pressure_accel(D)<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"g["<<D<<"]="<<g[D]<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"F["<<D<<"]["<<i<<"]"<<F[D][i]<<" ";
	cout<<endl;
	for(int D=0;D<3;D++) cout<<"F["<<D<<"]["<<i<<"]/mass="<<F[D][i]/mass<<" ";
	cout<<endl;
}

//密度の修正
void calc_modified_density(vector<mpselastic> &PART, elastic &ELAST)
{
	//体積と粒子数を用いたものはfreeon()で実行
	for(int i=0;i<PART.size();i++)
	{
//		cout<<"initial density: "<<PART[i].get_density()<<endl;
		double ratio=PART[i].PND/PART[i].PND0;
		PART[i].set_density(ELAST.get_density()*ratio);
//		cout<<"current density: "<<PART[i].get_density()<<endl;
	}
}

void rigid_calc(vector<mpselastic> &PART, elastic &ELAST)
{
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2) //剛体粒子
		{

		}
	}
}

//数値微分によるヤコビ行列
//中心差分近似。最適な刻み幅の存在に注意！！！
void numerical_jacobian(double **J, double *r, double *r0, double *q)
{
	//関数ポインタテーブル
	double (* const rotation[])(double *r, double *r0, double *q)={rotx, roty, rotz};

	const double h=1.0e-3;//刻み幅。アダプティブに変更すべき
	double result0=0.0;

	//ゼロになる場合は除く
	result0=((rotation[0](r, r0, q))-(rotation[0](r, r0, q)))/(2*h);

}

void calc_nonlinear_elastic(vector<mpselastic> &PART, elastic &ELAST, int t, double **F)
{
	vector<vector<double>> residual_acceleration(3, vector<double>(PART.size(), 0.0)); //残差加速度

	//初期化
	for(int D=0;D<ELAST.get_dimension();D++)
	{
		for(int i=0;i<PART.size();i++)
		{
			residual_acceleration[D][i]=F[D][i];
		}
	}

	while(true)
	{	
	//	calc_quaternion(ELAST, PART); //クォータニオンの計算//OK
		calc_quaternion_using_vectorSTL(ELAST, PART);
	//	calc_angular_velocity(ELAST, PART); //角速度の計算//OK
		calc_angular_velocity_using_vectorSTL(ELAST, PART);
		calc_r0_ij(PART); //初期位置ベクトルの回転
		calc_nonlinear_accel_for_3D(PART, ELAST);
		calc_pressure_and_contact(PART, ELAST);	//圧力と接触反力による加速度を計算
		PART[2].check_acceleration(); //粒子の加速度をチェック

		calc_residual_acceleration(PART, ELAST, residual_acceleration, t, F);

		calc_nonlinear_velocity_and_position(PART, ELAST, residual_acceleration);
		calc_hamiltonian(ELAST, t);
	}
}

void calc_nonlinear_accel_for_3D(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //これは既にre*le
	double vis=ELAST.get_nensei(); //粘度（動粘度ではない）
	double g=ELAST.get_g();

	double KE=0.0;	//運動エネルギー
	double EE1=0.0;	//ひずみエネルギー
	double EE2=0.0;	//体積ひずみによる
	double PE=0.0;	//ポテンシャル
	double ground=ELAST.get_ground_position(); //床のｚ座標の取得

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
		double lambda_average=0.0; //lambdaの平均値
		double particle_density=0.0; //粒子数密度・・・接触判定で使う 

		//WALLは弾性変形しないので除外
//		if(PART[i].type==ELASTIC || PART[i].type==INELASTIC || PART[i].type==BOELASTIC)
		{
			double EE1_temp=0.0;
			double EE2_temp=0.0;
			double volumetric_strain=0.0; //体積ひずみ

			//運動エネルギーの更新
			for(int D=0;D<3;D++) KE+=PART[i].u[D]*PART[i].u[D];

			//位置ポテンシャルの更新
			PE+=(PART[i].r[2]-ground);

			int neighboursN0=static_cast<int>(PART[i].get_initial_neighboursID().size());

			for(int k=0;k<neighboursN0;k++) //pressureと合わせないと狂う？
			{

				int j=PART[i].get_initial_neighboursID()[k];

				//初期配置でNONELASTICなものが近くにある場合は条件を分ける
				//これがないとエネルギー収支もおかしくなるはず・・・
				if(PART[j].type==ELASTIC)
				{				
					//iからみた初期配置での相対座標
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++){
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					//現在配置での相対座標
					double r_ij[3], r_ji[3];
					for(int D=0;D<3;D++)
					{
						r_ij[D]=PART[j].r[D]-PART[i].r[D];
						r_ji[D]=-r_ij[D];
					}

					double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離
					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離
			
					//回転行列Riによるr_ij_Initの回転・・・以降はr_ij_zeroを使う
					//qiは時々刻々変わるので使い回しできない・・・PARTのパラメータとして与えるのが良い
					double r_ij_zero[3], r_ji_zero[3];

					for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
					for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

				//重み関数
					double w=kernel(re, dis0);

				//粒子数密度の更新
					particle_density+=w;

					//回転後の初期配置相対位置の単位ベクトル
					double n_ij[3], n_ji[3];
					for(int D=0;D<3;D++) 
					{
						n_ij[D]=r_ij_zero[D]/dis0;
						n_ji[D]=r_ji_zero[D]/dis0;
					}

					//変位ベクトル
					double U_ij[3], U_ji[3]; 
					for(int D=0;D<3;D++)
					{
						U_ij[D]=r_ij[D]-r_ij_zero[D];
						U_ji[D]=r_ji[D]-r_ji_zero[D];
					}

				//ひずみベクトル
					double E_ij[3], E_ji[3];
					for(int D=0;D<3;D++)
					{
						E_ij[D]=U_ij[D]/dis0;
						E_ji[D]=U_ji[D]/dis0;
					}

				//ヤング率とポアソン比の計算
					double En_ij[3], En_ji[3];//縦ひずみベクトル
					double naiseki=0.0;
					for(int D=0;D<3;D++) naiseki+=E_ij[D]*r_ij[D];
					for(int D=0;D<3;D++) En_ij[D]=naiseki*r_ij[D]/dis/dis;	//ひずみベクトルをr_ij方向へ正射影
					double ns=sqrt(En_ij[0]*En_ij[0]+En_ij[1]*En_ij[1]+En_ij[2]*En_ij[2]); //縦ひずみ

					//シリコーンの場合
//					youngs_modulus=-6.11619E-11*x^9+1.66522E-8*x^8-1.93253E-6*x^7+1.2459E-4*x^6-4.87141E-3*x^5+0.118275*x^4-1.75252*x^3+14.9237*x^2-64.3273*x+123.411;
					double youngs_modulus=(((((((((-6.11619E-11*ns+1.66522E-8)*ns-1.93253E-6)*ns+1.2459E-4)*ns-4.87141E-3)*ns+0.118275)*ns-1.75252)*ns+14.9237)*ns-64.3273)*ns+123.411)*100;
//					-9.6723E-12*x^10+1.22758E-9*x^9-6.73685E-8*x^8+2.09375E-7*x^7-4.05552E-5*x^6+5.08152E-4*x^5-0.00414923*x^4+0.0218461*x^3-0.072922*x^2+0.151892*x^1+0.17713
					double poisson_ratio=(((((((((-9.6723E-12*ns+1.22758E-9)*ns-6.73685E-8)*ns+2.09375E-7)*ns-4.05552E-5)*ns+5.08152E-4)*ns-0.00414923)*ns+0.0218461)*ns-0.072922)*ns+0.151892)*ns+0.17713;

//					cout<<"youngs_modulus: "<<youngs_modulus<<", poisson_ratio: "<<poisson_ratio<<endl;

					//ラメ定数mu
					double shear_modulus=youngs_modulus/(2.0*(1.0+poisson_ratio));
					//ラメ定数lambda
					double lambda=(poisson_ratio*youngs_modulus)/((1.0+poisson_ratio)*(1.0-2.0*poisson_ratio));

//					cout<<"shear_modulus: "<<shear_modulus<<", lambda: "<<lambda<<endl;

				//圧力計算・・・e_vol_iのΣ[]（変位の発散）を計算
					for(int D=0;D<3;D++) volumetric_strain+=E_ij[D]*n_ij[D]*w; //体積ひずみ計算の準備
					PART[i].P+=volumetric_strain*lambda;
					EE2_temp+=PART[i].P*volumetric_strain;
					//PART[i].P*=w;・・・何故かこれだとうまくいかない・・・wが累乗で足されるので当たり前((()*w+smt)*w+smt)*w...

					lambda_average+=lambda;

				//ひずみエネルギー
					for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w*shear_modulus;

				//応力ベクトル
					double sigma_ij[3], sigma_ji[3];
					for(int D=0;D<3;D++)
					{
						sigma_ij[D]=2*shear_modulus*w*E_ij[D]/PART[i].PND/dis0;
						sigma_ji[D]=2*shear_modulus*w*E_ji[D]/PART[j].PND/dis0;

						sigma_ij[D]-=sigma_ji[D];
					}
					PART[i].add_stress_accel(sigma_ij);

	//				if(fabs((dis-dis0)/dis0)>1.04) w=0;//破壊条件

					//ひずみ速度の計算
					double strain_vi[3], strain_vj[3];

					//iからみたひずみ速度
					strain_vi[0]=((PART[j].u[0]-PART[i].u[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
					strain_vi[1]=((PART[j].u[1]-PART[i].u[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y軸周りのひずみ速度
					strain_vi[2]=((PART[j].u[2]-PART[i].u[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z軸周りのひずみ速度
					//anglar_u1..3は配列に書き換え

					//jからみたひずみ速度
					strain_vj[0]=((PART[i].u[0]-PART[j].u[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
					strain_vj[1]=((PART[i].u[1]-PART[j].u[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y軸周りのひずみ速度
					strain_vj[2]=((PART[i].u[2]-PART[j].u[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z軸周りのひずみ速度
			
					double sigma_v_ij[3], sigma_v_ji[3];
					for(int D=0;D<3;D++)
					{
						sigma_v_ij[D]=2*vis*w*strain_vi[D]/PART[i].PND/dis;
						sigma_v_ji[D]=2*vis*w*strain_vj[D]/PART[j].PND/dis;

						sigma_v_ij[D]-=sigma_v_ji[D];
					}

					PART[i].add_stress_visco_accel(sigma_v_ij);
				}
			}//for(int k=0;k<neighboursN0;k++)ループ終了

			PART[i].P*=-dimension/PART[i].PND;//体積ひずみによる圧力が求められた（これはエネルギー計算する前に求めておくこと）
			volumetric_strain*=-dimension/PART[i].PND;//体積ひずみが求められた（これはエネルギー計算する前に求めておくこと）
			PART[i].volumetric_strain=volumetric_strain;
			//接触判定（下壁も含む）・・・ここで減衰を入れないとエラー？
			//加速度計算はELASTICだけ考えればよいのでif(.type==ELASTIC)のループの中に入れても変わらない
//			double n00=PART[i].PND0;
//			if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=(PART[i].PND-n00)/n00; //粒子数密度の増加が見られれば置換する

			double coef=dimension/PART[i].get_density();

			PART[i].mul_stress_accel(coef);
			PART[i].mul_stress_visco_accel(coef);

			//ひずみエネルギー
//			EE1+=(EE1_temp)*(shear_modulus*mass*dimension/PART[i].PND/density);	//第二項のΣ
//			EE2+=0.5*lambda*mass*(PART[i].P*PART[i].P)/density;//第三項のΣ
			EE1+=EE1_temp/PART[i].PND/PART[i].get_density();//第二項のΣ 密度はそれぞれ違う
			EE2+=EE2_temp/PART[i].get_density();//第三項のΣ 密度はそれぞれ違う

			/****エネルギー計算にもlambda_averageを使うべき・・・？****/

//			PART[i].P*=lambda;//体積ひずみによる圧力が求められた
		}//if(PART[i].type==ELASTIC || PART[i].type==INELASTIC ||PART[i].type==BOELASTIC)終了

		//それ自体は動かないがWALLも加速度（反力）を有している
		//粒子数の増減を考慮しなければならないのでこれはifの外に置く
		//・・・PART[i].Pはすべての粒子を考慮する必要があるが、エネルギー計算に入れるべきではない
		double n00=PART[i].PND0;
//2012-11-26 圧力と粒子数密度を比較していた・・・！
//		if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=lambda*(PART[i].PND-n00)/n00; //粒子数密度の増加が見られれば置換する

		lambda_average/=static_cast<int>(PART[i].get_initial_neighboursID().size());

		//粒子数密度の増加が見られれば置換する
		//このPは弾性体も壁も受けている圧力であることに注意

		if(PART[i].P<lambda_average*((PART[i].PND-n00)/n00)) PART[i].P=lambda_average*(PART[i].PND-n00)/n00; //・・・数値的に何が起きているか調べる！！

	}//for(i=0;i<PART.size();i++)

	KE*=0.5*mass;
	PE*=-g*mass;
	EE1*=mass*dimension;
	EE2*=0.5*mass;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//全系の弾性エネルギー
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);

}

void calc_residual_acceleration(vector<mpselastic> &PART, elastic &ELAST, vector<vector<double>> &res_accel, int t, double **F)
{
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	g[A_Z]=ELAST.get_g();

	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			for(int D=0;D<DIMENSION;D++) res_accel[D][i]=PART[i].get_stress_accel(D)+PART[i].get_pressure_accel(D)+PART[i].get_stress_visco_accel(D)+g[D]+F[D][i]/mass;
		}
	}
}

void calc_nonlinear_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, vector<vector<double>> &residual_acceleration)
{
	//F[D][i]これは粒子一つ一つに対してx, y, z方向それぞれの成分を持つ(節点力法を利用すべき) //Fの単位は[N]!!!
	double dt=ELAST.get_dt();

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
//			check_velocity_and_position(PART, i, mass, g, F);

			for(int D=0;D<DIMENSION;D++){
				PART[i].u[D]+=dt*residual_acceleration[D][i];
				PART[i].r[D]+=dt*PART[i].u[D];
			}
		}
	}
}

/*********テスト領域***********/
/*********テスト領域***********/
/*********テスト領域***********/
/*********テスト領域***********/
/*********テスト領域***********/
/*********テスト領域***********/
/*********テスト領域***********/
/*********テスト領域***********/
/*********テスト領域***********/


//角速度の更新・・・スケーリングかける（というか手計算で直接解いても良いのでは？）
void calc_angular_velocity(elastic &ELAST, vector<mpselastic> &PART)
{
	double re=ELAST.get_r();
	bool pivot_check=ELAST.get_pivot_check();

	//対称係数行列
	double **A=new double*[3]; for(int D=0;D<3;D++) A[D]=new double[3];
	//右辺ベクトル
	double b[3];

	#pragma omp parallel
	{
		for(int i=0;i<PART.size();i++)
		{
			//行列の初期化
			for(int m=0;m<3;m++)
			{
				for(int n=0;n<3;n++) A[m][n]=0.0;
				b[m]=0.0;
			}

			size_t neighboursN0=PART[i].get_initial_neighboursID().size(); //この方が意味が通りやすい

			if(PART[i].type==(int)ELASTIC) //WALLの角速度は計算しない
			{
				for(int k=0;k<neighboursN0;k++)
				{
					int j=PART[i].get_initial_neighboursID()[k];//参照を取ってきたほうが早い

					double X=PART[j].r[A_X]-PART[i].r[A_X];
					double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
					double Z=PART[j].r[A_Z]-PART[i].r[A_Z];
			
					double vX=PART[j].u[A_X]-PART[i].u[A_X];
					double vY=PART[j].u[A_Y]-PART[i].u[A_Y];
					double vZ=PART[j].u[A_Z]-PART[i].u[A_Z];

	//				double r[3]={X, Y, Z};
	//				double v[3]={vX, vY, vZ};

					double dis0=PART[i].get_initial_distancebps()[k];	//初期粒子間距離
	//				double dis=PART[i].get_current_distancebps()[k]; //NG！影響半径内にある粒子がjとは限らない！・・・この関数必要？
					double dis=sqrt(X*X+Y*Y+Z*Z);//これはOK。iとjとの関係を考えているので

					double w=kernel(re, dis0);

					//係数行列
					A[0][0]+=(Y*Y+Z*Z)*w/dis/dis;
					A[1][1]+=(Z*Z+X*X)*w/dis/dis;
					A[2][2]+=(X*X+Y*Y)*w/dis/dis;			//A[2][2]=(X*X+Y*Y)*w/dis/dis;＋＝が正解・・・2012-09-13
					A[0][1]-=(X*Y)*w/dis/dis;
					A[0][2]-=(Z*Y)*w/dis/dis;
					A[1][2]-=(Y*Z)*w/dis/dis;
	//				A[0][1]=A[1][0]=-(X*Y)*w/dis/dis;		//これは=-ではない・・・
	//				A[0][2]=A[2][0]=-(Z*Y)*w/dis/dis;
	//				A[1][2]=A[2][1]=-(Y*Z)*w/dis/dis;

					//右辺ベクトル
					b[0]+=(Y*vZ-Z*vY)*w/dis/dis;
					b[1]+=(Z*vX-X*vZ)*w/dis/dis;
					b[2]+=(X*vY-Y*vX)*w/dis/dis;
	//				vector_product(r, v, b);				//この方式だと+=が合計されない（bがその都度リセットされる）
	//				for(int D=0;D<3;D++) b[D]*=w/dis/dis;

				}

				//要素が対称なので数え上げが終わってから代入
				A[1][0]=A[0][1];
				A[2][0]=A[0][2];
				A[2][1]=A[1][2];

				//LU分解法でAω=bを解く。解はbに格納される
				lu_decomposition(A, b, 3, pivot_check);
	//			cout<<"angular velocity: "<<endl;
				for(int D=0;D<3;D++)
				{
					PART[i].ang_u[D]=b[D];
	//				cout<<"PART["<<i<<"].ang_u["<<D<<"]="<<PART[i].ang_u[D]<<" ";
				}
	//			cout<<endl;
			}
		}
	}

	for(int D=0;D<3;D++) delete [] A[D];
	delete [] A;
}

//クォータニオンの計算・・・規格化条件はΣしてはいけない
//|q[i]|=<1なので、値がおかしくなったらクォータニオンの定義域を考えて強制的にリセットすることが有効・・・？
void calc_quaternion(elastic &ELAST, vector<mpselastic> &PART)
{
	double re=ELAST.get_r();
	double r0[3];
	bool pivot_check=ELAST.get_pivot_check();

	const int MAX_ITERATION=1000;		//最大繰り返し回数
	const double EPSILON=1.0e-6;//pow(10.0, -12);	//トレランス。マシンイプシロンはDBL_EPSILONで得る 約1.0e-16
	int particles_over_max_iteration=0;

	//クォータニオン用の係数行列用関数ポインタテーブル
	double (* const rotation[])(double *r, double *r0, double *q)={rotx, roty, rotz};
	double (* const calc_jacobi_matrix[3][4])(double *r, double *r0, double *q)={
		{rotx_x, rotx_y, rotx_z, rotx_s},
		{roty_x, roty_y, roty_z, roty_s},
		{rotz_x, rotz_y, rotz_z, rotz_s},
	};
	double (* const jacobi_norm[4])(double *q)={rotn_x, rotn_y, rotn_z, rotn_s};
	double d[4];//右辺ベクトル
	double** jacobi_matrix=new double*[4];//ヤコビ行列

	for(int D=0;D<4;D++) jacobi_matrix[D]=new double[4];

	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST) //WALLの角速度は計算しない
		{
			int itr=0; //反復回数のリセット
			size_t neighboursN0=PART[i].get_initial_neighboursID().size();//re内に含まれる初期粒子数・・・PART[i].Nは初期ではない！
			double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};

//			for(int m=0;m<4;m++) cout<<"itr: "<<itr<<", q["<<m<<"]="<<qi[m]<<endl;
//			cout<<"PART["<<i<<"].N="<<PART[i].N<<", neighboursN0="<<neighboursN0<<endl;

//			cout<<"q初期値: "<<endl;
			for(int m=0;m<4;m++){
				d[m]=0.0;
				for(int n=0;n<4;n++) jacobi_matrix[m][n]=0.0;//jacobi行列の初期化
//				cout<<"q["<<m<<"]="<<qi[m]<<" ";
			}
//			cout<<endl;

			do{
//				for(int j=0;j<neighbours0;j++)//re内で和を取る。これはIDと照合していなかった単純なエラー！！！
				for(int k=0;k<neighboursN0;k++)//re内で和を取る
				{
					//j: ID, k: 配列の添字！！ IDと粒子の添字の対照には気をつけること！！
					int j=PART[i].get_initial_neighboursID()[k];

					if(j!=i)//このifはいらない
					{
						//iからみた相対座標・・・現在の影響半径内の粒子ではなく、初期に近傍にあった粒子のIDを用いる。従ってPART[i].get_current_position();は無意味！！！
						double X=PART[j].r_temp[A_X]-PART[i].r_temp[A_X];
						double Y=PART[j].r_temp[A_Y]-PART[i].r_temp[A_Y];
						double Z=PART[j].r_temp[A_Z]-PART[i].r_temp[A_Z];

						double dis0=PART[i].get_initial_distancebps()[k];//r_init//これは[i][j]ではなく[i][k]！！
						PART[i].get_initial_neighbours_position(k, r0);

						double w=kernel(re, dis0);

//						cout<<"dis0="<<dis0<<", r0[0]="<<r0[0]<<", r0[1]="<<r0[1]<<", r0[2]="<<r0[2]<<endl;

						//右辺ベクトルの作成
						double r[3]={X, Y, Z};//r_ij
						for(int D=0;D<3;D++) d[D]-=rotation[D](r, r0, qi)*w/dis0/dis0;//Σに注意！(-1)を忘れないように.和を取る前に初期化！

						//ヤコビ行列の作成(改善→係数の和を渡すようにする！)
						for(int row=0;row<3;row++)
							for(int col=0;col<4;col++)
								jacobi_matrix[row][col]+=calc_jacobi_matrix[row][col](r, r0, qi)*w/dis0/dis0;
					}
				}

				//規格化条件
				d[3]=-rotn(qi);
				for(int col=0;col<4;col++) jacobi_matrix[3][col]=jacobi_norm[col](qi); //規格化条件の偏微分
				
				lu_decomposition(jacobi_matrix, d, 4, pivot_check);
//				pivot_gauss(jacobi_matrix, d, 4, pivot_check); //ヤコビアンは必ずJ(1, 1)>0なのでpivot選択なしの方が高速→改良
				for(int D=0;D<4;D++) qi[D]+=d[D]; //qi[]の更新
		
				//反復回数のインクリメント

				itr++;
			}while(check_q_norm(qi, EPSILON) && itr<MAX_ITERATION);
//			}while(vector_norm2(d, 0, 4)>EPSILON && itr<MAX_ITERATION);//ここは本来なら残差を用いるべき

			if(itr==MAX_ITERATION){
				particles_over_max_iteration++;
			}else{
				for(int D=0;D<4;D++) PART[i].ang[D]=qi[D];
//				std::cout<<"PART["<<i<<"], finished, iteration: "<<itr<<", matrix converged"<<std::endl;
//				cout<<"q convergence: "<<boolalpha<<check_q_residue(qi, EPSILON)<<"; ";
//				for(int m=0;m<4;m++) cout<<"q["<<m<<"]="<<qi[m]<<" ";
//				cout<<"iteration: "<<itr<<endl;
			}
		}
	}

	cout<<"particles_over_max_iteration: "<<particles_over_max_iteration<<endl;
	for(int D=0;D<4;D++) delete [] jacobi_matrix[D];
	delete [] jacobi_matrix;
}

void calc_velocity_and_position(vector<mpselastic> &PART, elastic &ELAST, double **F)
{
	//F[D][i]これは粒子一つ一つに対してx, y, z方向それぞれの成分を持つ(節点力法を利用すべき) //Fの単位は[N]!!!
	int dimension=ELAST.get_dimension();
	double dt=ELAST.get_dt();
	double density=ELAST.get_density();
	double mass=ELAST.get_mass();
	double g[3]={0.0, 0.0, 0.0};
	int symplectic=ELAST.get_symp_flag();
	int symplectic_order=ELAST.get_symp_order();
	
	if(dimension==2) g[A_Y]=ELAST.get_g(); 
	if(dimension==3) g[A_Z]=ELAST.get_g();

//	cout<<"mass= "<<mass<<endl;
//	cout<<"density= "<<density<<endl;
//	cout<<"PART[2].density= "<<PART[2].get_density()<<endl;

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			if(dimension==3)
			{
				//3Dではシンプレクティック条件を満たさないのに適用していた！オイラー法で計算すると・・・
				if(symplectic==OFF)
				{
					for(int D=0;D<dimension;D++){
//							double u_temp=PART[i].u[D];
//							PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+ELAST.get_P_visco_stress(D, i)+g[D]+F[D][i]/mass);
//							PART[i].u[D]+=dt*(PART[i].get_normal(D)+PART[i].get_shear(D)+PART[i].get_pressure(D)+PART[i].get_normal_visco(D)+PART[i].get_shear_visco(D)+g[D]+F[D][i]/mass);
					}
				}
				else
				{
//					check_velocity_and_position(PART, i, mass, g, F);

					for(int D=0;D<dimension;D++){
						PART[i].u[D]+=dt*(PART[i].get_stress_accel(D)+PART[i].get_pressure_accel(D)+PART[i].get_stress_visco_accel(D)+g[D]+F[D][i]/mass);
						PART[i].r[D]+=dt*PART[i].u[D];
					}
				}
			}
		}
	}
}

void calc_contact(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double le=ELAST.get_le();
	double re=ELAST.get_r();
	double density=ELAST.get_density();

	vector<int>::iterator ip;
	map<int, double>::iterator mp; //map探索用のiterator

	//mainのreload_INDEX2()実行しておかないとPART[i]の情報がめちゃくちゃになる

	//膨張・収縮加速度と接触反力の計算
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
//		countOK=countNOT=0;
		size_t neighbourN=PART[i].get_current_neighboursID().size();//現在位置での周辺粒子数を取得
		double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};//現在位置でのiのクォータニオンを取得

		//現在re内に存在しなければw=0なので結局0。findで探せば良い
		for(int k=0;k<neighbourN;k++)
		{
			int j=PART[i].get_current_neighboursID()[k];//現在周辺にある粒子のIDを取得

			//現在配置での相対座標
			double r_ij[3], r_ji[3];
			for(int D=0;D<3;D++)
			{
				r_ij[D]=PART[j].r[D]-PART[i].r[D];
				r_ji[D]=-r_ij[D];
			}

			double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離

			//初期配置でre内にあった粒子のIDを探す
			ip=find(PART[i].get_initial_neighboursID().begin(), PART[i].get_initial_neighboursID().end(), j);
			
			//粒子jが初期配置でre内にある場合
			if(ip!=PART[i].get_initial_neighboursID().end())//見つかった場合（現在周辺には初期配置のIDがある）
			{
				//iからみた初期配置での相対座標
				//iがELASTICかそうでないかで場合分け
				if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
				{
					if(PART[j].type==ELASTIC)//最初から材質で分けると弾性体どうしの衝突に対応できなくなる
					{
						//弾性体の内力では重み付けにdis0を使うべき（初期配置からの変形が重要なので）
						double r_ij_Init[3], r_ji_Init[3]; 
						for(int D=0;D<3;D++)
						{
							r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
							r_ji_Init[D]=-r_ij_Init[D];
						}

						//回転行列Riによるr_ij_Initの回転・・・以降はr_ij_zeroを使う
						//qiは時々刻々変わるので使い回しできない・・・PARTのパラメータとして与えるのが良い
						double r_ij_zero[3], r_ji_zero[3];
						double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

						rotate_r0(r_ij_Init, qi, r_ij_zero);
						rotate_r0(r_ji_Init, qj, r_ji_zero);

						double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離

						double w=kernel(re, dis0);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//重み関数の加重平均のとり方を変更。iの粒子数密度で計算する
							press_accel[D]=(PART[i].P*r_ij_zero[D]-PART[j].P*r_ji_zero[D])*w/dis0/dis0;
							//重みは初期の距離を使うが勾配は現在配置の物を使う
	//						press_accel[D]=((PART[i].P*r_ij[D]/PART[i].PND)-(PART[j].P*r_ji[D]/PART[j].PND))*w/dis/dis;
						//dis0^2で除しているということは，内積を考えている可能性がある．

						PART[i].add_pressure(press_accel);

					}else if(PART[j].type==WALL){
					
						//変形しない物体に対しては重み付けにdisを使うべき（初期配置からの変形は生じず抵抗力の向きだけが重要なので）
						//境界条件はステップごとに変わる
						//接触する場合距離が近づくので

						double w=kernel(re, dis);

						double press_accel[3];
						for(int D=0;D<3;D++)
							press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;

						PART[i].add_pressure(press_accel);//圧力の加速度を加える
					}				
				}
				else if(PART[i].type==WALL)
				{
					if(PART[j].type==ELASTIC)//最初から材質で分けると弾性体どうしの衝突に対応できなくなる
					{
						double w=kernel(re, dis);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//重み関数の加重平均のとり方を変更。iの粒子数密度で計算する
							press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/PART[i].PND/dis/dis;
							//重みは初期の距離を使うが勾配は現在配置の物を使う
	//						press_accel[D]=((PART[i].P*r_ij[D]/PART[i].PND)-(PART[j].P*r_ji[D]/PART[j].PND))*w/dis/dis;
						//dis0^2で除しているということは，内積を考えている可能性がある．

						PART[i].add_pressure(press_accel);

					}		
				}
			}
			else//初期配置でre内に粒子jがない場合（周辺に粒子が接近してきたor離れていった）
			{
				//粒子iがELASTICかWALLかで分ける
				if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
				{
					if(PART[j].type==ELASTIC)//最初から材質で分けると弾性体どうしの衝突に対応できなくなる
					{
						double r_ij_Init[3], r_ji_Init[3]; 
						for(int D=0;D<3;D++)
						{
							r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
							r_ji_Init[D]=-r_ij_Init[D];
						}

						double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離

						//回転行列Riによるr_ij_Initの回転・・・以降はr_ij_zeroを使う
						//qiは時々刻々変わるので使い回しできない・・・PARTのパラメータとして与えるのが良い
						double r_ij_zero[3], r_ji_zero[3];
						double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

						rotate_r0(r_ij_Init, qi, r_ij_zero);
						rotate_r0(r_ji_Init, qj, r_ji_zero);
			
						double w=kernel(re, dis0);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//圧力勾配を加速度として計算していることに相当する（粒子方向の単位ベクトル向きに勾配ベクトルがかかる）
	//						press_accel[D]=((PART[i].P*r_ij_zero[D]/PART[i].PND)-(PART[j].P*r_ji_zero[D]/PART[j].PND))*w/dis0/dis0;
							press_accel[D]=(PART[i].P*r_ij_zero[D]-PART[j].P*r_ji_zero[D])*w/dis0/dis0;

						PART[i].add_pressure(press_accel);

					}else if(PART[j].type==WALL){//上と同様

						double w=kernel(re, dis);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//反力は膨張と逆向きにかかるはず・・・
							press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;

						PART[i].add_pressure(press_accel);
					}
				}
				else if(PART[i].type==WALL)
				{
					if(PART[j].type==ELASTIC)//最初から材質で分けると弾性体どうしの衝突に対応できなくなる
					{
						double w=kernel(re, dis);

						double press_accel[3];
						for(int D=0;D<3;D++)
							//圧力勾配を加速度として計算していることに相当する（粒子方向の単位ベクトル向きに勾配ベクトルがかかる）
	//						press_accel[D]=((PART[i].P*r_ij_zero[D]/PART[i].PND)-(PART[j].P*r_ji_zero[D]/PART[j].PND))*w/dis0/dis0;
							press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;

						PART[i].add_pressure(press_accel);

					}
				}
			}

		}//for(int k=0;k<neighbourN;k++)・・・kループ終了

		double coef=-dimension/PART[i].get_density()/PART[i].PND;
		PART[i].mul_pressure(coef);
//		ELAST.set_P_visco_stress(D, i, (ELAST.get_P_visco_stress(D, i)*dimension*(-1)/ePND[i]/density));//n0でなくn[i]を使っている・・・
	}
}


void calc_accel_for_3D_ver_3(vector<mpselastic> &PART, elastic &ELAST)
{
//	cout<<"3D analysis start"<<endl;
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double mag_lambda=ELAST.get_mag_lambda();
	double elas_lambda=ELAST.get_elas_lambda();
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //これは既にre*le
	double vis=ELAST.get_nensei(); //粘度（動粘度ではない）
	double g=ELAST.get_g();

	double KE=0.0;	//運動エネルギー
	double EE1=0.0;	//ひずみエネルギー
	double EE2=0.0;	//体積ひずみによる
	double PE=0.0;	//ポテンシャル
	double ground=ELAST.get_ground_position(); //床のｚ座標の取得

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
		//WALLは弾性変形しないので除外
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			double EE1_temp=0.0;

			//運動エネルギーの更新
			for(int D=0;D<3;D++) KE+=PART[i].u[D]*PART[i].u[D];

			//位置ポテンシャルの更新
			PE+=(PART[i].r[2]-ground);

			int neighboursN0=static_cast<int>(PART[i].get_initial_neighboursID().size());

			for(int k=0;k<neighboursN0;k++) //pressureと合わせないと狂う？
			{
				int j=PART[i].get_initial_neighboursID()[k];

				//初期配置でNONELASTICなものが近くにある場合は条件を分ける
				//周辺粒子はWALLでも良いのでは・・・このifいらない・・・(2012-11-29)
//				if(PART[j].type==ELASTIC)
				{				
					//iからみた初期配置での相対座標
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++){
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					//現在配置での相対座標
					double r_ij[3], r_ji[3];
					for(int D=0;D<3;D++)
					{
						r_ij[D]=PART[j].r[D]-PART[i].r[D];
						r_ji[D]=-r_ij[D];
					}

					double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離
					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離
			
					double r_ij_zero[3], r_ji_zero[3];

					//calc_r0_ij(PART)で"現在位置のクォータニオンを使って"回転させた初期位置ベクトルを取得
					for(int D=0;D<3;D++) r_ij_zero[D]=PART[i].get_r0_ij()[k].get_comp()[D];
					for(int D=0;D<3;D++) r_ji_zero[D]=PART[i].get_r0_ji()[k].get_comp()[D];

					double w=kernel(re, dis0);

					//回転後の初期配置相対位置の単位ベクトル
					double n_ij[3], n_ji[3];
					for(int D=0;D<3;D++) 
					{
						n_ij[D]=r_ij_zero[D]/dis0;
						n_ji[D]=r_ji_zero[D]/dis0;
					}

					//変位ベクトル
					double U_ij[3], U_ji[3]; 
					for(int D=0;D<3;D++)
					{
						U_ij[D]=r_ij[D]-r_ij_zero[D];
						U_ji[D]=r_ji[D]-r_ji_zero[D];
					}

					//ひずみベクトル
					double E_ij[3], E_ji[3];
					for(int D=0;D<3;D++)
					{
						E_ij[D]=U_ij[D]/dis0;
						E_ji[D]=U_ji[D]/dis0;
					}
				
					//圧力計算の準備・・・e_vol_iのΣ[]（変位の発散）を計算
					//垂直ひずみの和を取る（応力テンソルの対角成分の和に重み付け）
					for(int D=0;D<3;D++) PART[i].P+=E_ij[D]*n_ij[D]*w; 

					//ひずみエネルギー
					for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;

					//応力ベクトル
					double sigma_ij[3], sigma_ji[3];
					for(int D=0;D<3;D++)
					{
						if(PART[i].type==MAGELAST){
						sigma_ij[D]=2*mag_shear_modulus*w*E_ij[D]/PART[i].PND/dis0;
						sigma_ji[D]=2*mag_shear_modulus*w*E_ji[D]/PART[j].PND/dis0;
						}
						else if(PART[i].type==ELASTIC){
						sigma_ij[D]=2*elas_shear_modulus*w*E_ij[D]/PART[i].PND/dis0;
						sigma_ji[D]=2*elas_shear_modulus*w*E_ji[D]/PART[j].PND/dis0;
						}

						sigma_ij[D]-=sigma_ji[D];
					}
					PART[i].add_stress_accel(sigma_ij);

	//				if(fabs((dis-dis0)/dis0)>1.04) w=0;//破壊条件

					//ひずみ速度の計算
					double strain_vi[3], strain_vj[3];

					//iからみたひずみ速度
					strain_vi[0]=((PART[j].u[0]-PART[i].u[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
					strain_vi[1]=((PART[j].u[1]-PART[i].u[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y軸周りのひずみ速度
					strain_vi[2]=((PART[j].u[2]-PART[i].u[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z軸周りのひずみ速度
					//anglar_u1..3は配列に書き換え

					//jからみたひずみ速度
					strain_vj[0]=((PART[i].u[0]-PART[j].u[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
					strain_vj[1]=((PART[i].u[1]-PART[j].u[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y軸周りのひずみ速度
					strain_vj[2]=((PART[i].u[2]-PART[j].u[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z軸周りのひずみ速度
			
					double sigma_v_ij[3], sigma_v_ji[3];
					for(int D=0;D<3;D++)
					{
						sigma_v_ij[D]=2*vis*w*strain_vi[D]/PART[i].PND/dis;
						sigma_v_ji[D]=2*vis*w*strain_vj[D]/PART[j].PND/dis;

						sigma_v_ij[D]-=sigma_v_ji[D];
					}

					PART[i].add_stress_visco_accel(sigma_v_ij);
				}//if(PART[j].type==ELASTIC)終了・・・ここでWALLからのせん断
			}//for(int k=0;k<neighboursN0;k++)ループ終了

			PART[i].P*=-dimension/PART[i].PND;//体積ひずみが求められた（これはエネルギー計算する前に求めておくこと）

			//接触判定（下壁も含む）・・・ここで減衰を入れないとエラー？
			//加速度計算はELASTICだけ考えればよいのでif(.type==ELASTIC)のループの中に入れても変わらない
			double n00=PART[i].PND0;
//			if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=(PART[i].PND-n00)/n00; //粒子数密度の増加が見られれば置換する

			double coef=dimension/PART[i].get_density();

			PART[i].mul_stress_accel(coef);
			PART[i].mul_stress_visco_accel(coef);

			//ひずみエネルギー
//			EE1+=(EE1_temp)*(shear_modulus*mass*dimension/PART[i].PND/density);	//第二項のΣ
//			EE2+=0.5*lambda*mass*(PART[i].P*PART[i].P)/density;//第三項のΣ
			if(PART[i].type==MAGELAST){
			EE1+=(EE1_temp/PART[i].PND/PART[i].get_density())*mag_shear_modulus*mass*dimension;//第二項のΣ 密度はそれぞれ違う
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*mag_lambda;//第三項のΣ 密度はそれぞれ違う

			PART[i].P*=mag_lambda;//体積ひずみによる圧力が求められた
			}
			else if(PART[i].type==ELASTIC){
			EE1+=(EE1_temp/PART[i].PND/PART[i].get_density())*elas_shear_modulus*mass*dimension;//第二項のΣ 密度はそれぞれ違う
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*elas_lambda;//第三項のΣ 密度はそれぞれ違う

			PART[i].P*=elas_lambda;//体積ひずみによる圧力が求められた
			}
		}//if(PART[i].type==ELASTIC)終了

		//それ自体は動かないがWALLも加速度（反力）を有している
		//・・・PART[i].Pはすべての粒子を考慮する必要があるが、エネルギー計算に入れるべきではない
		double n00=PART[i].PND0;
//		if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=lambda*(PART[i].PND-n00)/n00; //粒子数密度の増加が見られれば置換する
		//2012-11-26 圧力と粒子数密度を比較していた・・・！
		if(PART[i].type==MAGELAST){
			if(PART[i].P<mag_lambda*((PART[i].PND-n00)/n00)) PART[i].P=mag_lambda*(PART[i].PND-n00)/n00; //粒子数密度の増加が見られれば置換する
		}
		else if(PART[i].type==ELASTIC){
			if(PART[i].P<elas_lambda*((PART[i].PND-n00)/n00)) PART[i].P=elas_lambda*(PART[i].PND-n00)/n00;
		}

	}//for(i=0;i<PART.size();i++)ループ終了

	KE*=0.5*mass;
	PE*=-g*mass;
//	EE1*=shear_modulus*mass*dimension;
//	EE2*=0.5*mass*lambda;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//全系の弾性エネルギー
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);
}

void calc_accel_for_3D_ver_2(vector<mpselastic> &PART, elastic &ELAST)
{
//	cout<<"3D analysis start"<<endl;
	double dimension=static_cast<double>(ELAST.get_dimension());
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double mag_lambda=ELAST.get_mag_lambda();
	double elas_lambda=ELAST.get_elas_lambda();
	double mass=ELAST.get_mass();
	double density=ELAST.get_density();
	double re=ELAST.get_r(); //これは既にre*le
	double vis=ELAST.get_nensei(); //粘度（動粘度ではない）
	double g=ELAST.get_g();

	double KE=0.0;	//運動エネルギー
	double EE1=0.0;	//ひずみエネルギー
	double EE2=0.0;	//体積ひずみによる
	double PE=0.0;	//ポテンシャル
	double ground=ELAST.get_ground_position(); //床のｚ座標の取得

	for(int i=0;i<PART.size();i++)
	{   
		//WALLは弾性変形しないので除外
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{
			double EE1_temp=0.0;

			//運動エネルギーの更新
			for(int D=0;D<3;D++) KE+=PART[i].u[D]*PART[i].u[D];

			//位置ポテンシャルの更新
			PE+=(PART[i].r[2]-ground);

			double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]}; //iのクォータニオン
			size_t neighboursN0=PART[i].get_initial_neighboursID().size();

			for(int k=0;k<neighboursN0;k++) //pressureと合わせないと狂う？
			{
				int j=PART[i].get_initial_neighboursID()[k];

				double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]}; //jのクォータニオン

//				PART[i].get_initial_neighbours_position(k, X0, Y0, Z0);　ムダもいいところ
//				double dis0=PART[i].get_initial_distancebps()[k];//初期粒子間距離

				//iからみた初期配置での相対座標
				double r_ij_Init[3], r_ji_Init[3]; 
				for(int D=0;D<3;D++){
					r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
					r_ji_Init[D]=-r_ij_Init[D];
				}

				//現在配置での相対座標
				double r_ij[3], r_ji[3];
				for(int D=0;D<3;D++)
				{
					r_ij[D]=PART[j].r[D]-PART[i].r[D];
					r_ji[D]=-r_ij[D];
				}

				double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離
				double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離
			
				//回転行列Riによるr_ij_Initの回転・・・以降はr_ij_zeroを使う
				//qiは時々刻々変わるので使い回しできない・・・PARTのパラメータとして与えるのが良い
				double r_ij_zero[3], r_ji_zero[3];
				rotate_r0(r_ij_Init, qi, r_ij_zero);
				rotate_r0(r_ji_Init, qj, r_ji_zero);

				double w=kernel(re, dis0);

				//回転後の初期配置相対位置の単位ベクトル
				double n_ij[3], n_ji[3];
				for(int D=0;D<3;D++) 
				{
					n_ij[D]=r_ij_zero[D]/dis0;
					n_ji[D]=r_ji_zero[D]/dis0;
				}

				//変位ベクトル
				double U_ij[3], U_ji[3]; 
				for(int D=0;D<3;D++)
				{
					U_ij[D]=r_ij[D]-r_ij_zero[D];
					U_ji[D]=r_ji[D]-r_ji_zero[D];
				}

				//ひずみベクトル
				double E_ij[3], E_ji[3];
				for(int D=0;D<3;D++)
				{
					E_ij[D]=U_ij[D]/dis0;
					E_ji[D]=U_ji[D]/dis0;
				}
				
				//圧力計算の準備・・・e_vol_iのΣ[]（変位の発散）を計算
				for(int D=0;D<3;D++) PART[i].P+=E_ij[D]*n_ij[D]*w; 
//				PART[i].P*=w;

				//ひずみエネルギー
				for(int D=0;D<3;D++) EE1_temp+=E_ij[D]*E_ij[D]*w;
//				EE1_temp*=w;

				//応力ベクトル
				double sigma_ij[3], sigma_ji[3];
				for(int D=0;D<3;D++)
				{
					if(PART[i].type==MAGELAST){
					sigma_ij[D]=2*mag_shear_modulus*w*E_ij[D]/PART[i].PND/dis0;
					sigma_ji[D]=2*mag_shear_modulus*w*E_ji[D]/PART[i].PND/dis0;
					}
					else if(PART[i].type==ELASTIC){
					sigma_ji[D]=2*elas_shear_modulus*w*E_ji[D]/PART[j].PND/dis0;
					sigma_ij[D]=2*elas_shear_modulus*w*E_ij[D]/PART[j].PND/dis0;
					}

					sigma_ij[D]-=sigma_ji[D];
				}
				PART[i].add_stress_accel(sigma_ij);

//				if(fabs((dis-dis0)/dis0)>1.04) w=0;//破壊条件

				//ひずみ速度の計算
				double strain_vi[3], strain_vj[3];

				//iからみたひずみ速度
				strain_vi[0]=((PART[j].u[0]-PART[i].u[0])-(PART[i].ang_u[1]*r_ij[2]-PART[i].ang_u[2]*r_ij[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
				strain_vi[1]=((PART[j].u[1]-PART[i].u[1])-(PART[i].ang_u[2]*r_ij[0]-PART[i].ang_u[0]*r_ij[2]))/dis; //Y軸周りのひずみ速度
				strain_vi[2]=((PART[j].u[2]-PART[i].u[2])-(PART[i].ang_u[0]*r_ij[1]-PART[i].ang_u[1]*r_ij[0]))/dis; //Z軸周りのひずみ速度
				//anglar_u1..3は配列に書き換え

				//jからみたひずみ速度
				strain_vj[0]=((PART[i].u[0]-PART[j].u[0])-(PART[j].ang_u[1]*r_ji[2]-PART[j].ang_u[2]*r_ji[1]))/dis; //X軸周りのひずみ速度 これも正射影する必要あり？？
				strain_vj[1]=((PART[i].u[1]-PART[j].u[1])-(PART[j].ang_u[2]*r_ji[0]-PART[j].ang_u[0]*r_ji[2]))/dis; //Y軸周りのひずみ速度
				strain_vj[2]=((PART[i].u[2]-PART[j].u[2])-(PART[j].ang_u[0]*r_ji[1]-PART[j].ang_u[1]*r_ji[0]))/dis; //Z軸周りのひずみ速度
			
				double sigma_v_ij[3], sigma_v_ji[3];
				for(int D=0;D<3;D++)
				{
					sigma_v_ij[D]=2*vis*w*strain_vi[D]/PART[i].PND/dis;
					sigma_v_ji[D]=2*vis*w*strain_vj[D]/PART[j].PND/dis;

					sigma_v_ij[D]-=sigma_v_ji[D];
				}

				PART[i].add_stress_visco_accel(sigma_v_ij);
			}//for(int k=0;k<neighboursN0;k++)ループ終了

			PART[i].P*=-dimension/PART[i].PND;//体積ひずみが求められた（これはエネルギー計算する前に求めておくこと）

			//接触判定（下壁も含む）・・・ここで減衰を入れないとエラー？
			//加速度計算はELASTICだけ考えればよいのでif(.type==ELASTIC)のループの中に入れても変わらない
			double n00=PART[i].PND0;
			if(PART[i].P<(PART[i].PND-n00)/n00) PART[i].P=(PART[i].PND-n00)/n00; //粒子数密度の増加が見られれば置換する

			//ELASTICでないものも計算していた
			double coef=dimension/PART[i].get_density();

			PART[i].mul_stress_accel(coef);
			PART[i].mul_stress_visco_accel(coef);

			//ひずみエネルギー
//			EE1+=(EE1_temp)*(shear_modulus*mass*dimension/PART[i].PND/density);	//第二項のΣ
//			EE2+=0.5*lambda*mass*(PART[i].P*PART[i].P)/density;//第三項のΣ
			
			
			if(PART[i].type==MAGELAST){
				EE1+=(EE1_temp/PART[i].PND/PART[i].get_density())*mag_shear_modulus*mass*dimension;//第二項のΣ 密度はそれぞれ違う
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*mag_lambda;//第三項のΣ 密度はそれぞれ違う
			PART[i].P*=mag_lambda;//体積ひずみによる圧力が求められた
			}
			else if(PART[i].type==ELASTIC){
				EE1+=(EE1_temp/PART[i].PND/PART[i].get_density())*elas_shear_modulus*mass*dimension;//第二項のΣ 密度はそれぞれ違う
			EE2+=((PART[i].P*PART[i].P)/PART[i].get_density())*0.5*mass*elas_lambda;//第三項のΣ 密度はそれぞれ違う
			PART[i].P*=elas_lambda;
			}
		}//if(PART[i].type==ELASTIC || PART[i].type==INELASTIC ||PART[i].type==BOELASTIC)終了

		//それ自体は動かないがWALLも加速度（反力）を有している
		//・・・PART[i].Pはすべての粒子を考慮する必要があるが、エネルギー計算に入れるるべきではない
		//

	}//for(i=0;i<PART.size();i++)ループ終了

	KE*=0.5*mass;
	PE*=-g*mass;
//	EE1*=shear_modulus*mass*dimension;
//	EE2*=0.5*mass*lambda;
	ELAST.set_kinetic(KE);
	ELAST.set_elastic_energy(EE1+EE2);//全系の弾性エネルギー
	ELAST.set_elastic_energy1(EE1);
	ELAST.set_elastic_energy2(EE2);
	ELAST.set_potential(PE);
	ELAST.set_hamiltonian(KE+EE1+EE2+PE);
}


void calc_pressure_and_contact_ver_3(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double le=ELAST.get_le();
	double re=ELAST.get_r();
	double density=ELAST.get_density();

	//mainのreload_INDEX2()実行しておかないとPART[i]の情報がめちゃくちゃになる

	//膨張・収縮加速度と接触反力の計算
	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{
		//初期配置のもの＋
		//(PART[i].PND>PART[i].PND0)となった粒子を計算する
		size_t neighbourN0=PART[i].get_initial_neighboursID().size();//現在位置での周辺粒子数を取得
		double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};//現在位置でのiのクォータニオンを取得

		//初期配置にあるかどうかは気にしない
		for(int k=0;k<neighbourN0;k++)
		{
			int j=PART[i].get_current_neighboursID()[k];//現在周辺にある粒子のIDを取得

			//現在配置での相対座標
			double r_ij[3], r_ji[3];
			for(int D=0;D<3;D++)
			{
				r_ij[D]=PART[j].r[D]-PART[i].r[D];
				r_ji[D]=-r_ij[D];
			}

			double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離

			//iからみた初期配置での相対座標
			//iがELASTICかそうでないかで場合分け
			if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
			{
				if(PART[j].type==ELASTIC)//最初から材質で分けると弾性体どうしの衝突に対応できなくなる
				{
					//弾性体の内力では重み付けにdis0を使うべき（初期配置からの変形が重要なので）
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++)
					{
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離

					double w=kernel(re, dis0);

					double press_accel[3];
					for(int D=0;D<3;D++){
						//重み関数の加重平均のとり方を変更。iの粒子数密度で計算する
						//重みは初期の距離を使うが勾配は現在配置の物を使う
						//press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis0/dis0;・・・勾配にdis0を使ってはいけない！！！
						press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;
					}

					PART[i].add_pressure(press_accel);

				}else if(PART[j].type==WALL){
					
					//変形しない物体に対しては重み付けにdisを使うべき（初期配置からの変形は生じず抵抗力の向きだけが重要なので）
					//接触する場合距離が近づくので境界条件はステップごとに変わる

					double w=kernel(re, dis);

					double press_accel[3];
					for(int D=0;D<3;D++){
						press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;
					}

					PART[i].add_pressure(press_accel);//圧力の加速度を加える
				}				
			}
			else if(PART[i].type==WALL)
			{
				if(PART[j].type==ELASTIC)//最初から材質で分けると弾性体どうしの衝突に対応できなくなる
				{
					double w=kernel(re, dis);

					double press_accel[3];
					for(int D=0;D<3;D++){
						//重み関数の加重平均のとり方を変更。iの粒子数密度で計算する
						//重みは初期の距離を使うが勾配は現在配置の物を使う
						press_accel[D]=(PART[i].P*r_ij[D]-PART[j].P*r_ji[D])*w/dis/dis;
					}
					PART[i].add_pressure(press_accel);
				}		
			}
		}//for(int k=0;k<neighbourN;k++)・・・kループ終了



		double coef=-dimension/PART[i].get_density()/PART[i].PND;
		PART[i].mul_pressure(coef);
//		ELAST.set_P_visco_stress(D, i, (ELAST.get_P_visco_stress(D, i)*dimension*(-1)/ePND[i]/density));//n0でなくn[i]を使っている・・・
	}//	for(int i=0;i<PART.size();i++)
}

void calc_pressure_and_contact(vector<mpselastic> &PART, elastic &ELAST)
{
	double dimension=static_cast<double>(ELAST.get_dimension());
	double le=ELAST.get_le();
	double re=ELAST.get_r();
	double density=ELAST.get_density();

	vector<int>::iterator ip;
	map<int, double>::iterator mp; //map探索用のiterator

	//mainのreload_INDEX2()実行しておかないとPART[i]の情報がめちゃくちゃになる

	#pragma omp parallel for
	for(int i=0;i<PART.size();i++)
	{   
//		countOK=countNOT=0;
		size_t neighbourN=PART[i].get_current_neighboursID().size();//現在位置での周辺粒子数を取得
		double qi[4]={PART[i].ang[0], PART[i].ang[1], PART[i].ang[2], PART[i].ang[3]};//現在位置でのiのクォータニオンを取得

		//現在re内に存在しなければw=0なので結局0。findで探せば良い
		for(int k=0;k<neighbourN;k++)
		{
			int j=PART[i].get_current_neighboursID()[k];//現在周辺にある粒子のIDを取得

			//現在配置での相対座標
			double r_ij[3], r_ji[3];
			for(int D=0;D<3;D++)
			{
				r_ij[D]=PART[j].r[D]-PART[i].r[D];
				r_ji[D]=-r_ij[D];
			}

			double dis=0.0; for(int D=0;D<3;D++){dis+=r_ij[D]*r_ij[D];} dis=sqrt(dis); //現在粒子間距離

			//初期配置でre内にあった粒子のIDを探す
			ip=find(PART[i].get_initial_neighboursID().begin(), PART[i].get_initial_neighboursID().end(), j);
			
			//初期配置でre内に粒子jがある場合
			if(ip!=PART[i].get_initial_neighboursID().end())//見つかった場合（現在周辺には初期配置のIDがある）
			{
				//iからみた初期配置での相対座標
				if(PART[j].type==ELASTIC)//最初から材質で分けると弾性体どうしの衝突に対応できなくなる
				{
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++)
					{
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					//回転行列Riによるr_ij_Initの回転・・・以降はr_ij_zeroを使う
					//qiは時々刻々変わるので使い回しできない・・・PARTのパラメータとして与えるのが良い
					double r_ij_zero[3], r_ji_zero[3];
					double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

					rotate_r0(r_ij_Init, qi, r_ij_zero);
					rotate_r0(r_ji_Init, qj, r_ji_zero);

					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離
			
					double w=kernel(re, dis0);

					double press_accel[3];
					for(int D=0;D<3;D++)
						press_accel[D]=((PART[i].P*r_ij_zero[D]/PART[i].PND)-(PART[j].P*r_ji_zero[D]/PART[j].PND))*w/dis0/dis0;
					//dis0^2で除しているということは，内積を考えている可能性がある．

					PART[i].add_pressure(press_accel);

				}else if(PART[j].type==WALL){

					double w=kernel(re, dis);

					double press_accel[3];
					for(int D=0;D<3;D++)
						press_accel[D]=((PART[i].P*r_ij[D]/PART[i].PND)-(PART[j].P*r_ji[D]/PART[j].PND))*w/dis/dis;

					PART[i].add_pressure(press_accel);//圧力の加速度を加える
				}
			}
			else//初期配置でre内に粒子jがない場合
			{
				if(PART[j].type==ELASTIC)//最初から材質で分けると弾性体どうしの衝突に対応できなくなる
				{
					double r_ij_Init[3], r_ji_Init[3]; 
					for(int D=0;D<3;D++)
					{
						r_ij_Init[D]=PART[j].initial_r[D]-PART[i].initial_r[D];
						r_ji_Init[D]=-r_ij_Init[D];
					}

					double dis0=0.0; for(int D=0;D<3;D++){dis0+=r_ij_Init[D]*r_ij_Init[D];} dis0=sqrt(dis0); //初期粒子間距離

					//回転行列Riによるr_ij_Initの回転・・・以降はr_ij_zeroを使う
					//qiは時々刻々変わるので使い回しできない・・・PARTのパラメータとして与えるのが良い
					double r_ij_zero[3], r_ji_zero[3];
					double qj[4]={PART[j].ang[0], PART[j].ang[1], PART[j].ang[2], PART[j].ang[3]};

					rotate_r0(r_ij_Init, qi, r_ij_zero);
					rotate_r0(r_ji_Init, qj, r_ji_zero);
			
					double w=kernel(re, dis0);

					double press_accel[3];
					for(int D=0;D<3;D++)
						press_accel[D]=((PART[i].P*r_ij_zero[D]/PART[i].PND)-(PART[j].P*r_ji_zero[D]/PART[j].PND))*w/dis0/dis0;

					PART[i].add_pressure(press_accel);

				}else if(PART[j].type==WALL){

					double w=kernel(re, dis);

					double press_accel[3];
					for(int D=0;D<3;D++)
						press_accel[D]=((PART[i].P*r_ij[D]/PART[i].PND)-(PART[j].P*r_ji[D]/PART[j].PND))*w/dis/dis;

					PART[i].add_pressure(press_accel);
				}
			}

		}//for(int k=0;k<neighbourN;k++)・・・kループ終了

		double coef=-dimension/PART[i].get_density();
		PART[i].mul_pressure(coef);
	}
}

//角度計算 これを更新しないとおかしなことが起こる・・・
void calc_angular_velocity_and_angle(vector<mpselastic> &PART, elastic &ELAST)
{
	unsigned dimension=ELAST.get_dimension();
	double dt=ELAST.get_dt();

	//角度
	if(dimension==2)
	{
		for(int i=0;i<PART.size();i++) PART[i].ang[0]+=dt*PART[i].ang_u[0];
	}
	if(dimension==3)
	{	
		for(int i=0;i<PART.size();i++){
			double omega_norm=sqrt(PART[i].ang_u[0]*PART[i].ang_u[0]+PART[i].ang_u[1]*PART[i].ang_u[1]+PART[i].ang_u[2]*PART[i].ang_u[2]);

			double v[3]={0.0, 0.0, 0.0};//回転軸ベクトル

			if(omega_norm!=0){
				for(int D=0;D<DIMENSION; D++) v[D]=PART[i].ang_u[D]/omega_norm;
			}
			//クォータニオン修正量
			double theta=dt*omega_norm;//(4.33)
			
			double q[4];
			for(int D=0;D<DIMENSION;D++) q[D]=v[D]*sin(theta/2.0);
			q[3]=cos(theta/2.0);
		    
			///角度修正
			double qi[4];
			for(int D=0;D<4;D++) qi[D]=PART[i].ang[D];//粒子iのクォータニオン
			
			//OK??
			PART[i].ang[0]=q[3]*qi[0]+q[0]*qi[3]+q[1]*qi[2]-q[2]*qi[1];
			PART[i].ang[1]=q[3]*qi[1]+q[1]*qi[3]+q[2]*qi[0]-q[0]*qi[2];
			PART[i].ang[2]=q[3]*qi[2]+q[2]*qi[3]+q[0]*qi[1]-q[1]*qi[0];
			PART[i].ang[3]=q[3]*qi[3]-q[0]*qi[0]-q[1]*qi[1]-q[2]*qi[2];
		}
	}
}

//角速度の更新
void calc_torque_for_2D(vector<mpselastic> &PART, elastic &ELAST)
{
	double mass=ELAST.get_mass();
	double mag_shear_modulus=ELAST.get_mag_shear_modulus();
	double elas_shear_modulus=ELAST.get_elas_shear_modulus();
	double inertia=ELAST.get_inertia();
	double r=ELAST.get_r();
	double dt=ELAST.get_dt();

	double density=ELAST.get_density();
	double dimension=static_cast<double>(ELAST.get_dimension());
	double force[2];

	for(int i=0;i<PART.size();i++)
	{   
		size_t neighboursN0=PART[i].get_initial_neighboursID().size();
		if(PART[i].type==ELASTIC || PART[i].type==MAGELAST)
		{   
			//重要！重み関数は影響半径内の粒子に作用させる！NUMのループは「影響半径内の演算」という意味がある！iとjの相互作用を計算
			for(int k=0;k<neighboursN0;k++) //粒子iと粒子iの影響半径に含まれる粒子の相互作用を計算 NUM[i]=PART[i].N;(初期配置における各粒子の周辺粒子数) i!=jは満たされている！！
			{
				int j=PART[i].get_initial_neighboursID()[k];
				double X0[2];
				double dis0=0.0;
				PART[i].get_initial_neighbours_position(k, X0); X0[A_Z]=0.0;
				dis0=PART[i].get_initial_distancebps()[k];
				
				double X=PART[j].r[A_X]-PART[i].r[A_X];
				double Y=PART[j].r[A_Y]-PART[i].r[A_Y];
				double dis=sqrt(X*X+Y*Y);

				if(dis==0) cout<<"dis=0 "<<i<<" "<<j<<endl;

				double theta=(PART[i].ang[0]+PART[j].ang[0])/2;//相対角度 (4.8)
				double RX0=X0[0]*cos(theta)-X0[1]*sin(theta);//回転行列×r0 初期粒子を回転
				double RY0=X0[0]*sin(theta)+X0[1]*cos(theta);//回転行列×r0

				//相対変位の計算: u_ij=r_ij-R*r0_ij
				double U_X=X-RX0;//変位のX成分 (4.6)
				double U_Y=Y-RY0;//変位のY成分
			
				double orth_project=(U_X*X+U_Y*Y)/(dis*dis);

				double U_n_X=orth_project*X; //変位のr(ij)に平行な"X"成分
				double U_n_Y=orth_project*Y; //変位のr(ij)に平行な"Y"成分

				double U_s_X=U_X-U_n_X; //変位のr(ij)に垂直な"X"成分
				double U_s_Y=U_Y-U_n_Y; //変位のr(ij)に垂直な"Y"成分

				double EsX=U_s_X/dis0; //ひずみのr(ij)に垂直な"X"成分
				double EsY=U_s_Y/dis0; //ひずみのr(ij)に垂直な"Y"成分
			
				double r1=ELAST.get_r();
				double w=kernel(r1,dis0);
				double constantF=0;
				if(PART[i].type==MAGELAST){
				constantF=(4*mass*dimension*mag_shear_modulus*w)/(PART[i].get_density()*PART[i].PND*dis0);
				}
				else if(PART[i].type==ELASTIC){
				constantF=(4*mass*dimension*elas_shear_modulus*w)/(PART[i].get_density()*PART[i].PND*dis0);
				}
				force[0]=constantF*EsX;
				force[1]=constantF*EsY;
				double torque=Y*force[0]-X*force[1];

				PART[i].ang_u[0]+=dt*(-0.5)*torque/inertia;
				PART[j].ang_u[0]+=dt*(-0.5)*torque/inertia;
			}
		}
	}
}

//GNUplot用応力ファイル
void output_stress_for_GNU(vector<mpselastic> &PART, elastic &ELAST, int t)
{
	int dimension=ELAST.get_dimension();
	double le=ELAST.get_distancebp();
	double magni=ELAST.get_times()*le*le;

	if(t==1 || !(t%ELAST.get_interval())){
			char FileName1[128], FileName2[128];

		//ファイル名の準備
			sprintf_s(FileName1, "./stress/normal%04d.dat", t);
			sprintf_s(FileName2, "./stress/shear%04d.dat", t);

			ofstream fout1(FileName1);	//normal
			ofstream fout2(FileName2);	//shear

//			if(dimension==2) for(int i=0;i<PART.size();i++) fout1<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].get_normal(A_X, i)*magni<<" "<<PART[i].get_normal(A_Y, i)*magni<<endl;
//			else if(dimension==3) for(int i=0;i<PART.size();i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fout1<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].get_normal(A_X, i)*magni<<"\t"<<PART[i].get_normal(A_Z, i)*magni<<endl;

//			if(dimension==2) for(int i=0;i<PART.size();i++) fout2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].get_shear(A_X, i)*magni<<" "<<PART[i].get_shear(A_Y, i)*magni<<endl;
//			else if(dimension==3) for(int i=0;i<PART.size();i++) if(PART[i].r[A_Y]<le && PART[i].r[A_Y]>-le) fout2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Z]<<"\t"<<PART[i].get_shear(A_X, i)*magni<<"\t"<<PART[i].get_shear(A_Z, i)*magni<<endl;

			fout1.close();
			fout2.close();
	}
}

void quaternion(double P[], double axis[3], double angle){
	double Q[4]={cos(angle/2),axis[0]*sin(angle/2),axis[1]*sin(angle/2),axis[2]*sin(angle/2)};
	double R[4]={cos(angle/2),-axis[0]*sin(angle/2),-axis[1]*sin(angle/2),-axis[2]*sin(angle/2)};
	double QP[4]={-(Q[1]*R[0]+Q[2]*R[1]+Q[3]*R[2])};
}