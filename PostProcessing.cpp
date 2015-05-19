#include "stdafx.h"
#include "Micro_AVS.h"

void post_processing(mpsconfig &CON, vector<mpselastic> &PART, elastic &ELAST, int fluid_number,int particle_number,double dt,double Umax, int t,double TIME, double **F)
{
	double le=CON.get_distancebp();
	unsigned int timeA=GetTickCount();	//計算開始時間

	cout<<"各物理量開始: ";
	
	//AVSに粒子データ出力
	if(t==1 || t%CON.get_interval()==0) particle_movie_AVS(CON,PART,ELAST,fluid_number,particle_number,t,TIME);
		cout<<"AVS出力完了"<<endl;
		
	//AVS2に磁束密度、ローレンツ力出力
		if((t==1 || t%CON.get_EM_interval()==0) && CON.get_FEM_flag()==ON) particle_movie_AVS2(CON,t,PART,TIME,fluid_number,particle_number,F);
		cout<<"AVS2出力完了\n";

	//速度をプロット
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_speed(CON ,PART,particle_number,fluid_number);
		cout<<"速度プロット完了"<<endl;

	//圧力による加速度をプロット
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_pressure_acceleration(CON, PART);
		cout<<"加速度プロット完了"<<endl;

	//応力分布をプロット
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_stress_distribution(CON, PART);
		cout<<"応力分布プロット完了"<<endl;

	//せん断力をプロット
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_shear_force(CON, PART);
		cout<<"せん断加速度プロット完了"<<endl;

	//ひずみ速度による加速度をプロット
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_strainrate_acceleration(CON, PART);
		cout<<"ひずみ加速度プロット完了"<<endl;

	//ひずみ率をプロット
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_distortion_rate(CON, PART);
		cout<<"ひずみ率プロット完了"<<endl;

	//残差加速度をプロット
	if(CON.get_current_step()==1 || CON.get_current_step()%CON.get_interval()==0) plot_residual_acceleration(CON, PART, F);
		cout<<"残差加速度プロット完了"<<endl;

	//座標プロット
	ofstream gnu1("0.dat");//解析終了後の全粒子座標をプロット
	ofstream gnu2("suf.dat");//表面粒子だけをプロット
	
	if(CON.get_dimension()==2)
	{
		for(int i=0;i<particle_number;i++)
		{ 
			gnu1<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			if(PART[i].surface==ON)
			{
				gnu2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			}
		}
	}
	else if(CON.get_dimension()==3)
	{
		for(int i=0;i<particle_number;i++)
		{ 
			if(PART[i].surface==ON && (PART[i].type==FLUID || PART[i].type==ELASTIC || PART[i].type==MAGELAST||PART[i].type==HYPERELAST))//流体だけ表示したいとき
			{
				gnu2<<PART[i].r[A_X]<<"\t"<<PART[i].r[A_Y]<<"\t"<<PART[i].r[A_Z]<<endl;
			}
		}
	}
	gnu1.close();
	gnu2.close();
	//*/
		
	//平均粒子密度&圧力を表示
	double ave_n0=0;
	double ave_P=0;
	int count=0;
	for(int i=0;i<fluid_number;i++) 
	{
	    if(PART[i].surface==OFF)
	    {
	        ave_n0+=PART[i].PND; //平均粒子数密度
			ave_P+=fabs(PART[i].P); //平均圧力
			count++;
	    }
	}
	if(count!=0){ 
		ave_n0/=count;
		ave_P/=count;
	}

	cout<<"average n0="<<ave_n0<<" average P="<<ave_P<<"/"<<CON.get_ave_P_for_FEM_flag()<<" time="<<(GetTickCount()-timeA)*0.001<<"[sec]"<<endl;

	check_FEM_flag(CON, ELAST, ave_P);
	
	CON.set_current_step(t);
}

void plot_pressure_acceleration(mpsconfig &CON, vector<mpselastic> &PART)
{
	double le=CON.get_distancebp()*0.5;
	double times=1;//CON.get_pressure_times();0.001
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D解析時のpressure.datの出力面 0=YZ平面 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D解析時のpressure.datの出力面の座標
	int d1,d2,d3;							//3D解析時の出力に必要な次元
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	
	stringstream ss;
	ss<<"./Pressure"<<"/pressure_accel"<<CON.get_current_step()<<".dat";
	string filename=ss.str();
	
	ofstream vec(filename);

	if(vec.fail()){
		system("mkdir Pressure");
		ofstream vec(filename);//再試行
		if(vec.fail()){
			cout<<"Pressureディレクトリを作成できませんでした"<<endl;
			exit(1);
		}
	}
	stringstream ss2;
	ss2<<"./Pressure"<<"/pressure3D_"<<CON.get_current_step()<<".dat";
	string filename2=ss2.str();
	
	ofstream ve(filename2);

	if(ve.fail()){
		system("mkdir Pressure");
		ofstream ve(filename2);//再試行
		if(ve.fail()){
			cout<<"Pressureディレクトリを作成できませんでした"<<endl;
			exit(1);
		}
	}
	
	if(dimension==3)
	{
		//int d1,d2;				//出力に必要な次元
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}

		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
			{
				double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].get_pressure_accel(d1);
				double w=PART[i].get_pressure_accel(d2);
				vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				ve<<x<<" "<<z<<" "<<PART[i].P<<endl;
				if(x>xmax) xmax=x;
				if(z>ymax) ymax=z;
			}
		}
	}
	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	if(CON.get_legend_pressure()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_pressure()*times<<" "<<0*times<<endl;//最後に凡例出力
	vec.close();
	ve.close();
}

void plot_stress_distribution(mpsconfig &CON, vector<mpselastic> &PART)
{
	Micro_AVS Micro_avs;
	//構造格子データファイル
	int timestep=CON.get_current_step();
	double den=CON.get_density();

	stringstream sstre;
	sstre<<"./Stress/Stress"<<timestep<<".fld";
	string stForce=sstre.str();

	ofstream fout3(stForce);
	if(fout3.fail()){
		system("mkdir Stress");
		ofstream fout3(stForce);
		if(fout3.fail()){
			cout<<"./Stressフォルダを開けませんでした"<<endl;
			exit(1);
		}
	}

	fout3 << "# AVS field file" << endl;
	fout3 << "ndim=1" << endl;
	fout3 << "dim1=" << PART.size() <<endl;
	fout3 << "nspace=3" << endl;
	fout3 << "veclen=3" << endl;
	fout3 << "data=float" << endl;
	fout3 << "field=irregular" << endl;
	fout3 << "label=e-x e-y e-z" << endl << endl;
	fout3 << "variable 1 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=0 stride=6" << endl;
	fout3 << "variable 2 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=1 stride=6" << endl;
	fout3 << "variable 3 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=2 stride=6" << endl;
	fout3 << "coord    1 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=3 stride=6" << endl;
	fout3 << "coord    2 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=4 stride=6" << endl;
	fout3 << "coord    3 file=./Stress"<<timestep<<" filetype=ascii skip=1 offset=5 stride=6" << endl;
	fout3.close();

	//データファイル
	stringstream ss;
	ss<<"./Stress/Stress"<<timestep;
	string filename=ss.str();

	ofstream fout(filename);
	if(fout.fail()){
		cout<<"./Stressフォルダを開けませんでした"<<endl;
		exit(1);
	}

	fout<<"e-x e-y e-z x y z"<<endl;
	for(size_t i=0;i<PART.size();i++)
    {
		if(PART[i].type!=WALL){
//			Micro_avs.make_list(PART[i].r[A_X], PART[i].r[A_Y], PART[i].r[A_Z], PART[i].get_stress_accel(A_X)*den, PART[i].get_stress_accel(A_Y)*den, PART[i].get_stress_accel(A_Z)*den);
//		fout<<F[A_X][i]*times<<" "<<F[A_Y][i]*times<<" "<<F[A_Z][i]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		fout<<PART[i].get_stress_accel(A_X)*den<<" "<<PART[i].get_stress_accel(A_Y)*den<<" "<<PART[i].get_stress_accel(A_Z)*den<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
		}
	}
//	Micro_avs.Output_vector_MicroAVS("Stress",CON.get_current_step());
	fout.close();
}


void plot_shear_force(mpsconfig &CON, vector<mpselastic> &PART)
{
	Micro_AVS Micro_avs;
	double le=CON.get_distancebp()*0.5;
	double mass=CON.get_particle_mass();
	double times=CON.get_pressure_times();
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D解析時のpressure.datの出力面 0=YZ平面 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D解析時のpressure.datの出力面の座標
	int d1,d2,d3;							//3D解析時の出力に必要な次元
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	
/*	stringstream ss;
	ss<<"./Shear"<<"/shear"<<CON.get_current_step()<<".dat";
	string filename=ss.str();
	
	ofstream vec(filename);

	if(vec.fail()){
		system("mkdir Shear");
		ofstream vec(filename);//再試行
		if(vec.fail()){
			cout<<"Shearディレクトリを作成できませんでした"<<endl;
			exit(1);
		}
	}*/
	
	if(dimension==3)
	{
		//int d1,d2;				//出力に必要な次元
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}
		double u_max=0;
		double w_max=0;
		for(int i=0;i<PART.size();i++)
		{
	//		if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
	//		{
				Micro_avs.make_list(PART[i].r[A_X],PART[i].r[A_Y],PART[i].r[A_Z],PART[i].get_stress_accel(A_X)*mass,PART[i].get_stress_accel(A_Y)*mass,PART[i].get_stress_accel(A_Z)*mass);
	/*			double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].get_stress_accel(d1)*mass;
				double w=PART[i].get_stress_accel(d2)*mass;
				
				if(u_max<u)u_max=u;
				if(w_max<w)w_max=w;*/

	//			vec<<x<<" "<<z<<" "<<u<<" "<<w<<endl;
				
		//	}

		}
	}
	Micro_avs.Output_vector_MicroAVS("Shear",CON.get_current_step());
/*	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	if(CON.get_legend_pressure()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_pressure()*times<<" "<<0*times<<endl;//最後に凡例出力
	vec.close();*/
}

void plot_strainrate_acceleration(mpsconfig &CON, vector<mpselastic> &PART)
{
	Micro_AVS Micro_avs;
	double le=CON.get_distancebp()*0.5;
	double times=CON.get_pressure_times();
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D解析時のpressure.datの出力面 0=YZ平面 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D解析時のpressure.datの出力面の座標
	int d1,d2,d3;							//3D解析時の出力に必要な次元
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標

	double g[3]={0.0, 0.0, 0.0};
	if(dimension==2) g[A_Y]=CON.get_g(); 
	if(dimension==3) g[A_Z]=CON.get_g();

	double mass=CON.get_particle_mass();

//	stringstream ss;
//	ss<<"./StrainRate"<<"/strainrate"<<CON.get_current_step()<<".dat";
//	string filename=ss.str();
	
//	ofstream vec(filename);

//	if(vec.fail()){
//		system("mkdir StrainRate");
//		ofstream vec(filename);//再試行
//		if(vec.fail()){
//			cout<<"StrainRateディレクトリを作成できませんでした"<<endl;
//			exit(1);
//		}
//	}
	
	if(dimension==3)
	{
		//int d1,d2;				//出力に必要な次元
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}

		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
			{
				Micro_avs.make_list(PART[i].r[d1],PART[i].r[d2],PART[i].get_stress_visco_accel(d1),PART[i].get_stress_visco_accel(d2));
/*				double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].get_stress_visco_accel(d1);
				double w=PART[i].get_stress_visco_accel(d2);
				vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				if(x>xmax) xmax=x;
				if(z>ymax) ymax=z;*/
			}
		}
	}
	Micro_avs.Output_vector_MicroAVS("StrainRate",CON.get_current_step());
//	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
//	ymax+=4*le;
//	if(CON.get_legend_pressure()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_pressure()*times<<" "<<0*times<<endl;//最後に凡例出力
//	vec.close();
}

void plot_distortion_rate(mpsconfig &CON, vector<mpselastic> &PART)
{
	Micro_AVS Micro_avs;
	double le=CON.get_distancebp()*0.5;
	double times=CON.get_pressure_times();
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D解析時のpressure.datの出力面 0=YZ平面 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D解析時のpressure.datの出力面の座標
	int d1,d2,d3;							//3D解析時の出力に必要な次元
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標

	double g[3]={0.0, 0.0, 0.0};
	if(dimension==2) g[A_Y]=CON.get_g(); 
	if(dimension==3) g[A_Z]=CON.get_g();

	double mass=CON.get_particle_mass();

/*	stringstream ss;
	ss<<"./distortionRate"<<"/distortionRate"<<CON.get_current_step()<<".dat";
	string filename=ss.str();
	
	ofstream vec(filename);

	if(vec.fail()){
		system("mkdir distortionRate");
		ofstream vec(filename);//再試行
		if(vec.fail()){
			cout<<"distortionRateディレクトリを作成できませんでした"<<endl;
			exit(1);
		}
	}*/
	
	if(dimension==3)
	{
		//int d1,d2;				//出力に必要な次元
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}

		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
			{
				Micro_avs.make_list(PART[i].r[d1],PART[i].r[d2],PART[i].volumetric_strain,PART[i].get_stress_visco_accel(d2));
/*				double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].volumetric_strain;
				double w=PART[i].get_stress_visco_accel(d2);
	//			vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				vec<<x<<" "<<z<<" "<<u<<endl;
				if(x>xmax) xmax=x;
				if(z>ymax) ymax=z;*/
			}
		}
	}
	Micro_avs.Output_vector_MicroAVS("distortionRate",CON.get_current_step());
//	vec.close();
}

void plot_residual_acceleration(mpsconfig &CON, vector<mpselastic> &PART, double **F)
{
	double le=CON.get_distancebp()*0.5;
	double times=CON.get_pressure_times();
	int dimension=CON.get_dimension();
	int face=CON.get_speed_face();			//3D解析時のpressure.datの出力面 0=YZ平面 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D解析時のpressure.datの出力面の座標
	int d1,d2,d3;							//3D解析時の出力に必要な次元
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標

	double g[3]={0.0, 0.0, 0.0};
	if(dimension==2) g[A_Y]=CON.get_g(); 
	if(dimension==3) g[A_Z]=CON.get_g();

	double mass=CON.get_particle_mass();

	stringstream ss;
	ss<<"./Residual"<<"/residual"<<CON.get_current_step()<<".dat";
	string filename=ss.str();
	
	ofstream vec(filename);

	if(vec.fail()){
		system("mkdir Residual");
		ofstream vec(filename);//再試行
		if(vec.fail()){
			cout<<"Residualディレクトリを作成できませんでした"<<endl;
			exit(1);
		}
	}
	
	if(dimension==3)
	{
		//int d1,d2;				//出力に必要な次元
		if(face==0){
			d1=A_Y; d2=A_Z; d3=A_X;
		}
		else if(face==1)
		{
			d1=A_X; d2=A_Z; d3=A_Y;
		}

		for(int i=0;i<PART.size();i++)
		{
			if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
			{
				double x=PART[i].r[d1];
				double z=PART[i].r[d2];
				double u=PART[i].get_stress_accel(d1)+PART[i].get_pressure_accel(d1)+PART[i].get_stress_visco_accel(d1)+g[d1]+F[d1][i]/mass;
				double w=PART[i].get_stress_accel(d2)+PART[i].get_pressure_accel(d2)+PART[i].get_stress_visco_accel(d2)+g[d2]+F[d2][i]/mass;
				vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				if(x>xmax) xmax=x;
				if(z>ymax) ymax=z;
			}
		}
	}
	xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
	ymax+=4*le;
	if(CON.get_legend_pressure()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_pressure()*times<<" "<<0*times<<endl;//最後に凡例出力
	vec.close();
}

void check_FEM_flag(mpsconfig &CON, elastic &ELAST, double ave_P)
{
	if(CON.get_current_time()<CON.get_dt()){//timeはdoubleなので==は使わない
		ofstream fout("aveave_P_history.txt", ios::out);
		fout<<CON.get_current_time()<<"\t"<<ave_P<<endl;
		fout.close();
	}else{
		ofstream fout("aveave_P_history.txt", ios::app);
		fout<<CON.get_current_time()<<"\t"<<ave_P<<endl;
		fout.close();
	}

	//圧力が指定した値になるとFEMスイッチON
	if(CON.get_FEM_flag()==false){
		if(ave_P>CON.get_ave_P_for_FEM_flag()){
			CON.set_FEM_flag(true);//どっちかにしとけ・・・
			ELAST.set_FEM_switch(true);
		}
	}
}

//quickMPS用ポスト処理関数（36行）
void post_processing3(mpsconfig &CON,vector<mpselastic> &PART,int fluid_number,int particle_number,int t,double TIME)
{
	
	//restart用ファイル出力
	if(t==CON.get_step() || t%CON.get_autosave()==0 )
	{
		//restart用に粒子数と粒子データを記録
		ofstream hoge1("number.dat");
		hoge1<<particle_number<<endl;
		hoge1<<TIME<<endl;
		hoge1.close();
		
		FILE *hoge6;
		hoge6=fopen("restart_input.dat","w");//mps_input.datとは区別しないと、restartが失敗したとき困る
		if(hoge6==NULL){
			cout<<"ファイルオープンエラー"<<endl; 
			exit(1);
		}
		for(int i=0;i<particle_number;i++)//流体解析粒子から先に記述
		{
			fprintf( hoge6, "%d\t",i);
			fprintf( hoge6, "%5.10f\t",PART[i].r[A_X]);
			fprintf( hoge6, "%5.10f\t",PART[i].r[A_Y]);
			fprintf( hoge6, "%5.10f\t",PART[i].r[A_Z]);
			fprintf( hoge6, "%5.10f\t",PART[i].u[A_X]);  //速度x成分
			fprintf( hoge6, "%5.10f\t",PART[i].u[A_Y]);  //速度y成分
			fprintf( hoge6, "%5.10f\t",PART[i].u[A_Z]);  //速度z成分
			fprintf( hoge6, "%5.10f\t",PART[i].P); //圧力
			fprintf( hoge6, "%5.10f\t",PART[i].h); //エンタルピー
			fprintf( hoge6, "%5.10f\t",PART[i].val);
			fprintf( hoge6, "%d\n",PART[i].type);
			fprintf( hoge6, "%d\n",PART[i].materialID);
			fprintf( hoge6, "%d\n",PART[i].surface);
			fprintf( hoge6, "%d\n",PART[i].toFEM);
		}
		fclose(hoge6);
		//////////////////////////////////////////////
	}
}

//速度プロット関数（200行）
void plot_speed(mpsconfig &CON ,vector<mpselastic> &PART,int particle_number,int fluid_number)
{
	Micro_AVS Micro_avs;
	double le=CON.get_distancebp()*0.5;
	double times=CON.get_speedtimes();
	int d=CON.get_dimension();
	int NUM;								//AVSに出力する粒子数
	int startID=0;							//最初に出力する粒子のid
	int num=0;								//数えあげ変数
	int face=CON.get_speed_face();			//3D解析時のspeed.datの出力面 0=YZ平面 1=XZ 2=XY
	double face_p=CON.get_speed_face_p();	//3D解析時のspeed.datの出力面の座標
	int d1,d2,d3;								//3D解析時の出力に必要な次元
	double xmax=-100;						//出力粒子の最大横座標
	double ymax=-100;						//出力粒子の最大縦座標
	
	//AVS出力粒子数NUM計算
	if(CON.get_speed_plot_particle()==1) NUM=particle_number;	//全粒子出力
	else if(CON.get_speed_plot_particle()==2) NUM=fluid_number;//流体粒子のみ出力
	else if(CON.get_speed_plot_particle()==3)					//壁粒子のみ出力
	{
		NUM=particle_number; 
		startID=fluid_number;
	}

	{
/*		stringstream ss;
		ss<<"./Speed"<<"/speed"<<CON.get_current_step()<<".dat";
		string filename=ss.str();
	
		ofstream vec(filename);//絶対速度

		if(vec.fail()){
			system("mkdir Speed");
			ofstream vec(filename);//再試行
			if(vec.fail()){
				cout<<"Speedディレクトリを作成できませんでした"<<endl;
				exit(1);
			}
		}*/

		if(d==2)
		{
			for(int i=startID;i<NUM;i++)
			{
				Micro_avs.make_list(PART[i].r[A_X],PART[i].r[A_Y],PART[i].u[A_X],PART[i].u[A_Y]);
	/*			vec<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
				if(PART[i].r[A_X]>xmax) xmax=PART[i].r[A_X];
				if(PART[i].r[A_Y]>ymax) ymax=PART[i].r[A_Y];*/
			}
		}
		else if(d==3)
		{
			//int d1,d2;				//出力に必要な次元
			if(face==0) {d1=A_Y; d2=A_Z; d3=A_X;}
			else if(face==1) {d1=A_X; d2=A_Z; d3=A_Y;}
			if(CON.get_ax_sym_modify()==OFF)
			{
				for(int i=startID;i<NUM;i++)
				{
					if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
					{
						Micro_avs.make_list(PART[i].r[d1],PART[i].r[d2],PART[i].u[d1],PART[i].u[d2]);
					/*	double x=PART[i].r[d1];
						double z=PART[i].r[d2];
						double u=PART[i].u[d1];
						double w=PART[i].u[d2];
						vec<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
						if(x>xmax) xmax=x;
						if(z>ymax) ymax=z;*/
					}
				}
			}
			else if(CON.get_ax_sym_modify()==ON)
			{
				for(int i=startID;i<NUM;i++)
				{
					if(PART[i].type!=WALL){
						if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
						{
							double x=PART[i].r[d1];//出力に関与する軸 便宜上、変数名はxとなっているが、そうとは限らないことに注意
							double z=PART[i].r[d2];//出力に関与する軸
							double u=PART[i].u[d1];//出力に関与する軸
							double w=PART[i].u[d2];//出力に関与する軸

							double y=PART[i].r[d3];//出力に関与しない軸
							double v=PART[i].u[d3];//出力に関与しない軸

							double r=sqrt(x*x+y*y);//原点からの距離
			
							if(r>0)
							{
								double SIN,COS;
								if(x>0)
								{
									SIN=-y/r;
									COS=x/r;
								}
								if(x<0)
								{
									SIN=y/r;
									COS=-x/r;
								}
								double x2=COS*x-SIN*y;//回転後の座標　欲しいのはx2のみ。y2はいらない
								double u2=COS*u-SIN*v;//回転後の速度　欲しいのはu2のみ。v2はいらない
								Micro_avs.make_list(x2,z,u2,w);
						//		vec<<x2<<" "<<z<<" "<<u2*times<<" "<<w*times<<endl;
								//if(SIN*x+COS*y!=0) cout<<SIN*x+COS*y<<endl;
								if(x2>xmax) xmax=x2;
								if(z>ymax) ymax=z;
							}
							else Micro_avs.make_list(x,z,u,w);
						}
					}
				}
			}
		}
		Micro_avs.Output_vector_MicroAVS("Speed",CON.get_current_step());
/*		xmax+=4*le;//凡例を出す位置を保険をかけて少し斜めに移動
		ymax+=4*le;
//		if(CON.get_legend_speed()>0) vec<<xmax<<" "<<ymax<<" "<<CON.get_legend_speed()*times<<" "<<0*times<<endl;//最後に凡例出力
		vec.close();*/
	}
	/////////////////////////////

	/////////重心に対する相対速度出力
	if(CON.get_relative_speed()==ON)
	{
		double U=0;								//平均速度
		double V=0;
		ofstream vec2("relative_speed.dat");

		///平均速度計算
		if(d==2) for(int i=0;i<fluid_number;i++) {U+=PART[i].u[A_X]; V+=PART[i].u[A_Y];}
		else if(d==3) for(int i=0;i<fluid_number;i++) {U+=PART[i].u[d1]; V+=PART[i].u[d2];}//d1,d2はすでにもとまっている
		if(fluid_number>0) {U/=fluid_number; V/=fluid_number;}
		/////////////////////

		if(d==2)
		{
			for(int i=0;i<fluid_number;i++)
			{
				double x=PART[i].r[A_X];
				double y=PART[i].r[A_Y];
				double u=PART[i].u[A_X]-U;
				double v=PART[i].u[A_Y]-V;
				vec2<<x<<" "<<y<<" "<<u*times<<" "<<v*times<<endl;
			}
		}
		else if(d==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].r[face]>face_p-le && PART[i].r[face]<face_p+le)
				{
					double x=PART[i].r[d1];
					double z=PART[i].r[d2];
					double u=PART[i].u[d1]-U;
					double w=PART[i].u[d2]-V;
					vec2<<x<<" "<<z<<" "<<u*times<<" "<<w*times<<endl;
				}
			}
		}
		vec2.close();
	}/////////////////////


	if(CON.get_flat_speed_plot()==ON && d==3) //水平方向の速度出力
	{
		ofstream flat("flat_speed.dat");		
		double face_p2=CON.get_flat_speed_p();
		//for(int i=startID;i<NUM;i++)
		for(int i=0;i<fluid_number;i++)
		{
			if(PART[i].r[A_Z]>face_p2-le && PART[i].r[A_Z]<face_p2+le) flat<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<endl;
		}
		flat.close();
	}

	if(CON.get_speed_AVS()==ON && d==3)
	{
		num=0;
		for(int i=startID;i<NUM;i++) if(PART[i].r[face]<face_p) num++;
		
		ofstream fout2("speed_dist.fld");
		fout2 << "# AVS field file" << endl;
		fout2 << "ndim=1" << endl;
		fout2 << "dim1=" << num <<endl;
		fout2 << "nspace=3" << endl;
		fout2 << "veclen=3" << endl;
		fout2 << "data=float" << endl;
		fout2 << "field=irregular" << endl;
		fout2 << "label=e-x e-y e-z" << endl << endl;
		fout2 << "variable 1 file=./speedvec filetype=ascii skip=1 offset=0 stride=6" << endl;
		fout2 << "variable 2 file=./speedvec filetype=ascii skip=1 offset=1 stride=6" << endl;
		fout2 << "variable 3 file=./speedvec filetype=ascii skip=1 offset=2 stride=6" << endl;
		fout2 << "coord    1 file=./speedvec filetype=ascii skip=1 offset=3 stride=6" << endl;
		fout2 << "coord    2 file=./speedvec filetype=ascii skip=1 offset=4 stride=6" << endl;
		fout2 << "coord    3 file=./speedvec filetype=ascii skip=1 offset=5 stride=6" << endl;
		fout2.close();

		ofstream fout("speedvec");
		fout<<"e-x e-y e-z x y z"<<endl;
		for(int i=startID;i<NUM;i++)
		{
			if(PART[i].r[face]<face_p)
			{
				fout<<PART[i].u[A_X]*times<<" "<<PART[i].u[A_Y]*times<<" "<<PART[i].u[A_Z]*times<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			}
		}
		fout.close();
	}
}

//microAVS用粒子動画出力関数（400行）
void particle_movie_AVS(mpsconfig &CON,vector<mpselastic> &PART, elastic &ELAST, int fluid_number,int particle_number,int t,double T)
{
	//t:タイムステップ　T:総合時間
	double le=CON.get_distancebp();
	double times=CON.get_P_size_AVS();					//出力する粒子のサイズ(leの何倍か)
	bool cut=OFF;	//断面表示フラグ

	if(t==1)
	{
		ofstream fout("particle_movie.mgf");

		fout<<"# Micro AVS Geom:2.00"<<endl;
		fout<<CON.get_step()/CON.get_interval()+1<<endl;//microAVSに出力する総ステップ数。ファイル出力はCON.get_interval()回に1回と最初に行う。
		
		fout.close();
	}

	ofstream avs("particle_movie.mgf",ios :: app);

	if(CON.get_interval()==1) 
	{
		avs<<"step"<<t<<endl;
	}

	else
		avs<<"step"<<t/CON.get_interval()+1<<endl;
		avs<<"sphere"<<endl;
		avs<<"time="<<T<<endl;
		avs<<"color"<<endl;


	double red,green,blue;	//粒子の色を表現する3原色
	
	////////////////////////////粒子の動きだけを表示//////////////////////////				
	if(CON.get_AVS()==0)	
	{	
		double le=CON.get_distancebp();
		double mass=CON.get_particle_mass();


	////////////////////////////表示する粒子数を計算/////////////////////////
		int num=0;//表示する粒子数
		int num2=0;

		if(CON.get_dimension()==2) num=particle_number;//2次元では全粒子を表示
		else if(CON.get_dimension()==3){
			for(int i=0;i<PART.size();i++){
		//		if(PART[i].surface==ON) num++;
				num++;
			}
		}
		avs<<num<<endl;
		/////////////////////////////////

		if(CON.get_dimension()==2)
		{
			for(int i=0;i<PART.size();i++)
			{	
				double width=CON.get_max_pressure()-CON.get_min_pressure();//温度の幅
				double P=0.0;
				for(int D=0;D<3;D++) P+=PART[i].get_pressure_accel(D)*PART[i].get_pressure_accel(D);

				P=sqrt(P);

				double level=(P-CON.get_min_pressure())/width;
				red=level*level; //なぜ二乗している？
				green=-4*level*(level-1);//真ん中で1の放物線
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
		else if(CON.get_dimension()==3)
		{

			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2)
				{
					red=0;
					green=1;//真ん中で1の放物線
					blue=0;
				}
				else if(PART[i].type==ELASTIC)
				{
					red=0;
					green=0;
					blue=1;//真ん中で1の放物線
				}
				else if(PART[i].type==WALL)
				{
					red=1;//1
					green=0;//0
					blue=0;//真ん中で1の放物線 0
				}//*/
				else if(PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2)
				{
					red=1;
					green=1;
					blue=0;//真ん中で1の放物線
				}
				else if(PART[i].type==HYPERELAST)
				{
					red=0;
					green=1;
					blue=1;
				}

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////粒子の圧力をコンター表示//////////////////////////////////
	else if(CON.get_AVS()==1)
	{
		int num=0;//表示する粒子数
		num=particle_number;
		
		avs<<num<<endl;

		///////////////圧力の最大値と最小値をもとめる//////////
		double maxP=0;//
		double minP=0;
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].P>maxP) maxP=PART[i].P;
			if(PART[i].P<minP) minP=PART[i].P;
		}
		/////maxP,minPが求まった。
		//////////////////////////////////////////////////////

		double width=maxP-minP;
		if(width<0.0001) width=0.0001;
			
		if(CON.get_dimension()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				double level=(PART[i].P-minP)/width;
				red=level*level;
				green=-4*level*(level-1);//真ん中で1の放物線
				blue=1-red;
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
		else if(CON.get_dimension()==3)
		{
			for(int i=0;i<particle_number;i++)//2013-02-26 fluid_number→particle_number
			{
				double level=(PART[i].P-minP)/width;
				red=level*level;
				green=-4*level*(level-1);//真ん中(level==1/2)で1の放物線
				blue=1-red;
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
		
	}
	else if(CON.get_AVS()==2)	//粒子の温度をコンター表示
	{
		//////////////////////////////////表示する粒子数を計算
		int num=0;//表示する粒子数
		if(CON.get_dimension()==2) num=particle_number;//2次元では全粒子を表示
		else if(CON.get_dimension()==3) num=fluid_number;//3次元では流体粒子のみを表示
		
		avs<<num<<endl;
		/////////////////////////////////

		double le=CON.get_distancebp();
		double mass=CON.get_particle_mass();
		double T;//粒子iの温度
		double width=CON.get_maxT()-CON.get_minT();//温度の幅

		if(CON.get_dimension()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				double hs0=mass*CON.get_Cp()*CON.get_MP();//融解開始点のエンタルピー
				double hs1=hs0+CON.get_latent_H()*mass;//融解終了点のエンタルピー
				if(PART[i].h<hs0) T=PART[i].h/mass/CON.get_Cp();
				else if(hs0<=PART[i].h && PART[i].h<=hs1) T=CON.get_MP();
				else if(hs1<PART[i].h) T=CON.get_MP()+(PART[i].h-hs1)/mass/CON.get_Cp();
				double level=(T-CON.get_minT())/width;
				red=level*level;
				green=-4*level*(level-1);//真ん中で1の放物線
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
		else if(CON.get_dimension()==3)
		{
			for(int i=0;i<particle_number;i++)
			{
				double hs0=mass*CON.get_Cp()*CON.get_MP();//融解開始点のエンタルピー
				double hs1=hs0+CON.get_latent_H()*mass;//融解終了点のエンタルピー
				if(PART[i].h<hs0) T=PART[i].h/mass/CON.get_Cp();
				else if(hs0<=PART[i].h && PART[i].h<=hs1) T=CON.get_MP();
				else if(hs1<PART[i].h) T=CON.get_MP()+(PART[i].h-hs1)/mass/CON.get_Cp();
				double level=(T-CON.get_minT())/width;
				red=level*level;
				green=-4*level*(level-1);//真ん中で1の放物線
				blue=1-red;

				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
				
				avs<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力
	
				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
		}
	}
	else if(CON.get_AVS()==3)//壁粒子は非表示
	{
		//////////////////////////////////表示する粒子数を計算
		int num=fluid_number;//表示する粒子数
		avs<<num<<endl;
		/////////////////////////////////

		////粒子出力
		if(CON.get_dimension()==2) cout<<"error in AVS() 2Dは非対応"<<endl;      
		else if(CON.get_dimension()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if((PART[i].type==FLUID || PART[i].type==ELASTIC ) && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//壁粒子
			    
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
				avs<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力

				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力        
			}
		}
	}
	else if(CON.get_AVS()==4)//流体表面粒子のみ
	{
		//////////////////////////////////表示する粒子数を計算
		int num=0;//表示する粒子数
		if(CON.get_dimension()==2) num=particle_number;//2次元では全粒子を表示
		else if(CON.get_dimension()==3) for(int i=0;i<particle_number;i++) if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==ON) num++;//3次元の場合、内部流体は表示しない	
		avs<<num<<endl;
		/////////////////////////////////

		////粒子出力
		if(CON.get_dimension()==2)
		{
			for(int i=0;i<particle_number;i++)
			{
				if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==OFF) {red=0;green=0;blue=1;}
				else if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        	else {red=0.5;green=0.5;blue=0;}//壁粒子
			   
				avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
				avs<<CON.get_distancebp()/2<<" ";//粒子の大きさ出力

				avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
			}
			
		}
		else if(CON.get_dimension()==3)
		{
			for(int i=0;i<fluid_number;i++)
			{
				if(PART[i].surface==ON) 
				{
					if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==OFF) {red=0;green=0;blue=1;}
					else if((PART[i].type==FLUID || PART[i].type==ELASTIC) && PART[i].surface==ON) {red=0;green=0.5;blue=0.5;}
	        		else {red=0.5;green=0.5;blue=0;}//壁粒子
			    
					avs<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
			
					avs<<CON.get_distancebp()/2<<" ";//粒子の大きさ出力

					avs<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}        
			}
		}
	}
	if(CON.get_AVS()==5)	//壁だけを表示.ただの確認用
	{

	}
	else if(CON.get_AVS()==6){}	//特定の粒子の動きだけを表示。粒子の選択はプログラムを直接変更するしかない
	


	avs.close();
		////////////////////

	if(CON.get_dimension()==3 && CON.get_AVS()==0)//断面表示フラグ
	{
		if(CON.get_cut_x()==ON)
		{	
			if(t==1)
			{
				ofstream fout("half_x_particle_movie.mgf");

				fout<<"# Micro AVS Geom:2.00"<<endl;
				fout<<CON.get_step()/CON.get_interval()+1<<endl;//microAVSに出力する総ステップ数。ファイル出力はCON.get_interval()回に1回と最初に行う。
				fout.close();
			}

			ofstream avs2("half_x_particle_movie.mgf",ios :: app);
	
			if(CON.get_interval()==1) avs2<<"step"<<t<<endl;

			else
				avs2<<"step"<<t/CON.get_interval()+1<<endl;
				avs2<<"sphere"<<endl;
				avs2<<"time="<<T<<endl;
				avs2<<"color"<<endl;
		

			double red,green,blue;	//粒子の色を表現する3原色
	
			////////////////////////////粒子の動きだけを表示//////////////////////////				

			double le=CON.get_distancebp();
			double mass=CON.get_particle_mass();
			double width=CON.get_max_pressure()-CON.get_min_pressure();//温度の幅

		////////////////////////////表示する粒子数を計算/////////////////////////
			int num=0;//表示する粒子数


			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].r[A_Y]>=-0.00025 && PART[i].surface==ON)
				{
					num++;
				}
			}

			avs2<<num<<endl;
			/////////////////////////////////

			double P=0.0;
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].r[A_Y]>=-0.00025 && PART[i].surface==ON){
				//圧力の計算

				for(int D=0;D<3;D++) P+=PART[i].get_pressure_accel(D)*PART[i].get_pressure_accel(D);

				P=sqrt(P);

				double level=(P-CON.get_min_pressure())/width;
				if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2){
				red=0;
				green=1;//真ん中で1の放物線
				blue=0;
				}
				else if(PART[i].type==ELASTIC){
				red=0;
				green=0;
				blue=1;//真ん中で1の放物線
				}
				else if(PART[i].type==WALL ){
				red=1;
				green=0;
				blue=0;//真ん中で1の放物線
				}
				else if(PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
				red=1;
				green=1;
				blue=0;//真ん中で1の放物線
				}
				else if(PART[i].type==HYPERELAST)
				{
					red=0;
					green=1;
					blue=1;
				}

				avs2<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
		
				avs2<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力

				avs2<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}
			}
			avs2.close();
		}

		if(CON.get_cut_y()==ON)  	//断面表示フラグ
		{	
			if(t==1)
			{
				ofstream fout("half_y_particle_movie.mgf");

				fout<<"# Micro AVS Geom:2.00"<<endl;
				fout<<CON.get_step()/CON.get_interval()+1<<endl;//microAVSに出力する総ステップ数。ファイル出力はCON.get_interval()回に1回と最初に行う。
				fout.close();
			}

			ofstream avs3("half_y_particle_movie.mgf",ios :: app);
	
			if(CON.get_interval()==1) avs3<<"step"<<t<<endl;

			else
				avs3<<"step"<<t/CON.get_interval()+1<<endl;
				avs3<<"sphere"<<endl;
				avs3<<"time="<<T<<endl;
				avs3<<"color"<<endl;
		

			double red,green,blue;	//粒子の色を表現する3原色
	
			////////////////////////////粒子の動きだけを表示//////////////////////////				

			double le=CON.get_distancebp();
			double mass=CON.get_particle_mass();
			double width=CON.get_max_pressure()-CON.get_min_pressure();//温度の幅

		////////////////////////////表示する粒子数を計算/////////////////////////
			int num=0;//表示する粒子数


			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].r[A_X]>=-0.00025 && PART[i].surface==ON)
				{
					num++;
				}
			}

			avs3<<num<<endl;
			/////////////////////////////////

			double P=0.0;
			for(int i=0;i<PART.size();i++)
			{
				if(PART[i].r[A_X]>=-0.00025 && PART[i].surface==ON){
				//圧力の計算

				for(int D=0;D<3;D++) P+=PART[i].get_pressure_accel(D)*PART[i].get_pressure_accel(D);

				P=sqrt(P);

				double level=(P-CON.get_min_pressure())/width;
				if(PART[i].type==MAGELAST || PART[i].type==MAGELAST2){
				red=0;
				green=1;//真ん中で1の放物線
				blue=0;
				}
				else if(PART[i].type==ELASTIC){
				red=0;
				green=0;
				blue=1;//真ん中で1の放物線
				}
				else if(PART[i].type==WALL ){
				red=1;
				green=0;
				blue=0;//真ん中で1の放物線
				}
				else if(PART[i].type==TERMINAL1 || PART[i].type==TERMINAL2){
				red=1;
				green=1;
				blue=0;//真ん中で1の放物線
				}
				else if(PART[i].type==HYPERELAST)
				{
					red=0;
					green=1;
					blue=1;
				}

				avs3<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<" ";//座標出力
		
				avs3<<CON.get_distancebp()*times<<" ";//粒子の大きさ出力

				avs3<<red<<" "<<green<<" "<<blue<<endl;//色出力
				}
			}
			avs3.close();
		}
	}
}





//圧力など粒子の持つ情報をコンター図で表示する関数
void particle_movie_AVS2(mpsconfig &CON,int t,vector<mpselastic> &PART,double TIME,int fluid_number, int particle_number, double **F)
{
	//参考にしている書式はmicroAVSのヘルプであなたのデータは？→「非構造格子型データ（アスキー）の書式」
	
	double le=CON.get_distancebp();
	int STEP=CON.get_step()/CON.get_EM_interval()+1;		//出力する総ステップ数
	int *input=new int[particle_number];			//input[]=ONなら出力

	if(t==1) 
	{
		ofstream fp("particle_movie2.inp");			
		fp<<STEP<<endl;						//microAVSに出力する総ｽﾃｯﾌﾟ数。ファイル出力はCON->get_interval()回に1回と最初に行う。
		fp<<"data_geom"<<endl;
		fp.close();
	}


	//mainファイル書き込み
	ofstream fp("particle_movie2.inp",ios :: app);
	fp<<"step"<<t/CON.get_EM_interval()+1<<" TIME="<<TIME<<endl;

	//出力粒子数算出
	int num=0;//表示する粒子数
	if(CON.get_dimension()==2)
	{
		num=particle_number;//2次元では全粒子を表示
		for(int i=0;i<particle_number;i++) input[i]=ON;
	}
	else if(CON.get_dimension()==3)
	{
		for(int i=0;i<particle_number;i++)
		{
			if(PART[i].type!=WALL || PART[i].type!=ELASTIC)
			{
				input[i]=ON;
				num++;//3次元の場合、内部流体は表示しない
			}
			else input[i]=OFF;
		}
	}

	fp<<num<<" "<<num<<endl;	//節点数と要素数出力 この場合は両方とも節点数をいれておく

	//節点番号とその座標の出力 
	int count=0;
	for(int i=0;i<particle_number;i++)
	{
//		if(input[i]==ON)
		{
			 fp<<count<<" "<<PART[i].r[A_X]<<" "<<PART[i].r[A_Y]<<" "<<PART[i].r[A_Z]<<endl;
			 count++;
		}
	}

	//要素番号と要素形状の種類、そして要素を構成する節点番号出力
	for(int i=0;i<num;i++)
	{
		fp<<i<<"  0 pt "<<i<<endl;
	}

	//fp<<"2 3"<<endl;//節点の情報量が2で、要素の情報量が3ということ。
	fp<<"5 0"<<endl;//節点の情報量が8で、要素の情報量が0ということ。
	fp<<"5 1 1 1 1 1"<<endl;	//この行の詳細はヘルプを参照
	fp<<"surface,"<<endl;
	fp<<"type,\n";
	fp<<"MF_x,"<<endl;
	fp<<"MF_y,"<<endl;
	fp<<"MF_z,"<<endl;


	//各節点の情報値入力
	count=0;
	for(int i=0;i<particle_number;i++)
	{
//		if(input[i]==ON)
		{
			fp<<count<<" "<<PART[i].surface<<" "<<PART[i].type<<" "<<F[A_X][i]<<" "<<F[A_Y][i]<<" "<<F[A_Z][i]<<"\n";
			count++;
		}
	}

	
	fp.close();

	delete [] input;

	
}