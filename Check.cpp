#include "stdafx.h"
#include "Check.h"
#include "CONFIG.h"
#include "PART.h"


Check::Check():checkf("./Check/"){}	//コンストラクタ
Check::~Check(){}					//デストラクタ

string Check::Set_y_m_d()
{
	time_t timer;
	struct tm *utc;	/*協定世界時(UTC)*/	
	stringstream ss;		//
	int  year, mon, day;
	timer=time(NULL);
	utc = gmtime(&timer);
	year=utc->tm_year+1900;
	mon=utc->tm_mon;
	day=utc->tm_mday;
	
	ss<<year<<"_"<<mon<<"_"<<day;	//"/"は使えない。
	return ss.str();

}

void Check::Out_put_config()
{
	mpsconfig CON; //解析条件格納クラス
	stringstream day;
	double Mh=CON.get_magnet_H();//磁石の高さ
	double Mc=CON.get_magnet_Z();//自作中心座標
	double MMd=(abs((Mh/2)+Mc)+0.006)*1000;	//0.006?なぜ決め打ちなのか15/5/24
	day<<checkf<<Set_y_m_d()<<".dat";
	system("mkdir Check");
	ofstream checkconfig(day.str());
	if(checkconfig.fail()){
		system("mkdir Check");
		ofstream checkconfig(day.str());//同じ階層でないと書き込めない？ここでファイルを作ると空のファイルができる。
		if(!checkconfig){
			cout<<"can't open file(checkconfig) "<<endl;	
			getchar();
		}
	}
	checkconfig<<"///////解析条件////////"<<endl;
	checkconfig<<"時間刻み="<<CON.get_dt()<<endl;
	checkconfig<<"step="<<CON.get_step()<<endl;
	checkconfig<<"interval="<<CON.get_interval()<<endl;
	if(CON.get_FEM_flag()==0) checkconfig<<"FEM=OFF"<<endl;
	else checkconfig<<"FEM=ON"<<endl;
	if(CON.get_nonlinear_elastic_flag()==0) checkconfig<<"elastic nonlinearity=OFF"<<endl;
	else checkconfig<<"elastic nonlinearity=ON"<<endl;
	cout<<endl;
	checkconfig<<"///////モデル条件////////"<<endl;
	checkconfig<<"model_number="<<CON.get_model_number()<<endl;
	if(CON.get_model_number()!=21)
	{
		checkconfig<<"viscosity="<<CON.get_nensei()<<endl;
		checkconfig<<"MREと磁石の距離[mm]="<<MMd<<endl;
		checkconfig<<"MREの比透磁率="<<CON.get_RP()<<endl;
		checkconfig<<"MREのヤング率="<<CON.get_E_m()<<endl;
		checkconfig<<"MREのポアソン比="<<CON.get_v_m()<<endl;
		checkconfig<<"ICCG法FEMの解析度="<<CON.get_FEMCGep()<<"\n";
		checkconfig<<"球の半径＝"<<CON.get_R1()<<"\n";
/*		if(CON.get_avoid_step()!=0)
			checkconfig<<"避けているstep＝"<<CON.get_avoid_step()<<","<<CON.get_avoid_step2()<<","<<CON.get_avoid_step3()<<","<<CON.get_avoid_step4()<<","<<CON.get_avoid_step5()<<","<<CON.get_avoid_step6()<<","<<CON.get_avoid_step7()<<"\n";*/
	}
	if(CON.get_model_number()==21)
	{
		checkconfig<<"c10="<<CON.get_c10()<<endl;
		checkconfig<<"c01="<<CON.get_c01()<<endl;
		checkconfig<<"hyper_density="<<CON.get_hyper_density()<<endl;
	}


	checkconfig.close();
}

void Check::Courant_condition(vector<mpselastic> &PART)
{
	mpsconfig CON;
	double dt=CON.get_dt();
	double L=CON.get_distancebp();
	double kv=CON.get_nensei()/CON.get_density();//動粘性係数
	double CFL=0 ,CFLmax=0;
	double Cmax=0.2;	//クーラン数制限0.2
	double vv=0;	//速度ベクトルの大きさ
	double vvmax=0.00000000000001;
	double new_dt=0;
	double di=0;
	stringstream day2;
	day2<<checkf<<"Courant and diffusion number check.dat";
	ofstream cfl(day2.str(), ios::app);
	if(cfl.fail()){
		system("mkdir Check");
		ofstream cfl(day2.str(), ios::app);
		if(cfl.fail()){
			cout<<"can not opning check fail"<<endl;
			getchar();
		}
	}
	for(int i=0;i<PART.size();i++){
	vv=pow(pow(PART[i].u[A_X],2)+pow(PART[i].u[A_Y],2)+pow(PART[i].u[A_Z],2),0.5);
	CFL=(dt*vv)/L;
		if(CFL>CFLmax){ 
			CFLmax=CFL;
			vvmax=vv;
		}
	}

	if(CFLmax>Cmax){ 
		new_dt=(L*Cmax)/vvmax;
		cfl<<"CFL error: CFL="<<CFLmax<<" ,necessary dt="<<new_dt<<endl;
	}
	di=(kv*dt)/(L*L);
	new_dt=0.5*(L*L)/kv;
	if(di>=0.51) {
		cfl<<"diffusion error: di="<<di<<" ,necessary dt="<<new_dt<<endl;
		getchar();
	}
}