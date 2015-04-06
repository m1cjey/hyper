#include "Micro_AVS.h"
#include "CONFIG.h"
std::vector<std::string> Micro_AVS::Filelist(0);
bool Micro_AVS::fast_flag(true);
Micro_AVS::Micro_AVS()
{}

void Micro_AVS::make_list(double xr,double yr,double xdata,double ydata){
	Xr.push_back(xr);
	Yr.push_back(yr);
	Xdata.push_back(xdata);
	Ydata.push_back(ydata);
	Dimension=2;
}

void Micro_AVS::make_list(double xr,double yr,double zr,double xdata,double ydata,double zdata){
	Xr.push_back(xr);
	Yr.push_back(yr);
	Zr.push_back(zr);
	Xdata.push_back(xdata);
	Ydata.push_back(ydata);
	Zdata.push_back(zdata);
	Dimension=3;
}

void Micro_AVS::Output_vector_MicroAVS(std::string fname,int step){
	
	Filename=fname;
	Step=step;

	std::stringstream ss;
	ss<<"./"<<Filename<<"./"<<Filename<<".dat";
	std::ofstream check_file(ss.str());

	//フォルダーの確認
	if(check_file.fail()){
		std::stringstream makefolder;
		makefolder<<"mkdir "<<Filename;
		std::string a=makefolder.str();
		const char* folde=a.c_str();
		system(folde);
		std::ofstream check_file(ss.str());
		if(check_file.fail()){
			std::cout<<Filename<<"フォルダを開けませんでした"<<std::endl;
			exit(1);
		}
		Filelist.push_back(Filename); //フォルダが作られたので名前を記憶
	}
	else {
		_Foldername_storing();
	}
	check_file.close();

	if(Dimension==2){
		_Wright_vector_2d();
	}
	else if(Dimension==3){
		_Wright_vector_3d();
	}
	else std::cout<<"リストを作成してください"<<std::endl;

}

void Micro_AVS::Output_mgf_MicroAVS(std::string fname,int max_step){
	Filename=fname;
	Step=max_step;

	std::stringstream ss;
	ss<<"./"<<Filename<<"./"<<Filename<<".dat"; //フォルダーを作るための不要なデータファイル
	std::ofstream check_file(ss.str());

	//フォルダーの確認
	if(check_file.fail()){
		std::stringstream makefolder;
		makefolder<<"mkdir "<<Filename;
		std::string a=makefolder.str();
		const char* folde=a.c_str();
		system(folde);
		std::ofstream check_file(ss.str());
		if(check_file.fail()){
			std::cout<<Filename<<"フォルダを開けませんでした"<<std::endl;
			exit(1);
		}
		Filelist.push_back(Filename); //フォルダが作られたので名前を記憶
	}
	else {
		_Foldername_storing();
	}
	check_file.close();

	if(Dimension==2){
		//_Wright_mgf_2d();
	}
	else if(Dimension==3){
		_Wright_mgf_3d();
	}
	else std::cout<<"リストを作成してください"<<std::endl;
}

void Micro_AVS::Get_folder_name(){
	std::ofstream nlist("Folder_name.dat");
	if(Filename.size()!=0){
		for(int fol=0;fol<Filename.size();fol++){
			nlist<<Filename[fol]<<std::endl;
		}
	}
	else nlist<<"folder is nothing";
	nlist.close();
}

void Micro_AVS::_Wright_vector_2d(){

	std::stringstream ss2;
	ss2<<"./"<<Filename<<"./"<<Filename<<"_"<<Step<<".fld";
	std::ofstream check_file(ss2.str());

	//
	std::ofstream output_file(ss2.str());
	output_file << "# AVS field file" << std::endl;
	output_file << "ndim=1" << std::endl;
	output_file << "dim1=" << Xr.size() <<std::endl;
	output_file << "nspace=2" << std::endl;
	output_file << "veclen=2" << std::endl;
	output_file << "data=float" << std::endl;
	output_file << "field=irregular" << std::endl;
	output_file << "label=e-x e-y" << std::endl << std::endl;
	output_file << "variable 1 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=0 stride=4" << std::endl;
	output_file << "variable 2 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=1 stride=4" << std::endl;
	output_file << "coord    1 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=2 stride=4" << std::endl;
	output_file << "coord    2 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=3 stride=4" << std::endl;
	output_file.close();
	

	//データ部
	std::stringstream ss2_d;
	ss2_d<<"./"<<Filename<<"./"<<Filename<<"_"<<Step;
	std::ofstream dat_file(ss2_d.str());

	dat_file<<"e-x e-y x y "<<std::endl;
	for(int idata=0;idata<Xr.size();idata++){
		dat_file<<Xdata[idata]<<" "<<Ydata[idata]<<" "<<Xr[idata]<<" "<<Yr[idata]<<std::endl;
	}
	dat_file.close();

	//データリセット
	Dimension=0;
	Step=0;
	Xr.clear();
	Yr.clear();
	Xdata.clear();
	Ydata.clear();

}

void Micro_AVS::_Wright_vector_3d(){

	std::stringstream ss3;
	ss3<<"./"<<Filename<<"./"<<Filename<<"_"<<Step<<".fld";
	std::ofstream check_file(ss3.str());

	//
	std::ofstream output_file(ss3.str());
	output_file << "# AVS field file" << std::endl;
	output_file << "ndim=1" << std::endl;
	output_file << "dim1=" << Xr.size() <<std::endl;
	output_file << "nspace=3" << std::endl;
	output_file << "veclen=3" << std::endl;
	output_file << "data=float" << std::endl;
	output_file << "field=irregular" << std::endl;
	output_file << "label=e-x e-y e-z" << std::endl << std::endl;
	output_file << "variable 1 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=0 stride=6" << std::endl;
	output_file << "variable 2 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=1 stride=6" << std::endl;
	output_file << "variable 3 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=2 stride=6" << std::endl;
	output_file << "coord    1 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=3 stride=6" << std::endl;
	output_file << "coord    2 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=4 stride=6" << std::endl;
	output_file << "coord    3 file=./"<<Filename<<"_"<<Step<<" filetype=ascii skip=1 offset=5 stride=6" << std::endl;
	output_file.close();

	//データ部
	std::stringstream ss3_d;
	ss3_d<<"./"<<Filename<<"./"<<Filename<<"_"<<Step;
	std::ofstream dat_file(ss3_d.str());

	dat_file<<"e-x e-y e-z x y z"<<std::endl;
	for(int idata=0;idata<Xr.size();idata++){	
		dat_file<<Xdata[idata]<<" "<<Ydata[idata]<<" "<<Zdata[idata]<<" "<<Xr[idata]<<" "<<Yr[idata]<<" "<<Zr[idata]<<std::endl;
	}
	dat_file.close();

	//データリセット
	Dimension=0;
	Step=0;
	Xr.clear();
	Yr.clear();
	Zr.clear();
	Xdata.clear();
	Ydata.clear();
	Zdata.clear();

}

void Micro_AVS::_Wright_mgf_3d(){
	mpsconfig CON;
	
	std::stringstream ssm3;
	ssm3<<"./"<<Filename<<"./"<<Filename<<".mgf";
	std::ofstream check_file(ssm3.str());
	
	//
	if(!(check_file.fail())){
	std::ofstream output_file(ssm3.str());
	output_file<<"# Micro AVS Geom:2.00"<<std::endl;
	output_file.close();
//	fast_flag=false;
	}
	std::ofstream output_file2(ssm3.str(),std::ios::app);
	output_file2<<"sphere"<<endl;
	output_file2<<"particle"<<endl;
	output_file2<<"color"<<endl;
	output_file2<<Xr.size()<<endl;
	for(int i=0;i<Xr.size();i++){
		output_file2<<Xr[i]<<" "<<Yr[i]<<" "<<Zr[i]<<" ";//座標出力
		output_file2<<CON.get_distancebp()/2<<" ";//粒子の大きさ出力
		output_file2<<Xdata[i]<<" "<<Ydata[i]<<" "<<1.0<<endl;//座標出力
	}

	output_file2.close();
	//データリセット
	Dimension=0;
	Step=0;
	Xr.clear();
	Yr.clear();
	Zr.clear();
	Xdata.clear();
	Ydata.clear();
	Zdata.clear();
}

void Micro_AVS::_Foldername_storing(){
	bool storing=true;
	if(Filelist.size()!=0){
		for(int filenum=0;filenum<Filelist.size();filenum++){
				if(Filelist[filenum]==Filename); //フォルダ名はすでに記憶されている
				else storing=false;				//まだ記憶されていないフォルダ名
			}
		if(storing==false) Filelist.push_back(Filename); //フォルダ名を記憶
	}
	else Filelist.push_back(Filename); //最初のフォルダ名を記憶
}