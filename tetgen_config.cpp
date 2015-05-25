#include "stdafx.h"

tetgen_config::tetgen_config()
{
	//Γ“d–¶‰»‘•’u‚Μ΅–@
	radius_column=0.0005;		//‰~’“d‹Ι”Όa
	length_column=0.02;			//‰~’“d‹Ι’·‚³
	height_plate=0.0075;		//•½”Β“d‹Ι‚‚³
	length_plate=0.02;			//•½”Β“d‹Ικ•Σ’·‚³
	thickness_plate=0.0005;		//•½”Β“d‹Ιϊ‚³
	length_base=0.02;			//“y‘δκ•Σ’·‚³
	thickness_base=0.005;		//“y‘δϊ‚³
	
	//ƒƒbƒVƒ…‚Μ‘e‚³(‹«E‚Μ1—v‘f‚Μ•Σ‚Μ’·‚³‚π‚Η‚κ‚®‚η‚Ά‚Ι‚·‚ι‚©)
	fine_air=0.01;				//‹σ‹C—Μζ‹«E‚Μ‘e‚³		0.01
	fine_plate_t=0.00025;		//•½”Β“d‹Ιϊ‚έ•ϋό‚Μ‘e‚³	0.00025
	fine_plate_L=0.00050;		//•½”Β“d‹Ι•½–Κ•ϋό‚Μ‘e‚³	0.0005
	fine_column_L=0.00015;		//‰~’“d‹Ι‚Μ’·‚³•ϋό‚Μ‘e‚³	0.00015
	fine_base=0.0010;			//“y‘δ•\–Κ‚Μ‘e‚³			0.0010

	//…“H•\–Κ‚ΜƒƒbƒVƒ…‘w‚Μέ’θ
	num_layer_out=1;	//—¬‘ΜO‘¤ƒƒbƒVƒ…‘w”
	num_layer_in=0;		//—¬‘Μ“ΰ‘¤ƒƒbƒVƒ…‘w”
	thick_layer=0.3;	//‹«EƒƒbƒVƒ…1‘w‚Μϊ‚³(le‚Μ‰½”{‚©)

	//’·‚Ά—¬‘Μ—v‘f‚πν‚·‚ιθ‡’l 
	del_length=3.0;		//le‚Μ‰½”{Θγ‚Μ•Σ‚π‚Β—v‘f‚πΑ‚·‚©2.0
}

tetgen_config::tetgen_config(mpsconfig &CON)
{
	//¥Ξ΅–@ cf. MPSTOFEM_MRE()<-MPS_TO_FEM3D.cpp
	magnet_height=CON.get_magnet_H();
	magnet_radius=CON.get_magnet_r();

	//ƒƒbƒVƒ…‚Μ‘e‚³(‹«E‚Μ1—v‘f‚Μ•Σ‚Μ’·‚³‚π‚Η‚κ‚®‚η‚Ά‚Ι‚·‚ι‚©)
	fine_air=0.01;				//‹σ‹C—Μζ‹«E‚Μ‘e‚³		0.01
	fine_plate_t=0.00025;		//•½”Β“d‹Ιϊ‚έ•ϋό‚Μ‘e‚³	0.00025
	fine_plate_L=0.00050;		//•½”Β“d‹Ι•½–Κ•ϋό‚Μ‘e‚³	0.0005
	fine_column_L=0.00015;		//‰~’“d‹Ι‚Μ’·‚³•ϋό‚Μ‘e‚³	0.00015
	fine_base=0.0010;			//“y‘δ•\–Κ‚Μ‘e‚³			0.0010

	//…“H•\–Κ‚ΜƒƒbƒVƒ…‘w‚Μέ’θ
	num_layer_out=1;	//—¬‘ΜO‘¤ƒƒbƒVƒ…‘w”
	num_layer_in=0;		//—¬‘Μ“ΰ‘¤ƒƒbƒVƒ…‘w”
	thick_layer=0.3;	//‹«EƒƒbƒVƒ…1‘w‚Μϊ‚³(le‚Μ‰½”{‚©)

	//’·‚Ά—¬‘Μ—v‘f‚πν‚·‚ιθ‡’l 
	del_length=2.0;		//le‚Μ‰½”{Θγ‚Μ•Σ‚π‚Β—v‘f‚πΑ‚·‚©
}
