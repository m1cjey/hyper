class Post
{
	bool plot_B;	//磁束密度出力
	bool plot_NF;	//節点力出力
public:
	Post();
	bool get_plot_B(){return plot_B;}
	bool get_plot_NF(){return plot_NF;}
};

//表示用フラグの設定
Post::Post()
{
	plot_B=true;
	plot_NF=true;
}