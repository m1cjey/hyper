#include"stdafx.h"

/********************ベクトルノルムの計算********************/

//1ノルム
double vector_norm1(double *d, int start, int end)
{
	double norm=0.0;
	
	for(int i=start;i<end;i++) norm+=fabs(d[i]);
	//for(int i=m;i<n;i++) norm+=fabs(d[i]);
	//cout<<norm<<endl;
//	cout<<"finished, n:"<<n<<endl;

	return norm;
}

//2ノルム
double vector_norm2(double *d, int start, int end)
{
	double sum=0.0;

	for(int i=start;i<end;i++) sum+=d[i]*d[i];
//	cout<<sqrt(sum)<<endl;

	return sqrt(sum);
}

/********************連立一次方程式の計算********************/

//ガウスの消去法
void pivot_gauss(double **A, double *b, int size, bool pivot_check)
{
	int i, j, k, ip;
	double alpha, temp;
	double amax, eps=pow(2.0, -12.0);

	for(k=0;k<size-1;k++)
	{
		amax=fabs(A[k][k]);
		ip=k;

		//ピボット選択と行交換
		if(pivot_check==ON)
		{
			for(i=k+1;i<size;i++)//ピボット探索
			{
				if(fabs(A[i][k])>amax)
				{
					amax=fabs(A[i][k]);
					ip=i;
				}
			}

			if(amax<eps){cout<<"行列は正則ではない。計算中断"<<endl; return;}	//正則性の判定

			//行交換
			if(ip!=k)//ピボットの行が注目している行と同じ行なら交換の必要なし
			{
				for(j=k;j<size;j++)//各列ごとに交換作業（jは列番号）
				{
					temp=A[k][j];
					A[k][j]=A[ip][j];
					A[ip][j]=temp;
				}
				temp=b[k];//temp=x[k];
				b[k]=b[ip];//x[k]=x[ip];
				b[ip]=temp;//x[ip]=temp;
			}
		}

		//前進消去（第k列のk+1行以下は全て0になる...はず）
		for(i=k+1;i<size;i++)
		{
			alpha=-1*A[i][k]/A[k][k];
			for(j=k+1;j<size;j++)
			{
				A[i][j]+=alpha*A[k][j];
			}
			b[i]+=alpha*b[k];
		}
	}

	//後退代入
	b[size-1]/=A[size-1][size-1];
	for(k=size-2;k>=0;k--)//k=size-1は既に求めている
	{
		temp=b[k];
		for(j=k+1;j<size;j++)
		{
			temp=temp-A[k][j]*b[j];
		}
		b[k]=temp/A[k][k];
	}
}

//LU分解を用いる size:=ROW(正方行列を仮定)
void lu_decomposition(double **A, double *b, size_t size, bool pivot_check)
{
	int i, j, k, ip;
	double alpha, temp;
	double amax;
	int *p=new int[size];

	for(k=0;k<size-1;k++)
	{
		//ピボット選択
		amax=fabs(A[k][k]);//対角成分
		ip=k;

		for(i=k+1;i<size;i++)
		{
			if(fabs(A[i][k])>amax)
			{
				amax=fabs(A[i][k]);
				ip=i;
			}
		}

		//正則性の判定
		if(amax<EPS){
			cout<<"行列は正則ではない。計算終了"<<endl; 
			return;
		}

		p[k]=ip; //ipを配列p[i]に保存

		//行交換
		if(ip!=k)
		{
			for(j=k;j<size;j++)
			{
				temp=A[k][j];
				A[k][j]=A[ip][j];
				A[ip][j]=temp;
			}
		}

		//前進消去
		for(i=k+1;i<size;i++)
		{
			alpha=-A[i][k]/A[k][k];
			A[i][k]=alpha;
			for(j=k+1;j<size;j++)
			{
				A[i][j]+=alpha*A[k][j];
			}
		}
	}

	for(k=0;k<size-1;k++)
	{
		temp=b[k];
		b[k]=b[p[k]];//???
		b[p[k]]=temp;
		//前進代入
		for(i=k+1;i<size;i++)
		{
			b[i]+=A[i][k]*b[k];
		}
	}

	//後退代入・・・A[size-1][size-1]!=0は保証されている？？？
	b[size-1]/=A[size-1][size-1];
	for(k=size-2;k>=0;k--)
	{
		temp=b[k];
		for(j=k+1;j<size;j++)
		{
			temp-=A[k][j]*b[j];
		}
		b[k]=temp/A[k][k];
	}

	delete [] p;
}

//LU分解を用いる(vectorを使う)
void lu_decomposition(vector<vector<double>> &A, vector<double> &b)
{
	int i, j, k, ip;
	size_t size=b.size();
	double alpha, temp;
	double amax;
	vector<int> p(size); //宣言と同時にサイズを確保しておく

	for(k=0;k<size-1;k++)
	{
		//ピボット選択
		amax=fabs(A[k][k]);//対角成分
		ip=k;

		for(i=k+1;i<size;i++)
		{
			if(fabs(A[i][k])>amax)
			{
				amax=fabs(A[i][k]);
				ip=i;
			}
		}

		//正則性の判定
		if(amax<EPS){
			cout<<"行列は正則ではない。計算終了"<<endl; 
			exit(1);
		}

		p[k]=ip; //ipを配列p[i]に保存

		//行交換
		if(ip!=k)
		{
			for(j=k;j<size;j++)
			{
				temp=A[k][j];
				A[k][j]=A[ip][j];
				A[ip][j]=temp;
			}
		}

		//前進消去
		for(i=k+1;i<size;i++)
		{
			alpha=-A[i][k]/A[k][k];
			A[i][k]=alpha;
			for(j=k+1;j<size;j++)
			{
				A[i][j]+=alpha*A[k][j];
			}
		}
	}

	for(k=0;k<size-1;k++)
	{
		temp=b[k];
		b[k]=b[p[k]];//???
		b[p[k]]=temp;
		//前進代入
		for(i=k+1;i<size;i++)
		{
			b[i]+=A[i][k]*b[k];
		}
	}

	//後退代入・・・A[size-1][size-1]!=0は保証されている？？？
	b[size-1]/=A[size-1][size-1];
	for(k=size-2;k>=0;k--)
	{
		temp=b[k];
		for(j=k+1;j<size;j++)
		{
			temp-=A[k][j]*b[j];
		}
		b[k]=temp/A[k][k];
	}
}


//Cholesky分解法を用いる
//対称正定値行列にのみ適用可能
void modified_cholesky_decomposition(double **A, double *b, const int size)
{
	//修正コレスキー分解
	int i, j, k;
	double tmp;

	for(i=2;i<=size;i++)
	{
		for(j=1;j<=i-1;j++)
		{
			tmp=0.0;
			for(k=1;k<=j-1;k++)
			{
				tmp+=A[i][k]*A[k][k]*A[j][k];
			}
			A[i][j]=(A[i][j]-tmp)/A[i][j];
		}
		tmp=0.0;
		for(k=1;k<=j-1;k++)
		{
			tmp+=A[i][k]*A[i][k]*A[k][k];
		}
		A[i][i]=A[i][i]-tmp;
	}

	//コレスキー分解による行列解法
	//LDy=b
	b[1]=b[1]/A[1][1];
	for(i=2;i<=size;i++)
	{
		tmp=0.0;
		for(j=1;j<=i-1;j++)
		{
			tmp+=A[j][j]*A[i][j]*b[j];
		}
		b[i]=(b[i]-tmp)/A[i][i];
	}

	//L^t x=y
	for(i=size-1;i>=1;i--)
	{
		tmp=0.0;
		for(j=i+1;j<=size;j++)
		{
			tmp+=A[j][i]*b[j];
		}
		b[i]=b[i]-tmp;
	}
	//*bは解ベクトル
}

//CG法（共役勾配法）
void conjugate_gradient(double **A, double *b, double *x, const int size)
{
	const double EPSILON=1e-8;
	const int MAX_ITERATION=100;
	double *r=new double[size];
	double *p=new double[size];
	double *tmp=new double[size];
	double alpha, beta, work;

	int iteration=0;//反復回数

	//tmp=A*bを計算
	matrix_vector_product(A, b, tmp, size, size);

	for(int i=0;i<size;i++)
	{
		p[i]=b[i]-tmp[i];
		r[i]=p[i];
	}

	do{
		//alphaの計算
		matrix_vector_product(A, p, tmp, size, size);
		
		work=alpha=0.0;
		for(int i=0;i<size;i++) work+=p[i]*tmp[i];//内積
		for(int i=0;i<size;i++) alpha+=p[i]*r[i];
		alpha/=work;

		//x_(k+1)とr_(k+1)の計算
		for(int i=0;i<size;i++) x[i]+=alpha*p[i];
		for(int i=0;i<size;i++) r[i]-=alpha*tmp[i];

		//収束判定
		if(vector_norm1(r, 0, size)<EPSILON){
			break;//do-whileを抜ける
		}else{
			iteration++;
			beta=0.0;
			for(int i=0;i<size;i++) beta+=r[i]*tmp[i];
			beta*=-1;
			beta/=work;
			for(int i=0;i<size;i++) p[i]=r[i]+beta*p[i];
		}

	}while(iteration<MAX_ITERATION);

	if(iteration==MAX_ITERATION){
		std::cout<<"Iteration reached MAX_ITERATION: "<<MAX_ITERATION<<endl;
		//try-catch
	}else{
		//std::cout<<"Iteration: "<<iteration<<endl;
	}

	delete [] r;
	delete [] p;
	delete [] tmp;
}

//行列とベクトルの積
void matrix_vector_product(double **A, double *b, double *x, int row, int col)
{
	for(int i=0;i<row;i++)
	{
		double temp=0.0;
		for(int j=0;j<col;j++)
		{
			temp+=A[i][j]*b[j];
		}
		x[i]=temp;
	}
}