//連立一次方程式の解法
void pivot_gauss(double **A, double *b, const int size, bool pivot_check);//ガウスの消去法
void lu_decomposition(double **A, double *b, size_t size, bool pivot_check);//LU分解解法
void lu_decomposition(vector<vector<double>> &A, vector<double> &b);

//Choresky分解
void modified_cholesky_decomposition(double **A, double *b, int size);

//共役勾配法
void conjugate_gradient(double **A, double *b, double *x, int size);
void matrix_vector_product(double **A, double *b, double *x, int row, int col);

//ベクトルノルム
double vector_norm1(double*v, const int start, const int end);
double vector_norm2(double*v, const int start, const int end);
