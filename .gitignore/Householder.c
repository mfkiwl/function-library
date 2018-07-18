//需要的矩阵函数------------------------------------------------------------------------------------------
int matinv(double *A, int n);
int matmul(double *A, int n, int m, int k, double *B, double *C);
int kmat(double *A, int n, int m, double k);
int matps(double *A, const int n, const int m, double *B, double *C, const char a[2]);

//Householder变换---------------------------------------------------------------------------------------
int Housholder(double *A, const int m, const int n, double *b,int inter) {//A(mxn),b(mx1),inter<=n为迭代次数
	const int k = inter-1;//矩阵下标k
	double sigk, betak;
	double *u, *y,z;
	u = (double*)malloc(sizeof(double)*m);
	y = (double*)malloc(sizeof(double)*n);


	double sumA2 = 0;
	for (int i = k;i < m;++i)
		sumA2 += A[i*n + k] * A[i*n + k];
	sigk = sqrt(sumA2);
	if (A[k*n + k] < 0)
		sigk = -sigk;
	betak = 1 / (sigk*(sigk+A[k*n + k]));
	
	for (int i = 0;i < m;++i) {
		if (i < k) u[i] = 0;
		else if (i == k) u[i] = sigk + A[k*n + k];
		else u[i] = A[i*n + k];
	}
	double *Aj;
	Aj = (double*)malloc(sizeof(double)*m);
	double uAj[1] = { 0 };

	for (int j = 0;j < n;++j) {
		if (j < k) y[j] = 0;
		else if (j == k) y[j] =1;
		else {
			for (int i = 0;i < m;++i)Aj[i] = A[j+i*n];
			matmul(u, 1, m, 1, Aj, uAj);
			y[j] = betak*uAj[0];
		}
	}

	double ub[1] = {0};
	matmul(u, 1, m, 1, b, ub);
	z = betak*ub[0];
	double *uy;
	uy = (double*)malloc(sizeof(double)*n*m);
	matmul(u, m, 1, n, y, uy);
	matps(A, m, n, uy, A, "-");
	kmat(u, m, 1, z);
	matps(b, m, 1, u, b, "-");

	return ++inter;
}
