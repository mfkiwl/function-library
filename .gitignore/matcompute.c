//矩阵转置
int matTran(double*A, int n) {
	double *B;
	B = (double*)malloc(sizeof(double)*n*n);
	for (int i = 0;i < n;i++)
		for (int j = 0;j < n;++j) {
			B[j*n + i] = A[i*n + j];
		}
	for (int i = 0;i < n;i++)
		for (int j = 0;j < n;++j) {
			A[i*n + j] = B[i*n + j];
		}
	free(B);
	return 0;
}
//矩阵相乘
int matmul(double *A, int n, int m, int k, double *B, double *C) {//A(nxm) B(mxk) C(nxk)
	for (int i = 0;i < n;++i)
		for (int j = 0;j < k;++j) {
			double s = 0;
			for (int g = 0;g < m;g++) {
				s += A[i*m + g] * B[j + k*g];
			}
			C[i*k + j] = s;
		}
	return 0;
}
//矩阵乘以常数
int kmat(double *A,const int n,const int m,double k) {//A(nxm),常数k
	for (int i = 0;i < n;++i)
		for (int j = 0;j < m;++j)
			A[i*m + j] = k*A[i*m + j];
	return 0;
}

//矩阵加减
void matps(double *A, const int n,const int m, double *B, double *C,const char a[2]) {//C=A+B or C=A-B
	if (a == "+") {
		for (int i = 0;i < n;++i)
			for (int j = 0;j < m;++j)
				C[i*m + j] = A[i*m + j] + B[i*m + j];
	}

	if (a == "-") {
		for (int i = 0;i < n;++i)
			for (int j = 0;j < m;++j)
				C[i*m + j] = A[i*m + j] - B[i*m + j];
	}


}
