int CholeskyUppertri(double *D, int n, double *L) {//D为输出矩阵，n为矩阵的维度，L是输出的上三角矩阵
	double suml;int K = 0;
	for (int i = 0;i < n*n;++i)
		L[i] = 0;
	L[0] = sqrt(D[0]);
	for (int i = 0;i < n;++i) L[i * n] = D[i * n] / L[0];
	for (K = 1;K < n;K++) {
		Lkk(D, n, L, &K);
		Lik(D, n, L, &K);
	}
	matTran(L, n);//转置
	return 0;

}
int Lkk(double *D, int n, double *L, int *K) {
	double suml = 0;
	for (int i = 0;i <= *K - 1;++i) {
		suml += L[*K*n + i] * L[*K*n + i];
	}
	L[*K*n + *K] = sqrt(D[*K*n + *K] - suml);
	return 0;
}
int Lik(double *D, int n, double *L, int *K) {
	double suml;
	for (int i = *K + 1;i < n;++i) {
		suml = 0;
		for (int j = 0;j <= *K - 1;++j) suml += L[i*n + j] * L[*K*n + j];
		L[i*n + *K] = (D[i*n + *K] - suml) / L[*K*n + *K];
	}
	return 0;
}
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
