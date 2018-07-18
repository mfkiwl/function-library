#include <string.h>

/* new integer matrix /整数矩阵的内存分配 ----------------------------------------------------------
* allocate memory of integer matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern int *imat(int n, int m)
{
	int *p;

	if (n <= 0 || m <= 0) return NULL;
	if (!(p = (int *)malloc(sizeof(int)*n*m))) {
		printf("integer matrix memory allocation error: n=%d,m=%d\n", n, m);
	}
	return p;
}
/** new matrix  /矩阵内存的分配------------------------------------------------------------------
* allocate memory of matrix  /矩阵内存的分配
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *mat(int n, int m)
{
	double *p;

	if (n <= 0 || m <= 0) return NULL;
	if (!(p = (double *)malloc(sizeof(double)*n*m))) {
		printf("matrix memory allocation error: n=%d,m=%d\n", n, m);
	}
	return p;
}
/* copy matrix /复制矩阵----------------------------------------------------------------
* copy matrix
* args   : double *A        O   destination matrix A (n x m)
*          double *B        I   source matrix B (n x m)
*          int    n,m       I   number of rows and columns of matrix
* return : none
*-----------------------------------------------------------------------------*/
extern void matcpy(double *A, const double *B, int n, int m)
{
	memcpy(A, B, sizeof(double)*n*m);
}

/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double *A, int n, int *indx, double *d)
{
	double big, s, tmp, *vv = mat(n, 1);
	int i, imax = 0, j, k;

	*d = 1.0;
	for (i = 0;i<n;i++) {
		big = 0.0; for (j = 0;j<n;j++) if ((tmp = fabs(A[i + j*n]))>big) big = tmp;
		if (big>0.0) vv[i] = 1.0 / big; else { free(vv); return -1; }
	}
	for (j = 0;j<n;j++) {
		for (i = 0;i<j;i++) {
			s = A[i + j*n]; for (k = 0;k<i;k++) s -= A[i + k*n] * A[k + j*n]; A[i + j*n] = s;
		}
		big = 0.0;
		for (i = j;i<n;i++) {
			s = A[i + j*n]; for (k = 0;k<j;k++) s -= A[i + k*n] * A[k + j*n]; A[i + j*n] = s;
			if ((tmp = vv[i] * fabs(s)) >= big) { big = tmp; imax = i; }
		}
		if (j != imax) {
			for (k = 0;k<n;k++) {
				tmp = A[imax + k*n]; A[imax + k*n] = A[j + k*n]; A[j + k*n] = tmp;
			}
			*d = -(*d); vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (A[j + j*n] == 0.0) { free(vv); return -1; }
		if (j != n - 1) {
			tmp = 1.0 / A[j + j*n]; for (i = j + 1;i<n;i++) A[i + j*n] *= tmp;
		}
	}
	free(vv);
	return 0;
}

/* LU back-substitution ------------------------------------------------------*/
static void lubksb(const double *A, int n, const int *indx, double *b)
{
	double s;
	int i, ii = -1, ip, j;

	for (i = 0;i<n;i++) {
		ip = indx[i]; s = b[ip]; b[ip] = b[i];
		if (ii >= 0) for (j = ii;j<i;j++) s -= A[i + j*n] * b[j]; else if (s) ii = i;
		b[i] = s;
	}
	for (i = n - 1;i >= 0;i--) {
		s = b[i]; for (j = i + 1;j<n;j++) s -= A[i + j*n] * b[j]; b[i] = s / A[i + i*n];
	}
}
/* inverse of matrix ---------------------------------------------------------*/
extern int matinv(double *A, int n,double *C)//A为输入矩阵，C为输出矩阵
{
	
	double d, *B;
	int i, j, *indx;

	indx = imat(n, 1); B = mat(n, n); matcpy(B, A, n, n);

	if (ludcmp(B, n, indx, &d)) { free(indx); free(B); return -1; }
	for (j = 0;j<n;j++) {
		for (i = 0;i<n;i++) A[i + j*n] = 0.0; A[j + j*n] = 1.0;
		lubksb(B, n, indx, A + j*n);
	}

	
	free(indx); free(B);
	return 0;
}
