/* kf.c */

#include <string.h>
#include <math.h>
#include "kf.h"

static void mat_add(float *C, float *A, float *B, int count)
{
	int i;
	for (i = 0; i < count; i++)
		C[i] = A[i] + B[i];
}

static void mat_sub(float *C, float *A, float *B, int count)
{
	int i;
	for (i = 0; i < count; i++)
		C[i] = A[i] - B[i];
}

static void mat_mul(float *C, float *A, float *B, int ar, int ac, int bc)
{
	int i, j, k;
	for (i = 0; i < ar; i++)
		for (j = 0; j < bc; j++) {
			float s = 0;
			for (k = 0; k < ac; k++)
				s += A[i * ac + k] * B[k * bc + j];
			C[i * bc + j] = s;
		}
}

static void mat_mul_bt(float *C, float *A, float *B, int ar, int ac, int br)
{
	int i, j, k;
	for (i = 0; i < ar; i++)
		for (j = 0; j < br; j++) {
			float s = 0;
			for (k = 0; k < ac; k++)
				s += A[i * ac + k] * B[j * ac + k];
			C[i * br + j] = s;
		}
}

static void mat_transpose(float *At, float *A, int rows, int cols)
{
	int i, j;
	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
			At[j * rows + i] = A[i * cols + j];
}

static void mat_chol(float *S, int m)
{
	int i, j, k;
	for (j = 0; j < m; j++) {
		float d = S[j * m + j];
		for (k = 0; k < j; k++)
			d -= S[j * m + k] * S[j * m + k];
		d = sqrtf(d);
		S[j * m + j] = d;
		for (i = j + 1; i < m; i++) {
			float s = S[i * m + j];
			for (k = 0; k < j; k++)
				s -= S[i * m + k] * S[j * m + k];
			S[i * m + j] = s / d;
		}
	}
}

static void mat_chol_solve(float *L, float *B, int m, int n)
{
	int i, j, k;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++) {
			float s = B[i * n + j];
			for (k = 0; k < i; k++)
				s -= L[i * m + k] * B[k * n + j];
			B[i * n + j] = s / L[i * m + i];
		}
	for (i = m; i-- > 0; )
		for (j = 0; j < n; j++) {
			float s = B[i * n + j];
			for (k = i + 1; k < m; k++)
				s -= L[k * m + i] * B[k * n + j];
			B[i * n + j] = s / L[i * m + i];
		}
}

void kf_init(struct kf_filter *f)
{
	int i;
	memset(f->_I_nn, 0, sizeof(f->_I_nn));
	for (i = 0; i < KF_N; i++)
		f->_I_nn[i * KF_N + i] = 1;
}

void kf_reset(struct kf_filter *f, float *x, float *P)
{
	memcpy(f->x, x, sizeof(f->x));
	memcpy(f->P, P, sizeof(f->P));
}

void kf_predict(struct kf_filter *f, float *u, float *F, float *B, float *Q)
{
	mat_mul(f->_tmp_n1, F, f->x, KF_N, KF_N, 1);
	mat_mul(f->_tmp_n1_2, B, u, KF_N, KF_C, 1);
	mat_add(f->_tmp_n1, f->_tmp_n1, f->_tmp_n1_2, KF_N * 1);
	memcpy(f->x, f->_tmp_n1, sizeof(f->x));

	mat_mul(f->_tmp_nn, F, f->P, KF_N, KF_N, KF_N);
	mat_mul_bt(f->P, f->_tmp_nn, F, KF_N, KF_N, KF_N);
	mat_add(f->P, f->P, Q, KF_N * KF_N);
}

void kf_update(struct kf_filter *f, float *z, float *H, float *R)
{
	mat_mul(f->_tmp_mn, H, f->P, KF_M, KF_N, KF_N);
	mat_mul_bt(f->_S, f->_tmp_mn, H, KF_M, KF_N, KF_M);
	mat_add(f->_S, f->_S, R, KF_M * KF_M);
	mat_chol(f->_S, KF_M);
	mat_chol_solve(f->_S, f->_tmp_mn, KF_M, KF_N);
	mat_transpose(f->_K, f->_tmp_mn, KF_M, KF_N);

	mat_mul(f->_y, H, f->x, KF_M, KF_N, 1);
	mat_sub(f->_y, z, f->_y, KF_M * 1);
	mat_mul(f->_tmp_n1, f->_K, f->_y, KF_N, KF_M, 1);
	mat_add(f->x, f->x, f->_tmp_n1, KF_N * 1);

	mat_mul(f->_tmp_nn, f->_K, H, KF_N, KF_M, KF_N);
	mat_sub(f->_tmp_nn, f->_I_nn, f->_tmp_nn, KF_N * KF_N);
	mat_mul(f->_tmp_nn_2, f->_tmp_nn, f->P, KF_N, KF_N, KF_N);
	mat_mul_bt(f->P, f->_tmp_nn_2, f->_tmp_nn, KF_N, KF_N, KF_N);
	mat_mul(f->_tmp_nm, f->_K, R, KF_N, KF_M, KF_M);
	mat_mul_bt(f->_tmp_nn_2, f->_tmp_nm, f->_K, KF_N, KF_M, KF_N);
	mat_add(f->P, f->P, f->_tmp_nn_2, KF_N * KF_N);
}
