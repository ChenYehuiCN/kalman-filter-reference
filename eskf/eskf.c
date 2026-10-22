/* eskf.c */

#include <string.h>
#include <math.h>
#include "eskf.h"

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

void eskf_init(
	struct eskf_filter *f, 
	void (*compute_Fx)(float *Fx, float *x_nom, void *params), 
	void (*compute_Fi)(float *Fi, float *x_nom, void *params), 
	void (*nominal_propagate)(float *x_nom, void *params), 
	void (*compute_Hx)(float *Hx, float *x_nom, void *params), 
	void (*measure_model)(float *z_pred, float *x_nom, void *params), 
	void (*inject)(float *x_nom, float *dx, void *params), 
	void (*compute_G)(float *G, float *dx, void *params)
)
{
	int i;
	f->compute_Fx = compute_Fx;
	f->compute_Fi = compute_Fi;
	f->nominal_propagate = nominal_propagate;
	f->compute_Hx = compute_Hx;
	f->measure_model = measure_model;
	f->inject = inject;
	f->compute_G = compute_G;
	memset(f->dx, 0, sizeof(f->dx));
	memset(f->_I_nn, 0, sizeof(f->_I_nn));
	for (i = 0; i < ESKF_N; i++)
		f->_I_nn[i * ESKF_N + i] = 1;
}

void eskf_reset(struct eskf_filter *f, float *x_nom, float *P)
{
	memcpy(f->x_nom, x_nom, sizeof(f->x_nom));
	memcpy(f->P, P, sizeof(f->P));
	memset(f->dx, 0, sizeof(f->dx));
}

void eskf_predict(struct eskf_filter *f, float *Q, void *params)
{
	f->compute_Fx(f->Fx, f->x_nom, params);
	f->compute_Fi(f->Fi, f->x_nom, params);

	memset(f->dx, 0, sizeof(f->dx));
	f->nominal_propagate(f->x_nom, params);

	mat_mul(f->_tmp_nn, f->Fx, f->P, ESKF_N, ESKF_N, ESKF_N);
	mat_mul_bt(f->P, f->_tmp_nn, f->Fx, ESKF_N, ESKF_N, ESKF_N);
	mat_mul(f->_tmp_nw, f->Fi, Q, ESKF_N, ESKF_NW, ESKF_NW);
	mat_mul_bt(f->_tmp_nn, f->_tmp_nw, f->Fi, ESKF_N, ESKF_NW, ESKF_N);
	mat_add(f->P, f->P, f->_tmp_nn, ESKF_N * ESKF_N);
}

void eskf_update(struct eskf_filter *f, float *z, float *R, void *params)
{
	f->compute_Hx(f->Hx, f->x_nom, params);
	mat_mul(f->_tmp_mn, f->Hx, f->P, ESKF_M, ESKF_N, ESKF_N);
	mat_mul_bt(f->_S, f->_tmp_mn, f->Hx, ESKF_M, ESKF_N, ESKF_M);
	mat_add(f->_S, f->_S, R, ESKF_M * ESKF_M);
	mat_chol(f->_S, ESKF_M);
	mat_chol_solve(f->_S, f->_tmp_mn, ESKF_M, ESKF_N);
	mat_transpose(f->_K, f->_tmp_mn, ESKF_M, ESKF_N);

	f->measure_model(f->_y, f->x_nom, params);
	mat_sub(f->_y, z, f->_y, ESKF_M * 1);
	mat_mul(f->dx, f->_K, f->_y, ESKF_N, ESKF_M, 1);

	mat_mul(f->_tmp_nn, f->_K, f->Hx, ESKF_N, ESKF_M, ESKF_N);
	mat_sub(f->_tmp_nn, f->_I_nn, f->_tmp_nn, ESKF_N * ESKF_N);
	mat_mul(f->_tmp_nn_2, f->_tmp_nn, f->P, ESKF_N, ESKF_N, ESKF_N);
	mat_mul_bt(f->P, f->_tmp_nn_2, f->_tmp_nn, ESKF_N, ESKF_N, ESKF_N);
	mat_mul(f->_tmp_nm, f->_K, R, ESKF_N, ESKF_M, ESKF_M);
	mat_mul_bt(f->_tmp_nn_2, f->_tmp_nm, f->_K, ESKF_N, ESKF_M, ESKF_N);
	mat_add(f->P, f->P, f->_tmp_nn_2, ESKF_N * ESKF_N);

	f->inject(f->x_nom, f->dx, params);

	f->compute_G(f->G, f->dx, params);
	mat_mul(f->_tmp_nn, f->G, f->P, ESKF_N, ESKF_N, ESKF_N);
	mat_mul_bt(f->P, f->_tmp_nn, f->G, ESKF_N, ESKF_N, ESKF_N);

	memset(f->dx, 0, sizeof(f->dx));
}
