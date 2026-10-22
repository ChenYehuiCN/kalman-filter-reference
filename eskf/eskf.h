/* eskf.h */

#ifndef ESKF_H
#define ESKF_H

#ifndef ESKF_N
#define ESKF_N 1
#endif
#ifndef ESKF_N_NOM
#define ESKF_N_NOM 1
#endif
#ifndef ESKF_NW
#define ESKF_NW 1
#endif
#ifndef ESKF_M
#define ESKF_M 1
#endif

struct eskf_filter {
	void (*compute_Fx)(float *Fx, float *x_nom, void *params);
	void (*compute_Fi)(float *Fi, float *x_nom, void *params);
	void (*nominal_propagate)(float *x_nom, void *params);
	void (*compute_Hx)(float *Hx, float *x_nom, void *params);
	void (*measure_model)(float *z_pred, float *x_nom, void *params);
	void (*inject)(float *x_nom, float *dx, void *params);
	void (*compute_G)(float *G, float *dx, void *params);

	float Fx[ESKF_N * ESKF_N];
	float Fi[ESKF_N * ESKF_NW];
	float Hx[ESKF_M * ESKF_N];
	float G[ESKF_N * ESKF_N];

	float x_nom[ESKF_N_NOM * 1];
	float dx[ESKF_N * 1];
	float P[ESKF_N * ESKF_N];

	float _I_nn[ESKF_N * ESKF_N];
	float _tmp_nn[ESKF_N * ESKF_N];
	float _tmp_nn_2[ESKF_N * ESKF_N];
	float _tmp_nw[ESKF_N * ESKF_NW];
	float _tmp_mn[ESKF_M * ESKF_N];
	float _tmp_nm[ESKF_N * ESKF_M];
	float _S[ESKF_M * ESKF_M];
	float _K[ESKF_N * ESKF_M];
	float _y[ESKF_M * 1];
};

void eskf_init(
	struct eskf_filter *f, 
	void (*compute_Fx)(float *Fx, float *x_nom, void *params), 
	void (*compute_Fi)(float *Fi, float *x_nom, void *params), 
	void (*nominal_propagate)(float *x_nom, void *params), 
	void (*compute_Hx)(float *Hx, float *x_nom, void *params), 
	void (*measure_model)(float *z_pred, float *x_nom, void *params), 
	void (*inject)(float *x_nom, float *dx, void *params), 
	void (*compute_G)(float *G, float *dx, void *params)
);
void eskf_reset(struct eskf_filter *f, float *x_nom, float *P);
void eskf_predict(struct eskf_filter *f, float *Q, void *params);
void eskf_update(struct eskf_filter *f, float *z, float *R, void *params);

#endif
