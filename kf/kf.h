/* kf.h */

#ifndef KF_H
#define KF_H

#ifndef KF_N
#define KF_N 1
#endif
#ifndef KF_M
#define KF_M 1
#endif
#ifndef KF_C
#define KF_C 1
#endif

struct kf_filter {
	float x[KF_N * 1];
	float P[KF_N * KF_N];

	float _I_nn[KF_N * KF_N];
	float _tmp_nn[KF_N * KF_N];
	float _tmp_nn_2[KF_N * KF_N];
	float _tmp_mn[KF_M * KF_N];
	float _tmp_nm[KF_N * KF_M];
	float _S[KF_M * KF_M];
	float _K[KF_N * KF_M];
	float _y[KF_M * 1];
	float _tmp_n1[KF_N * 1];
	float _tmp_n1_2[KF_N * 1];
};

void kf_init(struct kf_filter *f);
void kf_reset(struct kf_filter *f, float *x, float *P);
void kf_predict(struct kf_filter *f, float *u, float *F, float *B, float *Q);
void kf_update(struct kf_filter *f, float *z, float *H, float *R);

#endif
