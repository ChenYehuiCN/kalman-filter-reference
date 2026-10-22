/*
 * kf_tracking -- linear Kalman filter example.
 *
 * A 2D target flies a figure-eight.  The filter uses a constant-acceleration
 * (CA) model driven by a known jerk input and observes the target position
 * with a noisy GPS-like sensor.  The sensor noise is intentionally raised
 * from t = 12 s to 16 s to show that R may change at every step.
 *
 * State (KF_N = 6):  x = [px, py, vx, vy, ax, ay]^T
 * Input (KF_C = 2):  u = [jx, jy]^T        (jerk)
 * Meas. (KF_M = 2):  z = [px, py]^T
 *
 * Build (dimensions are checked below):
 *     cc -I../../kf -DKF_N=6 -DKF_M=2 -DKF_C=2 main.c ../../kf/kf.c -lm
 *
 * The true state is simulated with the same discrete model plus a small
 * unmodelled acceleration increment w ~ N(0, ACCEL_NOISE^2).  The filter's Q
 * is built from exactly that disturbance model, so the covariance is
 * statistically meaningful.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kf.h"

#if KF_N != 6 || KF_M != 2 || KF_C != 2
#error "Compile this example with KF_N=6, KF_M=2 and KF_C=2"
#endif

#define PI_F 3.14159265358979323846f

#define T_END 30.0f
#define DT 0.02f
#define STEPS ((int)(T_END / DT + 0.5f))

/* Figure-eight: px = A sin(w t), py = B sin(2 w t). */
#define TRAJ_A 80.0f
#define TRAJ_B 45.0f
#define TRAJ_W 0.35f

/* Standard deviation of the acceleration increment per step [m/s^2]. */
#define ACCEL_NOISE 0.0005f

/* ------------------------------------------------------------------ */
/* Small deterministic RNG (xorshift64* + Box-Muller).                */
/* ------------------------------------------------------------------ */

struct rng {
	uint64_t state;
};

static uint64_t rng_u64(struct rng *r)
{
	uint64_t x = r->state;

	x ^= x >> 12;
	x ^= x << 25;
	x ^= x >> 27;
	r->state = x;
	return x * 2685821657736338717ULL;
}

static float rng_uniform(struct rng *r)
{
	/* 24 random bits -> [0, 1). */
	return (float)((double)(rng_u64(r) >> 40) * (1.0 / 16777216.0));
}

static float rng_normal(struct rng *r)
{
	float u1 = rng_uniform(r);
	float u2 = rng_uniform(r);

	if (u1 < 1.0e-7f)
		u1 = 1.0e-7f;
	return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * PI_F * u2);
}

/* ------------------------------------------------------------------ */
/* Truth and control input.                                           */
/* ------------------------------------------------------------------ */

static void truth_analytic(float t, float *x)
{
	float s1 = sinf(TRAJ_W * t);
	float c1 = cosf(TRAJ_W * t);
	float s2 = sinf(2.0f * TRAJ_W * t);
	float c2 = cosf(2.0f * TRAJ_W * t);

	x[0] = TRAJ_A * s1;
	x[1] = TRAJ_B * s2;
	x[2] = TRAJ_A * TRAJ_W * c1;
	x[3] = 2.0f * TRAJ_B * TRAJ_W * c2;
	x[4] = -TRAJ_A * TRAJ_W * TRAJ_W * s1;
	x[5] = -4.0f * TRAJ_B * TRAJ_W * TRAJ_W * s2;
}

static void jerk_analytic(float t, float *u)
{
	float w2 = TRAJ_W * TRAJ_W;
	float w3 = w2 * TRAJ_W;

	u[0] = -TRAJ_A * w3 * cosf(TRAJ_W * t);
	u[1] = -8.0f * TRAJ_B * w3 * cosf(2.0f * TRAJ_W * t);
}

/* ------------------------------------------------------------------ */
/* Model matrices.                                                    */
/* ------------------------------------------------------------------ */

/*
 * F = [ I  dt I  dt^2/2 I ]
 *     [ 0   I     dt I    ]
 *     [ 0   0      I      ]
 *
 * B = [ dt^3/6 I ; dt^2/2 I ; dt I ]
 */
static void build_FB(float *F, float *B, float dt)
{
	int i;

	memset(F, 0, sizeof(float) * KF_N * KF_N);
	memset(B, 0, sizeof(float) * KF_N * KF_C);

	for (i = 0; i < 2; i++) {
		F[(0 + i) * KF_N + (0 + i)] = 1.0f;
		F[(0 + i) * KF_N + (2 + i)] = dt;
		F[(0 + i) * KF_N + (4 + i)] = 0.5f * dt * dt;

		F[(2 + i) * KF_N + (2 + i)] = 1.0f;
		F[(2 + i) * KF_N + (4 + i)] = dt;

		F[(4 + i) * KF_N + (4 + i)] = 1.0f;

		B[(0 + i) * KF_C + i] = dt * dt * dt / 6.0f;
		B[(2 + i) * KF_C + i] = 0.5f * dt * dt;
		B[(4 + i) * KF_C + i] = dt;
	}
}

/*
 * G maps the acceleration disturbance w = [wx, wy]^T into the state:
 *     x_{k+1} = F x_k + B u_k + G w_k
 * A velocity jump w*dt and position jump w*dt^2/2 are included so that the
 * disturbance is applied before the acceleration state is carried forward.
 */
static void build_G(float *G, float dt)
{
	int i;

	memset(G, 0, sizeof(float) * KF_N * KF_C);
	for (i = 0; i < 2; i++) {
		G[(0 + i) * KF_C + i] = 0.5f * dt * dt;
		G[(2 + i) * KF_C + i] = dt;
		G[(4 + i) * KF_C + i] = 1.0f;
	}
}

/* Q = G * (ACCEL_NOISE^2 I) * G^T. */
static void build_Q(float *Q, float dt)
{
	float G[KF_N * KF_C];
	float var = ACCEL_NOISE * ACCEL_NOISE;
	int i, j, k;

	build_G(G, dt);
	for (i = 0; i < KF_N; i++) {
		for (j = 0; j < KF_N; j++) {
			float s = 0.0f;

			for (k = 0; k < KF_C; k++)
				s += G[i * KF_C + k] * G[j * KF_C + k];
			Q[i * KF_N + j] = var * s;
		}
	}
}

static void build_H(float *H)
{
	memset(H, 0, sizeof(float) * KF_M * KF_N);
	H[0 * KF_N + 0] = 1.0f;
	H[1 * KF_N + 1] = 1.0f;
}

static void build_R(float *R, float sigma)
{
	int i;

	memset(R, 0, sizeof(float) * KF_M * KF_M);
	for (i = 0; i < KF_M; i++)
		R[i * KF_M + i] = sigma * sigma;
}

static float measurement_sigma(float t)
{
	/* A short GPS degradation: the filter is told about it through R. */
	if (t >= 12.0f && t < 16.0f)
		return 3.0f;
	return 0.75f;
}

/* ------------------------------------------------------------------ */
/* Truth simulation helpers.                                          */
/* ------------------------------------------------------------------ */

/* x <- F x + B u (F and B are constant for this example). */
static void propagate_deterministic(float *x, const float *F, const float *B,
				    const float *u)
{
	float out[KF_N];
	int i, j;

	for (i = 0; i < KF_N; i++) {
		float s = 0.0f;

		for (j = 0; j < KF_N; j++)
			s += F[i * KF_N + j] * x[j];
		for (j = 0; j < KF_C; j++)
			s += B[i * KF_C + j] * u[j];
		out[i] = s;
	}
	memcpy(x, out, sizeof(out));
}

static void propagate_truth(float *x, const float *F, const float *B,
			    const float *G, const float *u, struct rng *rng)
{
	float w[KF_C];
	int i;

	propagate_deterministic(x, F, B, u);

	w[0] = ACCEL_NOISE * rng_normal(rng);
	w[1] = ACCEL_NOISE * rng_normal(rng);
	for (i = 0; i < KF_N; i++)
		x[i] += G[i * KF_C + 0] * w[0] + G[i * KF_C + 1] * w[1];
}

static void make_measurement(float *z, const float *x, float sigma,
			     struct rng *rng)
{
	int i;

	for (i = 0; i < KF_M; i++)
		z[i] = x[i] + sigma * rng_normal(rng);
}

/* ------------------------------------------------------------------ */
/* CSV output.                                                        */
/* ------------------------------------------------------------------ */

static void write_header(FILE *fp)
{
	fprintf(fp,
		"t,true_px,true_py,true_vx,true_vy,true_ax,true_ay,"
		"meas_px,meas_py,"
		"est_px,est_py,est_vx,est_vy,est_ax,est_ay,"
		"Pxx,Pxy,Pyy,Pvx,Pvy,Pax,Pay,"
		"jx,jy,sigma\n");
}

static void write_row(FILE *fp, float t, const float *xt, const float *z,
		      const struct kf_filter *f, const float *u, float sigma)
{
	fprintf(fp,
		"%.4f,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"
		"%.9g,%.9g,"
		"%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"
		"%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,"
		"%.9g,%.9g,%.9g\n",
		(double)t,
		(double)xt[0], (double)xt[1],
		(double)xt[2], (double)xt[3],
		(double)xt[4], (double)xt[5],
		(double)z[0], (double)z[1],
		(double)f->x[0], (double)f->x[1],
		(double)f->x[2], (double)f->x[3],
		(double)f->x[4], (double)f->x[5],
		(double)f->P[0 * KF_N + 0], (double)f->P[0 * KF_N + 1],
		(double)f->P[1 * KF_N + 1],
		(double)f->P[2 * KF_N + 2], (double)f->P[3 * KF_N + 3],
		(double)f->P[4 * KF_N + 4], (double)f->P[5 * KF_N + 5],
		(double)u[0], (double)u[1], (double)sigma);
}

/* ------------------------------------------------------------------ */
/* main                                                               */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
	const char *path = (argc > 1) ? argv[1] : "kf_tracking.csv";
	struct kf_filter f;
	struct rng rng = { 0x9e3779b97f4a7c15ULL };
	float F[KF_N * KF_N], B[KF_N * KF_C], G[KF_N * KF_C];
	float H[KF_M * KF_N], Q[KF_N * KF_N], R[KF_M * KF_M];
	float P0[KF_N * KF_N];
	float x_true[KF_N], x0[KF_N], z[KF_M], u[KF_C];
	float sigma;
	int k;
	FILE *fp;

	fp = fopen(path, "w");
	if (fp == NULL) {
		perror(path);
		return 1;
	}

	build_FB(F, B, DT);
	build_G(G, DT);
	build_Q(Q, DT);
	build_H(H);

	/* Ground-truth initial state. */
	truth_analytic(0.0f, x_true);

	/* Deliberately offset initial estimate: a rough prior. */
	x0[0] = x_true[0] + 3.0f;
	x0[1] = x_true[1] - 2.0f;
	x0[2] = x_true[2] + 2.0f;
	x0[3] = x_true[3] - 3.0f;
	x0[4] = x_true[4] + 0.5f;
	x0[5] = x_true[5] - 0.6f;

	memset(P0, 0, sizeof(P0));
	P0[0 * KF_N + 0] = 9.0f;
	P0[1 * KF_N + 1] = 9.0f;
	P0[2 * KF_N + 2] = 9.0f;
	P0[3 * KF_N + 3] = 9.0f;
	P0[4 * KF_N + 4] = 2.0f;
	P0[5 * KF_N + 5] = 2.0f;

	kf_init(&f);
	kf_reset(&f, x0, P0);

	write_header(fp);

	/* Row at t = 0: the measurement is shown but not yet consumed. */
	sigma = measurement_sigma(0.0f);
	build_R(R, sigma);
	make_measurement(z, x_true, sigma, &rng);
	jerk_analytic(0.0f, u);
	write_row(fp, 0.0f, x_true, z, &f, u, sigma);

	for (k = 1; k <= STEPS; k++) {
		float t = (float)k * DT;

		/* The control is held over [t-dt, t], so sample at the midpoint. */
		jerk_analytic(t - 0.5f * DT, u);
		propagate_truth(x_true, F, B, G, u, &rng);

		sigma = measurement_sigma(t);
		build_R(R, sigma);
		make_measurement(z, x_true, sigma, &rng);

		kf_predict(&f, u, F, B, Q);
		kf_update(&f, z, H, R);
		write_row(fp, t, x_true, z, &f, u, sigma);
	}

	fclose(fp);

	printf("KF tracking: %d steps written to %s\n", STEPS + 1, path);
	printf("final position error: %.3f m\n",
	       (double)hypotf(f.x[0] - x_true[0], f.x[1] - x_true[1]));

	return 0;
}
