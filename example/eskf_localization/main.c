/*
 * eskf_localization -- error-state Kalman filter example.
 *
 * A differential-drive robot drives through a field of three known landmarks.
 * The nonlinear nominal dynamics are those of a unicycle; a depth-camera-like
 * sensor observes the position of every landmark in the robot body frame.
 * The ESKF is corrected at 10 Hz while the nominal state is propagated at
 * 100 Hz, so predict and update are deliberately called at different rates.
 *
 * Nominal state (ESKF_N_NOM = 3):  x = [px, py, yaw]^T
 * Error state   (ESKF_N = 3):      dx = [dpx, dpy, dyaw]^T
 * Control                          u = [v, omega]^T
 * Measurement   (ESKF_M = 6):      z = [body_x_1, body_y_1, ..., body_x_3, body_y_3]^T
 *
 * Build (dimensions are checked below):
 *     cc -I../../eskf -DESKF_N=3 -DESKF_N_NOM=3 -DESKF_M=6 -DESKF_NW=2 \
 *        main.c ../../eskf/eskf.c -lm
 *
 * The control vector and time step are carried in the `params` structure,
 * because the callbacks and eskf_predict() no longer take u/dt directly.
 *
 * The filter callbacks demonstrate every hook expected by eskf_init():
 * nonlinear propagation, Fx, Fi, measurement model, Hx, error injection and
 * the injection Jacobian G (identity here because the injection is additive
 * apart from angle wrapping).
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "eskf.h"

#if ESKF_N != 3 || ESKF_N_NOM != 3 || ESKF_M != 6 || ESKF_NW != 2
#error "Compile this example with ESKF_N=3, ESKF_N_NOM=3, ESKF_M=6 and ESKF_NW=2"
#endif
#if (ESKF_M % 2) != 0
#error "ESKF_M must be twice the number of 2D landmarks"
#endif

#define NUM_LANDMARKS (ESKF_M / 2)
#define CONTROL_DIM 2

#define PI_F 3.14159265358979323846f

#define T_END 60.0f
#define DT 0.01f
#define STEPS ((int)(T_END / DT + 0.5f))
#define MEAS_EVERY 10 /* 100 Hz prediction, 10 Hz landmark update */

#define SIGMA_V 0.10f    /* wheel speed noise [m/s] */
#define SIGMA_W 0.04f    /* yaw rate noise [rad/s] */
#define SIGMA_MEAS 0.50f /* body-frame landmark noise [m] */

/*
 * Slight inflation of the filter's process noise.  The truth uses SIGMA_V/W
 * directly; the extra factor absorbs linearization error in Fi and makes the
 * covariance a little conservative rather than over-confident.
 */
#ifndef Q_SCALE
#define Q_SCALE 2.0f
#endif

struct model {
	float landmark[NUM_LANDMARKS][2];

	/* Control and time step are inputs to the nonlinear model callbacks. */
	float u[CONTROL_DIM];
	float dt;
};

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
/* Unicycle model and its Jacobians.                                  */
/* ------------------------------------------------------------------ */

static float wrap_pi(float angle)
{
	while (angle > PI_F)
		angle -= 2.0f * PI_F;
	while (angle < -PI_F)
		angle += 2.0f * PI_F;
	return angle;
}

/*
 * Exact constant-twist integration of
 *     px' = v cos(yaw),  py' = v sin(yaw),  yaw' = omega
 * over one interval.  The small-angle branch avoids 0/0 and, more importantly,
 * avoids catastrophic cancellation in single precision when |omega| * dt is
 * tiny (the difference of two nearby cosines loses most of its bits).
 *
 * Safe for in-place operation (propagate(x, u, dt, x)).
 */
static void propagate(const float *state, const float *u, float dt, float *out)
{
	float px = state[0];
	float py = state[1];
	float yaw = state[2];
	float v = u[0];
	float w = u[1];
	float yaw2 = yaw + w * dt;

	if (fabsf(w) * dt > 1.0e-3f) {
		out[0] = px + (v / w) * (sinf(yaw2) - sinf(yaw));
		out[1] = py + (v / w) * (cosf(yaw) - cosf(yaw2));
	} else {
		out[0] = px + v * cosf(yaw) * dt;
		out[1] = py + v * sinf(yaw) * dt;
	}
	out[2] = wrap_pi(yaw2);
}

/*
 * Fx = d f / d x  (3x3),  Fi = d f / d u  (3x2), evaluated at the current
 * nominal state and control.  The formulas are the analytic derivative of the
 * exact unicycle propagation above.  The same small-angle threshold is used
 * so that propagate() and its linearization switch branches together.
 */
static void motion_jacobians(const float *state, const float *u, float dt,
			     float *Fx, float *Fi)
{
	float v = u[0];
	float w = u[1];
	float yaw = state[2];
	float dt2 = dt * dt;

	memset(Fx, 0, sizeof(float) * ESKF_N * ESKF_N);
	memset(Fi, 0, sizeof(float) * ESKF_N * ESKF_NW);
	Fx[0] = 1.0f;
	Fx[4] = 1.0f;
	Fx[8] = 1.0f;

	if (fabsf(w) * dt > 1.0e-3f) {
		float yaw2 = yaw + w * dt;
		float s1 = sinf(yaw);
		float c1 = cosf(yaw);
		float s2 = sinf(yaw2);
		float c2 = cosf(yaw2);

		/* d(position)/d(yaw). */
		Fx[0 * ESKF_N + 2] = (v / w) * (c2 - c1);
		Fx[1 * ESKF_N + 2] = (v / w) * (s2 - s1);

		/* d(position, yaw)/d(v, w). */
		Fi[0 * ESKF_NW + 0] = (s2 - s1) / w;
		Fi[0 * ESKF_NW + 1] =
			v * (w * dt * c2 - (s2 - s1)) / (w * w);
		Fi[1 * ESKF_NW + 0] = (c1 - c2) / w;
		Fi[1 * ESKF_NW + 1] =
			v * (w * dt * s2 - (c1 - c2)) / (w * w);
		Fi[2 * ESKF_NW + 1] = dt;
	} else {
		/* First-order limit as omega -> 0. */
		float s1 = sinf(yaw);
		float c1 = cosf(yaw);

		Fx[0 * ESKF_N + 2] = -v * s1 * dt - 0.5f * v * c1 * w * dt2;
		Fx[1 * ESKF_N + 2] = v * c1 * dt - 0.5f * v * s1 * w * dt2;

		Fi[0 * ESKF_NW + 0] = c1 * dt - 0.5f * s1 * w * dt2;
		Fi[0 * ESKF_NW + 1] = -0.5f * v * s1 * dt2;
		Fi[1 * ESKF_NW + 0] = s1 * dt + 0.5f * c1 * w * dt2;
		Fi[1 * ESKF_NW + 1] = 0.5f * v * c1 * dt2;
		Fi[2 * ESKF_NW + 1] = dt;
	}
}

/*
 * Landmark position in the robot body frame:
 *     [ c  s ] [lx - px]
 *     [-s  c ] [ly - py]
 */
static void predict_landmark(float *obs, const float *state,
			     const float *landmark)
{
	float c = cosf(state[2]);
	float s = sinf(state[2]);
	float dx = landmark[0] - state[0];
	float dy = landmark[1] - state[1];

	obs[0] = c * dx + s * dy;
	obs[1] = -s * dx + c * dy;
}

/* ------------------------------------------------------------------ */
/* ESKF callbacks.                                                    */
/* ------------------------------------------------------------------ */

static void cb_compute_Fx(float *Fx, float *x_nom, void *params)
{
	const struct model *m = (const struct model *)params;
	float Fi[ESKF_N * ESKF_NW];

	motion_jacobians(x_nom, m->u, m->dt, Fx, Fi);
}

static void cb_compute_Fi(float *Fi, float *x_nom, void *params)
{
	const struct model *m = (const struct model *)params;
	float Fx[ESKF_N * ESKF_N];

	motion_jacobians(x_nom, m->u, m->dt, Fx, Fi);
}

static void cb_nominal_propagate(float *x_nom, void *params)
{
	const struct model *m = (const struct model *)params;

	propagate(x_nom, m->u, m->dt, x_nom);
}

static void cb_measure_model(float *z_pred, float *x_nom, void *params)
{
	const struct model *m = (const struct model *)params;
	int i;

	for (i = 0; i < NUM_LANDMARKS; i++)
		predict_landmark(&z_pred[2 * i], x_nom, m->landmark[i]);
}

static void cb_compute_Hx(float *Hx, float *x_nom, void *params)
{
	const struct model *m = (const struct model *)params;
	float c = cosf(x_nom[2]);
	float s = sinf(x_nom[2]);
	int i;

	memset(Hx, 0, sizeof(float) * ESKF_M * ESKF_N);
	for (i = 0; i < NUM_LANDMARKS; i++) {
		float dx = m->landmark[i][0] - x_nom[0];
		float dy = m->landmark[i][1] - x_nom[1];
		float hx = c * dx + s * dy;
		float hy = -s * dx + c * dy;
		int row = 2 * i;

		/* d hx / d(px, py, yaw) */
		Hx[(row + 0) * ESKF_N + 0] = -c;
		Hx[(row + 0) * ESKF_N + 1] = -s;
		Hx[(row + 0) * ESKF_N + 2] = hy;

		/* d hy / d(px, py, yaw) */
		Hx[(row + 1) * ESKF_N + 0] = s;
		Hx[(row + 1) * ESKF_N + 1] = -c;
		Hx[(row + 1) * ESKF_N + 2] = -hx;
	}
}

static void cb_inject(float *x_nom, float *dx, void *params)
{
	(void)params;
	x_nom[0] += dx[0];
	x_nom[1] += dx[1];
	x_nom[2] = wrap_pi(x_nom[2] + dx[2]);
}

/*
 * G = d(inject)/d(dx) at dx = 0.  Position and yaw errors are additive
 * (yaw is only wrapped), so G is the identity.  It is still supplied because
 * the library applies P <- G P G^T after every injection.
 */
static void cb_compute_G(float *G, float *dx, void *params)
{
	int i;

	(void)dx;
	(void)params;
	memset(G, 0, sizeof(float) * ESKF_N * ESKF_N);
	for (i = 0; i < ESKF_N; i++)
		G[i * ESKF_N + i] = 1.0f;
}

/* ------------------------------------------------------------------ */
/* Simulation helpers.                                                */
/* ------------------------------------------------------------------ */

static void control(float t, float *u)
{
	u[0] = 2.0f + 0.6f * sinf(0.17f * t);
	u[1] = 0.55f * sinf(0.23f * t) + 0.20f * sinf(0.61f * t);
}

static void make_measurement(float *z, const float *x_true,
			     const struct model *m, struct rng *rng)
{
	int i;

	for (i = 0; i < NUM_LANDMARKS; i++) {
		float obs[2];

		predict_landmark(obs, x_true, m->landmark[i]);
		z[2 * i + 0] = obs[0] + SIGMA_MEAS * rng_normal(rng);
		z[2 * i + 1] = obs[1] + SIGMA_MEAS * rng_normal(rng);
	}
}

/* ------------------------------------------------------------------ */
/* CSV output.                                                        */
/* ------------------------------------------------------------------ */

static void write_header(FILE *fp)
{
	fprintf(fp,
		"t,updated,"
		"true_x,true_y,true_yaw,"
		"est_x,est_y,est_yaw,"
		"dr_x,dr_y,dr_yaw,"
		"Pxx,Pxy,Pyy,Pyaw,"
		"v_cmd,w_cmd\n");
}

static void write_row(FILE *fp, float t, int updated, const float *x_true,
		      const struct eskf_filter *f, const float *x_dr,
		      const float *u)
{
	fprintf(fp,
		"%.4f,%d,"
		"%.9g,%.9g,%.9g,"
		"%.9g,%.9g,%.9g,"
		"%.9g,%.9g,%.9g,"
		"%.9g,%.9g,%.9g,%.9g,"
		"%.9g,%.9g\n",
		(double)t, updated,
		(double)x_true[0], (double)x_true[1], (double)x_true[2],
		(double)f->x_nom[0], (double)f->x_nom[1],
		(double)f->x_nom[2],
		(double)x_dr[0], (double)x_dr[1], (double)x_dr[2],
		(double)f->P[0 * ESKF_N + 0], (double)f->P[0 * ESKF_N + 1],
		(double)f->P[1 * ESKF_N + 1], (double)f->P[2 * ESKF_N + 2],
		(double)u[0], (double)u[1]);
}

/* ------------------------------------------------------------------ */
/* main                                                               */
/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
	const char *path = (argc > 1) ? argv[1] : "eskf_localization.csv";
	struct model model = {
		.landmark = {
			{ 14.0f, 9.0f },
			{ -11.0f, 16.0f },
			{ -15.0f, -18.0f },
		},
	};
	struct eskf_filter f;
	struct rng rng = { 0x243f6a8885a308d3ULL };
	float Q[ESKF_NW * ESKF_NW] = {
		Q_SCALE * SIGMA_V * SIGMA_V, 0.0f,
		0.0f, Q_SCALE * SIGMA_W * SIGMA_W,
	};
	float R[ESKF_M * ESKF_M];
	float P0[ESKF_N * ESKF_N] = {
		4.0f, 0.0f, 0.0f,
		0.0f, 4.0f, 0.0f,
		0.0f, 0.0f, 0.35f * 0.35f,
	};
	float x_nom0[ESKF_N_NOM] = { 1.2f, -0.8f, 0.20f };
	float x_true[ESKF_N_NOM] = { 0.0f, 0.0f, 0.0f };
	float x_dr[ESKF_N_NOM];
	float u[CONTROL_DIM], u_true[CONTROL_DIM];
	float z[ESKF_M];
	int k, i, updated;
	FILE *fp;

	fp = fopen(path, "w");
	if (fp == NULL) {
		perror(path);
		return 1;
	}

	memset(R, 0, sizeof(R));
	for (i = 0; i < ESKF_M; i++)
		R[i * ESKF_M + i] = SIGMA_MEAS * SIGMA_MEAS;

	eskf_init(&f, cb_compute_Fx, cb_compute_Fi, cb_nominal_propagate,
		  cb_compute_Hx, cb_measure_model, cb_inject, cb_compute_G);
	eskf_reset(&f, x_nom0, P0);

	/* Dead reckoning starts from the same (offset) initial estimate. */
	memcpy(x_dr, x_nom0, sizeof(x_dr));

	write_header(fp);

	control(0.0f, u);
	model.u[0] = u[0];
	model.u[1] = u[1];
	model.dt = DT;
	write_row(fp, 0.0f, 0, x_true, &f, x_dr, u);

	for (k = 1; k <= STEPS; k++) {
		float t = (float)k * DT;

		/* Control for the interval [t-dt, t], plus true execution noise. */
		control(t - 0.5f * DT, u);
		model.u[0] = u[0];
		model.u[1] = u[1];
		model.dt = DT;

		u_true[0] = u[0] + SIGMA_V * rng_normal(&rng);
		u_true[1] = u[1] + SIGMA_W * rng_normal(&rng);

		propagate(x_true, u_true, DT, x_true);
		propagate(x_dr, u, DT, x_dr);

		/*
		 * Predict at 100 Hz.  The filter reads the commanded control and
		 * dt from `params`; execution noise is modelled by Fi * Q * Fi^T.
		 */
		eskf_predict(&f, Q, &model);

		updated = 0;
		if ((k % MEAS_EVERY) == 0) {
			/* Correct at 10 Hz with all three body-frame landmarks. */
			make_measurement(z, x_true, &model, &rng);
			eskf_update(&f, z, R, &model);
			updated = 1;
		}

		write_row(fp, t, updated, x_true, &f, x_dr, u);
	}

	fclose(fp);

	{
		float ex = f.x_nom[0] - x_true[0];
		float ey = f.x_nom[1] - x_true[1];
		float eyaw = wrap_pi(f.x_nom[2] - x_true[2]);
		float edr = hypotf(x_dr[0] - x_true[0], x_dr[1] - x_true[1]);

		printf("ESKF localization: %d steps written to %s\n", STEPS + 1,
		       path);
		printf("final ESKF position error: %.3f m, yaw error: %.3f deg\n",
		       (double)hypotf(ex, ey),
		       (double)(eyaw * 180.0f / PI_F));
		printf("final dead-reckoning position error: %.3f m\n",
		       (double)edr);
	}

	return 0;
}
