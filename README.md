[简体中文](README_zh-CN.md)

# Kalman Filter Reference (KF / ESKF)

A compact, fixed-size C implementation of the linear Kalman filter (KF) and the
error-state Kalman filter (ESKF).

This is **reference code**, not a turnkey library. The priority is to make the
math, data layout and call flow easy to read and adapt. The implementation is
deliberately general: it supports control inputs, time-varying matrices,
nonlinear ESKF callbacks and a separate process-noise Jacobian, even when a
particular application needs none of them. Users are expected to copy the
relevant parts and trim or specialize them for their own model.

The only external dependency is `libm` (`sqrtf`, `sinf`, `cosf`, ...). There is
no heap allocation and no hidden global state.

## Contents

- [At a glance](#at-a-glance)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
- [Dimensions and data layout](#dimensions-and-data-layout)
- [Linear Kalman filter (KF)](#linear-kalman-filter-kf)
- [Error-state Kalman filter (ESKF)](#error-state-kalman-filter-eskf)
- [Trimming and optimizing](#trimming-and-optimizing)
- [Examples](#examples)
- [Validation](#validation)
- [License](#license)

## At a glance

| Item | KF (`kf/`) | ESKF (`eskf/`) |
| --- | --- | --- |
| Model | linear state space | nonlinear nominal model + error state |
| State dimension | `KF_N` | `ESKF_N_NOM` (nominal), `ESKF_N` (error) |
| Measurement dimension | `KF_M` | `ESKF_M` |
| Input/noise dimension | `KF_C` (control) | `ESKF_NW` (noise) |
| Main call | `kf_predict()`, `kf_update()` | `eskf_predict()`, `eskf_update()` |
| Covariance update | Joseph form | Joseph form |
| Linear solve | Cholesky + forward/back substitution | same |
| Storage | fixed-size arrays in the struct | fixed-size arrays in the struct |
| Data type | `float` | `float` |

The headers default every dimension to `1`. This only exists so the headers can
be included in a scratch file; real use must define the dimensions and compile
the library `.c` file with the **same** macros (see
[Dimensions and data layout](#dimensions-and-data-layout)).

## Repository layout

```text
.
├── kf/
│   ├── kf.h              # KF API and struct
│   └── kf.c              # KF implementation
├── eskf/
│   ├── eskf.h            # ESKF API, callbacks and struct
│   └── eskf.c            # ESKF implementation
├── example/
│   ├── README.md         # example guide (English)
│   ├── README_zh-CN.md   # example guide (Chinese)
│   ├── kf_tracking/      # linear KF example
│   └── eskf_localization/# nonlinear ESKF example
├── LICENSE               # MIT
├── README.md
└── README_zh-CN.md
```

## Quick start

The repository has no top-level build system because the library is only two
source files. The examples show the intended build flags and a complete usage
pattern:

```bash
cd example
make demo          # build, run both simulations, generate all figures
```

Requirements: a C11 compiler, `make`, Python 3, `numpy` and `matplotlib`.
See [`example/README.md`](example/README.md) for model details, output columns
and tuning notes. A Chinese version is available at
[`example/README_zh-CN.md`](example/README_zh-CN.md).

### Minimal KF build

```bash
cc -I kf -DKF_N=6 -DKF_M=2 -DKF_C=2 \
   user.c kf/kf.c -lm -o user_kf
```

### Minimal ESKF build

```bash
cc -I eskf -DESKF_N=3 -DESKF_N_NOM=3 -DESKF_M=6 -DESKF_NW=2 \
   user.c eskf/eskf.c -lm -o user_eskf
```

The `-D` flags are not optional in practice: the struct sizes and the array
dimensions in the library `.c` file are compiled with these macros.

## Dimensions and data layout

### Compile-time configuration

| Macro | Meaning | Used by |
| --- | --- | --- |
| `KF_N` | number of states | `kf.h`, `kf.c` |
| `KF_M` | number of measurements | `kf.h`, `kf.c` |
| `KF_C` | number of control inputs | `kf.h`, `kf.c` |
| `ESKF_N` | error-state dimension (`dx`, `P`, `Fx`, `G`) | `eskf.h`, `eskf.c` |
| `ESKF_N_NOM` | nominal-state dimension (`x_nom`) | `eskf.h`, `eskf.c` |
| `ESKF_M` | measurement dimension | `eskf.h`, `eskf.c` |
| `ESKF_NW` | process-noise dimension (`Fi` columns, `Q` size) | `eskf.h`, `eskf.c` |

`ESKF_N_NOM` and `ESKF_N` are independent. They are often equal, but they may
differ (for example, a 4-parameter quaternion nominal state with a 3-dimensional
attitude error). `x_nom` uses `ESKF_N_NOM`; `dx`, `P`, `Fx`, `Hx` columns and
`G` use `ESKF_N`. Defining `ESKF_N` does not set `ESKF_N_NOM`; both default to
`1` independently, so define both for an ESKF build.

The library and the user translation unit must see identical values. Mixing
any dimension macro between `main.c` and the library `.c` file is undefined
behavior because the compiler lays out the filter struct differently in each
translation unit.

### Matrix storage

All matrices are stored **row-major**, with no leading dimension argument:

```c
A[i * cols + j]   // row i, column j
```

Vectors are plain arrays. For example, the KF state is `f.x[KF_N]`, the
measurement is `z[KF_M]`, and the covariance is `f.P[KF_N * KF_N]` with
`P[i * KF_N + j]`.

All matrices are `float`. If you need `double`, change the type consistently in
the structs, APIs and helper functions (there is no typedef to switch).

### Ownership and reentrancy

All working memory is inside the filter struct. There is no `malloc`, no static
mutable data in the library, and no hidden state beyond the struct. Two filters
can run independently in two structs. User callbacks may of course access their
own global state through `params`.

Fields whose names start with `_` are internal workspace (`_tmp_*`, `_S`, `_K`,
`_y`, ...). They are visible only because C has no private members; do not treat
them as part of the stable API.

## Linear Kalman filter (KF)

### Model

The library implements the standard linear Gaussian model

$$
x_k = F_{k-1} x_{k-1} + B_{k-1} u_{k-1} + w_{k-1},
\qquad w_{k-1} \sim \mathcal{N}(0, Q_{k-1})
$$

$$
z_k = H_k x_k + v_k,
\qquad v_k \sim \mathcal{N}(0, R_k)
$$

The API separates the time update from the measurement update:

- `kf_predict(f, u, F, B, Q)` propagates the state and covariance.
- `kf_update(f, z, H, R)` applies one measurement.

The matrices are arguments, not stored fields, so `F`, `B`, `H`, `Q` and `R`
may change from step to step. Calling `kf_update()` is optional when no
measurement is available.

### Predict

$$
\hat{x}^{-} = F\hat{x} + Bu
$$

$$
P^{-} = FPF^{\top} + Q
$$

`kf_predict()` performs these two operations.

### Update

$$
S = HP^{-}H^{\top} + R
$$

$$
K = P^{-}H^{\top}S^{-1}
$$

$$
\tilde{y} = z - H\hat{x}^{-}
$$

$$
\hat{x} = \hat{x}^{-} + K\tilde{y}
$$

$$
P = (I-KH)P^{-}(I-KH)^{\top} + KRK^{\top}
$$

The last equation is the **Joseph form**. It costs more than the simpler
$(I-KH)P^{-}$ but keeps the covariance symmetric and positive semidefinite in
the presence of rounding errors. The implementation does not explicitly
symmetrize `P`. `kf_update()` performs the update equations.

### Code map: predict

| Math | Code in `kf_predict()` |
| --- | --- |
| $\hat{x}^{-}=F\hat{x}+Bu$ | `_tmp_n1 = F x`; `_tmp_n1_2 = B u`; add; copy to `x` |
| $P^{-}=FPF^{\top}+Q$ | `_tmp_nn = F P`; `P = _tmp_nn F^T`; `P += Q` |

### Code map: update

| Math | Code in `kf_update()` |
| --- | --- |
| $S=HP^{-}H^{\top}+R$ | `_tmp_mn = H P`; `_S = _tmp_mn H^T`; `_S += R` |
| $S=LL^{\top}$ | `mat_chol(_S, KF_M)` overwrites `_S` with `L` |
| $K=P^{-}H^{\top}S^{-1}$ | `mat_chol_solve(_S, _tmp_mn, M, N)` solves $SX=HP^{-}$; `_K = X^{\top}` |
| $\tilde{y}=z-H\hat{x}^{-}$ | `_y = H x`; `_y = z - _y` |
| $\hat{x}=\hat{x}^{-}+K\tilde{y}$ | `_tmp_n1 = K _y`; `x += _tmp_n1` |
| Joseph form | `_tmp_nn = K H`; `I - K H`; two `mat_mul_bt` calls; add $KRK^{\top}$ |

The gain is computed by solving $SX=HP^{-}$ with a Cholesky factorization of
$S$, then transposing $X$. This avoids forming $S^{-1}$ explicitly. The Cholesky
routine assumes $S$ is symmetric positive definite; it does not pivot.

### API

```c
struct kf_filter f;

kf_init(&f);                 // once: initialize the identity workspace
kf_reset(&f, x0, P0);        // copy initial state and covariance

for (;;) {
    kf_predict(&f, u, F, B, Q);        // time update
    if (measurement_available)
        kf_update(&f, z, H, R);        // measurement update
}
```

| Function | Behavior |
| --- | --- |
| `void kf_init(struct kf_filter *f)` | Zeroes `_I_nn` and sets its diagonal to one. Does **not** initialize `x` or `P`. |
| `void kf_reset(struct kf_filter *f, float *x, float *P)` | Copies `x` and `P` into the struct. Does not modify the workspace matrices. |
| `void kf_predict(struct kf_filter *f, float *u, float *F, float *B, float *Q)` | Predicts `x` and `P` with the control, transition and process-noise matrices. |
| `void kf_update(struct kf_filter *f, float *z, float *H, float *R)` | Applies one measurement update to `x` and `P`. |

Argument shapes:

| Argument | Shape | Used by |
| --- | --- | --- |
| `u` | `KF_C` | `kf_predict()` |
| `F` | `KF_N × KF_N` | `kf_predict()` |
| `B` | `KF_N × KF_C` | `kf_predict()` |
| `Q` | `KF_N × KF_N` | `kf_predict()` |
| `z` | `KF_M` | `kf_update()` |
| `H` | `KF_M × KF_N` | `kf_update()` |
| `R` | `KF_M × KF_M` | `kf_update()` |

All pointers passed to a function must be valid. If the application has no
control input, pass a zero `B` and a zero `u`; do not pass `NULL`, because
`kf_predict()` always evaluates $Bu$.

### KF limitations

- `kf_update()` has no missing-measurement, gating or outlier-rejection logic.
  Skip the call or handle those cases outside the filter.
- `kf_update()` may be called multiple times (for example, one sensor at a
  time). This performs a sequence of updates; if the measurements are
  correlated, stack them into one `z/H/R` instead.
- The dimensions are compile-time constants; runtime-varying state sizes are
  not supported.
- `S` must be positive definite for the Cholesky solve. A positive-definite
  `R` and a positive-semidefinite `P` normally guarantee this, but the library
  does not check it.

## Error-state Kalman filter (ESKF)

### Idea

The ESKF keeps two things:

- a **nominal state** $x \in \mathbb{R}^{N_{\mathrm{nom}}}$, which is
  propagated with the full nonlinear model;
- an **error state** $\delta x \in \mathbb{R}^{N}$ and its covariance $P$,
  which are assumed small and propagated with a linearized model.

Between calls the error state is zero. The Kalman correction is computed in
the error coordinates, injected into the nominal state, and the error state is
reset to zero.

The nominal state has dimension `ESKF_N_NOM`; the error state has dimension
`ESKF_N`. They are independent compile-time constants and may differ. A common
example is a 4-parameter quaternion nominal state with a 3-dimensional attitude
error: `x_nom` has 4 entries, while `dx`, `P`, `Fx`, `G` are 3-dimensional and
`Hx` has 3 columns. The user-provided `inject()` callback maps the nominal state
and error state back to the nominal dimension. If the dimensions are equal and
the injection is additive, this reduces to the usual ESKF form.

### Model

Let the nominal state be $x \in \mathbb{R}^{N_{\mathrm{nom}}}$ and the error
state be $\delta x \in \mathbb{R}^{N}$.

Nominal propagation:

$$
x_k = f(x_{k-1}, u_{k-1}, \Delta t)
$$

Error-state propagation:

$$
\delta x_k = F_x\,\delta x_{k-1} + F_i\,w_{k-1},
\qquad w_{k-1} \sim \mathcal{N}(0, Q_{k-1})
$$

where

$$
F_x \in \mathbb{R}^{N \times N},
\qquad
F_i \in \mathbb{R}^{N \times N_W}.
$$

When $N_{\mathrm{nom}} = N$, $F_x$ is the usual state Jacobian
$\partial f/\partial x$. When the dimensions differ, $F_x$ is the Jacobian of
the chosen local error-state propagation; the user supplies it through
`compute_Fx()`.

$F_i$ maps the process-noise vector (`ESKF_NW` entries) into the error state.
If noise is additive in the error state, set $F_i = I$ and
`ESKF_NW = ESKF_N`; then $Q$ is the usual state process-noise covariance.

The control $u$ and time step $\Delta t$ above are part of the user's model.
The library does not pass them to the callbacks; the callbacks read them from
`params` or from their own context. `eskf_predict()` only receives `Q` and
`params`.

Measurement model:

$$
z = h(x) + v,
\qquad v \sim \mathcal{N}(0,R)
$$

### Predict

$$
x^{-} = f(x,u,\Delta t)
$$

$$
P^{-} = F_x P F_x^{\top} + F_i Q F_i^{\top}
$$

The public entry point is `eskf_predict(f, Q, params)`. The implementation
calls `compute_Fx()` and `compute_Fi()` **before** `nominal_propagate()`, so
both Jacobians see the pre-propagation nominal state. It then propagates the
nominal state and `P`.

### Update

The measurement Jacobian is evaluated at $\delta x = 0$:

$$
H = \left.\frac{\partial h}{\partial \delta x}\right|_{\delta x=0}
$$

$$
S = HP^{-}H^{\top} + R
$$

$$
K = P^{-}H^{\top}S^{-1}
$$

$$
\delta\hat{x} = K\left(z - h(x^{-})\right)
$$

$$
P = (I-KH)P^{-}(I-KH)^{\top} + KRK^{\top}
$$

Here $\boxplus$ denotes the model-specific injection operation. Inject the
correction and reset the error state:

$$
x \leftarrow x^{-} \boxplus \delta\hat{x}
$$

$$
G = \frac{\partial(x \boxplus \delta x)}{\partial \delta x},
\qquad
P \leftarrow GPG^{\top}
$$

$$
\delta x \leftarrow 0
$$

For a purely additive injection $x \boxplus \delta x = x + \delta x$, the
injection Jacobian is $G=I$.

### Code map

| Math | Code in `eskf_predict()` / `eskf_update()` |
| --- | --- |
| $F_x$, $F_i$ | `compute_Fx()`, `compute_Fi()` before propagation; `Fx` is `ESKF_N × ESKF_N`, `Fi` is `ESKF_N × ESKF_NW` |
| $x^{-}=f(x,u,\Delta t)$ | `nominal_propagate()`; `x_nom` is `ESKF_N_NOM` long |
| $P^{-}=F_xPF_x^{\top}+F_iQF_i^{\top}$ | `_tmp_nn = Fx P`; `P = _tmp_nn Fx^T`; `_tmp_nw = Fi Q`; `_tmp_nn = _tmp_nw Fi^T`; add |
| $H$ | `compute_Hx()` at the pre-injection nominal state; `Hx` is `ESKF_M × ESKF_N` |
| $S=HP^{-}H^{\top}+R$ | `_tmp_mn = Hx P`; `_S = _tmp_mn Hx^T`; `_S += R` |
| $S=LL^{\top}$ | `mat_chol(_S, ESKF_M)` |
| $K=P^{-}H^{\top}S^{-1}$ | `mat_chol_solve()` then transpose into `_K` |
| $h(x^{-})$ | `measure_model()` before injection |
| $\delta\hat{x}=K(z-h)$ | `_y = z - _y`; `dx = K _y` |
| Joseph form | same sequence as the KF |
| $x\leftarrow x\boxplus\delta x$ | `inject(x_nom, dx, params)` |
| $P\leftarrow GPG^{\top}$ | `compute_G(G, dx, params)`; two `mat_mul_bt` calls |
| reset | `memset(dx, 0, ...)` |

`Fx`, `dx`, `P` and `G` live in the error-state dimension `ESKF_N`. `x_nom`
lives in the nominal dimension `ESKF_N_NOM`. `Hx` is `ESKF_M × ESKF_N`.

### API and call order

```c
struct eskf_filter f;

eskf_init(&f, compute_Fx, compute_Fi, nominal_propagate,
          compute_Hx, measure_model, inject, compute_G);
eskf_reset(&f, x_nom0, P0);

for (;;) {
    /* Put the current control and dt where the callbacks can read them. */
    eskf_predict(&f, Q, params);
    if (measurement_available)
        eskf_update(&f, z, R, params);
}
```

The callback order in `eskf_init()` is `compute_Fx`, `compute_Fi`,
`nominal_propagate`, `compute_Hx`, `measure_model`, `inject`, `compute_G`.

| Function | Behavior |
| --- | --- |
| `eskf_init(...)` | Stores the seven callbacks, zeros `dx` and fills the identity workspace. Does **not** initialize `x_nom` or `P`. |
| `eskf_reset(&f, x_nom, P)` | Copies the nominal state (`ESKF_N_NOM` entries) and covariance (`ESKF_N × ESKF_N`), zeros `dx`. |
| `eskf_predict(&f, Q, params)` | Computes `Fx`/`Fi`, zeros `dx`, propagates the nominal state, and propagates `P`. |
| `eskf_update(&f, z, R, params)` | Computes `Hx`, the gain, the innovation and `dx`; applies the Joseph covariance update, injects `dx`, applies `G`, and zeros `dx`. |

| Argument | Shape | Notes |
| --- | --- | --- |
| `Q` | `ESKF_NW × ESKF_NW` | process-noise covariance in the noise space |
| `z` | `ESKF_M` | measurement |
| `R` | `ESKF_M × ESKF_M` | measurement-noise covariance |
| `params` | `void *` | forwarded to all callbacks |

The controls and the time step are **not** library arguments. If
`compute_Fx()`, `compute_Fi()` or `nominal_propagate()` need them, store them
in the `params` object (or any user context) before calling
`eskf_predict()`. The same `params` pointer is passed to `eskf_update()`.

Unlike the KF, the ESKF has no compile-time control dimension. `ESKF_NW` sizes
the noise vector and `Q`, but not the user's control vector; the callbacks
interpret their own control data. The example happens to use two controls and
two noise inputs, but that is not a library constraint.

### Callback contract

```c
void compute_Fx(float *Fx, float *x_nom, void *params);
void compute_Fi(float *Fi, float *x_nom, void *params);
void nominal_propagate(float *x_nom, void *params);
void compute_Hx(float *Hx, float *x_nom, void *params);
void measure_model(float *z_pred, float *x_nom, void *params);
void inject(float *x_nom, float *dx, void *params);
void compute_G(float *G, float *dx, void *params);
```

| Callback | When it is called | Output shape | Notes |
| --- | --- | --- | --- |
| `compute_Fx` | in `eskf_predict()`, before nominal propagation | `ESKF_N × ESKF_N` | error-state transition Jacobian at the old nominal state |
| `compute_Fi` | in `eskf_predict()`, before nominal propagation | `ESKF_N × ESKF_NW` | noise Jacobian, not necessarily the control matrix |
| `nominal_propagate` | in `eskf_predict()`, after the Jacobians | modifies `x_nom` (`ESKF_N_NOM`) in place | full nonlinear propagation |
| `compute_Hx` | at the start of `eskf_update()` | `ESKF_M × ESKF_N` | measurement Jacobian at the pre-injection state |
| `measure_model` | in `eskf_update()`, before injection | `ESKF_M` | predicted measurement $h(x_{\text{nom}})$ |
| `inject` | in `eskf_update()`, after the covariance update | modifies `x_nom` (`ESKF_N_NOM`) in place | applies $x \leftarrow x\boxplus\delta x$; receives `dx` and `params` |
| `compute_G` | after `inject()`, before resetting `dx` | `ESKF_N × ESKF_N` | injection Jacobian; receives the applied `dx` |

All seven callbacks are required. The library calls them through function
pointers without checking for `NULL`. If a callback is trivial (for example
$G=I$), still provide a function that fills the output array.

`params` is forwarded to all seven callbacks, including `inject`. The library
does not interpret it. Put model data, control inputs and the current time step
there when the callbacks need them.

### Measurement residual and angle wrapping

The update uses the raw residual

$$
\tilde{y} = z - h(x_{\text{nom}})
$$

The library does **not** wrap angles, normalize quaternions or otherwise
transform the residual. If a measurement is an angle, make sure `measure_model`
returns a prediction on the same branch as `z`, or preprocess `z`/the
measurement model accordingly. The example avoids this by using Cartesian
body-frame landmark coordinates rather than range/bearing angles.

### ESKF limitations

- `ESKF_N_NOM` and `ESKF_N` are independent compile-time dimensions. The user
  must keep `inject` (maps `ESKF_N_NOM` + `ESKF_N` to `ESKF_N_NOM`) and `G`
  (`ESKF_N × ESKF_N`) consistent with them.
- Controls and the time step are not passed to the callbacks. The callbacks
  must obtain them from `params` or their own context, and must use the same
  values in `compute_Fx`, `compute_Fi` and `nominal_propagate`.
- The residual is not wrapped or gated; handle angle wrap, outliers and
  missing measurements outside the core update.
- `compute_G` receives only `dx` and `params`, not `x_nom`. If your injection
  Jacobian depends on the post-injection nominal state, stash that state in
  `params` during `inject`, or specialize the library.
- The covariance update is Joseph form; `P` is not explicitly symmetrized.
- `S` must be positive definite for Cholesky; no pivoting or regularization is
  performed.
- One update per call. Multiple sensors can be fused by stacking measurements
  or by calling `eskf_update()` multiple times, but independent noise
  assumptions must hold.
- All dimensions are compile-time constants and all data are `float`.

## Trimming and optimizing

The code intentionally keeps a general matrix path. A real product can usually
remove a large fraction of it. Typical changes:

- **No control input:** remove `B`/`u`, the control term in the KF predict, and
  the `_tmp_n1_2` workspace if you specialize the function.
- **No process-noise Jacobian:** if noise is additive in the error state, set
  `Fi = I`, `ESKF_NW = ESKF_N` and remove the `Fi Q Fi^T` multiplication.
- **Equal nominal/error dimensions:** if `ESKF_N_NOM == ESKF_N` and the
  injection is additive, the nominal state can use the same layout as the error
  state and `G = I`.
- **No injection Jacobian:** return `G = I` and remove the two `mat_mul_bt`
  calls after injection.
- **Constant matrices:** build `F`, `B`, `H`, `Q` and `R` once outside the
  loop instead of passing them every step.
- **Sparse structure:** if `H` selects states directly, replace the generic
  matrix multiplies with indexing. If `Q` or `R` is diagonal, exploit it.
- **Small fixed dimensions:** unroll the loops, store only the upper triangle
  of `P`, and avoid generic `mat_mul` calls.
- **Avoid forming `K`:** solve directly for the correction, or fuse the
  Cholesky solve with the state update.
- **Double precision or fixed point:** replace `float` consistently, not just
  in one file.
- **Multiple sensors:** stack measurements or write specialized update
  functions; the generic `H`/`R`/`z` path handles any `ESKF_M`/`KF_M` but may
  do unnecessary work.
- **Check the assumptions:** specialize only after deciding how you will handle
  positive definiteness, angle wrap, missing data and outliers.

The example directories show two concrete specialization patterns but keep
the library calls visible, so they are a good starting point for a trimmed
version.

## Examples

| Example | Filter | Dimensions | Highlights |
| --- | --- | --- | --- |
| [`example/kf_tracking`](example/kf_tracking) | linear KF | `KF_N=6`, `KF_M=2`, `KF_C=2` | 2D constant-acceleration tracking, known jerk input, time-varying GPS noise, Joseph-form covariance |
| [`example/eskf_localization`](example/eskf_localization) | nonlinear ESKF | `ESKF_N=3`, `ESKF_N_NOM=3`, `ESKF_M=6`, `ESKF_NW=2` | unicycle motion, three known landmarks, body-frame measurements, 100 Hz predict / 10 Hz update, dead-reckoning comparison |

```bash
cd example
make demo
```

The examples write CSV files and Python plots under `example/*/output/`.
See [`example/README.md`](example/README.md) for model equations, output
columns, numerical notes and tuning parameters.

## Validation

The examples compile with:

```bash
cc -std=c11 -Wall -Wextra -Wpedantic ...
```

and have also been run with AddressSanitizer and UndefinedBehaviorSanitizer.
The example output is deterministic for a fixed random seed. The project does
not ship a unit-test framework; the two examples are the integration tests.

## License

MIT. See [LICENSE](LICENSE).

Copyright (c) 2026 Chen Yehui.
