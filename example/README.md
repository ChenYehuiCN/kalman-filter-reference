[简体中文](README_zh-CN.md)

[Back to project README](../README.md)

# Examples

Two complete examples use the library as-is. The C program runs the simulation
and calls the library; Python only reads the generated CSV and plots the
results. The library sources are **not modified**.

| Example | Filter | Dimensions | Demonstrates |
| --- | --- | --- | --- |
| `kf_tracking/` | linear KF | `KF_N=6`, `KF_M=2`, `KF_C=2` | 2D constant-acceleration tracking, known jerk input, time-varying GPS noise |
| `eskf_localization/` | nonlinear ESKF | `ESKF_N=3`, `ESKF_N_NOM=3`, `ESKF_M=6`, `ESKF_NW=2` | unicycle motion, three known landmarks, body-frame observations, 100 Hz predict / 10 Hz update |

```text
example/
├── Makefile
├── README.md             # this file (English)
├── README_zh-CN.md       # Chinese version
├── kf_tracking/
│   ├── main.c            # simulation + kf_init/kf_reset/kf_predict/kf_update
│   ├── Makefile          # passes KF_N/KF_M/KF_C to main.c and ../../kf/kf.c
│   ├── plot.py           # CSV -> 3 PNG figures
│   └── output/           # generated
└── eskf_localization/
    ├── main.c            # simulation + 7 ESKF callbacks + predict/update
    ├── Makefile          # passes ESKF_N/ESKF_N_NOM/ESKF_M/ESKF_NW to both .c files
    ├── plot.py           # CSV -> 3 PNG figures
    └── output/           # generated
```

## Quick start

Requirements: a C11 compiler, `make`, Python 3, `numpy` and `matplotlib`.

```bash
cd example

make            # build both examples
make run        # run both simulations and write CSV files
make plot       # run both and generate all PNG figures
make demo       # alias for make plot
make clean      # remove binaries, CSV files and PNG figures
```

Run one example on its own:

```bash
make -C kf_tracking demo
make -C eskf_localization demo

# interactive plots (requires a GUI backend)
python3 kf_tracking/plot.py kf_tracking/output/kf_tracking.csv --show
python3 eskf_localization/plot.py eskf_localization/output/eskf_localization.csv --show
```

Generated figures (in each example's `output/` directory):

| File | Content |
| --- | --- |
| `kf_tracking/output/kf_trajectory.png` | 2D path, GPS samples, 2σ position ellipses |
| `kf_tracking/output/kf_states.png` | truth, estimate and ±2σ for position/velocity/acceleration |
| `kf_tracking/output/kf_errors.png` | steady-state error (`t ≥ 2 s`) and ±2σ, including the degraded GPS window |
| `eskf_localization/output/eskf_trajectory.png` | truth, ESKF, dead reckoning, landmarks, 2σ ellipses |
| `eskf_localization/output/eskf_errors.png` | ESKF error and ±2σ; right column zooms `t ≥ 1 s` |
| `eskf_localization/output/eskf_vs_dr.png` | ESKF versus dead reckoning (log error scale) |

The default CSV path in each `plot.py` is resolved relative to the script, so
the scripts can also be invoked from the repository root:

```bash
python3 example/kf_tracking/plot.py
python3 example/eskf_localization/plot.py
```

---

## 1. Linear KF: 2D constant-acceleration tracking

### State, input and measurement

```text
x = [px, py, vx, vy, ax, ay]^T       (KF_N = 6)
u = [jx, jy]^T                       (KF_C = 2, jerk)
z = [px, py]^T                       (KF_M = 2)
```

The discrete model is a constant-acceleration (CA) model driven by a known
jerk input:

$$
F = \begin{bmatrix}
I & \Delta t I & \frac{\Delta t^2}{2} I \\
0 & I & \Delta t I \\
0 & 0 & I
\end{bmatrix},
\qquad
B = \begin{bmatrix}
\frac{\Delta t^3}{6} I \\
\frac{\Delta t^2}{2} I \\
\Delta t I
\end{bmatrix}
$$

$$
H = \begin{bmatrix} I & 0 & 0 \end{bmatrix}
$$

with `dt = 0.02 s`, `T_END = 30 s` and 1500 steps.

### Truth simulation

The ground truth follows the same discrete model plus a small unmodelled
acceleration increment:

$$
x_{k+1} = F x_k + B u_k + G w_k,
\qquad w_k \sim \mathcal{N}(0, \sigma_a^2 I)
$$

$$
G = \begin{bmatrix}
\frac{\Delta t^2}{2} I \\
\Delta t I \\
I
\end{bmatrix},
\qquad
Q = G(\sigma_a^2 I)G^{\top}
$$

with `ACCEL_NOISE = 0.0005 m/s²`. The filter's `Q` is built from exactly the
same disturbance model, so the covariance has a clear statistical meaning.

The true trajectory starts as a figure-eight:

```text
px(t) = 80 sin(0.35 t)
py(t) = 45 sin(0.70 t)
```

`main.c` differentiates the trajectory analytically to obtain the jerk `u`.
The control is sampled at the midpoint of each step before the truth state is
propagated. The figure-eight is only the initial condition/shape; the process
noise perturbs the actual path.

### Time-varying GPS noise

```text
0 ≤ t < 12 s       σ_gps = 0.75 m
12 ≤ t < 16 s      σ_gps = 3.0 m
16 ≤ t ≤ 30 s      σ_gps = 0.75 m

z = [px, py]^T + N(0, σ_gps² I)
R = σ_gps² I        (rebuilt at every step)
```

`kf_update()` receives a different `R` during the degraded interval. The
plotting script detects samples whose `sigma` is more than 1.5× the minimum and
shades them automatically.

### Initial condition

```text
x0 = truth(0) + [3.0, -2.0, 2.0, -3.0, 0.5, -0.6]
P0 = diag(9, 9, 9, 9, 2, 2)
```

### Library calls

```c
kf_init(&f);
kf_reset(&f, x0, P0);

for (k = 1; k <= STEPS; k++) {
    /* build u, F, B, H, Q, R, z for this step */
    kf_predict(&f, u, F, B, Q);
    kf_update(&f, z, H, R);
}
```

- `u/F/B/Q` are consumed by `kf_predict()`; `z/H/R` are consumed by
  `kf_update()`.
- `F`, `B`, `H` and `Q` are constant in this example; `R` changes with the GPS
  noise level.
- Prediction and update are separate calls. A step without a measurement would
  simply skip `kf_update()`.
- The library performs the Joseph-form covariance update and the Cholesky-based
  gain calculation internally.

### CSV columns

```text
t,
true_px, true_py, true_vx, true_vy, true_ax, true_ay,
meas_px, meas_py,
est_px, est_py, est_vx, est_vy, est_ax, est_ay,
Pxx, Pxy, Pyy, Pvx, Pvy, Pax, Pay,
jx, jy, sigma
```

`P*` columns are variances except `Pxy`, which is the position covariance.

### Results (fixed random seed)

```text
samples              : 1501 (including the t = 0 row)
2D position RMSE     : 0.213 m
2D velocity RMSE     : 0.442 m/s
2D acceleration RMSE : 0.336 m/s²

fraction inside 2σ: position x/y 97.3% / 91.1%
```

The error plot starts at `t = 2 s`; the initial convergence is visible in
`kf_states.png`, and the console reports the full-run RMSE. This is why the
number printed in the figure can differ from the full-run console output.

### Tuning

| Parameter | Effect |
| --- | --- |
| `TRAJ_A`, `TRAJ_B`, `TRAJ_W` | figure-eight shape and speed |
| `ACCEL_NOISE` | process-noise level (used by both truth and `Q`) |
| `measurement_sigma()` | normal/degraded GPS noise and the degraded interval |
| `x0`, `P0` | initial estimate and uncertainty |
| `T_END`, `DT` | simulation length and step size (`STEPS` is derived) |

If you change the GPS noise levels, the plot still finds the degraded samples
as long as they are at least 1.5× the minimum `sigma`.

---

## 2. Nonlinear ESKF: unicycle with body-frame landmarks

### Nominal state and motion model

```text
x  = [px, py, yaw]^T      (nominal, ESKF_N_NOM = 3)
dx = [dpx, dpy, dyaw]^T   (error,   ESKF_N = 3)
u  = [v, omega]^T         (control, user-defined)
```

The library has no compile-time control dimension: `u` is not passed through
`eskf_predict()` at all. The example stores `u` and `dt` in `struct model`,
which is passed as `params`; the callbacks read them from there. The control has
two entries here, and the process-noise dimension happens to be
`ESKF_NW = 2` as well. The nominal and error dimensions are both 3 in this
example; the library allows them to differ (for example, a 4-entry quaternion
nominal state with a 3-entry attitude error).

Continuous-time unicycle model:

$$
\dot p_x = v\cos\theta,\qquad
\dot p_y = v\sin\theta,\qquad
\dot\theta = \omega
$$

`nominal_propagate()` uses the exact constant-twist integration. For
$|\omega|\Delta t > 10^{-3}$:

$$
p_x^{+} = p_x + \frac{v}{\omega}
\left(\sin(\theta+\omega\Delta t)-\sin\theta\right)
$$

$$
p_y^{+} = p_y + \frac{v}{\omega}
\left(\cos\theta-\cos(\theta+\omega\Delta t)\right)
$$

$$
\theta^{+} = \theta + \omega\Delta t
$$

For small $|\omega|\Delta t$, it falls back to the straight-line update:

$$
p_x^{+} = p_x + v\cos\theta\,\Delta t,\qquad
p_y^{+} = p_y + v\sin\theta\,\Delta t
$$

The threshold is on `|omega| * dt`, not just `|omega|`. In single precision the
exact formulas subtract two nearby sines/cosines and lose most of their
significant bits when the angle increment is tiny. Using a too-small threshold
can make the analytic $F_x$ produce a huge spurious value and make the
covariance jump when the robot drives straight. Both `propagate()` and
`motion_jacobians()` use the same `|omega| * dt > 1e-3` branch.

### Error state and Jacobians

```text
dx = [δpx, δpy, δyaw]^T
```

The error-state transition matrix is the Jacobian of the exact propagation:

$$
F_x =
\begin{bmatrix}
1 & 0 & \partial p_x / \partial \theta \\
0 & 1 & \partial p_y / \partial \theta \\
0 & 0 & 1
\end{bmatrix},
\qquad
F_i =
\begin{bmatrix}
\partial p_x/\partial v & \partial p_x/\partial \omega \\
\partial p_y/\partial v & \partial p_y/\partial \omega \\
0 & \partial \theta/\partial \omega
\end{bmatrix}
$$

For the arc branch:

$$
\frac{\partial p_x}{\partial \theta}
= \frac{v}{\omega}\left(\cos(\theta+\omega\Delta t)-\cos\theta\right),
\qquad
\frac{\partial p_y}{\partial \theta}
= \frac{v}{\omega}\left(\sin(\theta+\omega\Delta t)-\sin\theta\right)
$$

$$
\frac{\partial p_x}{\partial v}
= \frac{\sin(\theta+\omega\Delta t)-\sin\theta}{\omega},
\qquad
\frac{\partial p_y}{\partial v}
= \frac{\cos\theta-\cos(\theta+\omega\Delta t)}{\omega}
$$

$$
\frac{\partial p_x}{\partial \omega}
= \frac{v\left(\omega\Delta t\cos(\theta+\omega\Delta t)
-(\sin(\theta+\omega\Delta t)-\sin\theta)\right)}{\omega^2}
$$

$$
\frac{\partial p_y}{\partial \omega}
= \frac{v\left(\omega\Delta t\sin(\theta+\omega\Delta t)
-(\cos\theta-\cos(\theta+\omega\Delta t))\right)}{\omega^2}
$$

$$
\frac{\partial \theta}{\partial \omega} = \Delta t
$$

The small-angle branch uses the corresponding first-order limits. The code
contains both cases; `motion_jacobians()` is the single source of truth for
`compute_Fx()` and `compute_Fi()`.

### Measurement model

Three known landmarks:

```text
L1 = (14, 9)
L2 = (-11, 16)
L3 = (-15, -18)
```

The sensor reports each landmark position in the robot body frame:

$$
h_i(x) = R(\theta)^{\top}(L_i-p)
= \begin{bmatrix}
c & s \\
-s & c
\end{bmatrix}
\begin{bmatrix}
L_{x,i}-p_x \\
L_{y,i}-p_y
\end{bmatrix},
\qquad c=\cos\theta,\ s=\sin\theta
$$

All three landmarks are stacked into a 6-dimensional measurement:

$$
z = \begin{bmatrix}
h_1 \\
h_2 \\
h_3
\end{bmatrix} + v,
\qquad v\sim\mathcal{N}(0,\sigma_m^2 I),
\qquad \sigma_m = 0.5\ \text{m}
$$

The measurement Jacobian for one landmark is

$$
H_i =
\begin{bmatrix}
-c & -s & h_y \\
s & -c & -h_x
\end{bmatrix}
$$

where $h_x,h_y$ are the components of $h_i$. Python and C use the same
landmark coordinates; if you change the map, update both `main.c` and
`plot.py`.

### Error injection and G

```text
px  += δpx
py  += δpy
yaw  = wrap_pi(yaw + δyaw)
```

The injection is additive apart from angle wrapping, so

$$
G = \frac{\partial(x\boxplus\delta x)}{\partial\delta x} = I
$$

`compute_G()` is still implemented and returns the identity; the library
applies $P\leftarrow GPG^{\top}$ after every injection.

### Noise model and multi-rate loop

Truth controls are perturbed by execution noise:

$$
v_{\text{true}} = v_{\text{cmd}} + \mathcal{N}(0,\sigma_v^2),
\qquad
\omega_{\text{true}} = \omega_{\text{cmd}} + \mathcal{N}(0,\sigma_\omega^2)
$$

with `SIGMA_V = 0.10 m/s`, `SIGMA_W = 0.04 rad/s`. The filter receives the
commanded controls; the execution noise enters the covariance through

$$
P^{-} = F_xPF_x^{\top}+F_iQF_i^{\top}.
$$

The filter uses an inflated noise covariance to absorb linearization error:

$$
Q = Q_{\text{SCALE}}\cdot\mathrm{diag}(\sigma_v^2,\sigma_\omega^2),
\qquad Q_{\text{SCALE}}=2.
$$

Schedule:

```text
dt = 0.01 s              -> 100 Hz predict
MEAS_EVERY = 10          -> 10 Hz measurement update
T_END = 60 s             -> 6000 steps, 600 updates
```

The loop calls `eskf_predict()` every step and `eskf_update()` only when
`k % MEAS_EVERY == 0`. Both filters expose prediction and update separately;
the ESKF example uses that to run them at different rates.

### Control profile and initial condition

```text
v(t)     = 2.0 + 0.6 sin(0.17 t)                  [m/s]
omega(t) = 0.55 sin(0.23 t) + 0.20 sin(0.61 t)    [rad/s]

initial estimate: x_nom0 = (1.2, -0.8, 0.20)
initial covariance: P0 = diag(4, 4, 0.35²)
```

The true state starts at `(0,0,0)`. Dead reckoning starts from the same offset
initial estimate and integrates only the commanded controls, without any
measurement updates.

### Library calls and callbacks

```c
eskf_init(&f, compute_Fx, compute_Fi, nominal_propagate,
          compute_Hx, measure_model, inject, compute_G);
eskf_reset(&f, x_nom0, P0);

for (k = 1; k <= STEPS; k++) {
    /* u, dt, Q, params */
    model.u[0] = u[0];
    model.u[1] = u[1];
    model.dt = DT;
    eskf_predict(&f, Q, &model);
    if (k % MEAS_EVERY == 0) {
        /* build z, R */
        eskf_update(&f, z, R, &model);
    }
}
```

All seven callbacks are implemented. `params` points to `struct model`, which
holds the three landmarks, the current control `u` and `dt`. The library
sources are untouched. Note the callback order in `eskf_init()`: `compute_Hx`
comes before `measure_model`.

### CSV columns

```text
t, updated,
true_x, true_y, true_yaw,
est_x, est_y, est_yaw,
dr_x, dr_y, dr_yaw,
Pxx, Pxy, Pyy, Pyaw,
v_cmd, w_cmd
```

The `yaw` columns are wrapped to $[-\pi,\pi]$. The Python script computes yaw
errors with a wrapped difference. `Pxx/Pyy/Pyaw` are variances; `Pxy` is the
position covariance.

### Results (fixed random seed)

```text
ESKF position RMSE : 0.078 m     final error: 0.015 m
ESKF yaw RMSE      : 0.518°      final error: 0.390°
dead-reckoning pos : 3.120 m     final error: 1.738 m
dead-reckoning yaw : 11.754°     final error: 13.011°

position accuracy improvement over DR: about 40×
fraction inside 2σ: x/y/yaw 97.5% / 93.7% / 99.4%
```

`eskf_errors.png` shows the full run in the left column and a `t ≥ 1 s` zoom in
the right column. `eskf_vs_dr.png` uses a log scale so both the initial
transient and the steady-state gap between ESKF and dead reckoning are visible.

### Tuning

| Parameter | Effect |
| --- | --- |
| `SIGMA_V`, `SIGMA_W` | true control noise (also sets the base filter `Q`) |
| `Q_SCALE` | inflation of the filter process noise |
| `SIGMA_MEAS`, `MEAS_EVERY` | measurement noise and update rate |
| `control()` | robot trajectory |
| `x_nom0`, `P0` | initial estimate and uncertainty |
| `model.landmark` | map; update `plot.py`'s `LANDMARKS` too |
| `ESKF_M` in the Makefile | two times the number of landmarks |
| `T_END`, `DT` | simulation length and step size (`STEPS` is derived) |

If you change the number of landmarks, update the Makefile's `-DESKF_M`
(`2 × landmark count`) and the `model.landmark` initializer. `NUM_LANDMARKS` in
`main.c` is derived from `ESKF_M`; `plot.py` still needs its own landmark list.

---

## Compile-time dimensions

The library structs contain fixed-size arrays, so these macros must be passed
to **both** the library `.c` file and the user `.c` file:

| Example | Makefile flags |
| --- | --- |
| KF | `-I../../kf -DKF_N=6 -DKF_M=2 -DKF_C=2` |
| ESKF | `-I../../eskf -DESKF_N=3 -DESKF_N_NOM=3 -DESKF_M=6 -DESKF_NW=2` |

Each `main.c` also has a compile-time check:

```c
#if KF_N != 6 || KF_M != 2 || KF_C != 2
#error "Compile this example with KF_N=6, KF_M=2 and KF_C=2"
#endif

#if ESKF_N != 3 || ESKF_N_NOM != 3 || ESKF_M != 6 || ESKF_NW != 2
#error "Compile this example with ESKF_N=3, ESKF_N_NOM=3, ESKF_M=6 and ESKF_NW=2"
#endif
```

If you compile without the `-D` flags, the headers fall back to dimension 1
and the check stops the build with a clear error instead of producing a
mismatched struct layout. `ESKF_N_NOM` is the nominal-state size; `ESKF_N` is
the error-state size. They are equal in this example but may differ in models
such as quaternion attitude estimation.

## Notes

- Both examples use a deterministic xorshift64* + Box-Muller RNG. Change the
  seed in `struct rng rng = {...}` to get different noise realizations.
- The C programs accept an optional output path argument. Running `./kf_demo`
  or `./eskf_demo` writes the CSV into the current directory.
- The plotting scripts default to the `Agg` backend and only save PNGs. Pass
  `--show` to open them with the available GUI backend.
- The generated CSV files and figures are ignored by Git; `make demo`
  regenerates them.
- The examples compile cleanly with `-std=c11 -Wall -Wextra -Wpedantic` and
  have been run with AddressSanitizer/UndefinedBehaviorSanitizer.
