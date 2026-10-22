[English](README.md)

[返回项目 README](../README_zh-CN.md)

# 示例程序

这里有两个把库当成黑盒使用的完整示例。C 程序负责仿真和调用库，Python 只
读取生成的 CSV 并画图。库源码**没有做任何修改**。

| 示例 | 滤波器 | 维数 | 演示内容 |
| --- | --- | --- | --- |
| `kf_tracking/` | 线性 KF | `KF_N=6`、`KF_M=2`、`KF_C=2` | 二维常加速度跟踪、已知 jerk 输入、时变 GPS 噪声 |
| `eskf_localization/` | 非线性 ESKF | `ESKF_N=3`、`ESKF_N_NOM=3`、`ESKF_M=6`、`ESKF_NW=2` | 单轨运动、三个已知路标、本体坐标系观测、100 Hz 预测 / 10 Hz 更新 |

```text
example/
├── Makefile
├── README.md             # 英文版
├── README_zh-CN.md       # 本文件
├── kf_tracking/
│   ├── main.c            # 仿真 + kf_init/kf_reset/kf_predict/kf_update
│   ├── Makefile          # 把 KF_N/KF_M/KF_C 同时传给 main.c 和 ../../kf/kf.c
│   ├── plot.py           # CSV -> 3 张 PNG
│   └── output/           # 运行后生成
└── eskf_localization/
    ├── main.c            # 仿真 + 7 个 ESKF 回调 + predict/update
    ├── Makefile          # 把 ESKF_N/ESKF_N_NOM/ESKF_M/ESKF_NW 同时传给两个 .c 文件
    ├── plot.py           # CSV -> 3 张 PNG
    └── output/           # 运行后生成
```

## 快速开始

依赖：C11 编译器、`make`、Python 3、`numpy`、`matplotlib`。

```bash
cd example

make            # 只编译两个示例
make run        # 运行两个仿真并写出 CSV
make plot       # 运行两个仿真并生成全部 PNG
make demo       # make plot 的别名
make clean      # 删除可执行文件、CSV 和 PNG
```

单独运行某个示例：

```bash
make -C kf_tracking demo
make -C eskf_localization demo

# 交互式查看（需要图形界面）
python3 kf_tracking/plot.py kf_tracking/output/kf_tracking.csv --show
python3 eskf_localization/plot.py eskf_localization/output/eskf_localization.csv --show
```

生成的图（在各自 `output/` 目录下）：

| 文件 | 内容 |
| --- | --- |
| `kf_tracking/output/kf_trajectory.png` | 二维轨迹、GPS 采样点、2σ 位置椭圆 |
| `kf_tracking/output/kf_states.png` | 位置/速度/加速度的真值、估计值和 ±2σ |
| `kf_tracking/output/kf_errors.png` | `t ≥ 2 s` 的稳态误差和 ±2σ，包含 GPS 变差区间 |
| `eskf_localization/output/eskf_trajectory.png` | 真值、ESKF、航位推算、路标、2σ 椭圆 |
| `eskf_localization/output/eskf_errors.png` | ESKF 误差和 ±2σ；右列放大 `t ≥ 1 s` |
| `eskf_localization/output/eskf_vs_dr.png` | ESKF 与航位推算的误差对比（对数坐标） |

每个 `plot.py` 的默认 CSV 路径是相对脚本目录解析的，所以也可以从仓库根目录
直接运行：

```bash
python3 example/kf_tracking/plot.py
python3 example/eskf_localization/plot.py
```

---

## 1. 线性 KF：二维常加速度目标跟踪

### 状态、输入与量测

```text
x = [px, py, vx, vy, ax, ay]^T       (KF_N = 6)
u = [jx, jy]^T                       (KF_C = 2, jerk)
z = [px, py]^T                       (KF_M = 2)
```

离散模型是常加速度（CA）模型，由已知 jerk 输入驱动：

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

其中 `dt = 0.02 s`、`T_END = 30 s`、1500 步。

### 真值仿真

真值使用相同的离散模型，外加一个小的未建模加速度增量：

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

其中 `ACCEL_NOISE = 0.0005 m/s²`。滤波器的 `Q` 使用完全相同的扰动模型，
因此协方差有明确的统计含义。

真值初始轨迹是“8”字形：

```text
px(t) = 80 sin(0.35 t)
py(t) = 45 sin(0.70 t)
```

`main.c` 对该轨迹解析求导得到 jerk `u`。控制量在每一步的中点采样，然后再
传播真值。8 字形只是初始形状；过程噪声会扰动实际路径。

### 时变 GPS 噪声

```text
0 ≤ t < 12 s       σ_gps = 0.75 m
12 ≤ t < 16 s      σ_gps = 3.0 m
16 ≤ t ≤ 30 s      σ_gps = 0.75 m

z = [px, py]^T + N(0, σ_gps² I)
R = σ_gps² I        （每一步重新构造）
```

在 GPS 变差区间，`kf_update()` 收到不同的 `R`。绘图脚本会把 `sigma` 大于
最小值 1.5 倍的样本自动识别为“坏 GPS”并加阴影。

### 初始条件

```text
x0 = truth(0) + [3.0, -2.0, 2.0, -3.0, 0.5, -0.6]
P0 = diag(9, 9, 9, 9, 2, 2)
```

### 库调用

```c
kf_init(&f);
kf_reset(&f, x0, P0);

for (k = 1; k <= STEPS; k++) {
    /* 构造这一步的 u, F, B, H, Q, R, z */
    kf_predict(&f, u, F, B, Q);
    kf_update(&f, z, H, R);
}
```

- `u/F/B/Q` 传给 `kf_predict()`，`z/H/R` 传给 `kf_update()`。
- 本示例中 `F`、`B`、`H`、`Q` 是常量；`R` 随 GPS 噪声等级变化。
- 预测和更新是两个独立调用；没有量测的步可以直接跳过 `kf_update()`。
- Joseph 形式协方差更新和基于 Cholesky 的增益计算都在库内部完成。

### CSV 列

```text
t,
true_px, true_py, true_vx, true_vy, true_ax, true_ay,
meas_px, meas_py,
est_px, est_py, est_vx, est_vy, est_ax, est_ay,
Pxx, Pxy, Pyy, Pvx, Pvy, Pax, Pay,
jx, jy, sigma
```

`P*` 列除了 `Pxy`（位置协方差）之外都是方差。

### 结果（固定随机种子）

```text
采样点            : 1501（含 t = 0 初始行）
二维位置 RMSE     : 0.213 m
二维速度 RMSE     : 0.442 m/s
二维加速度 RMSE   : 0.336 m/s²

位置 x/y 落在 2σ 内的比例: 97.3% / 91.1%
```

误差图从 `t = 2 s` 开始画；初始收敛过程在 `kf_states.png` 里，控制台输出
的是全程 RMSE。因此图里显示的 RMSE 可能和全程控制台输出不同，这是正常的。

### 调参

| 参数 | 作用 |
| --- | --- |
| `TRAJ_A`、`TRAJ_B`、`TRAJ_W` | 8 字形轨迹的形状和速度 |
| `ACCEL_NOISE` | 过程噪声强度（真值和 `Q` 同时使用） |
| `measurement_sigma()` | 正常/变差 GPS 噪声及变差区间 |
| `x0`、`P0` | 初始估计和初始不确定性 |
| `T_END`、`DT` | 仿真长度和步长（`STEPS` 自动计算） |

如果改变 GPS 噪声等级，只要“变差”样本的 `sigma` 仍大于最小值的 1.5 倍，
绘图脚本仍然能自动识别。

---

## 2. 非线性 ESKF：单轨机器人 + 本体坐标系路标

### 名义状态与运动模型

```text
x  = [px, py, yaw]^T      (名义，ESKF_N_NOM = 3)
dx = [dpx, dpy, dyaw]^T   (误差，ESKF_N = 3)
u  = [v, omega]^T         (控制量，用户定义)
```

库没有编译期的控制维数：`u` 根本不会通过 `eskf_predict()` 传入。本示例把
`u` 和 `dt` 放在 `struct model` 里，通过 `params` 传给回调。这里控制量有两个
分量，过程噪声维数恰好也是 `ESKF_NW = 2`。本示例中名义状态和误差状态都是
3 维；库允许两者不同（例如 4 元名义四元数配 3 维姿态误差）。

连续时间单轨（unicycle）模型：

$$
\dot p_x = v\cos\theta,\qquad
\dot p_y = v\sin\theta,\qquad
\dot\theta = \omega
$$

`nominal_propagate()` 使用恒定 `(v, omega)` 下的精确圆弧积分。当
$|\omega|\Delta t > 10^{-3}$ 时：

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

当 $|\omega|\Delta t$ 很小时，退化为直线积分：

$$
p_x^{+} = p_x + v\cos\theta\,\Delta t,\qquad
p_y^{+} = p_y + v\sin\theta\,\Delta t
$$

阈值判断的是 `|omega| * dt`，不是只看 `|omega|`。单精度下，精确公式里的
两个相近正弦/余弦相减会丢失大量有效位；阈值太小时，解析 $F_x$ 会产生
虚假的巨大值，机器人直行时协方差会突然跳变。`propagate()` 和
`motion_jacobians()` 使用完全相同的 `|omega| * dt > 1e-3` 分支。

### 误差状态与雅可比

```text
dx = [δpx, δpy, δyaw]^T
```

误差状态转移矩阵是精确传播的雅可比：

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

圆弧分支下：

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

小角度分支使用对应的一阶极限。代码里两个分支都有；`motion_jacobians()`
是 `compute_Fx()` 和 `compute_Fi()` 的唯一实现来源。

### 量测模型

三个已知路标：

```text
L1 = (14, 9)
L2 = (-11, 16)
L3 = (-15, -18)
```

传感器给出每个路标在机器人本体系下的位置：

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

三个路标拼成 6 维量测：

$$
z = \begin{bmatrix}
h_1 \\
h_2 \\
h_3
\end{bmatrix} + v,
\qquad v\sim\mathcal{N}(0,\sigma_m^2 I),
\qquad \sigma_m = 0.5\ \text{m}
$$

单个路标的量测雅可比为

$$
H_i =
\begin{bmatrix}
-c & -s & h_y \\
s & -c & -h_x
\end{bmatrix}
$$

其中 $h_x,h_y$ 是 $h_i$ 的两个分量。C 和 Python 使用相同的路标坐标；
如果修改地图，需要同时更新 `main.c` 和 `plot.py`。

### 误差注入与 G

```text
px  += δpx
py  += δpy
yaw  = wrap_pi(yaw + δyaw)
```

除了角度 wrap 之外，注入是加性的，因此

$$
G = \frac{\partial(x\boxplus\delta x)}{\partial\delta x} = I
$$

`compute_G()` 仍然实现了，并返回单位阵；库在每次注入后执行
$P\leftarrow GPG^{\top}$。

### 噪声模型与多速率循环

真值控制叠加执行噪声：

$$
v_{\text{true}} = v_{\text{cmd}} + \mathcal{N}(0,\sigma_v^2),
\qquad
\omega_{\text{true}} = \omega_{\text{cmd}} + \mathcal{N}(0,\sigma_\omega^2)
$$

其中 `SIGMA_V = 0.10 m/s`、`SIGMA_W = 0.04 rad/s`。滤波器只拿到命令控制量；
执行噪声通过下式进入协方差：

$$
P^{-} = F_xPF_x^{\top}+F_iQF_i^{\top}.
$$

滤波器使用的过程噪声在真值噪声基础上放大，以吸收线性化误差：

$$
Q = Q_{\text{SCALE}}\cdot\mathrm{diag}(\sigma_v^2,\sigma_\omega^2),
\qquad Q_{\text{SCALE}}=2.
$$

时间调度：

```text
dt = 0.01 s              -> 100 Hz 预测
MEAS_EVERY = 10          -> 10 Hz 量测更新
T_END = 60 s             -> 6000 步，600 次更新
```

循环每一步调用 `eskf_predict()`，仅当 `k % MEAS_EVERY == 0` 时调用
`eskf_update()`。两个滤波器都把预测和更新拆开了；ESKF 示例利用这一点让
二者以不同频率运行。

### 控制曲线与初始条件

```text
v(t)     = 2.0 + 0.6 sin(0.17 t)                  [m/s]
omega(t) = 0.55 sin(0.23 t) + 0.20 sin(0.61 t)    [rad/s]

初始估计: x_nom0 = (1.2, -0.8, 0.20)
初始协方差: P0 = diag(4, 4, 0.35²)
```

真值从 `(0,0,0)` 出发。航位推算（dead reckoning）从同一个偏移初值出发，
只积分命令控制量，不做任何量测更新。

### 库调用与回调

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
        /* 构造 z, R */
        eskf_update(&f, z, R, &model);
    }
}
```

7 个回调全部实现。`params` 指向 `struct model`，里面保存三个路标、当前控制
量 `u` 和 `dt`。库源码没有修改。注意 `eskf_init()` 的回调顺序：
`compute_Hx` 在 `measure_model` 前面。

### CSV 列

```text
t, updated,
true_x, true_y, true_yaw,
est_x, est_y, est_yaw,
dr_x, dr_y, dr_yaw,
Pxx, Pxy, Pyy, Pyaw,
v_cmd, w_cmd
```

`yaw` 列已经 wrap 到 $[-\pi,\pi]$。Python 脚本用 wrap 后的差值计算航向误差。
`Pxx/Pyy/Pyaw` 是方差，`Pxy` 是位置协方差。

### 结果（固定随机种子）

```text
ESKF 位置 RMSE : 0.078 m     最终误差: 0.015 m
ESKF 航向 RMSE : 0.518°      最终误差: 0.390°
航位推算位置   : 3.120 m     最终误差: 1.738 m
航位推算航向   : 11.754°     最终误差: 13.011°

位置精度相对航位推算提升约 40 倍
x/y/yaw 落在 2σ 内的比例: 97.5% / 93.7% / 99.4%
```

`eskf_errors.png` 左列是全程，右列是 `t ≥ 1 s` 的稳态放大。
`eskf_vs_dr.png` 使用对数纵轴，可以同时看清初始收敛和 ESKF 与航位推算的
稳态差距。

### 调参

| 参数 | 作用 |
| --- | --- |
| `SIGMA_V`、`SIGMA_W` | 真值控制噪声，同时决定滤波器 `Q` 的基础值 |
| `Q_SCALE` | 滤波器过程噪声的放大系数 |
| `SIGMA_MEAS`、`MEAS_EVERY` | 量测噪声和更新频率 |
| `control()` | 机器人轨迹 |
| `x_nom0`、`P0` | 初始估计和初始不确定性 |
| `model.landmark` | 地图；同时修改 `plot.py` 的 `LANDMARKS` |
| Makefile 里的 `ESKF_M` | 等于路标数的两倍 |
| `T_END`、`DT` | 仿真长度和步长（`STEPS` 自动计算） |

如果改变路标数量，要同时修改 Makefile 里的 `-DESKF_M`（`2 × 路标数`）和
`model.landmark` 初始化。`main.c` 的 `NUM_LANDMARKS` 会自动从 `ESKF_M`
推导；`plot.py` 仍然需要手动同步自己的路标列表。

---

## 编译期维数

库结构体里是固定长度数组，因此下面这些宏必须同时传给**库的 `.c` 文件**和
**用户的 `.c` 文件**：

| 示例 | Makefile 参数 |
| --- | --- |
| KF | `-I../../kf -DKF_N=6 -DKF_M=2 -DKF_C=2` |
| ESKF | `-I../../eskf -DESKF_N=3 -DESKF_N_NOM=3 -DESKF_M=6 -DESKF_NW=2` |

每个 `main.c` 还有编译期检查：

```c
#if KF_N != 6 || KF_M != 2 || KF_C != 2
#error "Compile this example with KF_N=6, KF_M=2 and KF_C=2"
#endif

#if ESKF_N != 3 || ESKF_N_NOM != 3 || ESKF_M != 6 || ESKF_NW != 2
#error "Compile this example with ESKF_N=3, ESKF_N_NOM=3, ESKF_M=6 and ESKF_NW=2"
#endif
```

如果不加 `-D` 直接编译，头文件会退回维数 1，检查会给出明确错误，而不是
产生结构体布局不一致的未定义行为。`ESKF_N_NOM` 是名义状态维数，
`ESKF_N` 是误差状态维数；本示例中两者相等，但在四元数姿态估计等模型里
可以不同。

## 说明

- 两个示例都使用固定种子的 xorshift64* + Box-Muller 随机数。修改
  `struct rng rng = {...}` 的种子即可改变噪声实现。
- C 程序接受一个可选的输出路径参数；直接运行 `./kf_demo` 或 `./eskf_demo`
  会写到当前目录。
- 绘图脚本默认使用 `Agg` 后端，只保存 PNG；加 `--show` 会用可用的 GUI 后端
  打开图片。
- 生成的 CSV 和图已被 Git 忽略；`make demo` 可以重新生成。
- 示例可以用 `-std=c11 -Wall -Wextra -Wpedantic` 干净编译，并且用
  AddressSanitizer/UndefinedBehaviorSanitizer 运行过。
