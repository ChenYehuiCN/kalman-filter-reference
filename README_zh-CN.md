[English](README.md)

# 卡尔曼滤波参考实现（KF / ESKF）

一个紧凑的、固定尺寸的 C 语言线性卡尔曼滤波（KF）和误差状态卡尔曼滤波
（ESKF）实现。

这是**参考代码**，不是开箱即用的库。它的首要目标是让公式、数据布局和调用
流程容易阅读和修改。实现刻意写得比较通用：支持控制输入、时变矩阵、ESKF
非线性回调和单独的过程噪声雅可比，即使某个具体应用根本用不到这些功能。
使用者应当把需要的部分复制出去，按自己的模型裁剪、特化。

唯一的外部依赖是 `libm`（`sqrtf`、`sinf`、`cosf` 等）。没有堆分配，也没有
库内部的全局可变状态。

## 目录

- [一览](#一览)
- [仓库结构](#仓库结构)
- [快速开始](#快速开始)
- [维数与数据布局](#维数与数据布局)
- [线性卡尔曼滤波（KF）](#线性卡尔曼滤波kf)
- [误差状态卡尔曼滤波（ESKF）](#误差状态卡尔曼滤波eskf)
- [裁剪与优化](#裁剪与优化)
- [示例](#示例)
- [验证](#验证)
- [许可证](#许可证)

## 一览

| 项目 | KF（`kf/`） | ESKF（`eskf/`） |
| --- | --- | --- |
| 模型 | 线性状态空间 | 非线性名义模型 + 误差状态 |
| 状态维数 | `KF_N` | `ESKF_N_NOM`（名义）、`ESKF_N`（误差） |
| 量测维数 | `KF_M` | `ESKF_M` |
| 输入/噪声维数 | `KF_C`（控制） | `ESKF_NW`（噪声） |
| 主要调用 | `kf_predict()`、`kf_update()` | `eskf_predict()`、`eskf_update()` |
| 协方差更新 | Joseph 形式 | Joseph 形式 |
| 线性求解 | Cholesky 分解 + 前代/回代 | 同左 |
| 存储 | 结构体内的固定长度数组 | 结构体内的固定长度数组 |
| 数据类型 | `float` | `float` |

头文件里的所有维数默认都是 `1`。这只是为了让头文件能单独编译；实际使用时
必须定义维数，并用**同一组宏**编译库的 `.c` 文件（见
[维数与数据布局](#维数与数据布局)）。

## 仓库结构

```text
.
├── kf/
│   ├── kf.h              # KF API 与结构体
│   └── kf.c              # KF 实现
├── eskf/
│   ├── eskf.h            # ESKF API、回调与结构体
│   └── eskf.c            # ESKF 实现
├── example/
│   ├── README.md         # 示例说明（英文）
│   ├── README_zh-CN.md   # 示例说明（中文）
│   ├── kf_tracking/      # 线性 KF 示例
│   └── eskf_localization/# 非线性 ESKF 示例
├── LICENSE               # MIT
├── README.md
└── README_zh-CN.md
```

## 快速开始

仓库没有顶层构建系统，因为库本身只有两个源文件。示例展示了推荐的编译参数
和完整用法：

```bash
cd example
make demo          # 编译、运行两个仿真并生成全部图
```

依赖：C11 编译器、`make`、Python 3、`numpy`、`matplotlib`。
模型细节、输出列和调参说明见 [`example/README_zh-CN.md`](example/README_zh-CN.md)。
英文版在 [`example/README.md`](example/README.md)。

### 最小 KF 编译命令

```bash
cc -I kf -DKF_N=6 -DKF_M=2 -DKF_C=2 \
   user.c kf/kf.c -lm -o user_kf
```

### 最小 ESKF 编译命令

```bash
cc -I eskf -DESKF_N=3 -DESKF_N_NOM=3 -DESKF_M=6 -DESKF_NW=2 \
   user.c eskf/eskf.c -lm -o user_eskf
```

实际使用时 `-D` 参数基本是必须的：结构体大小和库 `.c` 文件里的数组维数都
由这些宏决定。

## 维数与数据布局

### 编译期配置

| 宏 | 含义 | 使用者 |
| --- | --- | --- |
| `KF_N` | 状态维数 | `kf.h`、`kf.c` |
| `KF_M` | 量测维数 | `kf.h`、`kf.c` |
| `KF_C` | 控制输入维数 | `kf.h`、`kf.c` |
| `ESKF_N` | 误差状态维数（`dx`、`P`、`Fx`、`G`） | `eskf.h`、`eskf.c` |
| `ESKF_N_NOM` | 名义状态维数（`x_nom`） | `eskf.h`、`eskf.c` |
| `ESKF_M` | 量测维数 | `eskf.h`、`eskf.c` |
| `ESKF_NW` | 过程噪声维数（`Fi` 列数、`Q` 大小） | `eskf.h`、`eskf.c` |

`ESKF_N_NOM` 和 `ESKF_N` 是两个独立的维数，通常相等，也可以不相等（例如
4 参数四元数名义状态配 3 维姿态误差）。`x_nom` 使用 `ESKF_N_NOM`；`dx`、
`P`、`Fx`、`Hx` 的列数和 `G` 使用 `ESKF_N`。定义 `ESKF_N` 不会自动设置
`ESKF_N_NOM`；两者默认值都是 `1`，实际编译 ESKF 时要同时定义。

库和用户编译单元必须看到完全相同的值。只要库和用户代码看到的任何一个维数
宏不一致，两个编译单元对滤波器结构体的布局就会不同，这是未定义行为。

### 矩阵存储

所有矩阵都是**行优先**存储，没有 leading dimension 参数：

```c
A[i * cols + j]   // 第 i 行、第 j 列
```

向量就是普通数组。例如 KF 状态是 `f.x[KF_N]`，量测是 `z[KF_M]`，协方差是
`f.P[KF_N * KF_N]`，其中 `P[i * KF_N + j]` 对应第 i 行第 j 列。

所有矩阵元素都是 `float`。如果需要 `double`，要在结构体、API 和辅助函数里
统一修改（没有可切换的 typedef）。

### 所有权与可重入性

所有工作内存都在滤波器结构体里：没有 `malloc`，库内部没有静态可变数据，
除了结构体之外没有隐藏状态。两个滤波器可以用两个结构体独立运行。用户的
回调当然可以通过 `params` 访问自己的全局数据。

以下划线开头的字段是内部工作区（`_tmp_*`、`_S`、`_K`、`_y` 等）。它们只是
因为 C 没有私有成员才暴露出来；不要把它们的布局当成稳定 API。

## 线性卡尔曼滤波（KF）

### 模型

库实现标准的线性高斯模型：

$$
x_k = F_{k-1} x_{k-1} + B_{k-1} u_{k-1} + w_{k-1},
\qquad w_{k-1} \sim \mathcal{N}(0, Q_{k-1})
$$

$$
z_k = H_k x_k + v_k,
\qquad v_k \sim \mathcal{N}(0, R_k)
$$

API 把时间更新和量测更新拆成了两个函数：

- `kf_predict(f, u, F, B, Q)` 传播状态和协方差；
- `kf_update(f, z, H, R)` 应用一次量测。

矩阵都是函数参数，不作为结构体字段保存，因此 `F`、`B`、`H`、`Q`、`R` 可以
每一步变化。没有量测时可以不调用 `kf_update()`。

### 预测

$$
\hat{x}^{-} = F\hat{x} + Bu
$$

$$
P^{-} = FPF^{\top} + Q
$$

`kf_predict()` 完成这两步。

### 更新

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

最后一行是 **Joseph 形式**。它比简单的 $(I-KH)P^{-}$ 计算量大，但在有舍入
误差时能保持协方差对称、半正定。实现没有额外对 `P` 做对称化。
`kf_update()` 完成更新公式。

### 代码对应关系：预测

| 公式 | `kf_predict()` 中的代码 |
| --- | --- |
| $\hat{x}^{-}=F\hat{x}+Bu$ | `_tmp_n1 = F x`；`_tmp_n1_2 = B u`；相加；拷贝到 `x` |
| $P^{-}=FPF^{\top}+Q$ | `_tmp_nn = F P`；`P = _tmp_nn F^T`；`P += Q` |

### 代码对应关系：更新

| 公式 | `kf_update()` 中的代码 |
| --- | --- |
| $S=HP^{-}H^{\top}+R$ | `_tmp_mn = H P`；`_S = _tmp_mn H^T`；`_S += R` |
| $S=LL^{\top}$ | `mat_chol(_S, KF_M)`，原地把 `_S` 变成下三角 `L` |
| $K=P^{-}H^{\top}S^{-1}$ | `mat_chol_solve(_S, _tmp_mn, M, N)` 解 $SX=HP^{-}$；`_K = X^{\top}` |
| $\tilde{y}=z-H\hat{x}^{-}$ | `_y = H x`；`_y = z - _y` |
| $\hat{x}=\hat{x}^{-}+K\tilde{y}$ | `_tmp_n1 = K _y`；`x += _tmp_n1` |
| Joseph 形式 | `_tmp_nn = K H`；`I - K H`；两次 `mat_mul_bt`；再加 $KRK^{\top}$ |

增益是通过求解 $SX=HP^{-}$（对 $S$ 做 Cholesky 分解）再转置得到的，没有显式
求 $S^{-1}$。Cholesky 例程假设 $S$ 对称正定，不做主元选取。

### API

```c
struct kf_filter f;

kf_init(&f);                 // 只需一次：初始化单位阵工作区
kf_reset(&f, x0, P0);        // 拷贝初始状态和协方差

for (;;) {
    kf_predict(&f, u, F, B, Q);        // 时间更新
    if (measurement_available)
        kf_update(&f, z, H, R);        // 量测更新
}
```

| 函数 | 行为 |
| --- | --- |
| `void kf_init(struct kf_filter *f)` | 把 `_I_nn` 清零并对角置一。**不**初始化 `x` 和 `P`。 |
| `void kf_reset(struct kf_filter *f, float *x, float *P)` | 把 `x`、`P` 拷贝进结构体，不修改工作区矩阵。 |
| `void kf_predict(struct kf_filter *f, float *u, float *F, float *B, float *Q)` | 用控制、状态转移和过程噪声矩阵预测 `x`、`P`。 |
| `void kf_update(struct kf_filter *f, float *z, float *H, float *R)` | 对 `x`、`P` 应用一次量测更新。 |

参数形状：

| 参数 | 形状 | 使用者 |
| --- | --- | --- |
| `u` | `KF_C` | `kf_predict()` |
| `F` | `KF_N × KF_N` | `kf_predict()` |
| `B` | `KF_N × KF_C` | `kf_predict()` |
| `Q` | `KF_N × KF_N` | `kf_predict()` |
| `z` | `KF_M` | `kf_update()` |
| `H` | `KF_M × KF_N` | `kf_update()` |
| `R` | `KF_M × KF_M` | `kf_update()` |

传给某个函数的所有指针都必须有效。如果应用没有控制输入，传一个全零的 `B`
和全零的 `u`，不要传 `NULL`，因为 `kf_predict()` 总会计算 $Bu$。

### KF 的限制

- `kf_update()` 没有丢量测、门限（gating）或异常值剔除逻辑。需要跳过更新或
  在滤波器之外处理这些情况。
- `kf_update()` 可以被调用多次（例如逐个传感器更新），这等价于连续做多次
  量测更新；如果量测相关，应把它们堆叠成一个 `z/H/R`。
- 维数是编译期常量，不支持运行时改变状态维数。
- Cholesky 求解要求 `S` 正定。正定的 `R` 和半正定的 `P` 通常能保证这一点，
  但库不做检查。

## 误差状态卡尔曼滤波（ESKF）

### 基本思想

ESKF 维护两个东西：

- **名义状态** $x \in \mathbb{R}^{N_{\mathrm{nom}}}$：用完整非线性模型传播；
- **误差状态** $\delta x \in \mathbb{R}^{N}$ 及其协方差 $P$：假设很小，用线性
  化模型传播。

两次调用之间，误差状态为零。卡尔曼修正是在误差坐标里计算的，注入到名义
状态之后，误差状态再清零。

名义状态维数是 `ESKF_N_NOM`，误差状态维数是 `ESKF_N`。它们是两个独立的
编译期常量，可以相等也可以不等。典型例子是 4 参数四元数名义状态配 3 维姿态
误差：`x_nom` 有 4 个元素，而 `dx`、`P`、`Fx`、`G` 是 3 维，`Hx` 有 3 列。
用户提供的 `inject()` 负责把名义状态和误差状态映射回 `ESKF_N_NOM` 维。
如果两个维数相等且注入是加性的，就退化为常见的 ESKF 形式。

### 模型

令名义状态 $x \in \mathbb{R}^{N_{\mathrm{nom}}}$，误差状态
$\delta x \in \mathbb{R}^{N}$。

名义状态传播：

$$
x_k = f(x_{k-1}, u_{k-1}, \Delta t)
$$

误差状态传播：

$$
\delta x_k = F_x\,\delta x_{k-1} + F_i\,w_{k-1},
\qquad w_{k-1} \sim \mathcal{N}(0, Q_{k-1})
$$

其中

$$
F_x \in \mathbb{R}^{N \times N},
\qquad
F_i \in \mathbb{R}^{N \times N_W}.
$$

当 $N_{\mathrm{nom}} = N$ 时，$F_x$ 就是通常的状态雅可比
$\partial f/\partial x$。当两者不同时，$F_x$ 是所选局部误差状态传播的
雅可比，由 `compute_Fx()` 提供。

$F_i$ 把 `ESKF_NW` 维噪声向量映射到误差状态。如果噪声是误差状态加性的，令
$F_i = I$ 且 `ESKF_NW = ESKF_N`，此时 $Q$ 就是通常的状态过程噪声协方差。

上面的控制 $u$ 和时间步 $\Delta t$ 属于用户模型。库不会把它们传给回调；
回调需要从 `params` 或自己的上下文中读取。`eskf_predict()` 只接收 `Q` 和
`params`。

量测模型：

$$
z = h(x) + v,
\qquad v \sim \mathcal{N}(0,R)
$$

### 预测

$$
x^{-} = f(x,u,\Delta t)
$$

$$
P^{-} = F_x P F_x^{\top} + F_i Q F_i^{\top}
$$

公开入口是 `eskf_predict(f, Q, params)`。实现会在 `nominal_propagate()`
**之前**调用 `compute_Fx()` 和 `compute_Fi()`，所以两个雅可比看到的都是
传播前的名义状态；之后传播名义状态和 `P`。

### 更新

量测雅可比在 $\delta x = 0$ 处取值：

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

这里 $\boxplus$ 表示与模型相关的注入操作。注入修正并清零误差状态：

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

对于纯加性注入 $x \boxplus \delta x = x + \delta x$，注入雅可比就是 $G=I$。

### 代码对应关系

| 公式 | `eskf_predict()` / `eskf_update()` 中的代码 |
| --- | --- |
| $F_x$、$F_i$ | 传播前的 `compute_Fx()`、`compute_Fi()`；`Fx` 是 `ESKF_N × ESKF_N`，`Fi` 是 `ESKF_N × ESKF_NW` |
| $x^{-}=f(x,u,\Delta t)$ | `nominal_propagate()`；`x_nom` 长度为 `ESKF_N_NOM` |
| $P^{-}=F_xPF_x^{\top}+F_iQF_i^{\top}$ | `_tmp_nn = Fx P`；`P = _tmp_nn Fx^T`；`_tmp_nw = Fi Q`；`_tmp_nn = _tmp_nw Fi^T`；相加 |
| $H$ | 注入前的 `compute_Hx()`；`Hx` 是 `ESKF_M × ESKF_N` |
| $S=HP^{-}H^{\top}+R$ | `_tmp_mn = Hx P`；`_S = _tmp_mn Hx^T`；`_S += R` |
| $S=LL^{\top}$ | `mat_chol(_S, ESKF_M)` |
| $K=P^{-}H^{\top}S^{-1}$ | `mat_chol_solve()` 后转置到 `_K` |
| $h(x^{-})$ | 注入前的 `measure_model()` |
| $\delta\hat{x}=K(z-h)$ | `_y = z - _y`；`dx = K _y` |
| Joseph 形式 | 与 KF 相同的步骤 |
| $x\leftarrow x\boxplus\delta x$ | `inject(x_nom, dx, params)` |
| $P\leftarrow GPG^{\top}$ | `compute_G(G, dx, params)`；两次 `mat_mul_bt` |
| 清零 | `memset(dx, 0, ...)` |

`Fx`、`dx`、`P`、`G` 使用误差状态维数 `ESKF_N`；`x_nom` 使用名义状态维数
`ESKF_N_NOM`；`Hx` 的形状是 `ESKF_M × ESKF_N`。

### API 与调用顺序

```c
struct eskf_filter f;

eskf_init(&f, compute_Fx, compute_Fi, nominal_propagate,
          compute_Hx, measure_model, inject, compute_G);
eskf_reset(&f, x_nom0, P0);

for (;;) {
    /* 把当前控制量和 dt 放到回调能读到的地方。 */
    eskf_predict(&f, Q, params);
    if (measurement_available)
        eskf_update(&f, z, R, params);
}
```

`eskf_init()` 的回调顺序是 `compute_Fx`、`compute_Fi`、
`nominal_propagate`、`compute_Hx`、`measure_model`、`inject`、
`compute_G`。

| 函数 | 行为 |
| --- | --- |
| `eskf_init(...)` | 保存 7 个回调，清零 `dx`，初始化单位阵工作区。**不**初始化 `x_nom` 和 `P`。 |
| `eskf_reset(&f, x_nom, P)` | 拷贝名义状态（`ESKF_N_NOM` 个元素）和协方差（`ESKF_N × ESKF_N`），清零 `dx`。 |
| `eskf_predict(&f, Q, params)` | 计算 `Fx`/`Fi`，清零 `dx`，传播名义状态，传播 `P`。 |
| `eskf_update(&f, z, R, params)` | 计算 `Hx`、增益、残差和 `dx`；做 Joseph 协方差更新，注入 `dx`，应用 `G`，最后清零 `dx`。 |

| 参数 | 形状 | 说明 |
| --- | --- | --- |
| `Q` | `ESKF_NW × ESKF_NW` | 噪声空间中的过程噪声协方差 |
| `z` | `ESKF_M` | 量测向量 |
| `R` | `ESKF_M × ESKF_M` | 量测噪声协方差 |
| `params` | `void *` | 传给所有回调 |

控制量和时间步**不是**库的参数。如果 `compute_Fx()`、`compute_Fi()` 或
`nominal_propagate()` 需要它们，请在调用 `eskf_predict()` 之前把它们放进
`params` 对象（或用户自己的上下文）。同一个 `params` 指针也会传给
`eskf_update()`。

和 KF 不同，ESKF 没有编译期的控制维数。`ESKF_NW` 只决定噪声向量和 `Q` 的
大小，不决定用户控制数据的长度；回调自己解释控制数据。示例中恰好使用两个
控制量和两个噪声输入，但这不是库的约束。

### 回调契约

```c
void compute_Fx(float *Fx, float *x_nom, void *params);
void compute_Fi(float *Fi, float *x_nom, void *params);
void nominal_propagate(float *x_nom, void *params);
void compute_Hx(float *Hx, float *x_nom, void *params);
void measure_model(float *z_pred, float *x_nom, void *params);
void inject(float *x_nom, float *dx, void *params);
void compute_G(float *G, float *dx, void *params);
```

| 回调 | 调用时机 | 输出形状 | 说明 |
| --- | --- | --- | --- |
| `compute_Fx` | `eskf_predict()` 中，名义传播之前 | `ESKF_N × ESKF_N` | 旧名义状态处的误差状态转移雅可比 |
| `compute_Fi` | `eskf_predict()` 中，名义传播之前 | `ESKF_N × ESKF_NW` | 噪声雅可比，不一定是控制矩阵 |
| `nominal_propagate` | `eskf_predict()` 中，雅可比之后 | 原地修改 `x_nom`（`ESKF_N_NOM`） | 完整非线性传播 |
| `compute_Hx` | `eskf_update()` 开头 | `ESKF_M × ESKF_N` | 注入前状态处的量测雅可比 |
| `measure_model` | `eskf_update()` 中，注入之前 | `ESKF_M` | 预测量测 $h(x_{\text{nom}})$ |
| `inject` | `eskf_update()` 中，协方差更新之后 | 原地修改 `x_nom`（`ESKF_N_NOM`） | 执行 $x \leftarrow x\boxplus\delta x$；接收 `dx` 和 `params` |
| `compute_G` | `inject()` 之后、清零 `dx` 之前 | `ESKF_N × ESKF_N` | 注入雅可比；会收到已经应用的 `dx` |

七个回调都是必需的。库通过函数指针调用，不检查 `NULL`。如果某个回调很简单
（例如 $G=I$），也要提供一个填好输出数组的函数。

`params` 会传给全部 7 个回调（包括 `inject`）。库本身不解释 `params` 的
内容。回调需要模型数据、控制量或当前时间步时，把它们放在 `params` 里。

### 量测残差与角度 wrap

更新使用原始残差

$$
\tilde{y} = z - h(x_{\text{nom}})
$$

库**不会**自动 wrap 角度、归一化四元数或做其他残差变换。如果量测是角度，
需要保证 `measure_model` 返回的预测和 `z` 在同一分支上，或者在外部预处理
`z`/量测模型。示例通过使用本体坐标系下的笛卡尔路标坐标避开了这个问题。

### ESKF 的限制

- `ESKF_N_NOM` 和 `ESKF_N` 是两个独立的编译期维数。用户必须保证 `inject`
  （把 `ESKF_N_NOM` 名义状态和 `ESKF_N` 误差状态映射回 `ESKF_N_NOM`）和
  `G`（`ESKF_N × ESKF_N`）与它们一致。
- 控制量和时间步不会传给回调。回调必须从 `params` 或自己的上下文中读取，
  并且 `compute_Fx`、`compute_Fi`、`nominal_propagate` 要使用同一组值。
- 残差不会自动 wrap，也没有门限/异常值处理；角度分支、丢量测、异常值需要
  在核心更新之外处理。
- `compute_G` 只收到 `dx` 和 `params`，不收到 `x_nom`。如果你的注入雅可比
  依赖注入后的名义状态，可以在 `inject` 中把该状态存进 `params`，或特化库
  代码。
- 协方差更新是 Joseph 形式；不会显式对称化 `P`。
- Cholesky 要求 `S` 正定；没有主元选取或正则化。
- 每次调用只做一次更新。多传感器可以堆叠量测，也可以连续调用多次
  `eskf_update()`，但要保证噪声独立假设成立。
- 所有维数都是编译期常量，所有数据都是 `float`。

## 裁剪与优化

代码刻意保留了通用的矩阵路径。实际产品里通常可以删掉很大一部分。常见做法：

- **没有控制输入：** 去掉 `B`/`u`、KF 预测里的控制项，特化函数时还可以去掉
  `_tmp_n1_2` 工作区。
- **不需要过程噪声雅可比：** 如果噪声在误差状态里是加性的，令 `Fi = I`、
  `ESKF_NW = ESKF_N`，并删掉 `Fi Q Fi^T` 的乘法。
- **名义/误差同维：** 如果 `ESKF_N_NOM == ESKF_N` 且注入是加性的，名义状态
  可以和误差状态使用同一套布局，且 `G = I`。
- **不需要注入雅可比：** 令 `G = I`，删掉注入后的两次 `mat_mul_bt`。
- **常量矩阵：** 在循环外构造一次 `F`、`B`、`H`、`Q`、`R`，不要每步重新传。
- **稀疏结构：** 如果 `H` 只是选择状态，用索引代替通用矩阵乘法；如果 `Q`
  或 `R` 是对角阵，利用对角性。
- **小固定维数：** 展开循环，只存 `P` 的上三角，避免通用 `mat_mul`。
- **不显式构造 `K`：** 直接求解修正量，或把 Cholesky 求解和状态更新融合。
- **双精度/定点：** 统一替换 `float`，不要只改某一个文件。
- **多传感器：** 堆叠量测，或写专用的 update 函数；通用的 `H`/`R`/`z` 路径
  能处理任意 `ESKF_M`/`KF_M`，但可能做很多无用功。
- **先确定假设：** 在特化之前想清楚正定性、角度 wrap、丢数据和异常值怎么
  处理。

`example/` 下的两个目录展示了两种特化方式，同时保留了库调用，适合作为裁剪
的起点。

## 示例

| 示例 | 滤波器 | 维数 | 内容 |
| --- | --- | --- | --- |
| [`example/kf_tracking`](example/kf_tracking) | 线性 KF | `KF_N=6`、`KF_M=2`、`KF_C=2` | 二维常加速度跟踪、已知 jerk 输入、时变 GPS 噪声、Joseph 形式协方差 |
| [`example/eskf_localization`](example/eskf_localization) | 非线性 ESKF | `ESKF_N=3`、`ESKF_N_NOM=3`、`ESKF_M=6`、`ESKF_NW=2` | 单轨运动、三个已知路标、本体坐标系量测、100 Hz 预测 / 10 Hz 更新、航位推算对比 |

```bash
cd example
make demo
```

示例会在 `example/*/output/` 下生成 CSV 和 Python 图。模型公式、输出列、
数值细节和调参说明见 [`example/README_zh-CN.md`](example/README_zh-CN.md)，
英文版见 [`example/README.md`](example/README.md)。

## 验证

示例使用以下选项编译：

```bash
cc -std=c11 -Wall -Wextra -Wpedantic ...
```

并且用 AddressSanitizer 和 UndefinedBehaviorSanitizer 运行过。固定随机种子
时示例输出可复现。项目没有单独的单测框架；两个示例就是集成测试。

## 许可证

MIT，见 [LICENSE](LICENSE)。

Copyright (c) 2026 Chen Yehui.
