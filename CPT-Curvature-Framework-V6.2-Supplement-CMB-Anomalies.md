# 曲率相依 CPT 破壞框架 V6.2 補充：CMB 冷斑與大角度異常

> **版本定位**：本文件為 V6.1 的獨立補充稿（V6.2），不取代 V6.1 主體。  
> **任務**：在既有 V1–V6.1 框架內，自行推演 CMB 冷斑與三類大角度異常的**候選機制**與**可檢驗量**。  
> **約束**：不引用外部文獻或網路資料；實際觀測數值**暫留空**，僅給結構與量級口徑。  
> **文件版本**：V6.2（2026-05-09，對齊 V6.1-patch）

---

## 〇、V6.1-patch 依賴與可擴展性確認

### 0.1 修補後依賴關係

V6.2 沿用 V6.1 的 $f(R)$、$\xi$、W 投影，**不新增基本參數**。  
自 V6.1-patch（2026-05-09）起，擴展前提如下：

| V6.1 物件 | V6.2 用法 | 狀態 |
|-----------|-----------|------|
| $f(R)=\varepsilon_0\tanh(R/R_\ast)$ | $\mathcal{S}_\ell,\mathcal{A},K(\theta)$ | ✓ 可擴展 |
| $\rho_{\rm DE}^\star=4\pi f^2 M_{\rm Pl}^2 H^2$ | §七 量級相容 | ✓ 量級可對；$\Omega$ 份額見 V6.1 附錄 B3 |
| $3-4\pi f^2>0$ | **背景** $z\lesssim 1$ 才良定；LSS $f\to\varepsilon_0$ 时 $4\pi f^2\gg 3$ | ⚠ **已确认**（见 C7、B10） |
| $\xi(x,t)$ | 冷斑、相關缺失 | ✓ 可接數值 |

### 0.2 可擴展 vs 待閉合

**可擴展**：§二–§六 公式 → `compute_cmb_ir_v62.py`（待寫）；与 $w(z)$ 共用 $R_\ast^{\rm late}$；$TB/EB$ 接 V5 手征项。

**待閉合**：V6.1 **B3**（$\Omega$ 双份额）；**B10**（LSS 背景不可套用 §6.2 代数 Friedmann）；$\mathcal{I}_\ell,K(\theta)$ 无数值；IR 边界层未定理化。

---

## 目錄

1. [V6.2 在框架中的位置](#一v62-在框架中的位置)
2. [概念框架：IR CPT 邊界層](#二概念框架ir-cpt-邊界層)
3. [宇宙微波背景冷斑](#三宇宙微波背景冷斑)
4. [大角度異常 I：四極矩低谷](#四大角度異常-i四極矩低谷)
5. [大角度異常 II：對稱軸排列](#五大角度異常-ii對稱軸排列)
6. [大角度異常 III：大尺度溫度相關性缺失](#六大角度異常-iii大尺度溫度相關性缺失)
7. [與暗能量量級的定性相容性](#七與暗能量量級的定性相容性)
8. [可證偽預測（數據暫空）](#八可證偽預測數據暫空)
9. [限制與待補計算](#九限制與待補計算)
10. [小結](#十小結)

---

## 一、V6.2 在框架中的位置

### 1.1 承接 V6.1 的核心物件

V6.1 已建立：

- **統一耦合**：$f(R)=\varepsilon_0\tanh(R/R_\ast)$，進入 $S_{\rm int}\propto (c_5/M_\phi)f(R)\,\partial_\mu\phi\,\bar\psi\gamma^5\gamma^\mu\psi$。
- **三投影**：B（對湮滅不對稱）、C（可見–$\phi$ 能量交換）、W（白洞虛對長時標殘餘 → 暗能量量級）。
- **時間尺度分工**：短時標 CPT 鏡像**主動相消**；長時標壽命變異**殘餘**，由 de Sitter 視界封頂為 $\rho_{\rm DE}^\star\sim f^2(R)M_{\rm Pl}^2 H^2$。
- **隨機源**：$\xi(x,t)$，強度 $D(R)\sim\kappa f^2(R)M_{\rm Pl}^9$，相干體積 $V_{\rm coh}\sim H^{-3}$。

V6.2 **不引入新的基本參數**，只把上述物件投影到**最後散射面（LSS）的大角度（低 $\ell$）溫度場**。

### 1.2 V6.2 的新斷言（一句話）

> **CMB 大角度異常不是「另起爐灶的新物理」，而是 W 投影在 IR（超視界／大角度）上的幾何印記：  
> CPT 對稱部分在大尺度上被更有效抵消，殘餘由 $\phi$ 時間箭頭定向，並留下局域冷斑與相關性缺口。**

---

## 二、概念框架：IR CPT 邊界層

### 2.1 尺度分區

在 comoving 波數 $k$ 上，定義三类尺度（以 LSS 时刻 $t_\ast$、$H_\ast=H(t_\ast)$ 为基准）：

| 尺度 | 条件 | 框架中的主导机制 |
|------|------|------------------|
| **UV**（小角度，高 $\ell$） | $k \gg a_\ast H_\ast$ | B/C 微扰，$\Delta G[f(R)]$ 修正 $\lesssim f^2$，接近标准 SW/acoustic |
| **中间** | $k \sim a_\ast H_\ast$ | W 投影开始参与，但尚未饱和视界封顶 |
| **IR**（大角度，低 $\ell$） | $k \lesssim k_{\rm IR}$ | **IR CPT 边界层**：CPT 对称部分被系统性抵消 |

定义 IR 截断波数：
$$
k_{\rm IR} \;\equiv\; \eta_\ast\, a_\ast H_\ast,
\qquad \eta_\ast = \mathcal{O}(1)
$$
对应多极矩 $\ell_{\rm IR}\sim \eta_\ast^{-1} D_A H_\ast \sim \mathcal{O}(2\text{–}10)$（$D_A$ 为角直径距离）。

### 2.2 修正后的原初扰动谱（结构式）

在 V5/V6.1 的标准标量谱 $\mathcal{P}_\mathcal{R}^{(0)}(k)$ 上，V6.2 采用**最小修正**：
$$
\boxed{
\mathcal{P}_\mathcal{R}(k,t_\ast)
=
\mathcal{P}_\mathcal{R}^{(0)}(k)
\left[
1
- \mathcal{S}_\ell(k;\,f_\ast)
+ \mathcal{A}(k;\,\hat n_\phi)
\right]
}
$$

其中：

- $\mathcal{S}_\ell$：**各向同性**大尺度抑制因子（四极低谷、相关缺失的共同来源）。
- $\mathcal{A}$：**各向异性**因子，由 $\phi$ 时间箭头方向 $\hat n_\phi$ 控制（轴排列）。
- $f_\ast \equiv f(R_\ast)$，$R_\ast=R(t_\ast)$ 由 V6.1 背景方程给出。

**关键**：$\mathcal{S}_\ell$ 与 $\mathcal{A}$ 均 **正比于 $f^2$**（与 W 投影、暗能量封顶同源），这是 §七 相容性的核心。

### 2.3 从 V6.1 真空期望值到 CMB 温度

CMB 温度扰动在 SW 极限下：
$$
\frac{\delta T}{T}(\hat\Omega)
=
\frac{1}{3}\,\Psi(\hat\Omega,t_\ast)
+ \frac{\delta T}{T}\Big|_{\rm W}(\hat\Omega)
$$

第二项为 **W 投影 IR 残差** 的直接贡献（非 adiabatic 成分）：
$$
\frac{\delta T}{T}\Big|_{\rm W}(\hat\Omega)
=
-\frac{1}{4}\,
\frac{\delta\rho_{\rm DE}^{\rm loc}(\hat\Omega,t_\ast)}{\rho_\gamma(t_\ast)}
$$

其中 $\delta\rho_{\rm DE}^{\rm loc}$ 来自 $\xi(x,t)$ 在 LSS 时刻、角尺度 $\gtrsim \ell_{\rm IR}^{-1}$ 上的局域涨落。

---

## 三、宇宙微波背景冷斑

### 3.1 物理图像

**冷斑**在本框架中不是「空洞吞噬光子」，而是：

> 在 LSS 之前某一 comoving 区域内，白洞虚对寿命方差 $\sigma_\tau^2=f^2(R)\tau_0^2+\tau_{\rm th}^2$ 使 **W 投影局域地多吸收了有效能量到 $\phi$ 场**（C 机制局域版），  
> 等价于该 patch 在再组合时 **有效辐射温度低于** 周围平均。

### 3.2 局域 Langevin 涨落

由 V6.1 §3.1，局域未配对能量流 $\xi(x,t)$ 满足：
$$
\langle\xi(x,t)\rangle=0,\qquad
\langle\xi(x,t)\xi(x',t')\rangle
= D(R)\,\delta(t-t')\,\delta^3(x-x')
$$

在 LSS 时刻 $t_\ast$，对 comoving 半径 $R_p$ 的球 patch 积分：
$$
\delta\rho_{\rm DE}^{\rm loc}
\sim
\frac{1}{V_p}\int_{V_p} d^3x\int_{t_\ast-\Delta t}^{t_\ast} dt\,\xi(x,t)
$$

其中 $\Delta t \sim \tau_0 \sim M_{\rm Pl}^{-1}$ 在宇宙学单位下极短；**有效涨落**来自 $t \lesssim t_\ast$ 上累积的 **Ornstein–Uhlenbeck 型** 方差（V6.1 §3.3）：
$$
\langle(\delta\rho_{\rm DE}^{\rm loc})^2\rangle
\sim
\frac{D_{\rm eff}(R_\ast)}{2\gamma_\ast}\,
\frac{V_{\rm coh,\ast}}{V_p}
$$

- $D_{\rm eff}=D(R_\ast)/V_{\rm coh,\ast}$，$V_{\rm coh,\ast}\sim H_\ast^{-3}$。
- $\gamma_\ast=3H_\ast(1+w_{\rm DE,\ast})$，$w_{\rm DE,\ast}\approx -1$ ⇒ $\gamma_\ast$ 很小 ⇒ 方差**大**（但受 §3.3 封顶约束）。

### 3.3 冷斑温度偏移（推导）

温度偏移：
$$
\boxed{
\frac{\delta T}{T}\Big|_{\rm cold}
\;\approx\;
-\,\frac{1}{4}\,
\frac{\sqrt{\langle(\delta\rho_{\rm DE}^{\rm loc})^2\rangle}}{\rho_\gamma(t_\ast)}
\;\sim\;
-\,\frac{1}{4}\,
\frac{f(R_\ast)\,M_{\rm Pl}\,H_\ast}{\rho_\gamma(t_\ast)}
\cdot
\mathcal{G}(R_p/H_\ast^{-1})
}
$$

其中 $\mathcal{G}$ 为几何因子：

- $R_p \ll H_\ast^{-1}$：$\mathcal{G}\sim (R_p H_\ast)^{3/2}$（小 patch，涨落弱）。
- $R_p \sim H_\ast^{-1}$：$\mathcal{G}\sim \mathcal{O}(1)$（**冷斑典型角尺度** $\theta_c \sim H_\ast D_A \sim$ 几度）。
- $R_p \gg H_\ast^{-1}$：受视界封顶，$\mathcal{G}$ 饱和。

**符号**：负号表示 **能量转入 $\phi$ / DE 残差 ⇒ 光子温度降低**。

### 3.4 冷斑的空间数密度（结构预测）

虚对生成率 $\Gamma_\mathcal{P}\sim \rho_{\rm Pl}[1-f^2]^{-1/2}$。  
在 $f\ll 1$ 的 LSS 背景上，**极端冷斑**（$|\delta T/T|$ 超过 $\sigma$ 若干倍）来自 $\sigma_\tau$ 尾部分布，数密度粗略为：
$$
n_{\rm cold}(>\nu\sigma)
\;\sim\;
\Gamma_\mathcal{P}\,V_{\rm survey}\,
\mathrm{erfc}(\nu/\sqrt{2})
\cdot
f^2(R_\ast)
$$

⇒ **冷斑 abundance 与 $f^2$ 成正比**——与今日 $\rho_{\rm DE}$ 使用**同一耦合**。

### 3.5 与 V6.1 的衔接

- 冷斑 **不是** 第二套自由参数；强度由 $f(R_\ast)$ 与 $H_\ast$ 锁定。
- 若 $f(R_{\rm now})\approx 0.09$（V6.1 §4.4），则 $f(R_\ast)$ 由 $R_\ast/R_{\rm now}$ 及 $\tanh$ 形状决定（**数值留空**，待背景 ODE 代入）。

---

## 四、大角度异常 I：四极矩低谷

### 4.1 概念

**四极矩低谷**：$C_2^{TT}$（及相关的 $C_2^{EE}$）低于标准 $\Lambda$CDM 预期。

在本框架中，这是 **IR CPT 边界层** 的最典型表现：

> 整球天区在 $\ell=2$ 上近似**单一超视界模式**；  
> 此时 CPT 镜像对称的「正反成对」扰动分量在统计上被 **最大程度抵消**；  
> 观测到的 $C_2$ 主要来自 **$f^2$ 加权的破缺残差**。

### 4.2 抑制因子 $\mathcal{S}_\ell$ 的构造

定义 CPT 对称部分的投影算符（在 $k$ 空间）：
$$
\hat{\mathcal{P}}_{\rm CPT-symm}:\;
\mathcal{R}_{\rm sym}(k)=
\frac{1}{2}\big[\mathcal{R}(k)+\mathcal{R}^{\rm CPT}(k)\big]
$$

V6.1 §5.1 断言 $\langle \mathcal{R}_{\rm sym}\rangle=0$ 在 IR 极限 **严格成立**；  
对 **涨落谱** 而言，等价于 power 被减去一份对称分量：
$$
\mathcal{P}_\mathcal{R}(k)
=
\mathcal{P}_\mathcal{R}^{(0)}(k)
-
f^2(R_\ast)\,\mathcal{P}_\mathcal{R}^{(0)}(k)\,
\Theta\!\left(\frac{k_{\rm IR}-k}{k_{\rm IR}}\right)
$$

定义各向同性抑制：
$$
\mathcal{S}_\ell(k)
=
f^2(R_\ast)\,
\Theta\!\left(\frac{k_{\rm IR}-k}{k_{\rm IR}}\right)
$$

### 4.3 从 $\mathcal{P}_\mathcal{R}$ 到 $C_\ell^{TT}$

在 SW + 无 ISW 近似（仅 LSS 大角度）：
$$
C_\ell^{TT}
\approx
\frac{4\pi}{2\ell+1}
\int dk\,k^2
\left|\Delta_\ell^{\rm SW}(k)\right|^2
\mathcal{P}_\mathcal{R}(k)
$$

对 $\ell=2$，$\Delta_2^{\rm SW}(k)\approx \tfrac{1}{3} j_2(k D_A)$ 在 $k D_A \lesssim 1$ 非零。

因此：
$$
\boxed{
\frac{C_2^{TT}}{C_{2,\Lambda{\rm CDM}}^{TT}}
\;\approx\;
1 - f^2(R_\ast)\,\mathcal{I}_2
}
$$

其中
$$
\mathcal{I}_2
=
\frac{\displaystyle\int_0^{k_{\rm IR}} dk\,k^2 |j_2(kD_A)|^2 \mathcal{P}_\mathcal{R}^{(0)}(k)}
{\displaystyle\int_0^\infty dk\,k^2 |j_2(kD_A)|^2 \mathcal{P}_\mathcal{R}^{(0)}(k)}
\;\in\;(0,1)
$$

**性质**：

- $f\to 0$ ⇒ $C_2\to C_{2,\Lambda{\rm CDM}}$（恢复标准）。
- $f^2$ 越大 ⇒ 四极越低谷越深。
- **不需要额外参数**：深度由 $f(R_\ast)$ 与 $\mathcal{I}_2$（纯几何积分）决定。

### 4.4 向更高 $\ell$ 的过渡

对 general $\ell$：
$$
\frac{C_\ell^{TT}}{C_{\ell,\Lambda{\rm CDM}}^{TT}}
\approx
1 - f^2(R_\ast)\,\mathcal{I}_\ell,
\qquad
\mathcal{I}_\ell \downarrow \text{ as } \ell \uparrow
$$

⇒ **大角度异常集中在 $\ell \lesssim \ell_{\rm IR}$**，与观测上「低 $\ell$ 异常、高 $\ell$ 正常」**定性一致**。

---

## 五、大角度异常 II：对称轴排列

### 5.1 概念

**对称轴排列（axis alignment）**：低 $\ell$ 模式的偶极/四极/八极主轴彼此对齐，或与某个物理方向（如 $\phi$ 梯度、运动方向）对齐。

在本框架中，**唯一天然的优选方向** 是 **$\phi$ 时间箭头**：

$$
\hat n_\phi \;\equiv\; \frac{\nabla\phi}{|\nabla\phi|}
$$

在 FRW 背景 comoving 系中，$\phi=\phi(t)$ ⇒ $\hat n_\phi$ 指向 **时间演化方向在天空图的投影**（与观测者共动系中的「特殊方向」耦合）。

### 5.2 各向异性修正 $\mathcal{A}$

在 $\mathcal{P}_\mathcal{R}$ 上加入 **轴对称调制**（最小形式）：
$$
\mathcal{A}(k,\hat k)
=
\alpha_a\, f(R_\ast)\,
\frac{c_5\,\dot\phi_\ast}{M_\phi}\,
P_2(\hat k\cdot\hat n_\phi)
$$

- $P_2$ 为 Legendre 二阶多项式（四极型调制）。
- $\dot\phi_\ast=\dot\phi(t_\ast)$。
- $\alpha_a=\mathcal{O}(1)$，**不引入新耦合常数**（整体仍 $\propto f\cdot b_0$）。

于是：
$$
\mathcal{P}_\mathcal{R}(k,\hat k)
=
\mathcal{P}_\mathcal{R}^{(0)}(k)
\left[
1 - f^2(R_\ast)\mathcal{S}_\ell(k)
+ \alpha_a f(R_\ast)\frac{c_5\dot\phi_\ast}{M_\phi} P_2(\hat k\cdot\hat n_\phi)
\right]
$$

### 5.3 从功率到 $a_{\ell m}$ 的相关

球谐系数：
$$
a_{\ell m}=\int d\Omega\,Y_{\ell m}^*(\hat\Omega)\,\frac{\delta T}{T}(\hat\Omega)
$$

轴对称调制使 **不同 $\ell$ 的 $a_{\ell m}$ 通过共同的 $\hat n_\phi$ 关联。  
定义主轴对齐统计量（结构式，数值留空）：
$$
\mathcal{A}_{\rm align}
\;\equiv\;
\frac{
\sum_{\ell,\ell'} w_{\ell\ell'}\,
\hat n_\ell \cdot \hat n_{\ell'}
}{
\sum_{\ell,\ell'} w_{\ell\ell'}
}
$$

本框架预测：
$$
\boxed{
\mathcal{A}_{\rm align}
\;\propto\;
f(R_\ast)\,\frac{c_5|\dot\phi_\ast|}{M_\phi}
\;\sim\;
f(R_\ast)\,|b_{0,\ast}|
}
$$

**与四极低谷同源**：$\mathcal{A}$ 与 $\mathcal{S}_\ell$ 均 $\propto f$ 或 $f^2$，  
⇒ **轴排列与四极抑制应同时出现**，而非独立巧合。

### 5.4 与 parity-odd 的联系（V5/V6.1 延伸）

V5 已指出 parity-odd 的 $TB/EB$ 可能来自 B 机制。  
V6.2 补充：**轴排列方向 $\hat n_\phi$ 也是 parity-odd 模的优选轴**，  
⇒ 低 $\ell$ 的 $TB/EB$ 非零分量应与 $\hat n_\phi$ **共线或正交**（可证伪关联，数据暂空）。

---

## 六、大角度异常 III：大尺度温度相关性缺失

### 6.1 概念

**大尺度温度相关缺失**：在角度分离 $\theta \gtrsim \theta_{\rm IR}$ 的两点，观测到的 $C(\theta)=\langle \delta T(\hat n_1)\delta T(\hat n_2)\rangle$ **低于** 标准 $\Lambda$CDM 外推。

### 6.2 为什么相关性会缺失（核心论证）

标准相关来自 **同一原初扰动 $\mathcal{R}(k)$ 在两点间的相干**。

V6.1 的 W 投影在 IR 引入 **空间白噪声** $\xi(x,t)$：
$$
\langle\xi(x,t)\xi(x',t')\rangle \propto \delta^3(x-x')
$$

其效应是：

1. **减去** CPT 对称的相干部分（§4.2 的 $\mathcal{S}_\ell$）；
2. **加上** 不相干的随机残差 $\delta T|_{\rm W}$（§2.3）。

对两点相关函数：
$$
C(\theta)
=
C_{\Lambda{\rm CDM}}(\theta)\,[1-\epsilon_{\rm corr}(\theta)]
\;+\;
C_{\rm noise}(\theta)
$$

其中：

**相干抑制项**：
$$
\epsilon_{\rm corr}(\theta)
\;\approx\;
f^2(R_\ast)\,
\Theta\!\left(\theta - \theta_{\rm IR}\right)\,
\mathcal{J}(\theta)
$$
$\mathcal{J}(\theta)\in(0,1)$ 为角分离积分核（由 $j_\ell(kD_A)$ 重叠决定）。

**噪声项**（白噪声 ⇒ 仅 $\theta=0$ 奇异；有限角分辨率下展宽）：
$$
C_{\rm noise}(\theta)
\;\approx\;
\frac{f^2(R_\ast)\,M_{\rm Pl}^2 H_\ast^2}{\rho_\gamma(t_\ast)^2}\,
\delta_{\rm beam}(\theta)
$$

在 **有限 beam** 下，$\delta_{\rm beam}$ 展宽为 $\theta \lesssim \theta_{\rm beam}$ 的小角贡献；  
对 **$\theta \sim \theta_{\rm IR}$**，主导仍是 $\epsilon_{\rm corr}>0$ ⇒ **相关低于标准**。

### 6.3 形式推导（两点函数）

设
$$
\frac{\delta T}{T}(\hat n)
=
\frac{1}{3}\Psi(\hat n)
+ \delta T_W(\hat n)
$$

$\Psi$ 来自相干 $\mathcal{R}$；$\delta T_W$ 来自 $\xi$，且 $\langle\delta T_W(\hat n_1)\delta T_W(\hat n_2)\rangle \propto \delta(\hat n_1-\hat n_2)$。

两点相关：
$$
\langle \delta T_1 \delta T_2\rangle
=
\underbrace{\frac{1}{9}\langle\Psi_1\Psi_2\rangle}_{\text{相干}}
\;+\;
\underbrace{\langle\delta T_{W,1}\,\delta T_{W,2}\rangle}_{\text{仅 }\theta=0}
$$

对 **ensemble 平均** 且 **减去 CPT 对称相干** 后，有效相干振幅：
$$
\langle\Psi_1\Psi_2\rangle_{\rm eff}
=
\langle\Psi_1\Psi_2\rangle^{(0)}\,[1-f^2(R_\ast)\,K(\theta)]
$$

⇒
$$
\boxed{
\frac{C(\theta)}{C_{\Lambda{\rm CDM}}(\theta)}
\;\approx\;
1 - f^2(R_\ast)\,K(\theta),
\qquad \theta \gtrsim \theta_{\rm IR}
}
$$

$K(\theta)$ 单调从 $K(0)=0$ 增至 $K(\theta_{\rm IR})\sim \mathcal{O}(1)$。

**这就是「相关性缺失」的机制**：  
不是扰动「不存在」，而是 **IR 上 CPT 边界条件把对称相干部分系统性扣除**，  
留下低于 $\Lambda$CDM 的两点相关。

### 6.4 与 $C_\ell$ 的关系

对 full-sky：
$$
C_\ell = 2\pi\int_{-1}^{1} d(\cos\theta)\,P_\ell(\cos\theta)\,C(\theta)
$$

若 $C(\theta)<C_{\Lambda{\rm CDM}}(\theta)$ 在 $\theta \gtrsim \theta_{\rm IR}$，  
则 **所有 $\ell \lesssim \ell_{\rm IR}$ 的 $C_\ell$ 一并压低**——  
与 §4 四极低谷、§5 轴排列 **同一 IR 边界层** 的三个面向。

---

## 七、与暗能量量级的定性相容性

### 7.1 共同控制参数

| 现象 | 框架表达式 | 控制参数 |
|------|-----------|----------|
| 今日暗能量 | $\rho_{\rm DE,0}\sim 4\pi f^2(R_0)M_{\rm Pl}^2 H_0^2$ | $f(R_0)$ |
| 冷斑深度 | $|\delta T/T|_{\rm cold}\propto f(R_\ast)M_{\rm Pl}H_\ast/\rho_\gamma$ | $f(R_\ast)$ |
| 四极抑制 | $C_2/C_{2,0}\approx 1-f^2(R_\ast)\mathcal{I}_2$ | $f^2(R_\ast)$ |
| 轴排列 | $\mathcal{A}_{\rm align}\propto f(R_\ast)|b_{0,\ast}|$ | $f(R_\ast)$ |
| 相关缺失 | $C(\theta)/C_0(\theta)\approx 1-f^2(R_\ast)K(\theta)$ | $f^2(R_\ast)$ |

**全部依赖同一 $f(R)=\varepsilon_0\tanh(R/R_\ast)$**，无第二套「大角度专用」耦合。

### 7.2 量纲与数量级链条

1. V6.1 已用 **今日** $f(R_0)\approx 0.09$ 锁定 $\rho_{\rm DE,0}\sim 10^{-120}M_{\rm Pl}^4$ 量级（通过 $M_{\rm Pl}^2 H_0^2$ 封顶，而非 $M_{\rm Pl}^4$ 抵消）。
2. LSS 时刻 $R_\ast > R_0$（早期曲率更高）⇒ $f(R_\ast) > f(R_0)$（在 $\tanh$ 上升支）。
3. 因此 **同一套小 $f$ 机制** 在 IR 造成 **$\mathcal{O}(f^2)$–$\mathcal{O}(f)$ 的 CMB 大角度修正**，  
   而 **不会** 在 UV（高 $\ell$）引入 $\mathcal{O}(1)$ 破坏——与 V6.1「低曲率 $f\to 0$、高 $\ell$ 接近标准」一致。

### 7.3 定性相容陈述（patch 後修訂）

> **量級相容**：若 $f(R_0)$ 由 $\rho_{\rm DE,0}\sim 4\pi f^2 M_{\rm Pl}^2 H_0^2$ 鎖定，則 V6.2 大角度修正 $\propto f^2(R_\ast)$ **同階**，联合可证伪 **在量級上** 成立。

> **份额张力（V6.1 附录 B3）**：纯代数 §6.2 给 $\Omega_{\rm DE,0}=4\pi f^2/3\approx 0.053$（$f^2\sim 0.013$，$R$ 取 §4.4  convention），**不能**同时解释 $\Omega_{m,0}\approx 0.3$。  
> **背景域限制（B10）**：同一 $f(R)$ 在 LSS 使 $4\pi f^2\gg 3$，§6.2 代数背景 **不适用于** 再组合；V6.2 大角度修正取 **扰动层** $f(R_\ast)$，与晚期暗能量 **逻辑 decouple** 直至 V7 统一。

### 7.4 与 $w(z)$ 非单调的关联（V6.1 延伸）

V6.1 预测 $w(z)$ 在 $z\sim 1.5$ 附近有峰。  
V6.2 补充：$f(R(z))$ 随宇宙演化 ⇒ **大角度 CMB 异常幅度**（如 $\mathcal{I}_2 f^2$）  
与 **晚期 $w(z)$ 形状** 通过 **同一 $R_\ast^{\rm late}$** 关联——  
改变 $R_\ast^{\rm late}$ 以拟合 $w(z)$ 时，**$C_2$ 抑制深度应同步变化**（数值留空）。

---

## 八、可证伪预测（数据暂空）

### 8.1 预测清单

| 编号 | 预测 | 结构式 | 数据 |
|------|------|--------|------|
| V6.2-1 | 四极 TT 相对抑制 | $C_2^{TT}/C_{2,0}^{TT} \approx 1-f^2(R_\ast)\mathcal{I}_2$ | （留空） |
| V6.2-2 | 低 $\ell$ 轴对齐 | $\mathcal{A}_{\rm align} \propto f(R_\ast)|b_{0,\ast}|$ | （留空） |
| V6.2-3 | 大角相关缺失 | $C(\theta)/C_0(\theta) \approx 1-f^2(R_\ast)K(\theta)$ for $\theta\gtrsim\theta_{\rm IR}$ | （留空） |
| V6.2-4 | 冷斑深度–尺度 | $|\delta T/T|_{\rm cold} \propto f(R_\ast)\mathcal{G}(R_p H_\ast)$ | （留空） |
| V6.2-5 | 联合约束 | $f(R_0)$ 由 $\rho_{\rm DE,0}$ 定 ⇒ $f(R_\ast)$ 定 ⇒ V6.2-1–4 无独立自由度 | （留空） |
| V6.2-6 | parity-odd | 低 $\ell$ $TB/EB$ 主轴 $\parallel$ 或 $\perp$ $\hat n_\phi$ | （留空） |

### 8.2 否定路径

1. 若精确测量显示 **$C_2$ 无抑制** 且 **大角相关无缺失** 且 **无冷斑**，而 $\rho_{\rm DE,0}$ 仍维持 V6.1 所需量级 ⇒ **V6.2 IR 边界层被否证**。
2. 若四极抑制存在但 **与 $f(R_0)$ 反推的 $f(R_\ast)$ 不一致** ⇒ $\tanh$ 单参数形状不足，需 V7 修正 $f(R)$。
3. 若轴排列方向与 $\phi$ 梯度预测 **无关联** ⇒ 各向异性项 $\mathcal{A}$ 形式需修正。

---

## 九、限制与待补计算

1. **$\mathcal{I}_\ell$、$K(\theta)$ 的精确数值**：需指定 $\mathcal{P}_\mathcal{R}^{(0)}(k)$ 与 V6.1 背景 $H(z)$、$D_A(z)$ 联立积分（留待数值脚本）。
2. **再组合后 ISW 修正**：本文 SW 近似对大角度足够作定性，完整计算需加入 $\dot\Phi$ 项。
3. **$\Delta G[f(R)]$ 对 transfer function 的修正**：预期 $\lesssim f^2$，但应并入 Boltzmann 级计算。
4. **非高斯性**：冷斑尾部分布 $\Rightarrow$ $f_{\rm NL}^{\rm loc}$ 可能升高，需与 V5 $f_{\rm NL}$ 预测联立。
5. **本稿仍为候选机制**：IR CPT 边界层尚未从 $\Delta G_{\mu\nu}^{[f(R)]}$ 严格定理化（V7-A 路线）。

---

## 十、小结

V6.2 在 **不增加基本参数** 的前提下，把 V6.1 的 W 投影（白洞虚对、时间尺度残差、视界封顶）  
推广到 **CMB IR 波段**，给出：

1. **冷斑**：局域 $\xi$ 涨落 ⇒ $\delta T<0$ 的 patch。  
2. **四极低谷**：CPT 对称分量在 $\ell=2$ 上最大抵消 ⇒ $C_2$ 按 $f^2(R_\ast)$ 压低。  
3. **轴排列**：$\phi$ 时间箭头 ⇒ $\mathcal{P}_\mathcal{R}$ 轴对称调制 ⇒ 低 $\ell$ 主轴对齐。  
4. **相关缺失**：相干项被 $\mathcal{S}_\ell$ 扣除 + 白噪声 $\xi$ ⇒ $C(\theta)$ 低于标准。

四者 **与暗能量量级共用 $f(R)$** ⇒ 形成 V6.1+V6.2 的 **联合可证伪** 结构。  
实际数据对比 **全部留空**，待后续数值实现与观测填入。

---

## 附录 A：与 V6.1 的章节对照

| V6.1 章节 | V6.2 延伸 |
|-----------|-----------|
| §2 白洞虚对 | §3 冷斑；§6 相关缺失的噪声源 |
| §3 Langevin $\xi$ | §3.2–3.3；§6.2–6.3 |
| §4 视界封顶 | §7 与 $\rho_{\rm DE}$ 相容 |
| §5.3 W 投影 | §2 IR 边界层 |
| §6 $w(z)$ | §7.4 联合约束 |
| §7 V6-1–7 | §8 V6.2-1–6 新增 |

---

## 附錄 B：V6.2 可擴展性與待確認問題

| # | 項目 | 可擴展？ | 說明 |
|---|------|----------|------|
| C1 | 冷斑 $\delta T/T|_{\rm cold}$ | ✓ | 接 $\xi$ 局域方差；需 LSS $f(R_\ast)$ |
| C2 | 四极 $\mathcal{I}_2$ | ✓ | 需 $\mathcal{P}_\mathcal{R}^{(0)}(k)$ 数值积分 |
| C3 | 轴排列 $\mathcal{A}_{\rm align}$ | ✓ | 需 $\dot\phi_\ast$ 与 $\hat n_\phi$ |
| C4 | $K(\theta)$ 相关缺失 | ✓ | 需 full-sky $C(\theta)$ 模板 |
| C5 | 与 $\Omega_{m,0}$ 联合 | ✗（暂） | 依赖 V6.1 B3 / $\tilde Q$ 闭合 |
| C6 | Planck 数据对比 | ✗（暂） | 数据留空；待 C2–C4 算完后填入 |
| C7 | $4\pi f^2(R_\ast)<3$ at LSS | ✗（背景） | **已确认**：LSS $f^2\to 1$，$4\pi f^2\approx 12.6\gg 3$；§6.2 **背景** 在 LSS 奇异。**扰动**仍取 $f(R_\ast)\approx\varepsilon_0$ 作 IR 耦合幅度（与背景 decouple，V7 待统一） |
| C8 | 代数背景自洽域 $z\lesssim 1$ | ⚠ | 与 V6-2、BAO 高 $z$ 检验冲突；联合拟合须限 $z<1.2$ 或修正 $f(R)$ |

**结论**：V6.2 **结构上可扩展**至数值 CMB pipeline（C1–C4）；**物理上**依赖 V6.1-patch 的 $f(R)$ 链。**背景** §6.2 在 LSS 不可用（C7/B10），**扰动** IR 修正仍可用 $f(R_\ast)\to\varepsilon_0$。**不能**在 B3 未解前声称已闭合宇宙学参数。

---

*文件版本：V6.2 Supplement（2026-05-09，對齊 V6.1-patch）*  
*作者推论性笔记，仅供讨论，非经同行评审。*  
*完全为自行推演，不引用既有文献；观测数值暂空。*
