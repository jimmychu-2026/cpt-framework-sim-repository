# 曲率相依 CPT 破壞框架 V4.1：最小自洽修補版

> 本文件為 V4 的最小修補版（V4.1）。
> 目標不是新增現象，而是完成三件事：**量綱閉合、參數閉合、動力學閉合**。
> 文件版本：2026-05-07（V4.1）

---

## 目錄

1. 核心定位與 V4.1 修補目標  
2. 協變作用量（量綱修補）  
3. 場方程與守恆條件（條件式聲明）  
4. B/C 雙機制的最小一致定義  
5. 觀測預測的語氣修正  
6. 數值章節（ODE 閉合版）  
7. 參數表與可證偽量  
8. 已知限制與下一步  

---

## 一、核心定位與 V4.1 修補目標

V4.1 保留 V4 主軸：

- **機制 B**：對稱對湮滅/對創的不對稱率  
- **機制 C**：可見部門與 $\phi$ 部門的能量交換

但修補三個核心問題：

1. **量綱問題**：互動算符加入抑制尺度 $M_\phi$。  
2. **參數問題**：明確區分「基本參數」與「導出/擬合量」。  
3. **動力學問題**：互動暗能量交換項 $Q$ 改為可積分的封閉形式。  

---

## 二、協變作用量（量綱修補）

### 2.1 完整作用量

$$
S = S_{\rm grav} + S_\phi + S_\psi + S_{\rm int}
$$

$$
S_{\rm grav}=\frac{1}{16\pi G}\int d^4x\sqrt{-g}\,R
$$

$$
S_\phi=\int d^4x\sqrt{-g}\left[-\frac12 g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi - V(\phi)\right]
$$

$$
S_\psi=\int d^4x\sqrt{-g}\,\bar\psi(i\gamma^\mu\nabla_\mu-m)\psi
$$

$$
\boxed{
S_{\rm int}
=
\int d^4x\sqrt{-g}\;
\frac{c_5\,f(R)}{M_\phi}\,
\partial_\mu\phi\,
\bar\psi\gamma^5\gamma^\mu\psi
}
$$

其中：

- $c_5$：無因次耦合常數  
- $M_\phi$：有效場論抑制尺度（量綱修補關鍵）  
- $f(R)$：無因次曲率形狀函數

### 2.2 唯一耦合函數

$$
f(R)=\varepsilon_0\tanh\!\left(\frac{R}{R_\ast}\right)
$$

說明：

- $f(R)$ 僅負責曲率依賴形狀與飽和，不負責量綱。  
- 量綱由 $1/M_\phi$ 單獨承擔。  

---

## 三、場方程與守恆條件（條件式聲明）

### 3.1 修正 Dirac 方程

$$
\left(
i\gamma^\mu\nabla_\mu - m
+ \frac{c_5 f(R)}{M_\phi}\gamma^5\gamma^\mu\partial_\mu\phi
\right)\psi=0
$$

### 3.2 $\phi$ 場方程

$$
\Box\phi - V'(\phi)
=
\nabla_\mu\!\left[
\frac{c_5 f(R)}{M_\phi}\bar\psi\gamma^5\gamma^\mu\psi
\right]
$$

### 3.3 修正 Einstein 方程（形式）

$$
G_{\mu\nu}+\Delta G_{\mu\nu}^{[f(R)]}
=
8\pi G\,
\left(T_{\mu\nu}^{(\phi)}+T_{\mu\nu}^{(\psi)}+T_{\mu\nu}^{(\rm int)}\right)
$$

### 3.4 守恆聲明（V4.1 修正語氣）

本框架採用條件式表述：

- 在完整變分包含 $\Delta G_{\mu\nu}^{[f(R)]}$，且與 Bianchi 恆等式閉合時，總能動量守恆成立。  
- 目前文件提供的是可運算架構，完整張量展開列為待補附錄。  

---

## 四、B/C 雙機制的最小一致定義

### 4.1 機制 B（對湮滅/對創）

- 典型形式：$X+\bar X\leftrightarrow Y+\bar Y$  
- 粒子/反粒子能譜分裂：
$$
E_\pm(\vec p)=\sqrt{p^2+m^2}\mp \frac{c_5 f(R)}{M_\phi}\dot\phi\,s(\vec p)
$$
- 導致凍結時殘餘不對稱：
$$
\eta_B \sim \frac{c_5 f(R_{\rm freeze})}{M_\phi}\frac{\dot\phi_{\rm freeze}}{T_{\rm freeze}}
$$

### 4.2 機制 C（可見部門與 $\phi$ 部門能量交換）

V4.1 明確修正語意：

- 不再以「粒子數不守恆」作主敘述  
- 改為「**可見部門與 $\phi$ 部門之間的能量交換**」

連續性方程：

$$
\dot\rho_m+3H\rho_m=+Q,\qquad
\dot\rho_{\rm DE}+3H(1+w_{\rm DE})\rho_{\rm DE}=-Q
$$

採最小封閉參數化：

$$
\boxed{
Q=
\xi_c\,
\frac{c_5 f(R)}{M_\phi}\,
H\rho_m
}
$$

其中 $\xi_c$ 為封裝微觀截面與相空間平均的有效係數。

---

## 五、觀測預測的語氣修正（避免過度宣稱）

V4.1 對 parity-odd 與 GW 敘述改為：

- 原初 GW 手徵不對稱、CMB 的 $TB/EB$，在本框架下屬「可檢驗候選簽名」。  
- 是否為主導來源，需由完整張量微擾方程與轉移函數計算確認。  
- 因此目前不再宣稱「微觀根源已完全建立」，而是「提出可驗證機制」。  

---

## 六、數值章節（ODE 閉合版）

### 6.1 背景方程組

使用聯立 ODE：

$$
E^2(a)=\frac{H^2}{H_0^2}
=
\Omega_{m0}a^{-3}
+\Omega_{r0}a^{-4}
+\Omega_{\rm DE}(a)
$$

$$
\frac{d\rho_m}{dt}+3H\rho_m=+Q,\qquad
\frac{d\rho_{\rm DE}}{dt}+3H(1+w_{\rm DE})\rho_{\rm DE}=-Q
$$

$$
Q=\xi_c\frac{c_5 f(R)}{M_\phi}H\rho_m,\qquad
R=6(\dot H+2H^2)
$$

### 6.2 與舊版表格的關係

- V1/V4 的 $A,B,\beta$ 表格可保留作現象學對照。  
- 但在 V4.1 中它們明確標為「導出/擬合量」，不再列為基本參數。  
- 正式預測以 ODE 聯立解為主。  

---

## 七、參數表與可證偽量

### 7.1 基本參數（V4.1）

推薦最小集合：

$$
\{\varepsilon_0,\;R_\ast^{\rm late},\;R_\ast^{\rm infl},\;M_\phi,\;\xi_c\}
$$

其中：

- 前三者控制曲率形狀與轉折尺度  
- $M_\phi$ 控制高維算符抑制  
- $\xi_c$ 控制交換強度  

### 7.2 導出/擬合量

- $A, B, \beta$：背景擬合用參數，非基本微觀參數  
- $\Omega_{m0}, H_0$：宇宙學初始條件/觀測輸入  

### 7.3 真正可證偽關聯

1. $w(z)$ 非單調峰位與 $z_{\rm acc}$ 偏移的聯動  
2. 結構增長偏差 $\Delta f_g(z)$ 與峰深的比例  
3. parity-odd（若檢出）與 $\eta_B$ 的比例關聯  

---

## 八、已知限制與下一步

### 8.1 仍未解決

1. 暗能量量級問題（為何今日如此小）  
2. $f(R)$ 形狀仍屬現象學選擇  
3. $n_s,r$ 尚未由同一作用量完成精算  
4. 高曲率區仍缺完整量子重力一致性  

### 8.2 下一步計算清單

- 完整推導 $\Delta G_{\mu\nu}^{[f(R)]}$ 並檢查微擾穩定性  
- 聯立 ODE 跑全資料集（Planck+BAO+SNIa+WL）  
- 計算張量手徵譜並投影到 $TB/EB$ 可觀測量  
- 建立 $\eta_B$ 與 CMB/GW 之聯合似然  

---

## 一句話總結

> V4.1 保留「B/C 雙機制」敘事，但把最關鍵的理論漏洞先補到可計算層級：  
> **算符有量綱、參數有邊界、演化有閉環。**

---

*文件版本：V4.1（2026-05-07）*  
*作者推論性筆記，僅供討論，非經同行評審。*  
*完全為自行推演，不引用既有文獻。*  
