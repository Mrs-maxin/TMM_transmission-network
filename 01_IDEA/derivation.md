# 含PV节点输电网的TMM交流概率潮流计算方法

> 文件说明（不进入论文）：本文件按中国电机工程学报正文体例整理，只保留可进入小论文的方法内容，版本2026-10-21。符号与代码变量对应见`symbols.md`，数值证据见`evidence_map.md`中E-PF-PV-20261021、E-QLIM-IEEE14-20261021、E-ACC-L2-20261021。重写前的推导及风险评估相关公式见`derivation_archive_20261021.md`。文中【文献待核查】处需补充已发表的配电网GMM/TMM概率潮流前作、分块TMM前作及MATPOWER文献。

## 1 含PV节点输电网的交流潮流模型

### 1.1 节点类型与直角坐标潮流方程

设输电网共有$n$个节点，节点集合为$\mathcal N$，平衡节点为$s$，PV节点集合为$\mathcal N_{\rm PV}$，PQ节点集合为$\mathcal N_{\rm PQ}$，非平衡节点集合$\mathcal N^{*}=\mathcal N_{\rm PV}\cup\mathcal N_{\rm PQ}$。节点$i$的电压相量记为$\dot V_i=e_i+{\rm j}f_i$，节点导纳矩阵元素记为$Y_{ij}=G_{ij}+{\rm j}B_{ij}$。直角坐标下节点注入有功、无功及电压幅值平方为

$$
P_i=e_i\sum_{j\in\mathcal N}\left(G_{ij}e_j-B_{ij}f_j\right)+f_i\sum_{j\in\mathcal N}\left(G_{ij}f_j+B_{ij}e_j\right)
\tag{1}
$$

$$
Q_i=f_i\sum_{j\in\mathcal N}\left(G_{ij}e_j-B_{ij}f_j\right)-e_i\sum_{j\in\mathcal N}\left(G_{ij}f_j+B_{ij}e_j\right)
\tag{2}
$$

$$
U_i^{2}=e_i^{2}+f_i^{2}
\tag{3}
$$

平衡节点电压幅值与相角给定，承担系统有功、无功不平衡量及网损。PV节点给定有功注入和机端电压幅值，PQ节点给定有功、无功注入，即

$$
P_i^{\rm sp}=P_{{\rm G},i}-P_{{\rm L},i},\ i\in\mathcal N^{*};\qquad
Q_i^{\rm sp}=Q_{{\rm G},i}-Q_{{\rm L},i},\ i\in\mathcal N_{\rm PQ};\qquad
U_i^{\rm sp}=U_i^{\rm set},\ i\in\mathcal N_{\rm PV}
\tag{4}
$$

式中：$P_{{\rm G},i}$、$Q_{{\rm G},i}$为机组计划有功和无功出力，节点无机组时取零，机组无功到达限值转为PQ节点后$Q_{{\rm G},i}$取该限值；$P_{{\rm L},i}$、$Q_{{\rm L},i}$为负荷有功和无功；$U_i^{\rm set}$为机端电压设定值。

以全部非平衡节点电压的实部和虚部为未知量$\boldsymbol x=[\boldsymbol e_{\mathcal N^{*}}^{\rm T},\boldsymbol f_{\mathcal N^{*}}^{\rm T}]^{\rm T}\in\mathbb R^{2(n-1)}$，潮流方程写为

$$
\boldsymbol F(\boldsymbol x)=
\begin{bmatrix}\boldsymbol P_{\mathcal N^{*}}(\boldsymbol x)\\ \boldsymbol Q_{\mathcal N_{\rm PQ}}(\boldsymbol x)\\ \boldsymbol U^{2}_{\mathcal N_{\rm PV}}(\boldsymbol x)\end{bmatrix}
=\boldsymbol y=
\begin{bmatrix}\boldsymbol P^{\rm sp}_{\mathcal N^{*}}\\ \boldsymbol Q^{\rm sp}_{\mathcal N_{\rm PQ}}\\ (\boldsymbol U^{\rm set}_{\mathcal N_{\rm PV}})^{2}\end{bmatrix}
\tag{5}
$$

方程个数为$(n-1)+n_{\rm PQ}+n_{\rm PV}=2(n-1)$，与未知量个数相等。与配电网中非平衡节点均为PQ节点的情形相比，PV节点以电压幅值平方方程(3)替代无功平衡方程。电压求得后，平衡节点和仍在调压的PV节点的机组无功出力为

$$
Q_{{\rm G},i}=Q_i(\boldsymbol x)+Q_{{\rm L},i},\qquad i\in\mathcal N_{\rm PV}\cup\{s\}
\tag{6}
$$

### 1.2 潮流方程的二次结构

平衡节点电压给定后，式(5)中各方程均为$\boldsymbol x$的二次函数。在任一运行点$\boldsymbol x_0$处有

$$
\boldsymbol F(\boldsymbol x_0+\Delta\boldsymbol x)=\boldsymbol F(\boldsymbol x_0)+\boldsymbol J_0\Delta\boldsymbol x+\boldsymbol h(\Delta\boldsymbol x)
\tag{7}
$$

式中：$\boldsymbol J_0=\partial\boldsymbol F/\partial\boldsymbol x|_{\boldsymbol x_0}$为雅可比矩阵，其有功、无功行由式(1)、(2)对$e_j$、$f_j$求导得到，PV节点行为$\partial U_i^{2}/\partial e_i=2e_i$、$\partial U_i^{2}/\partial f_i=2f_i$，其余元素为零；$\boldsymbol h(\Delta\boldsymbol x)$为$\Delta\boldsymbol x$的二次齐次项，式(7)不含更高阶项，且$\boldsymbol h(\Delta\boldsymbol x)=[\boldsymbol F(\boldsymbol x_0+\Delta\boldsymbol x)+\boldsymbol F(\boldsymbol x_0-\Delta\boldsymbol x)]/2-\boldsymbol F(\boldsymbol x_0)$。

式(7)表明，计入PV节点后潮流方程仍保持与配电网相同的二次结构，已发表的基于展开点的批量潮流计算方法【文献待核查】可直接用于含PV节点的输电网。

## 2 源荷场景下的批量潮流计算

### 2.1 展开点与二阶初值

以牛顿-拉夫逊法求得展开点潮流解$\boldsymbol x_0$，其给定量为$\boldsymbol y_0=\boldsymbol F(\boldsymbol x_0)$。对源荷场景$\omega$的给定量$\boldsymbol y_\omega$，由式(7)取电压增量的二阶近似作为初值：

$$
\Delta\boldsymbol x_\omega^{(0)}=\boldsymbol J_0^{-1}\Delta\boldsymbol y_\omega-\boldsymbol J_0^{-1}\boldsymbol h\!\left(\boldsymbol J_0^{-1}\Delta\boldsymbol y_\omega\right),\qquad
\Delta\boldsymbol y_\omega=\boldsymbol y_\omega-\boldsymbol y_0
\tag{8}
$$

全部场景共用同一$\boldsymbol J_0^{-1}$，式(8)只需两次潮流方程求值和矩阵乘法，可对场景矩阵批量计算。

### 2.2 固定雅可比修正与收敛判据

以式(8)为初值，保持$\boldsymbol J_0$不变进行修正：

$$
\boldsymbol x_\omega^{(m+1)}=\boldsymbol x_\omega^{(m)}+\alpha_m\boldsymbol J_0^{-1}\left[\boldsymbol y_\omega-\boldsymbol F(\boldsymbol x_\omega^{(m)})\right]
\tag{9}
$$

式中：$\alpha_m$在$1,1/2,\ldots,1/2^{8}$中取使残差2-范数下降的最大值，每个展开点最多修正10次。场景$\omega$的收敛判据为

$$
\|\boldsymbol r_\omega\|_\infty\le\varepsilon_y,\qquad
\|\boldsymbol J_0^{-1}\|_2\,\|\boldsymbol r_\omega\|_2\le\varepsilon_x,\qquad
\boldsymbol r_\omega=\boldsymbol y_\omega-\boldsymbol F(\boldsymbol x_\omega)
\tag{10}
$$

式中：$\varepsilon_y=10^{-8}$ p.u.，$\varepsilon_x=10^{-7}$ p.u.；同时要求节点电压幅值位于0.5~1.5 p.u.的数值求解范围内，该范围不是运行限值。

首个展开点取给定量最接近全部场景均值的场景。若仍有场景不满足式(10)，取残差估计$\|\boldsymbol J_0^{-1}\|_2\|\boldsymbol r_\omega\|_2$最大的未收敛场景作为下一个展开点，重复式(8)、(9)，展开点不超过12个；其余未收敛场景逐点用牛顿-拉夫逊法求解。由此得到的全部场景电压均满足式(5)至$\varepsilon_y$，与逐点牛顿-拉夫逊法求解的是同一组交流潮流方程。

### 2.3 机组无功限值与运行方式

运行方式由各机组状态$\sigma_g\in\{-1,0,1\}$表示：$\sigma_g=0$表示机组调压，$\sigma_g=1$、$-1$分别表示无功固定于上限$Q_{{\rm G},g}^{\max}$和下限$Q_{{\rm G},g}^{\min}$。固定于限值的机组所在节点转为PQ节点，其电压幅值不再给定。对每个场景按以下规则计算：

1）全部PV机组调压，按2.1、2.2节求解；

2）按式(6)计算仍调压机组的无功出力，若$Q_{{\rm G},g}>Q_{{\rm G},g}^{\max}+\varepsilon_Q$则置$\sigma_g=1$，若$Q_{{\rm G},g}<Q_{{\rm G},g}^{\min}-\varepsilon_Q$则置$\sigma_g=-1$，同一次求解中的全部越限机组同时转换；

3）若出现新的越限，按新的运行方式重新求解；已转换的机组不再恢复调压，直至不再出现新的越限。

式中：$\varepsilon_Q=5\times10^{-6}$ Mvar。平衡节点不设无功限值。该规则与MATPOWER牛顿-拉夫逊潮流计入机组无功限值时的处理一致【文献待核查】。运行方式相同的场景分为一组批量计算，每种运行方式按2.1、2.2节建立各自的展开点。

## 3 节点电压的联合概率模型

### 3.1 电压坐标与标准化

全部非平衡节点（含PV节点）电压的实部和虚部构成电压坐标$\boldsymbol x$。由用于建立模型的源荷场景（下称建模场景）的电压样本计算均值$\bar{\boldsymbol x}$及协方差矩阵的特征分解，标准化坐标为

$$
\boldsymbol z=\boldsymbol D^{-1}\boldsymbol E_d^{\rm T}(\boldsymbol x-\bar{\boldsymbol x}),\qquad
\boldsymbol x=\bar{\boldsymbol x}+\boldsymbol E_d\boldsymbol D\boldsymbol z
\tag{11}
$$

式中：$\boldsymbol E_d$为前$d$个特征向量组成的矩阵；$\boldsymbol D={\rm diag}(\sigma_1,\ldots,\sigma_d)$为对应方向的电压标准差。$d$取累计方差损失不超过$10^{-10}$且数值可分辨的方向数，PV节点电压幅值固定等确定性关系对应的方向不保留。

### 3.2 多元t混合模型

标准化电压的概率密度采用多元t混合模型（TMM）：

$$
f(\boldsymbol z)=\sum_{k=1}^{K}\pi_k\,t_{\nu_k}(\boldsymbol z;\boldsymbol\mu_k,\boldsymbol\Sigma_k)
\tag{12}
$$

$$
t_{\nu}(\boldsymbol z;\boldsymbol\mu,\boldsymbol\Sigma)=
\frac{\Gamma\!\left(\frac{\nu+d}{2}\right)}{\Gamma\!\left(\frac{\nu}{2}\right)(\nu\pi)^{d/2}|\boldsymbol\Sigma|^{1/2}}
\left(1+\frac{\delta}{\nu}\right)^{-\frac{\nu+d}{2}},\qquad
\delta=(\boldsymbol z-\boldsymbol\mu)^{\rm T}\boldsymbol\Sigma^{-1}(\boldsymbol z-\boldsymbol\mu)
\tag{13}
$$

式中：$K$为分量数；$\pi_k$、$\boldsymbol\mu_k$、$\boldsymbol\Sigma_k$、$\nu_k$分别为第$k$个分量的权重、位置向量、散布矩阵和自由度，分量协方差为$\nu_k\boldsymbol\Sigma_k/(\nu_k-2)$。式(12)可等价表示为

$$
\boldsymbol z\mid C=k,S_k=s\sim\mathcal N(\boldsymbol\mu_k,s\boldsymbol\Sigma_k),\qquad
S_k\sim{\rm InvGamma}\!\left(\frac{\nu_k}{2},\frac{\nu_k}{2}\right)
\tag{14}
$$

式中：$C$为场景所属分量；$S_k$为第$k$个分量的共同波动强度，同一场景中所有电压方向按同一$\sqrt{S_k}$缩放。$\nu_k\to\infty$时$S_k\equiv1$，式(12)退化为高斯混合模型（GMM）。

### 3.3 共同波动强度作用方向与分块模型

式(14)中共同波动强度作用于全部电压方向。输电网中，只有由同一波动过程驱动的注入功率所引起的电压变化需要随该强度放大，机组调压约束及其他独立波动对应的电压方向不必共享该强度。为此采用分块TMM【文献待核查】，由潮流灵敏度确定强度的作用方向。

在运行点$\boldsymbol x^{(0)}$处，共用同一波动强度的注入功率所驱动的电压变化方向张成

$$
\boldsymbol H=\boldsymbol J_{0}^{-1}\boldsymbol L_{\mathcal W}
\tag{15}
$$

式中：$\boldsymbol J_{0}$为式(7)的雅可比矩阵，运行点取建模场景中出现最多的运行方式在平均负荷下的潮流解；$\boldsymbol L_{\mathcal W}$为0-1选择矩阵，选取共用同一波动强度的注入所在的方程行。本文算例中全部负荷共用同一波动强度，故选取全部非平衡负荷节点的有功行及PQ负荷节点的无功行；仍调压的PV节点上负荷无功由机组承担，不改变节点电压，不予选取。

将$\boldsymbol H$映射到标准化坐标并作奇异值分解：

$$
\tilde{\boldsymbol H}=\boldsymbol D^{-1}\boldsymbol E_d^{\rm T}\boldsymbol H=\boldsymbol U_H\boldsymbol\Lambda_H\boldsymbol R_H^{\rm T},\qquad
r=\#\{\lambda_{H,j}>10^{-10}\lambda_{H,1}\},\qquad
\boldsymbol Q_u=\boldsymbol U_H(:,1{:}r)
\tag{16}
$$

式中：$\lambda_{H,j}$为奇异值；$\boldsymbol Q_w$取$\boldsymbol Q_u$的正交补。记$\boldsymbol u=\boldsymbol Q_u^{\rm T}\boldsymbol z\in\mathbb R^{r}$，$\boldsymbol w=\boldsymbol Q_w^{\rm T}\boldsymbol z\in\mathbb R^{d-r}$，第$k$个分量为

$$
\boldsymbol u\mid C=k,S_k=s\sim\mathcal N(\boldsymbol m_k,s\boldsymbol\Sigma_k^{u}),\qquad
S_k\sim{\rm InvGamma}\!\left(\frac{\nu_k}{2},\frac{\nu_k}{2}\right)
\tag{17}
$$

$$
\boldsymbol w\mid\boldsymbol u,C=k\sim\mathcal N(\boldsymbol a_k+\boldsymbol G_k\boldsymbol u,\boldsymbol\Omega_k)
\tag{18}
$$

式中：$\boldsymbol m_k$、$\boldsymbol\Sigma_k^{u}$为强度作用方向上的位置向量与散布矩阵；$\boldsymbol a_k$、$\boldsymbol G_k$、$\boldsymbol\Omega_k$为其余方向的截距、对$\boldsymbol u$的回归系数矩阵和残差协方差矩阵。$r=d$时式(17)、(18)退化为式(14)；$\nu_k\to\infty$时退化为GMM。同一分量数下，分块TMM与标准TMM的参数个数相同。

### 3.4 参数估计与分量数确定

标准TMM采用期望最大化方法估计参数。第$k$个分量对场景$i$的后验概率记为$\gamma_{ik}$，强度权重为$\tau_{ik}=(\nu_k+d)/(\nu_k+\delta_{ik})$，$N_k=\sum_i\gamma_{ik}$，参数更新为

$$
\pi_k=\frac{N_k}{N},\qquad
\boldsymbol\mu_k=\frac{\sum_i\gamma_{ik}\tau_{ik}\boldsymbol z_i}{\sum_i\gamma_{ik}\tau_{ik}},\qquad
\boldsymbol\Sigma_k=\frac{1}{N_k}\sum_i\gamma_{ik}\tau_{ik}(\boldsymbol z_i-\boldsymbol\mu_k)(\boldsymbol z_i-\boldsymbol\mu_k)^{\rm T}
\tag{19}
$$

自由度由下式求根：

$$
\ln\frac{\nu_k}{2}-\psi\!\left(\frac{\nu_k}{2}\right)+1+\frac{1}{N_k}\sum_i\gamma_{ik}\left[\psi\!\left(\frac{\nu_k'+d}{2}\right)-\ln\frac{\nu_k'+\delta_{ik}}{2}-\tau_{ik}\right]=0
\tag{20}
$$

式中：$N$为建模场景数；$\psi$为双伽马函数；$\nu_k'$为上一次迭代的自由度。分块TMM采用期望条件最大化方法：$\boldsymbol u$块按式(19)、(20)更新，其中$d$换为$r$；$\boldsymbol w$块按加权最小二乘更新：

$$
[\boldsymbol a_k,\boldsymbol G_k]=\arg\min\sum_i\gamma_{ik}\|\boldsymbol w_i-\boldsymbol a_k-\boldsymbol G_k\boldsymbol u_i\|_2^{2},\qquad
\boldsymbol\Omega_k=\frac{1}{N_k}\sum_i\gamma_{ik}\boldsymbol\epsilon_{ik}\boldsymbol\epsilon_{ik}^{\rm T}
\tag{21}
$$

式中：$\boldsymbol\epsilon_{ik}=\boldsymbol w_i-\boldsymbol a_k-\boldsymbol G_k\boldsymbol u_i$；分块TMM的$\gamma_{ik}$由式(17)的t密度与式(18)的条件高斯密度之积计算。

计算设置如下：散布矩阵与残差协方差矩阵的特征值下限为$10^{-7}$；自由度取值范围为$[4.1,10^{5}]$；平均对数似然变化小于$10^{-6}$时停止迭代。初值按前三个电压主方向上的投影将建模场景等分为$K$组形成，取3组不同投影方向；先估计GMM，再以GMM结果为初值估计TMM，自由度初值取12和200；分块TMM以同一分量数的GMM为初值。估计结果需满足迭代收敛、对数似然不下降且各分量有效场景数$N_k>d+2$。

分量数由另一组独立源荷场景（下称定阶场景）的平均对数似然确定：

$$
K^{*}=\arg\max_{K}\frac{1}{N_{\rm v}}\sum_{v=1}^{N_{\rm v}}\ln\hat f_K(\boldsymbol z_v)
\tag{22}
$$

式中：$N_{\rm v}$为定阶场景数；$\hat f_K$为分量数为$K$时的估计密度。分量数确定不使用交流Monte Carlo结果。

## 4 节点电压与线路功率的概率分布

### 4.1 节点电压幅值与相角

由式(11)将$\boldsymbol z$还原为$\boldsymbol x$，节点$i$的电压幅值和相角为

$$
U_i=\sqrt{e_i^{2}+f_i^{2}},\qquad\theta_i=\arctan\!2(f_i,e_i)
\tag{23}
$$

相角以平衡节点电压相角为基准。

### 4.2 线路功率

计及线路充电电容与变压器变比、相移，支路$ij$首端电流和功率为

$$
\dot I_{ij}=y_{ij}^{\rm ff}\dot V_i+y_{ij}^{\rm ft}\dot V_j,\qquad
P_{ij}+{\rm j}Q_{ij}=S_{\rm B}\,\dot V_i\dot I_{ij}^{*}
\tag{24}
$$

式中：$y_{ij}^{\rm ff}$、$y_{ij}^{\rm ft}$为支路$\pi$型等值导纳参数；$S_{\rm B}$为功率基准值；末端功率同理。由式(24)，$P_{ij}$、$Q_{ij}$及节点注入功率均为$\boldsymbol x$的二次函数，记为$q(\boldsymbol x)=\boldsymbol x^{\rm T}\boldsymbol A\boldsymbol x+\boldsymbol b^{\rm T}\boldsymbol x+c$，其中常数项和一次项来自平衡节点电压。给定分量$k$时，

$$
{\rm E}\left[q\mid C=k\right]=q(\boldsymbol c_k)+{\rm tr}\!\left(\boldsymbol A\boldsymbol M\boldsymbol V_k\boldsymbol M^{\rm T}\right),\qquad
\boldsymbol c_k=\bar{\boldsymbol x}+\boldsymbol M\boldsymbol\mu_k,\quad\boldsymbol M=\boldsymbol E_d\boldsymbol D
\tag{25}
$$

式中：$\boldsymbol V_k$为分量$k$中$\boldsymbol z$的协方差矩阵，标准TMM为$\nu_k\boldsymbol\Sigma_k/(\nu_k-2)$；分块TMM中$\boldsymbol u$、$\boldsymbol w$两块的协方差分别为$\boldsymbol V_k^{u}=\nu_k\boldsymbol\Sigma_k^{u}/(\nu_k-2)$和$\boldsymbol G_k\boldsymbol V_k^{u}\boldsymbol G_k^{\rm T}+\boldsymbol\Omega_k$，互协方差为$\boldsymbol G_k\boldsymbol V_k^{u}$，按$\boldsymbol Q_u$、$\boldsymbol Q_w$变换回$\boldsymbol z$坐标。电气量均值为各分量均值按$\pi_k$加权之和。

### 4.3 概率分布的数值计算

电气量的累积分布函数按分量分层的随机化Sobol序列计算。第$k$个分量分配$N_k^{\rm s}=N^{\rm s}\pi_k$个积分点，每点取$d+1$维Sobol点$\boldsymbol\eta$，由前$d$维生成标准正态向量，由最后一维生成共同波动强度：

$$
\boldsymbol z=\boldsymbol\mu_k+\frac{\boldsymbol L_k\Phi^{-1}(\boldsymbol\eta_{1:d})}{\sqrt{\tau}},\qquad
\tau=\frac{2}{\nu_k}\Gamma^{-1}\!\left(\eta_{d+1};\frac{\nu_k}{2}\right)
\tag{26}
$$

式中：$\boldsymbol L_k$为$\boldsymbol\Sigma_k$的Cholesky因子；$\Phi^{-1}$为标准正态分布的逆函数；$\Gamma^{-1}(\cdot;a)$为正则化不完全伽马函数的逆函数。分块TMM中$\boldsymbol u$按式(26)生成，$\boldsymbol w=\boldsymbol a_k+\boldsymbol G_k\boldsymbol u+\boldsymbol L_k^{\Omega}\Phi^{-1}(\boldsymbol\eta_{r+1:d})$。每个积分点的权重为$\pi_k/N_k^{\rm s}$，由式(11)、(23)、(24)计算电压幅值、相角与线路功率，加权得到累积分布函数。本文取$N^{\rm s}=50000$。

## 5 计算流程

步骤1：按源荷不确定性模型生成建模场景与定阶场景。

步骤2：按2.1~2.3节对两组场景进行批量潮流计算，得到计及机组无功限值的节点电压场景。

步骤3：由建模场景电压按式(11)确定标准化坐标。

步骤4：按式(15)、(16)确定共同波动强度的作用方向。

步骤5：对各候选分量数按3.4节估计GMM、标准TMM和分块TMM参数，按式(22)确定分量数。

步骤6：按4.3节计算节点电压幅值、相角及线路有功、无功的概率分布，均值可按式(25)解析计算。

全过程只在各运行方式的展开点及少量未收敛场景处调用牛顿-拉夫逊潮流计算，不需要逐场景交流潮流计算。

## 6 算例源荷模型与精度评价指标

### 6.1 算例源荷模型

算例中各负荷节点有功、无功负荷按基准值乘以随机倍率：

$$
P_{{\rm L},i}=P_{{\rm L},i}^{0}\left[1+c\sqrt3\left(2F_{t,\nu}(h_i^{P})-1\right)\right],\qquad
h_i^{P}=\sqrt{\frac{\nu}{\chi^{2}_{\nu}}}\left(0.4\xi_{\rm g}+0.3\xi_{\rm a}+\sqrt{0.75}\,\xi_i^{P}\right)
\tag{27}
$$

$$
h_i^{Q}=\sqrt{\frac{\nu}{\chi^{2}_{\nu}}}\left[0.4\xi_{\rm g}+0.3\xi_{\rm a}+\sqrt{0.75}\left(0.6\xi_i^{P}+0.8\xi_i^{Q}\right)\right]
\tag{28}
$$

式中：$\xi_{\rm g}$、$\xi_{\rm a}$、$\xi_i^{P}$、$\xi_i^{Q}$为相互独立的标准正态变量，分别表示全系统、区域及节点有功、无功的波动；$\chi^{2}_{\nu}$为与之独立、自由度$\nu=12$的卡方变量，同一场景中全部负荷共用，构成共同波动强度；$F_{t,\nu}$为t分布的累积分布函数；$c=5\%$为负荷波动的变异系数。无功负荷倍率以$h_i^{Q}$按式(27)同理计算。每个负荷倍率在$1\pm\sqrt3c$范围内均匀分布，基准负荷为零的节点保持为零。

### 6.2 精度评价指标

以独立的交流Monte Carlo场景（下称检验场景）逐点牛顿-拉夫逊潮流结果为参照，对每个电气量计算

$$
e_\mu=|\hat\mu-\mu_{\rm MC}|,\qquad
e_\sigma=\left|\frac{\hat\sigma}{\sigma_{\rm MC}}-1\right|,\qquad
e_p=\left|\hat y_p-y_{{\rm MC},p}\right|,\qquad
D_{\rm KS}=\sup_{y}\left|\hat F(y)-F_{\rm MC}(y)\right|
\tag{29}
$$

式中：$\hat\mu$、$\hat\sigma$、$\hat y_p$、$\hat F$分别为计算所得均值、标准差、$p$分位数和累积分布函数，下标MC表示交流Monte Carlo结果；$p$取1%和99%，$e_p$同时以$\sigma_{\rm MC}$为单位报告。$e_\mu$反映计算结果的准确度，$e_\sigma$、$e_p$和$D_{\rm KS}$反映分布形状与尾部的精度。建模场景本身与检验场景之间的$D_{\rm KS}$作为有限场景数下的抽样误差水平。
