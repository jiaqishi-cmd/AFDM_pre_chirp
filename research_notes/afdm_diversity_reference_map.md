# AFDM分集增益相关研究地图

## 检索目标

本记录围绕当前课题中的三个问题整理文献：

1. AFDM满分集应当用什么数学对象证明；
2. 逐子载波或动态预啁啾参数改变后，满分集条件是否仍成立；
3. ML分集保证、MMSE实际性能和PAPR优化之间如何建立联系。

检索结论截至2026年9月8日。未检出不等于不存在，专利提交前仍需进行专利数据库和引文追踪检索。

## 第一优先级文献

### 1 Bemani等 AFDM基础论文

- 题目：Affine Frequency Division Multiplexing for Next Generation Wireless Communications
- 期刊：IEEE Transactions on Wireless Communications, 2023, 22(11), 8214-8229
- DOI：10.1109/TWC.2023.3260906
- 链接：https://arxiv.org/abs/2204.12798

可直接借鉴：

- 将接收模型写成y=Phi(x)h+w；
- 对错误事件delta=x-x_hat构造Phi(delta)；
- PEP高SNR指数由rank(Phi(delta))决定；
- 系统分集阶数是所有非零delta上的最小秩；
- 非零奇异值乘积决定满秩后的编码增益。

限制：证明针对ML检测和公共标量c2。它不能直接解释线性MMSE误码平台，也没有分析数据相关的分组c2模式。

### 2 Liu等 AFDM-PIM满分集论文

- 题目：Pre-Chirp-Domain Index Modulation for Full-Diversity Affine Frequency Division Multiplexing Toward 6G
- 期刊：IEEE Transactions on Wireless Communications, 2025, 24(9), 7331-7345
- DOI：10.1109/TWC.2025.3559997
- 链接：https://arxiv.org/abs/2410.00313

可直接借鉴：

- 允许不同子载波采用不同预啁啾参数；
- 使用Phi_hat(x_hat)-Phi(x)处理数据和预啁啾模式同时变化的错误事件；
- 通过满秩分析和PEP推导证明分集；
- 将预啁啾字母表设计为最大化不同模式间的最小欧氏距离；
- 给出“路径数量条件+预啁啾参数取无理数”的满分集条件。

对当前课题的影响：不能再把“逐子载波无理数c2保持满分集”作为新结论。我们的理想边信息模型中模式q固定，错误矩阵比AFDM-PIM更简单，但必须与其无理数证明和MED优化正面对照。

### 3 Li等 MMSE啁啾参数选择论文

- 题目：Chirp Parameter Selection for Affine Frequency Division Multiplexing With MMSE Equalization
- 期刊：IEEE Transactions on Communications, 2025, 73(7), 5079-5093
- DOI：10.1109/TCOMM.2024.3519521
- 链接：https://doi.org/10.1109/TCOMM.2024.3519521

可直接借鉴：

- 明确指出ML分集分析不能直接外推到线性MMSE；
- 推导AFDM在MMSE均衡下的理论BER下界；
- 使用各啁啾子载波SINR是否均衡作为参数优劣判据；
- 分析啁啾参数与MMSE性能损失之间的关系；
- 根据频率选择性和双选择性信道分别设计参数选择策略。

对当前课题的影响：交底书图4的固定信道平台应当对照该文的SINR和BER下界分析。仅用rank(Phi)不能解释MMSE结果。

### 4 Lu等 A2FDM论文

- 题目：Augmented Affine Frequency Division Multiplexing for Both Low PAPR Signaling and Diversity Gain Protection
- 期刊：IEEE Access, 2026, 14, 69426-69442
- DOI：10.1109/ACCESS.2026.3690709
- 链接：https://arxiv.org/abs/2512.17679

可直接借鉴：

- 显式考虑模N以后不同路径发生重叠的坏参数条件；
- 给出坏c1值c1=(dN+nu_i-nu_j)/(2N(l_j-l_i))；
- 先列出危险参数集合，再设计结构保证最低分集阶数；
- 用子块DFT扩频使可保证分集达到min(N_mu,L)；
- 同时讨论PAPR、复杂度、数据率和分集增益。

对当前课题的影响：宽泛的“低PAPR与分集保护”问题已被占据。我们应学习其“坏参数闭式条件+最低保证阶数”的论证方式，同时强调我们的对象是c2分组模式和有限星座循环，而不是错误c1导致的路径位置重叠。

## 第二优先级文献

### 5 Yin Coded AFDM误差性能

- 题目：Error Performance of Coded AFDM Systems in Doubly Selective Channels
- arXiv：https://arxiv.org/abs/2311.15595
- 正式期刊与DOI：[NOT FOUND]

可直接借鉴：定义Gram矩阵Phi(e)^H Phi(e)，用其非零特征值几何平均描述条件编码增益，并讨论分集阶数与编码增益的trade-off。

用途：我们的安全指标不能只有rank。应同时报告最小特征值、行列式或非零特征值乘积，否则“仍然满秩但非常接近退化”的候选无法被识别。

### 6 Yi等 ZP-AFDM误差率和接收机

- 题目：Error Rate Analysis and Low-Complexity Receiver Design for Zero-Padded AFDM
- 期刊：IEEE Transactions on Vehicular Technology, 2026, in press
- DOI：10.1109/TVT.2026.3689279
- 链接：https://arxiv.org/abs/2510.14507

可直接借鉴：分别推导ML的PEP/union bound和MMSE的逐子载波SINR近似BER。该结构适合我们建立“分集秩结论”和“线性接收平台结论”两条平行证据链。

### 7 Mao等 联合参数与低复杂度MMSE

- 题目：Joint Chirp Parameter Selection and Low-Complexity MMSE Receiver Design for AFDM Systems
- arXiv：https://arxiv.org/abs/2607.20227
- 状态：2026年7月预印本

可直接借鉴：用DAFT的循环-对角结构构造简化BER指标，联合优化c1和MMSE矩阵稀疏化，并在复杂度约束下分层搜索。

对当前课题的影响：简单地增加“基于MMSE BER选择啁啾参数”或“层次搜索”已经不安全。我们的区别必须来自GPS相位循环、有限星座兼容条件或有限精度鲁棒性。

### 8 AFDM-IM与信道估计误差

- 题目：Affine Frequency Division Multiplexing With Index Modulation
- 会议：IEEE WCNC 2024
- DOI：10.1109/WCNC57260.2024.10570960

可直接借鉴：在信道估计误差下推导渐近紧的BER上界，并比较满足和不满足满分集条件时不同类型比特的保护程度。

用途：核心理论完成后，可用于增加不完美CSI扩展；当前不应优先于理想CSI下的结构问题。

## 对GPS基线的关键修正

GPS论文先由最大候选去相关得到

c2_m = plus_or_minus 1/(4m^2),

随后为满足无理数要求，又给出含pi的替代形式。但其数值实验部分仍明确使用plus_or_minus 1/(4m^2)。当前MATLAB工程复现的也是该有理数简化形式。

因此必须区分至少三种GPS实现：

1. GPS-rational：论文仿真使用的plus_or_minus 1/(4m^2)；
2. GPS-irrational：论文公式给出的含pi无理数替代；
3. GPS-quantized：将无理数候选量化到有限比特后的实际实现。

当前发现的严格BPSK循环秩退化，已经确认适用于GPS-rational。它是否适用于GPS-irrational和不同精度GPS-quantized，尚未验证。

## 当前未检出的方向

[NOT FOUND] 专门研究AFDM预啁啾参数有限精度、无理数理论保证与实际分集/MMSE性能之间差异的论文。

这一点可以作为后续检索和实验方向，但在完成IEEE、Google Patents、WIPO和CNIPA的专利检索前，不能声称为文献空白。

## 对后续推导和仿真的直接调整

1. 沿用Bemani的Phi(delta)和PEP框架，不另造分集定义。
2. 对照AFDM-PIM，分析公共已知模式q下的简化错误矩阵，并检查无理数条件是否足够。
3. 用coded AFDM的特征值乘积补充rank指标。
4. 用Li和Yi的MMSE SINR/BER框架解释图4平台。
5. 学习A2FDM的方式，推导GPS循环危险集合和最低可保证安全指标。
6. 新增GPS-rational、GPS-irrational和GPS-quantized三套基线。
7. 联合报告PAPR、rank、最小奇异值、编码增益、MMSE SINR离散度和BER。

## 文献检索后的首轮验证结果

GPS无理修正版在当前双路径支撑扫描中消除了严格秩退化，这与AFDM-PIM关于无理预啁啾满分集的结论方向一致。但是其k=2参数距离有理退化相位仅7.967e-4 rad，最近有限星座错误事件的最小奇异值仅5.634e-4。在28 dB固定危险信道中，其MMSE BER仍与有理GPS基本相同。

此外，GPS公式中的k越大，相位越接近原始plus_or_minus pi/2结构。k从1增加到8时，最近错误事件的最小奇异值从1.49e-2下降到1.27e-9，而PAPR基本不变。这表明“取无理数即可保证满分集”只解决高SNR渐近秩问题，不能保证有限SNR编码增益或线性MMSE可靠性。

k=1并不能在常规SNR下直接消除问题：其零噪声判决裕量约为4.44e-4，在28 dB危险双径下BER仍约为3.83e-3，到80 dB才在当前样本中未观察到错误。k=2的判决裕量仅6.35e-7，在80 dB下仍表现出接近有理GPS的误码。这进一步说明需要使用有限SNR编码增益和MMSE判决裕量，而不是只判断参数是否为无理数。

因此后续应重点借鉴coded AFDM的非零特征值乘积和Li等人的MMSE SINR指标，而不是仅用AFDM-PIM的满秩条件作为安全判据。

## 建议阅读顺序

1. Bemani 2023：建立标准PEP和分集语言。
2. AFDM-PIM 2025：学习逐子载波c2的满秩证明和MED设计。
3. Li 2025：学习MMSE性能为什么不同于ML分集。
4. A2FDM 2026：学习坏参数闭式条件和最低分集保证。
5. Coded AFDM：补充编码增益指标。
6. Yi 2026与Mao 2026：补充MMSE公式和复杂度约束。
