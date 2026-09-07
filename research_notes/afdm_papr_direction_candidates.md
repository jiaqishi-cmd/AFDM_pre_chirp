# AFDM PAPR 研究方向候选与导师讨论提纲

更新日期：2026-09-07

## 1. 当前结论

本阶段不应提前锁定“新的 `c2` 分组分配方法”和“接收端适配不同
`c2`”为两个创新点。更合理的做法是先把研究问题拆成四个彼此关联、
但可以独立验证的层次：

1. 候选信号怎样构造，才能在有限候选数下产生足够不同的峰值形态；
2. 怎样在组合空间内搜索，才能以较少计算接近枚举结果；
3. 怎样约束候选，使 AFDM 的通信属性、接收复杂度和边信息代价可控；
4. 怎样评价低 PAPR 带来的真实收益，而不只展示一条 CCDF 曲线。

当前最值得优先探索的是“相关性和 AFDM 结构约束驱动的候选码本/分组
设计”以及“具有 AFDM 特定下界或节点评价的预算受限树搜索”。这两条
可以连接成一个完整问题，也可以根据实验结果只保留其中一条。

## 2. 新近文献带来的范围修正

截至 2026 年 9 月，AFDM PAPR 已经不是只有 Yuan GPS 的窄小方向。
已出现下列路线：

- 分组预啁啾选择：Yuan GPS，逐组贪婪改变 `c2`；论文明确指出只比较了
  adjacent/comb 等分组，最优 grouping 仍待研究。
- 标量或码本式 chirp selection：CSM、Choi 的 per-slot scalar chirp
  offset 和 intrinsic side information。
- 连续预啁啾优化：UQP、multi-start projected gradient、subcarrier-wise
  continuous `c2` perturbation 和 spectral projected gradient。
- 变换/扩频与非线性处理：WHT/DCT/ZC/DFT spreading、carrier
  interferometry、DCT 加 companding、SLM 加 companding。
- 帧结构资源利用：pilot/guard interval 上的 ACE 或低幅符号填充、
  chirp-subcarrier reservation。
- 联合指标：PAPR-WISL、PAPR-OOBE、频谱掩模、功放非线性下 BER。
- 联合收发机学习：星座、bit labeling、precoder 和 detector 的联合优化。

因此，下列表述已经不适合作为主要创新：

- 使用 `c2` 自由度降低 AFDM PAPR；
- 构造若干 `c2` 候选并选择最小 PAPR；
- 单纯把 PSO、MCTS、遗传算法或梯度下降换到现有候选空间；
- 单纯声称不发送 side information 或让接收端检测候选索引；
- 单纯把 DCT、SLM、TR、ACE 或 companding 从 OFDM 搬到 AFDM。

## 3. OFDM 思想迁移地图

| OFDM 方法族 | 可以迁移到 AFDM 的对象 | 已知拥挤程度 | 对当前课题的作用 |
| --- | --- | --- | --- |
| SLM | `c2` 码本、相位码本、候选选择 | 很高 | 理论和 baseline |
| PTS | 子载波分组、部分波形叠加、相位组合 | 很高 | 当前模型最接近的母体 |
| Tree-PTS / sphere search | 分层决定每组候选、剪枝或固定预算搜索 | 高 | 可作为搜索骨架 |
| MCTS-PTS | UCT、rollout、搜索预算分配 | 高 | 可做实验候选，不能只靠名称主张创新 |
| Conditional expectation | 逐层固定离散变量、使用可计算代理指标 | 中 | 适合追求可解释的贪婪规则 |
| SDR / randomization / cross entropy | 离散相位组合的连续松弛和随机恢复 | 高 | 适合作为性能上界或算法 baseline |
| TR / TI | 保留 DAFT/chirp 子载波生成峰值抵消信号 | 中高 | 接收端简单，但有速率/功率代价 |
| ACE | 在判决安全区域移动 QAM 符号 | 中高 | 可与 `c2` 选择形成混合方法 |
| Coding / Golay | 限定码字或导频使 PMEPR 有上界 | 中 | 更适合前导和导频，不适合当前数据块主线 |
| Clipping / companding | 直接修改时域幅度 | 很高 | 必做传统 baseline，需看 BER/OOBE |
| DFT/DCT/WHT/ZC spreading | 预编码降低相干叠加 | 很高 | AFDM 已有工作，适合 baseline |
| PA-aware selection | 以 EVM、ACLR、输出相关性或 BER 代替纯 PAPR | 中 | 可增强工程意义和评价完整性 |

这里最关键的认识是：AFDM 和 OFDM 的数学相似性使大量算法可以迁移，
也使“首次应用某个通用算法”很难形成稳定创新。创新需要落在 AFDM 的
二次相位、`c2` 约束、DAFT 域接收或高移动信道行为上。

## 4. 候选方向 A：相关性与结构约束驱动的 `c2` 码本设计

### 核心问题

在给定候选数和边信息开销下，怎样构造一组 `c2` pattern，使候选 AFDM
波形之间尽量不相关，同时维持可接受的结构风险与接收代价？

对第 `q` 个候选，写成

```text
s_q = Lambda_c1^H F^H D_q x,
D_q = diag(exp(j 2 pi c2_q[m] m^2)).
```

由于公共的 post-chirp 和 Fourier 变换是酉变换，有

```text
s_a^H s_b = x^H D_a^H D_b x
            = sum_m |x[m]|^2 exp(j 2 pi (c2_b[m]-c2_a[m]) m^2).
```

对于 PSK，归一化相关性只由两个 `c2` pattern 的二次相位差决定，可以
离线计算；对于 QAM，可研究平均相关性或最坏功率权重。这给出了一个比
“局部扰动看起来更温和”更明确的理论入口。

### 可选优化形式

若工程阈值有物理依据，优先采用约束式：

```text
minimize_C      mu(C)
subject to      R_struct(C) <= epsilon_struct
                B_SI(C) <= B_max
                P_index_error(C, SNR0) <= epsilon_rx
                c2_q[m] in Omega_m
```

其中 `mu(C)` 是候选 Gram 矩阵的最大非对角相关性。之后再以 PAPR CCDF
验证“低相关码本是否带来更好的选择增益”。

也可以把 CCDF 尾部或平滑峰值直接放入目标：

```text
minimize_C  E_x[softmax_n |s_q[n]|^2] + lambda * mu(C)
```

但权重 `lambda` 必须通过 Pareto 曲线或敏感性实验解释，不能只凭经验给
一个数。

### 可以形成的贡献

- 推导 `c2` pattern 与候选波形相关性的闭式关系；
- 在 full-diversity、phase-safety 或结构风险约束下设计低相关码本；
- 解释候选相关性、候选数与 PAPR 选择增益之间的关系；
- 设计固定的离线码本，避免逐帧传输 grouping 结构。

### 最小验证实验

固定 `N=64`、候选数和搜索预算，比较：GPS 候选、当前局部扰动候选、
随机候选、低相关候选。先画 `max/mean correlation` 与 CCDF `10^-3`
处 PAPR 的散点关系，再比较 BER 和结构指标。

### 放弃条件

若低相关性与 PAPR 尾部改善在不同调制和帧数下没有稳定关系，或优化码本
相对随机/当前码本只改善约 `0.2 dB` 以内，则不把它作为主创新。

### 当前判断

优先级：高。风险：中。它既回应 Yuan 对 grouping 的开放问题，也能为
候选构造提供定性分析和显式模型，最符合师兄关于盲审的经验。

## 5. 候选方向 B：相关性/峰值贡献驱动的分组设计

### 核心问题

GPS 默认等长 adjacent 或 comb grouping，但不同子载波对时域峰值和二次
相位扰动的敏感性不同。可以把 grouping 本身视为图划分或集合划分问题。

一种候选模型是：以 chirp 子载波为节点，以二次相位响应相似度、对主要
峰值样本的共同贡献，或候选切换后的波形差异为边权，求固定组数下的划分。

```text
minimize_P    alpha * within_group_coherence(P)
              + beta * candidate_coherence(P)
subject to    |G_v| in [M_min, M_max]
              number_of_groups = V
```

这里仍要先回答 grouping 是离线固定还是逐帧变化。为了控制边信息和专利
复杂度，第一版应优先做离线固定 grouping；逐帧 data-aware grouping 可
留作高风险扩展。

### Baseline

- equal adjacent grouping；
- comb/interleaved grouping；
- random grouping；
- unequal adjacent grouping；
- 所提相关性或峰值贡献 grouping。

### 最小验证实验与放弃条件

小规模下枚举全部或大量随机 partition，观察设计指标能否预测真正 PAPR。
若指标不能稳定筛出前 10% 的 partition，停止复杂聚类算法；如果有效，
再扩展到 `N=64/128/256`。

### 当前判断

优先级：高。风险：中高。它有明确的文献缺口，但单纯使用聚类算法并不
够，必须先建立 AFDM 相关的边权和机制解释。

## 6. 候选方向 C：具有 AFDM 节点下界的预算受限树搜索

### 为什么普通 MCTS 不够

OFDM PTS 已有 W-way tree、sphere decoding、MCTS 和多种 branch-and-bound；
CN117014280A 还直接保护了 MCTS 搜索 OFDM 分块相位因子的方案。因此，
“把 MCTS 用于 AFDM `c2` 搜索”可以做实验，但专利和论文新颖性都偏弱。

当前代码已经具备 partial waveform reuse 和 beam search。更值得研究的是
能否利用 AFDM 分组部分波形构造一个可解释的节点评价或下界。

若已选择组的部分和为 `p[n]`，剩余第 `g` 组候选部分波形为
`u_g,w[n]`，则一个保守的样本幅度下界为：

```text
L_n = max(|p[n]| - sum_g max_w |u_g,w[n]|, 0).
L_node = max_n L_n^2 / P_avg.
```

当 `L_node` 已不小于当前完整候选的最好 PAPR 时，可以剪掉该分支。
这个下界可能较松，但它提供了可验证的 AFDM/部分波形依据。还可研究：

- 优先处理对当前主峰最敏感的 group；
- 用少量 dominant samples 做快速节点排序，再用 4 倍过采样复核叶节点；
- 固定节点预算或时延预算，输出 anytime performance；
- 用较强的凸包距离或相位区间界替代简单三角不等式界。

### Baseline

- exhaustive search；
- Yuan/GPS sequential greedy；
- random search；
- 现有 beam search；
- 移植的 MCTS；
- 所提 bounded beam / branch-and-bound。

所有方法必须在相同“完整候选等效评价次数”或相同 wall-clock budget 下
比较，不能只比较迭代次数。

### 最小验证实验与放弃条件

先在 `V=6/8`、`W=3` 上统计剪枝率、找到最优解比例、PAPR gap 和时间。
若下界剪枝率长期低于约 30%，或维护树的开销使运行时间不优于 beam，
则退回“固定预算 beam + dominant-sample metric”，不再包装成精确搜索。

### 当前判断

优先级：高。风险：中。树搜索可以成为方法，但真正的贡献必须是节点界、
搜索顺序或预算机制，而不是 MCTS 名称。

## 7. 候选方向 D：连续优化引导的离散、可接收码本

连续 `c2` 优化通常能得到较低 PAPR，却可能带来逐块参数描述、接收失配
或实时迭代负担。可研究“连续解作为 teacher/initialization，最终投影到
小型结构化码本”，例如：

```text
continuous phase-safe solution
    -> quantization/codebook projection
    -> local tree refinement
    -> small candidate index
```

它的价值在于连接性能上界和可实现方案，但 2026 年已经有 UQP、projected
gradient 和 subcarrier-wise continuous perturbation，碰撞风险很高。

当前判断：优先级中，风险高。适合当性能上界或辅助算法，暂不建议作为
首选创新点。

## 8. 候选方向 E：发送端收益与接收可区分性的联合码本

将候选码本同时按两种目标设计：候选波形在发送端应低相关以提供 PAPR
选择增益；在接收端，错误候选应产生足够大的判决不一致或似然差距。

```text
minimize_C   J_tx(C) + lambda J_rx(C)
subject to   B_SI <= B_max, complexity <= Q_max
```

也可用约束形式：在 `P_index_error <= epsilon` 条件下最小化 PAPR 尾部。
这比“设计一个半盲检测器”更完整，因为码本和检测器共同决定可恢复性。

但 Choi 已经提出 well-separated scalar chirp-offset codebook 和 intrinsic
SI；affine-domain circular shift 也已有无 SI 检测。因此，只有在 group-wise
pattern、不同接收统计量、单次共享信道估计或更低复杂度方面形成明确区别，
才值得推进。

当前判断：优先级中，风险高。先深读 Choi，再决定是否保留。

## 9. 候选方向 F：从 PAPR 指标升级为功放感知的选择与评价

PAPR 是功放问题的代理指标。相同 PAPR 的两个候选，在实际 PA 下可能有
不同 EVM、ACLR、谱再生和 BER。可以引入 modified Rapp PA，比较：

- PAPR/CM；
- EVM；
- ACLR 或 OOBE；
- 不同 IBO 下的 BER/BLER；
- PA 输出与理想线性输出的相关性。

若形成方法，可在候选选择时使用 PA-aware objective；若时间不足，也至少
把这部分作为当前方法的完整评价，回答“降低 PAPR 到底换来了什么”。

由于已有 PAPR-OOBE、spectral mask、PA nonlinearity 相关 AFDM 工作，单纯
加入第二指标不够新。更现实的定位是主方法的工程验证和参数选择依据。

当前判断：作为评价层优先级高；作为独立创新优先级低。

## 10. 暂不优先的方向

- 纯深度学习或端到端收发机：已有 AFDM transceiver joint optimization、
  pilot-embedded GIS-ACE-Net 和 convolutional autoencoder，训练与解释成本
  也不适合当前毕业时间。
- 单纯 DCT/WHT/ZC/DFT/CI spreading：AFDM 已有较多工作，作为 baseline
  更合适。
- 单纯 companding、clipping 或 SLM：领域过于成熟。
- 完整信道估计、均衡或通用检测器：会使论文主线转向接收机，工作量和
  风险都明显扩大。
- 匆忙转向 ISAC：PAPR-WISL、chirp-subcarrier reservation 等近期工作已
  进入该方向，且需要新增 sensing 指标和系统模型。

## 11. 师兄经验转化成研究设计规则

### 11.1 先讲清逻辑，再增加算法

每一章都应能回答四句话：

1. 我要解决什么具体问题；
2. 参考方法在什么条件下存在什么不足；
3. 我改变了哪个变量、约束或机制；
4. 为什么这个改变可能有效，实验怎样验证这个因果解释。

如果只能写成“依次使用 A、B、C 算法”，说明仍在堆砌方法。

### 11.2 问题必须显式建模

至少明确：优化变量、目标函数、约束、已知量和复杂度预算。若同时考虑
多个指标，不必自动使用加权和：

- 有明确工程阈值时，用约束式更容易解释，例如最小化 PAPR，同时要求
  BER、结构风险、SI 和复杂度不超过阈值；
- 没有公认阈值时，给 Pareto 曲线，再说明选择某个工作点的原因；
- 使用加权和时，要给权重敏感性实验或归一化方法，避免权重只为得到一条
  好看的曲线。

### 11.3 参数必须有依据

建议把参数分成三类：

- 标准或文献参数：如 3GPP TR 38.901 的 TDL/CDL 类模型、Yuan 的
  `N/V/W` 设置；
- 理论约束参数：如 full-diversity 条件、phase-safety bound、SI bit 数；
- 实验选择参数：如 `delta`、beam width、tree budget，应通过 sweep、
  Pareto knee 或消融实验选取。

“沿用某文献参数”只能证明可比性，不能单独证明参数合理。最好同时给出
参数变化趋势。

### 11.4 Baseline 同时覆盖经典与近期

经典 baseline 建议至少包括 conventional AFDM、clipping/companding 中
一种，以及 PTS/SLM 思想的代表。近期 baseline 应优先覆盖 GPS、CSM 或
Choi 中与最终方向最接近者，并按需要加入 spreading 或 continuous `c2`
optimization。所有 baseline 使用相同调制、功率、过采样、候选预算、
信道 realization 和接收信息。

### 11.5 文献结构

大部分正文引用应来自近两三年 AFDM/OTFS/新波形研究；经典 OFDM PAPR、
PTS、SLM、TR 等奠基论文可以更老，但应明确它们承担的是理论源头或经典
baseline，而不是当前研究现状。

## 12. 建议的第一轮筛选实验

第一轮目的不是做正式论文图，而是用较低成本筛掉弱方向。

1. **候选相关性试验**：生成 100 至 500 个不同码本/partition，检验
   coherence 与 PAPR CCDF 分位点的相关性。
2. **grouping 预测试验**：小 `N` 枚举 partition，判断所设计 grouping
   指标能否预测优质分组。
3. **树搜索试验**：在相同 node budget 下比较 greedy、beam、MCTS 和
   bounded search，记录 PAPR gap、命中最优率、时间和方差。
4. **功放验证试验**：只对前三种有希望的方法加入 Rapp PA，检查 PAPR
   排名是否与 EVM/BER/ACLR 排名一致。

第一轮完成后，再从结果中组合论文方向。例如：

- 若 A 有效、C 有效：形成“结构约束码本设计 + 低复杂度搜索”；
- 若 A 有效、C 无效：形成“码本/分组设计 + 性能与接收代价分析”；
- 若 A 无效、C 有效：保留当前候选构造，主攻 AFDM 特定搜索；
- 若 A/C 均弱：回到 PA-aware 或 pilot/guard-assisted 方向重新评估。

## 13. 和导师讨论时建议提出的问题

1. 老师更认可“候选/码本构造的理论解释”，还是“低复杂度组合搜索”作为
   主创新？
2. 毕业论文是否接受一个较强主方法，加上一章完整接收/功放/复杂度分析，
   还是必须两个彼此更独立的方法？
3. 专利是否继续围绕当前局部扰动组合方案收窄，论文允许在其上扩展码本
   与搜索，而不要求二者完全同构？
4. 正式仿真应采用哪类场景作为主场景：文献常用离散 LTV 路径模型，还是
   增加 3GPP TDL/CDL/高速移动参数作为外部有效性验证？
5. 对近期 2026 年 AFDM PAPR 工作，老师认为需要实现哪些强 baseline，
   哪些只需文献对比？

## 14. 关键来源

- Yuan et al., PAPR Reduction with Pre-chirp Selection for AFDM:
  https://arxiv.org/abs/2406.14064
- Li et al., Chirp Parameter Selected Mapping for Low PAPR AFDM:
  https://www.eurecom.edu/en/publication/8431
- Choi, Low-PAPR AFDM With Intrinsic SI and Low Complexity:
  https://doi.org/10.1109/LWC.2026.3655748
- Cui et al., Adaptive c2-Perturbed AFDM Waveform Design:
  https://arxiv.org/abs/2606.04698
- Savaux et al., Joint PAPR and OOBE Reduction for AFDM:
  https://arxiv.org/abs/2609.01255
- Liu et al., AFDM Transceiver Optimization for PAPR Reduction:
  https://doi.org/10.1109/TCOMM.2026.3700535
- Karthiga and Deepa, UQP pre-chirp selection:
  https://doi.org/10.1016/j.phycom.2025.102731
- Karthiga and Deepa, multi-start projected gradient pre-chirp tuning:
  https://doi.org/10.1016/j.phycom.2026.103035
- Huang et al., pilot-embedded GIS-ACE-Net:
  https://doi.org/10.1109/LCOMM.2025.3643615
- Ye et al., MCTS for OFDM PTS:
  https://doi.org/10.1109/ICCC52777.2021.9580407
- Qi et al., W-way tree PTS:
  https://doi.org/10.1109/LCOMM.2012.072012.121228
- Alavi et al., sphere decoding for PTS:
  https://doi.org/10.1109/LCOMM.2005.11014
- Afrasiabi-Gorgani and Wunder, conditional expectations for PAPR/CM:
  https://arxiv.org/abs/1909.10639
- Tsuda and Umeno, SDR and randomization for PTS:
  https://doi.org/10.1587/transcom.2019EBP3243
- CN117014280A, MCTS-based OFDM PAPR suppression:
  https://patents.google.com/patent/CN117014280A/zh
- WO2017064448A1, tree-search tone reservation for OFDM:
  https://patents.google.com/patent/WO2017064448A1/en
- 3GPP TR 38.901 channel models:
  https://portal.3gpp.org/desktopmodules/Specifications/SpecificationDetails.aspx?specificationId=3173

部分 2026 论文目前只获得摘要或出版元数据，正式确定选题前仍需获取全文，
逐项核对问题模型、权利要求和实验设置。
