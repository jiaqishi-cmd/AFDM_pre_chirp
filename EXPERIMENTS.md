# AFDM / GPS / Proposed Experiment Index

本文档用于冻结当前工程的实验脉络，避免继续无序增加脚本和结果文件。当前阶段建议先维护这条证据链，再决定是否做目录级重构。

## 当前主线

### 1. GPS worst-case vulnerability

目标：说明 Yuan GPS-AFDM 的 per-subcarrier `c2_m` 虽然能改善 PAPR，但在某些结构化两径信道下更容易匹配 Bemani key-equation，从而出现 worst-case 性能风险。

主脚本：

- `matlab_afdm/experiments/bemani/main_focused_bemani_gps_search.m`
- `matlab_afdm/experiments/caseA/run_caseA_fixed_ber_deep.m`
- `matlab_afdm/experiments/caseA/run_caseA_theta_scan.m`
- `matlab_afdm/experiments/caseA/plot_caseA_theta_scan_enhanced.m`

核心结果：

- `results/focused_mean_error_20260508_112141.png`
- `results/focused_phi_sigma_ratio_20260508_112141.png`
- `results/focused_phi_col_corr_20260508_112141.png`
- `results/caseA_fixed_ber_deep_20260509_105956.png`
- `results/fig_caseA_theta_raw_smooth_20260508_170746.png`
- `results/fig_caseA_theta_ratio_20260508_170746.png`
- `results/fig_caseA_theta_statistics_20260508_170746.png`

当前结论：

- Case A: `N=64, V=4, M=16, L=16` 是结构上比较危险的两径 profile。
- GPS 在 Bemani key-equation mismatch 上明显小于 baseline。
- GPS 的 `Phi(delta)` 指标中出现更低 `sigma_ratio` 和更高 column correlation。
- 在 fixed phase `h=[1,-1]/sqrt(2)` 下，BPSK BER 出现明显 high-SNR floor。
- theta scan 显示 GPS 对路径增益相位更敏感，存在多个高 BER 区域。

注意：

- 这条线讨论的是 worst-case / stress-test，不代表 Rayleigh average channel 下必然出现 BER floor。
- 主要依据是 Bemani key-equation、`Phi(delta)` 和 finite-constellation BER，不再使用 total `H_eff` rank 作为主结论。

## 2. Proposed average-channel BER

目标：验证 proposed group-wise small `c2` perturbation 在普通随机 doubly selective channel 下没有明显平均 BER 损失。

主脚本：

- `matlab_afdm/experiments/random_channel/run_choi_style_random_channel_ber.m`

核心结果：

- `results/fig_choi_style_random_channel_ber_20260509_154723.png`
- `results/fig_choi_style_random_channel_ber_vs_delta_20260509_154723.png`
- `results/results_choi_style_random_channel_ber.mat`

当前结论：

- 在 Choi-style random doubly selective channel 下，baseline、GPS、proposed BER 曲线整体接近。
- 没有观察到 proposed 因 small perturbation 带来明显平均 BER 损失。
- 这与 Case A fixed-channel stress test 互补：random channel 看平均表现，Case A 看 worst-case vulnerability。

## 3. Proposed PAPR and delta sweep

目标：研究 proposed 的 `{c2-delta, c2, c2+delta}` 三点扰动对 PAPR、结构风险和 fixed-channel BER 的影响。

主脚本：

- `matlab_afdm/experiments/papr/run_papr_ccdf_comparison.m`
- `matlab_afdm/experiments/papr/run_delta_sweep_proposed.m`
- `matlab_afdm/experiments/papr/run_proposed_delta_ccdf_curves.m`
- `matlab_afdm/experiments/papr/run_proposed_delta_ber_curves.m`

核心结果：

- `results/papr_ccdf_comparison_20260508_165645.png`
- `results/papr_summary_comparison_20260508_165645.png`
- `results/proposed_delta_ccdf_curves_20260509_111701.png`
- `results/proposed_delta_ber_curves_20260509_112841.png`
- `results/delta_sweep_papr_20260509_104315.png`
- `results/delta_sweep_ber_20260509_104315.png`
- `results/delta_sweep_structural_20260509_104315.png`
- `results/delta_sweep_selection_20260509_104315.png`

当前结论：

- Proposed 可以获得 PAPR 改善。
- 较大的 `delta/c2` 往往改善 PAPR 更明显，但需要和结构风险、BER 保持平衡。
- 当前推荐继续重点观察 `delta/c2 = 0.1` 和 `0.2`，其中 `0.2` 在已有结果中表现较稳。

## 4. Channel-independent structural risk

目标：完全不依赖信道，只从 `c2_m` pattern 的相位离散度、星座差分对齐度等角度衡量结构风险。

主脚本与函数：

- `matlab_afdm/calc_c2_structural_metrics.m`
- `matlab_afdm/calc_phase_structural_metrics.m`
- `matlab_afdm/calc_c2_candidate_correlation.m`
- `matlab_afdm/experiments/structural/test_c2_structural_metrics.m`

核心结果：

- `results/c2_structural_phase_hist_20260508_163359.png`
- `results/c2_structural_risk_bar_20260508_163359.png`
- `results/results_c2_structural_metrics.csv`
- `results/results_c2_structural_metrics.mat`

当前结论：

- GPS 的粗相位集合会带来较高 structural risk。
- Proposed small perturbation 的结构风险更接近 baseline。
- `R_struct` 不是 full-diversity 证明，只是 channel-independent risk indicator，需要和 key-equation / BER 交叉验证。

## 5. Search complexity and partial waveform reuse

目标：降低 proposed PAPR search 的计算量，验证 partial waveform reuse 是否等价于 full recompute，并研究不同 `M` 和 `topK` 下的复杂度-性能折中。

主脚本：

- `matlab_afdm/experiments/complexity/run_papr_search_complexity_study.m`
- `matlab_afdm/experiments/complexity/run_partial_reuse_theory_and_timing.m`
- `matlab_afdm/experiments/complexity/run_partial_reuse_M_sweep.m`
- `matlab_afdm/experiments/complexity/run_partial_reuse_topK_sweep.m`

核心结果：

- `results/results_partial_reuse_theory_and_timing.mat`
- `results/results_partial_reuse_M_sweep.mat`
- `results/results_partial_reuse_topK_sweep.mat`
- `results/fig_partial_reuse_M_runtime.png`
- `results/fig_partial_reuse_M_speedup.png`
- `results/fig_partial_reuse_M_ifft_count.png`
- `results/fig_topK_papr_1e3.png`
- `results/fig_topK_papr_1e2.png`
- `results/fig_topK_runtime.png`
- `results/fig_topK_eval_count.png`
- `results/fig_topK_complexity_tradeoff.png`

当前结论：

- Partial waveform reuse 与 full recompute 数值一致：self-check relative error 约 `2e-16`，pattern agreement `100%`。
- `M=64,V=8,W=3,topK=8,searchOS=2,finalOS=4` 下，IFFT count 从 `164` 降到 `48`，runtime speedup 约 `1.7x`。
- M 越大，reuse 的 measured speedup 趋势上升。
- `topK=4/8` 是较实用折中；`topK=16` PAPR 更好但 runtime 增加明显。

## 6. Common-pilot dynamic c2 reconstruction

目标：验证“固定公共训练块估计物理信道，再为数据相关的候选 `c2` 模式重构有效信道并完成模式检测”是否具备继续研究的可行性。

主脚本：

- `matlab_afdm/experiments/rx_detection/run_common_pilot_c2_reconstruction_validation.m`
- `matlab_afdm/experiments/rx_detection/run_common_pilot_c2_spacing_sweep.m`
- `matlab_afdm/experiments/rx_detection/run_common_pilot_support_error_sweep.m`

代表结果：

- `results/common_pilot_c2_validation_20260909_095440.csv`
- `results/common_pilot_c2_validation_20260909_095440.mat`
- `results/common_pilot_c2_validation_20260909_095440.png`
- `results/common_pilot_c2_spacing_20260909_100556.csv`
- `results/common_pilot_c2_spacing_20260909_100556.png`
- `results/common_pilot_support_error_20260909_100942.csv`
- `results/common_pilot_support_error_20260909_100942.png`

当前结论：

- 在已知粗时延/多普勒支撑、公共训练块 LS 估计路径增益的受控条件下，100 帧/点的路径增益 NMSE 约从 15 dB 的 `2.5e-3~2.8e-3` 降至 25 dB 的 `2.6e-4~2.8e-4`。
- 对四个分组小扰动 `c2` 向量模式，20 dB 的模式识别率为 `90%~99%`，25 dB 为 `96%~100%`；分数多普勒增大时模式检测更敏感。
- 已知正确模式时，估计 CSI 相对完美 CSI 的 BER 损失较小；未知模式时的主要额外损失来自模式误判。
- 四模式 PAPR 平均收益约 `0.8~1.2 dB`，说明候选可分辨性与 PAPR 候选多样性需要联合设计。
- 在 20 dB、分数多普勒 `0.4` 下，`delta/c2` 从 `1/64` 增至 `1/16` 时 PAPR 收益由 `0.43 dB` 增至 `1.01 dB`，随后基本饱和；模式正确率没有随间隔单调增加。
- 较大的 `delta` 会提高部分几何距离，但一旦模式判错，错误模式与真实模式的失配更严重，因此 detected-mode BER 反而由约 `1.92e-2` 上升至 `2.75e-2`。
- 接收机多普勒支撑误差为 `0/0.05/0.1/0.2` bin 时，模式正确率约为 `94.7%/89.3%/70.3%/46.0%`，detected-mode BER 约为 `2.38e-2/4.52e-2/1.05e-1/2.17e-1`。该接收链对 off-grid 支撑误差较敏感。

注意：

- 当前使用一个完整已知训练块，尚未优化导频开销。
- 当前假设路径时延和多普勒支撑已经由粗估计获得，只对复路径增益做 LS；不能将结果表述为完整信道估计方案。
- 当前 `min_mode_separation` 只是有效信道 Frobenius 距离，尚不能准确预测模式错误率，需要发展 receiver-aware 距离或错误概率界。
- 新增的交叉模式距离正确比较了 `A_j H A_i^H` 与 `A_j H A_j^H`，但在严重分数多普勒下仍不能单独预测误判；码本评价还必须考虑误判概率和误判后的 BER 代价。

## 7. Fractional-Doppler waveform-level gate

目标：检查整数多普勒下构造的安全 `c2` 码本是否会被分数多普勒泄漏击穿，并在进入接收机研究前验证分数信道模型。

主脚本：

- `matlab_afdm/experiments/structural/run_fractional_model_crosscheck.m`
- `matlab_afdm/experiments/structural/run_fractional_doppler_c2_safety_gate.m`

核心结果：

- 当前分数多普勒实现与显式逐样本复指数模型的相对误差约 `4e-16`。
- 下载的公开参考包使用非整数矩阵幂 `W^nu`；其整数结果一致，但 `beta=0.5` 时相对显式模型误差约 `0.97`，不能原样作为分数多普勒真值。
- 在三个已知整数危险支撑 `(2,2)`、`(5,-3)`、`(8,0)` 上，GPS-rational 的精确秩亏只出现在整数点；分数偏移使其最小奇异值转为正值。
- local-proposed 与 integer-safe GPS 在当前候选事件集合和 `beta in [-0.5,0.5]` 下没有出现新秩亏。
- `|beta|=0.5` 时单路径最大一个元素只保留约 `40.5%` 能量，前三个元素约 `85.6%`，95%能量约需八个元素。
- 不同 `c2` 方案的路径矩阵幅度差异低于 `4.5e-16`，说明 `c2` 不控制泄漏幅度，主要影响相位和码字差几何。

当前判断：分数多普勒会显著影响稀疏性和接收机，但首轮证据没有表明它会使整数安全码本失效。暂不将其提升为第二核心章；作为第一章鲁棒性验证继续保留。

补充收尾：二维扫描确认风险仅由相对分数多普勒 `beta2-beta1` 决定，相同相对偏移下的最小奇异值差异不超过 `1.11e-15`。两径具有相同分数部分时，GPS 的秩亏和约 `4e-3` BER 平台继续存在；相对偏移为 `0.05` 后平台在当前采样中消失。

## 8. MIMO PAPR shallow gate

目标：判断多天线独立 `c2` 选择是否自然形成联合 PAPR 优化问题。

主脚本：

- `matlab_afdm/experiments/papr/run_mimo_c2_papr_separability_gate.m`

核心结果：

- 无跨天线约束时，联合 minimax PAPR 问题严格可分，每根天线独立选择就是全局最优；`Nt=2/4` 的穷举验证误差为零。
- 强制共用一个 `c2` 模式后才产生耦合，`Nt=2/4/8` 的平均最坏天线 PAPR 代价约为 `0.264/0.526/0.768 dB`。
- MIMO 独立选择所需模式信息为 `Nt*log2(Q)` 比特，共同模式只需 `log2(Q)` 比特，但会损失 PAPR 性能。

当前判断：简单 SISO 到 MIMO 推广不构成第二项工作。只有加入共同模式、STBC、预编码、空间码字距离或耦合 PA 约束后才出现非平凡问题，这将显著增加工程量，暂时保持低优先级。

## 9. Explicit c2 side-information baseline

目标：用最简单的显式边带建立发送端模式索引到接收端匹配 `c2` 解调的端到端接口，并量化模式错误的块级传播。

主函数与脚本：

- `matlab_afdm/+afdm/+tx/encode_mode_index.m`
- `matlab_afdm/+afdm/+rx/decode_mode_index.m`
- `matlab_afdm/experiments/rx_detection/run_explicit_c2_side_info_baseline.m`

核心结果：

- 四模式只需2个信息比特；未编码BPSK开销为16QAM、`N=64`数据块的 `0.781%`，三重复为 `2.344%`。
- 当数据/边带SNR为 `10/0 dB` 时，未编码模式错误率为 `16.7%`，三重复降至 `1.6%`；对应数据BER由oracle的 `0.1282` 分别变为 `0.1582` 和 `0.1310`。
- 当数据/边带SNR为 `15/5 dB` 时，未编码模式错误率为 `1.7%`，三重复在1000帧中未观察到模式错误。
- 无噪声下四种模式、重复次数 `1/3/5` 的编解码检查全部通过。

当前判断：显式边带能以很低开销可靠闭合系统，应作为工程基线而非论文核心创新。

## 10. Multipath safety-bound gate

目标：检查路径对安全码本是否会遗漏仅由三/四径联合相关引起的秩亏，并评估Gershgorin/互相关下界能否降低一般多径认证复杂度。

主脚本：

- `matlab_afdm/experiments/structural/run_multipath_safety_bound_gate.m`

核心结果：

- 1652个候选事件、三组P=3/4支撑中，Gershgorin下界与真实最小特征值相关系数约为 `0.98~0.99`，可直接给 `95.5%~99.9%` 的事件正安全证书。
- 14个精确秩亏全部来自GPS-rational的已知路径对循环，没有发现“两两安全、整体才秩亏”的事件。
- pairwise-safe GPS在当前P=3/4事件集合中没有秩亏，但并非所有支撑下的最小特征值都优于baseline/local-proposed。
- 所有方案的最坏Gershgorin下界仍为负，完整认证仍需对剩余事件计算精确特征值。
- 对P=3/4小矩阵，当前MATLAB下界实现比直接特征值分解慢约 `1.5~3.7` 倍，尚无实测复杂度优势。

当前判断：多径下界适合作为第一章的解释性诊断或混合预筛，但当前没有形成可独立支撑第二章的新失效机制与复杂度收益。

## 11. Collective-only multipath search

目标：定向寻找“所有路径对均满秩，但三径/四径联合响应秩亏或严重病态”的有限星座事件，检验路径对安全认证是否充分。

主脚本：

- `matlab_afdm/experiments/structural/run_collective_only_rank_search.m`
- `matlab_afdm/experiments/structural/run_n64_collective_rank_construction.m`

核心结果：

- 合法 `N=9,V=3` 系统穷举全部19682个BPSK差分、9种支撑、84组三径组合和44个模式，找到564个collective-only事件。
- 典型GPS三径事件整体奇异值为 `[1.2247,1.2247,1.46e-15]`，秩为2；任意两列奇异值为 `[1.2247,0.7071]`，相关系数仅0.5。
- GPS-rational的8个路径对安全模式全部出现collective-only秩亏，共384个事件；GPS-offset的8个模式全部出现，共180个事件；local-proposed的27个模式当前未发现精确事件。
- 合法二次幂 `N=8` 穷举未发现collective-only事件，表明精确退化依赖循环长度和系统尺寸。
- 当前 `N=64` 的定向四径构造未发现精确秩亏，但local-proposed模式 `[1,3,1,1]` 出现 `lambda_min=3.45e-10`、条件数约 `7.61e4` 的联合近退化；任意路径对仍有最小奇异值0.707和相关系数0.5。

更新判断：两径认证在一般多径下确实不充分。多径联合安全认证重新提升为第二核心章候选；下一门槛是验证N=64近退化对PEP/BER的影响和PAPR在线选择命中危险模式的概率。

## 推荐默认配置

除非专门做 sweep，后续正式图建议先固定：

```matlab
N_or_M = 64;
V = 8;
W = 3;
delta_ratio = 0.2;
search_os = 2;
final_os = 4;
topK = 8;   % 工程折中
% topK = 16; % 追求更低 PAPR 时使用
```

## 暂时不作为主结果的脚本

以下脚本主要用于探索、调试或早期搜索；保留但不优先放进论文主线：

- `matlab_afdm/experiments/bemani/main_bemani_gps_key_equation_search.m`
- `matlab_afdm/experiments/bemani/main_gps_unique_rank_loss_search.m`
- `matlab_afdm/experiments/bemani/main_gps_bestcase_ber_snr.m`
- `matlab_afdm/experiments/caseA/main_caseA_phase_scan_dmin.m`
- `matlab_afdm/experiments/caseA/main_caseA_besttheta_ber_snr.m`
- `matlab_afdm/experiments/bemani/main_focused_bestcase_ber_snr.m`
- `matlab_afdm/run_c2_channel_diagnostics.m`
- `matlab_afdm/run_c2_exposure_ber.m`
- `matlab_afdm/run_bestcase_ber_snr.m`
- `matlab_afdm/scan_theta_dmin_caseA.m`

## Results cleanup policy

- `results/` 根目录保留当前较有价值的最终图、`.mat` 和 `.csv`。
- `results/archive_intermediate/` 存放重复的 smoke test 图或早期不满意版本。
- 不直接删除中间结果；确认无用后再统一清理。

## 下一步建议

1. 暂停扩展新方向，先验证 `N=64` 联合近退化事件是否具有实际BER/PEP影响，以及PAPR选择器命中危险模式的概率。
2. 若影响成立，比较“路径对安全”和“多径联合安全”筛选后的码本规模、PAPR损失与可靠性收益；若不成立，则不把该方向提升为第二核心章。
3. 固定 proposed 推荐配置，并对最终需要进报告的图重新跑高统计量版本：
   - PAPR CCDF: `numFrames >= 1e4`
   - random channel BER: 每个 SNR 尽量保证足够 error count 或明确标注统计不稳定点
   - Case A fixed-channel BER: 保留 stress-test 定位，不和 average BER 混为一谈
4. 后续若要重构代码目录，先把被调用函数从脚本中拆成公共 `utilities`，再移动脚本，避免 MATLAB path 断裂。
