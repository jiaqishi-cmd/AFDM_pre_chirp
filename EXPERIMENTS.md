# AFDM / GPS / Proposed Experiment Index

本文档用于冻结当前工程的实验脉络，避免继续无序增加脚本和结果文件。当前阶段建议先维护这条证据链，再决定是否做目录级重构。

## 当前主线

### 1. GPS worst-case vulnerability

目标：说明 Yuan GPS-AFDM 的 per-subcarrier `c2_m` 虽然能改善 PAPR，但在某些结构化两径信道下更容易匹配 Bemani key-equation，从而出现 worst-case 性能风险。

主脚本：

- `matlab_afdm/main_focused_bemani_gps_search.m`
- `matlab_afdm/run_caseA_fixed_ber_deep.m`
- `matlab_afdm/run_caseA_theta_scan.m`
- `matlab_afdm/plot_caseA_theta_scan_enhanced.m`

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

- `matlab_afdm/run_choi_style_random_channel_ber.m`

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

- `matlab_afdm/run_papr_ccdf_comparison.m`
- `matlab_afdm/run_delta_sweep_proposed.m`
- `matlab_afdm/run_proposed_delta_ccdf_curves.m`
- `matlab_afdm/run_proposed_delta_ber_curves.m`

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
- `matlab_afdm/test_c2_structural_metrics.m`

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

- `matlab_afdm/run_papr_search_complexity_study.m`
- `matlab_afdm/run_partial_reuse_theory_and_timing.m`
- `matlab_afdm/run_partial_reuse_M_sweep.m`
- `matlab_afdm/run_partial_reuse_topK_sweep.m`

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

- `matlab_afdm/main_bemani_gps_key_equation_search.m`
- `matlab_afdm/main_gps_unique_rank_loss_search.m`
- `matlab_afdm/main_gps_bestcase_ber_snr.m`
- `matlab_afdm/main_caseA_phase_scan_dmin.m`
- `matlab_afdm/main_caseA_besttheta_ber_snr.m`
- `matlab_afdm/main_focused_bestcase_ber_snr.m`
- `matlab_afdm/run_c2_channel_diagnostics.m`
- `matlab_afdm/run_c2_exposure_ber.m`
- `matlab_afdm/run_bestcase_ber_snr.m`
- `matlab_afdm/scan_theta_dmin_caseA.m`

## Results cleanup policy

- `results/` 根目录保留当前较有价值的最终图、`.mat` 和 `.csv`。
- `results/archive_intermediate/` 存放重复的 smoke test 图或早期不满意版本。
- 不直接删除中间结果；确认无用后再统一清理。

## 下一步建议

1. 暂停新增实验类型。
2. 固定 proposed 推荐配置。
3. 对最终需要进报告的图，重新跑高统计量版本：
   - PAPR CCDF: `numFrames >= 1e4`
   - random channel BER: 每个 SNR 尽量保证足够 error count 或明确标注统计不稳定点
   - Case A fixed-channel BER: 保留 stress-test 定位，不和 average BER 混为一谈
4. 后续若要重构代码目录，先把被调用函数从脚本中拆成公共 `utilities`，再移动脚本，避免 MATLAB path 断裂。
