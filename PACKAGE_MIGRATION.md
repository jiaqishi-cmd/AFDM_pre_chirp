# MATLAB Package Migration Plan

目标是逐步从全局 `addpath` 调用迁移到 MATLAB package 命名空间调用，减少函数名冲突和路径污染。

## Long-Term Target

最终形态建议：

```text
matlab_afdm/
  +afdm/
    +config/
    +tx/
    +rx/
    +channel/
    +chirp/
    +metrics/
    +search/
    +experiments/
  experiments/
```

新代码优先使用 package 调用，例如：

```matlab
metrics = afdm.metrics.c2_structural(c2Vec, cfg);
corr = afdm.metrics.c2_candidate_correlation(candidateVecs, cfg);
```

## Current Migration Status

第一阶段已完成：

- 新增 `matlab_afdm/+afdm/+metrics/`
- 新增 package functions:
  - `afdm.metrics.c2_structural`
  - `afdm.metrics.phase_structural`
  - `afdm.metrics.c2_candidate_correlation`
- 原有函数保留为 compatibility wrappers:
  - `calc_c2_structural_metrics`
  - `calc_phase_structural_metrics`
  - `calc_c2_candidate_correlation`
- 已将结构风险实验和 delta sweep 中的 metrics 调用改为 `afdm.metrics.*`

## Important MATLAB Constraint

MATLAB package 并不是完全不需要路径。它要求 package 父目录在 MATLAB 当前目录或 path 中。

也就是说：

- 若当前目录是 `matlab_afdm/`，可以直接调用 `afdm.metrics.*`
- 若从别处运行，仍需要：
  - 打开 MATLAB Project，或
  - 将 `matlab_afdm/` 加入 path，或
  - 从 `matlab_afdm/` 目录启动实验

Package 的收益不是“绝对不需要任何 path”，而是：

- 不需要把每个子目录都加入 path
- 函数命名空间清楚
- 后续可以只暴露 `matlab_afdm/` 一个根目录

## Recommended Next Steps

### Stage 2: chirp/search package

迁移稳定的 pre-chirp 和 PAPR search helper：

```text
+afdm/+chirp/
+afdm/+search/
```

候选函数：

- `build_c2m_gps_pattern`
- `build_c2m_proposed_pattern`
- `build_group_index`
- `gps_candidate_set`
- `greedy_group_papr_selection`

### Stage 3: tx/rx/channel package

迁移核心链路：

```text
+afdm/+tx/
+afdm/+rx/
+afdm/+channel/
```

候选函数：

- `idaft_mod`
- `compute_papr`
- `afdm_tx_engine`
- `daft_demod`
- `mmse_equalize`
- `estimate_effective_channel`
- `multipath_channel`
- `add_awgn`

这一步需要系统性改调用名，风险较高，应在单独 commit 中做，并配套 smoke tests。

### Stage 4: scripts call package only

实验脚本只保留：

```matlab
rootDir = find_afdm_root(...);
cd(rootDir);
```

然后调用 package 函数，不再依赖 `genpath` 暴露所有子目录。

## Compatibility Policy

每迁移一个模块，旧函数先保留 wrapper 至少一个阶段：

```matlab
function y = old_function(varargin)
    y = afdm.module.new_function(varargin{:});
end
```

等所有脚本都切到 package 调用后，再删除旧 wrapper。
