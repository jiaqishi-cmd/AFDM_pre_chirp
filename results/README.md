# Results Directory

`results/` 根目录保留当前阶段较有价值的输出图、`.mat` 和 `.csv`。

已经归档的中间图放在：

- `results/archive_intermediate/`

归档原则：

- 保留较新的正式输出；
- 将明显重复的 smoke test 图、早期低统计量图移动到 `archive_intermediate/`；
- 暂时不删除结果文件，避免误删后续需要追溯的实验。

当前结果和脚本对应关系见仓库根目录：

- `EXPERIMENTS.md`
