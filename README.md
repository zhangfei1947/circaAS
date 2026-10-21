# circaAS 0.1.0

GitHub 项目名为 **circaAS**；Python 包及命令行名称保持 **circadian-as**。

第一版是可安装的 Python 核心包与命令行工具。BB 偏差修正为主引擎，连续 BLN 保留统一接口但尚未实现。当前版本用于探索性筛选，不承诺全基因组 FDR 控制。

先稳定计数输入、模型数值、结果语义和输出协议，再以同一核心封装 R 接口；现在直接移植成 R 会额外引入优化器与随机数差异。核心代码已按输入、覆盖诊断、引擎、运行与追溯分层，不是一次性分析脚本。

## 安装与运行

需要 Python 3.10+、NumPy、SciPy。可选图形依赖 matplotlib。

```bash
git clone https://github.com/zhangfei1947/circaAS.git
cd circaAS
python -m venv .venv
source .venv/bin/activate
python -m pip install .
circadian-as --version
circadian-as engines
```

不需要原始研究数据的最小示例（完全合成的演示计数，不用于性能评测）：

```bash
python examples/make_demo.py --out demo.jsonl
circadian-as run --input demo.jsonl --out demo-run --bootstrap 199 --workers 2
```

分析自己的 rMATS 原始计数：

```bash

circadian-as prepare-rmats \
  --counts-dir /path/to/rmats/output \
  --samples /path/to/samples.tsv \
  --bam-list-1 /path/to/all_bams.txt \
  --dataset my_dataset --mode JC --out prepared

circadian-as coverage --input prepared/records.jsonl --out coverage-run

circadian-as run --input prepared/records.jsonl \
  --out analysis-run --bootstrap 199 --workers 4 --seed 20261020

# event 路径由 results.tsv 的 artifact 列给出。
python -m pip install '.[plot]'
circadian-as plot --event analysis-run/events/EVENT_HASH.json --out event.png
```

samples.tsv 的字段与 BAM 对齐要求见[输入规范](docs/input-format.md)。所有输出目录必须尚不存在，避免覆盖既有结果。CLI 每个计算进程设置一个 BLAS 线程；`--workers` 控制事件并行度。

## 模型与覆盖

每个条件独立分析，使用原始 inclusion/skipping counts 与有效长度；不是对观测 PSI 直接作等权回归。固定周期为 24 小时，logit PSI 包含实验固定截距和 sin/cos 项；观测计数概率加上 `log(LI/LS)` offset。BB 分布处理抽样误差与过离散，零 reads 样本不贡献似然，不插补 PSI。

主检验沿用 V2 likelihood-ratio statistic；99 次 pilot 估计 null 离散度偏差，再用修正后的 null 进行参数 bootstrap。默认最终 199 次，失败模拟保守计入上尾，pilot/最终失败比例超过既定阈值则不发布主 p 值。数值代码只改相对 import，附原始 SHA256 和新旧同种子回归测试。

TempCompen r1/r2 按贯穿时间序列的两次实验处理。该版本不估计时间自相关，不搜索周期，不做菌株/温度差异节律检验，也不联合拟合 gene expression。条件化于实际事件深度 N 能反映 reads 少导致的信息减少，但不能恢复缺失相位，也不解决 PSI 依赖的非随机缺失。

覆盖门槛逐字保留 V2 实现，仅使用 N、时间和实验标签。另报告各绝对时间/相位深度、跨周期重复零覆盖相位、有效样本量、相位间隙和设计条件数。这些附加诊断不改变过滤条件。

## 四类结果

| status | 含义 |
|---|---|
| `rhythm_candidate` | 覆盖合格且检验数值可靠，原始 p ≤ alpha；探索性候选 |
| `not_detected` | 覆盖合格且检验可靠，未达到阈值；不能推断不振荡 |
| `coverage_insufficient` | 原覆盖规则不通过；p 缺失 |
| `model_unresolved` | 覆盖合格但拟合/重采样不能可靠完成；p 缺失 |

单一 splice support 等不可识别情形也属于 model_unresolved。覆盖足够不等于同时具备可靠的剪接比例推断条件。不提供“确定不振荡”标签。

199 次最终 bootstrap 的 p 分辨率为 0.005；提高次数减少 Monte Carlo 误差，不修复模型错设。事件 JSON 中的 Monte Carlo 区间仅描述条件于拟合 null 的最终模拟尾概率，排除 pilot、参数和模型不确定性，不能解释为完整节律 p 的置信区间。阈值附近事件标注 `threshold_ambiguous`，不自动改动分类。

## 可追溯输出

- `results.tsv`：每个输入事件一行，四类状态、原因、可用 p、效应量、相位与事件 JSON 路径。
- `coverage.tsv`：冻结覆盖规则的逐事件指标与失败原因。
- `events/*.json`：输入计数、时间、实验、有效长度、原注释、额外覆盖诊断、拟合、bootstrap、Monte Carlo 诊断和按实验预测曲线。
- `manifest.json`：输入及伴随来源 SHA256、参数、软件/依赖版本、代码与全部输出散列、状态计数及运行完成标记。

事件随机种子由 root seed 与无歧义 record_key 的 SHA256 生成，独立于事件输入次序和进程数量。输入错误直接拒绝运行；意外中断标记整个 run 为 failed，不伪装为完成。数值失败诊断中的历史 p=1 哨兵仅保留在内部对象，主输出不采用它。

## 验证与下一步

参见 [验证说明](validation/VALIDATION.md)、[输入规范](docs/input-format.md)、[引擎接口](docs/engine.md) 与 [本次运行结果](validation/RESULTS.md)。两套数据都是开发集，无认可真阳性；全量导入/覆盖扫描与固定面板拟合分开报告。

后续可先提供 R 的 `system2()` 包装与结果读入/绘图，让 R 和 CLI 共用相同引擎和 schema。连续 BLN 应在不改变 BB 主结果的条件下返回独立敏感性结果，并经过自己的校准验证；当前 `--engine bln_continuous` 会明确报未实现。

仓库包含源码、测试、模拟回归 fixture、验证脚本和汇总。原始 rMATS/BAM 数据、完整 332 项面板、逐事件分析输出和本地审计文件不随公开仓库分发；验证文档中这些路径描述的是原开发环境。汇总见 [published_summary.json](validation/published_summary.json)。

开发测试：`python -m unittest discover -s tests -v`。
