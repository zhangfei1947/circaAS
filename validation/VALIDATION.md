# 第一版真实数据验证方案

这里的验证分为完整导入正确性、历史固定面板运行、数值引擎回归和结果语义四部分。两套数据都属于开发集；没有预先认可的节律 AS 真阳性。真实数据中检出的候选数不能用来估计误报率或检出率，也不能证明内源性 circadian 调控。

## 完整原始输入导入

使用本地已经保存并留有 SHA256 的原始 rMATS JC.raw.input 与 fromGTF 快照，不从旧的 PSI 筛选结果重新构造输入。分别导入 Neurospora 和 TempCompen 所有五类事件。当前快照共有 111,920 个 event-condition：Neurospora 5,804、TempCompen 106,116。所有事件均应保留，包括全零、极低覆盖和单一 isoform 的事件。

`make_validation_panel.py` 不导入新软件或旧准备脚本。它独立按字面制表符读取原始 count 和 annotation，核对每个事件所有 I/S/N、LI/LS、样本 ID、时间、重复标签、原始完整 annotation 和 BAM 向量位置。Neurospora RI 的 ID140 有不成对的引号，必须保留这一物理行，不可用普通 CSV 引号解析使后续事件吞并。新软件保留原始 annotation；旧面板的 display gene 字符串曾去掉边缘引号，验证中只为该历史展示字段作等价比较。

原始 BAM 顺序与 samples.tsv 逐一核对，并进一步与 Neurospora plan.json 和 TempCompen Table_S4_sample_metadata.tsv 独立核对时间、条件及实验标签。TempCompen r1/r2 是贯穿时间序列的两次实验。新 records 的样本以 sample_id 排序，因此逐 count 比较要按 sample_id 对齐，不能误用位置索引。该排序不改变模型设计，但会改变 bootstrap 随机数与样本的配对；旧引擎回归测试必须使用同一排序和种子。

## 固定面板

只有完成全量导入核对后，再按原 V2 的 332 个 profile_id 选择新 records：Neurospora 147、TempCompen 185。原面板含按覆盖分层的固定 hash 选择及强覆盖节律 profile 补充；最低 total N=20、至少3个正覆盖绝对时间是其历史入选条件。它不代表随机抽取的生物学总体。新软件测试不按新 p 值、PSI 或振幅重新挑选事件。

运行 `make_validation_panel.py --help` 查看参数。输出 `records.jsonl` 与 `selection_audit.json`，后者记录原始文件、完整导入结果、历史参考面板和输出的 SHA256，以及全部所选身份。原始全量导入核对与固定面板运行结果分开报告。

## 结果状态与可追溯

运行输出必须对每个输入保留一条记录，只允许以下四类状态：

- `rhythm_candidate`：覆盖合格、数值检验可靠、固定周期下 p 达到配置的阈值。
- `not_detected`：覆盖合格且数值检验可靠，但未达到阈值；不等于确定不振荡。
- `coverage_insufficient`：覆盖规则不通过。显著性、振幅等推断不得伪装为零。
- `model_unresolved`：覆盖合格但数值计算/重采样可靠性不满足要求；不能计入阴性。

后两类没有可用的主检验 p 值。即使内部失败诊断保留某个中间数值，也不能进入候选判定或未检出分母。结果与 coverage 表应与每事件 JSON 一致，并核对 manifest 的输入散列、参数、代码散列与事件种子。用同一小子集、同一种子，顺序/多进程以及改变输入排列运行应得到相同事件结果（只允许时间等运行元数据不同）。

## 解释边界

第一版主要验证既有 BB 偏差修正的工程实现，不增加新的模拟性能主张。V6 中 BB 在部分高离散度/相同相位缺失场景表现相对较好，但并无统一 5% 误报保证。不得把 bootstrap 候选阈值直接称为已验证全基因组 FDR 控制；199 次 bootstrap 的 p 最小步长为 1/200。

连续 BLN 为预留敏感性接口，不自动与 BB 取交集，也不把未实现接口误标为运行成功。模型条件化于 N=I+S，未单独拟合基因表达，更不恢复完全缺失相位的信息。

## 已完成的验证（0.1.0）

- **完整导入独立重构通过**：111,920 条 event-condition 全部逐 count/metadata/完整 annotation/source_position 匹配原始快照；176 个样本与原 BAM 顺序和原始元数据一致。24 个关键原始输入文件 SHA256 与历史 source_manifest 匹配。见 [selection_audit.json](fixed_panel/selection_audit.json) 和 [source_snapshot_audit.json](fixed_panel/source_snapshot_audit.json)。
- **全量覆盖规则回归通过**：对 111,920 条结果逐一与 V2 冻结 coverage_catalog 对比，所有 gate 决定与 reason 完全相同，1,343,040 个数值字段比较通过（浮点舍入容差 relative 1e-10 / absolute 1e-8）。Neurospora 4,837 合格、967 不足；TempCompen 90,784 合格、15,332 不足。见 [full_coverage_parity.json](full_coverage_parity.json)。
- **固定 332 真实事件完整运行通过**：最终 `run_v1_bb199` 使用 BB pilot99＋最终199，根种子20261020，4 workers。19 个候选、132 个未检出、78 个覆盖不足、103 个模型未解决。103 个未解决包含 51 个 `insufficient_design_or_splice_support`，并非全部数值收敛失败。全部输入/输出 SHA、事件种子、四状态及主 p 缺失语义已独立检查。当前源码 SHA 与运行记录完全匹配。见 [final_run_audit.json](final_run_audit.json)。
- **报告增强未改变结果**：最终版相对预发布仅增加报告用 Monte Carlo 区间和追溯；移除新增 `monte_carlo` 字段后，332 个事件的整个 JSON 内容完全相同，包含统计量、p、状态、原因、拟合与 bootstrap 全部细节。见 [release_parity_audit.json](release_parity_audit.json)。
- **调度与输入顺序回归通过**：按 dataset×历史coverage状态的固定 hash 选4个真实事件，同根种子、最终bootstrap9，使用1 worker/原序和2 workers/逆序，4个事件全部JSON逐字段完全相同。该小预算仅测试工程可复现性；其1个未检出、1个未解决、2个覆盖不足不能作为新增性能评估。见 [repro_v1_audit.json](repro_v1_audit.json)。
- **验证器错误注入检查通过**：在临时副本中向未解决事件填入公共p=1、改变事件种子、篡改artifact，独立审计均拒绝。其中前两项同时更新文件散列，确保拒绝来自语义/种子检查而非仅散列不一致。见 [auditor_mutation_checks.json](auditor_mutation_checks.json)。

全量覆盖合格的 95,621 个 event-condition 没有全部运行 BB bootstrap；完整推断本次限于固定332面板。这里的软件回归与数据导入正确性不构成误报率/FDR或生物学真实性证明。候选与未检出仅覆盖151个数值检验成功事件，不能代表全部254个覆盖合格面板事件。
