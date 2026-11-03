# BIOAGE-REL-01 筛选合同 v1.0（冻结版；v1.1 为终局）

> 题目方向（未立题）：常规血液生化"生物年龄"（Levine PhenoAge、KDM 生物年龄、稳态失调 HD）的短间隔重测可靠性，及其对死亡关联估计的回归稀释。
> 数据：NHANES III 第二次检查子样本（发现）+ NHANES 2001–2002 第二天采血子样本（独立复现）+ NHANES III 公开死亡链接文件。
> 执行方：Codex。全部为 NCHS 公开数据，无需签字，无 AI 处理限制条款。

---

## 1. 执行纪律

1. 本合同在看任何效应之前冻结，记录 sha256。修改规则须升版本号；**v1.1 为终局**，只允许修正变量名/文件名写错，不改任何科学阈值。v1.1 仍停即停，不出 v1.2。
2. 独立单位 = 人（SEQN）。所有 N 按人报告。
3. 工程错误（文件格式、列名、包版本）在同一版本内修，不算救题；修不好写 `NOT_RUN`，不写成阴性。
4. 失败分层报告：工程 / 数据不可执行 / 测量不成立 / 科学结果 / 增量不足，五者分开。
5. 禁止看结果后更换生物年龄算法、换人群子集、换阈值以求"更好看"。
6. 不依赖作者或 NCHS 回信。
7. 生物年龄算法一律用已发表系数（BioAge R 包 `phenoage_calc` / `kdm_calc` / `hd_calc`，或 Levine 2018 原文系数），**不在第二次检查数据上重新训练**。

## 2. 冻结对象

| 项 | 冻结内容 |
|---|---|
| 主指标 | PhenoAgeAccel = PhenoAge 对实足年龄回归的残差（残差模型只在第一次检查数据上拟合，再应用于两次检查） |
| 次指标 | KDM-BA Accel（BioAge 默认 NHANES III 参数，按性别）；HD（BioAge 默认参考人群） |
| 人群 | NHANES III 第二次检查中年龄 ≥20 岁、两次检查 9 项 PhenoAge 指标均完整者；排除第二次检查时怀孕者 |
| 可靠性估计量 | ICC(2,1)（双向随机、绝对一致）；个体内 SD（年）；MDC95 = 1.96 × √2 × 个体内 SD（年）；Bland–Altman 偏差与一致性界限 |
| 分类稳定性 | 第一次检查"加速前 20%"在第二次检查中仍在前 20% 的比例 |
| 死亡分析 | Cox：每 1 SD PhenoAgeAccel 的 HR，分别用检查1、检查2、两次平均；回归稀释比 λ 与校正 HR；两次平均相对单次的 C-index 增益 |

## 3. 竞争解释 → 对应分析

| 竞争解释 | 区分它的分析 |
|---|---|
| 不稳定主要来自空腹/时段不一致（葡萄糖） | 按两次检查时段是否相同（MXPSESSR vs MXRSESSR）及空腹时长分层报告 ICC/MDC；另算"葡萄糖固定为两次平均"的变体 |
| 不稳定来自单一成分（CRP 或 RDW） | 成分方差分解：逐一把某成分固定为两次平均，看 MDC 降幅；按系数 × 个体内 SD 的 delta 法贡献排序 |
| 第二次检查子样本更健康、不代表总体 | 比较第二次检查参与者 vs 未参与者的第一次检查 PhenoAgeAccel、年龄、种族、性别分布 |
| 两周内的真实急性变化（感染等） | 敏感性：剔除任一次 CRP ≥1.0 mg/dL 者 |
| NHANES III CRP 检测下限高（0.3 mg/dL）导致低估噪声 | 报告两次检查 CRP 低于检测下限的比例；在 2001–2002 复现中（CRP 更灵敏）比较 |
| 检查1 属于 PhenoAge 训练样本，死亡关联是样本内 | 主死亡分析以检查2（从未用于训练）为暴露；检查1 结果并列报告 |

## 4. 数据清单与桌面预检结果

| 数据 | 来源 | 准入 | 已核 / 待核 |
|---|---|---|---|
| NHANES III 第二次检查实验室 LABSE（及 LAB2SE、EXAMSE） | CDC NCHS NHANES III Series 11 No. 3A | 公开下载 | **已核**：子样本 2,603 人（<12 岁 212，12–19 岁 231，20–39 岁 809，40–59 岁 578，≥60 岁 773）；同一流程、同一 MEC 人员、尽量同一时段；重复检查变量名第三位由 p 换为 r；LAB2SE 含 MXPRDAYS（两次检查间隔天数，1,726 人有值，869 人空白）与 MXRSESSR。**待核（A1）**：LABSE 中 9 项指标的重复变量是否齐全及实际完成 N |
| NHANES III 第一次检查 LAB / EXAM / 成人家庭问卷 | 同上 | 公开 | 9 项指标（白蛋白、ALP、肌酐、血糖、CRP、WBC、淋巴细胞%、MCV、RDW）均在主实验室文件 |
| NHANES III 公开死亡链接文件（随访至 2019-12-31） | NCHS Linked Mortality Files | 公开 | **待核（A1）**：第二次检查成人的 ELIGSTAT 合格比例与死亡数 |
| NHANES 2001–2002 第二天采血：L11_2_B（CRP，16–69 岁）、L25_2_B（血常规）、L40_2_B（生化） | NHANES 2001–2002 | 公开 | **已核**：L11_2_B 与 L25_2_B 为公开文件，含 LB2DAY（间隔天数）。**待核（A1）**：L25_2_B 是否含淋巴细胞%、MCV、RDW；三文件交集完整 N。注意 NHANES 2000 第二天数据仅限 RDC，**不得使用** |
| 校准参照 | BioAge R 包内置 NHANES3 数据中的 `phenoage0`；NHANES III 第二次检查肾功能变异研究（PMC3628297）报告的血清肌酐 CVw | 公开 | — |

查重（桌面，已完成，12 条检索式以内）：
- 未见用 NHANES III 第二次检查或 NHANES 1999–2002 第二天采血评估 PhenoAge / KDM / HD 重测可靠性的论文 → **NONE**。
- 相邻工作（ADJACENT，不撞）：Lacher 等 NCHS 报告给出 1999–2002 单项指标个体内变异；NHANES III 第二次检查已用于肾功能单项指标与 2001–2002 的 CRP 单项变异；Add Health 第 6 波重复测量子研究（143 人、受限数据）只报单项指标 ICC。
- 执行时 Codex 须补查近 12 个月 bioRxiv/medRxiv："phenotypic age" + reliability/repeatability/test-retest/within-person；命中同数据同比较 → `STOP_DUP`。

## 5. 分层任务

### A1 数据与语义门（预算 ≤半天）
- 列出 LABSE 中 9 项指标重复变量的确切名称、单位、检测下限编码、填充值（8888 等）处理。
- 计算：≥20 岁、两次检查 9 项均完整的人数 N_main；其中死亡链接合格者及死亡数 D；两次时段相同的人数。
- 计算 2001–2002 复现的完整 N_rep（9 项均在第二天文件中有值）。
- **PASS**：9 项重复变量全部存在；N_main ≥ 600；D ≥ 300。
- **STOP**：任一 PhenoAge 指标无重复测量，或 N_main < 600 → `STOP_DATA`。
- 复现层：N_rep ≥ 200 → 复现可做；否则复现降级为描述性，不影响主判定。
- KDM/HD 所需额外变量（如收缩压、BUN、尿酸、HbA1c、总胆固醇）若缺重复测量，KDM/HD 记 `NOT_RUN`，不影响 PhenoAge 主判定。

### A2 校准门（预算 ≤半天）
- 在 NHANES III 第一次检查全样本上复算 PhenoAge，与 BioAge 内置 `phenoage0` 比较：**PASS** = 相关 ≥0.99 且平均绝对差 ≤0.5 年。
- 复现第二次检查血清肌酐 CVw：落在 PMC3628297 报告值 ±1.5 个百分点内。
- 任一不过 → `CALIBRATION_FAIL`，先修流程，不出主结果。

### B 主分析（预算 ≤1 天）
- 第 2 节全部可靠性估计量，含 bootstrap 95% CI（按人重抽 2,000 次）。
- 第 3 节全部竞争解释分析。
- 按年龄组（20–39 / 40–59 / ≥60）、性别分层报告。
- 精度门：ICC 的 95% CI 半宽 ≤0.03。不满足 → 记"精度不足"，判定为 STOP。

### C 死亡与回归稀释（预算 ≤半天）
- 第 2 节死亡分析；校正年龄、性别、种族/族裔；PH 假设检验。
- 报告：λ、校正前后 HR、单次 vs 两次平均的 C-index 差（bootstrap CI）。

### D 独立复现（仅当 N_rep ≥200）
- 在 2001–2002 重复 B 的主估计量（PhenoAgeAccel 的 ICC、个体内 SD、MDC95）。

**终局判定（C 结束即终局，只分 GO / STOP）**
- **GO**：A1、A2 通过，且 B 精度门通过。（方向不预设：可靠性高或低都可成文，结果只决定叙事与期刊档次。）
- **STOP**：A1/A2 未过、精度不足、或查重命中 `EXACT`。

## 6. 输出

```
BIOAGE-REL-01/
  contract_v1.0.md  (+ sha256)
  status.json        # {"layer": "...", "verdict": "PASS|STOP|NOT_RUN", "reason": "..."}
  A1_codebook_check.csv
  A1_counts.csv
  A2_calibration.csv
  B_reliability.csv
  B_components.csv
  B_sensitivity.csv
  C_mortality.csv
  D_replication.csv
  REPORT.md
```
REPORT.md 段落顺序：工程｜数据｜科学结果｜是否继续投入｜与合同的偏离。

## 7. 已知风险（执行前告知）

- 第二次检查为非随机便利样本，NCHS 明确不建议套用总体权重 → 全部分析不加权，并做第 3 节代表性比较。
- NHANES III CRP 检测下限 0.3 mg/dL，低于下限值被替换为下限/√2，ln(CRP) 在下限处堆积，会让 CRP 的个体内变异显得偏小。
- MXPRDAYS 有 869 人空白；间隔分析只用有值者，主 ICC 不依赖间隔。
- 检查1 属于 PhenoAge 训练样本；死亡主分析用检查2。

## 8. 参考

- NHANES III Second Exam File Documentation, Series 11 No. 3A (NCHS, 1999): https://wwwn.cdc.gov/Nchs/data/nhanes3/3a/LAB2SE-acc.pdf
- Levine ME et al. An epigenetic biomarker of aging for lifespan and healthspan. Aging 2018.
- Liu Z et al. A new aging measure captures morbidity and mortality risk across diverse subpopulations from NHANES IV. PLoS Med 2018;15:e1002718.
- Kwon D, Belsky DW. BioAge toolkit. GeroScience 2021;43:2795–2808.
- Lacher DA et al. Biological variation of laboratory analytes based on the 1999–2002 NHANES. NHSR No. 21, 2010.
- NHANES 2001–2002 L11_2_B / L25_2_B / L40_2_B second-day documentation.
- NHANES III second-exam kidney-marker variability study: PMC3628297.