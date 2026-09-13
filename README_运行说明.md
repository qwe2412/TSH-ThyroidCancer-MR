# TSH 与甲状腺癌 · 全流程重跑脚本包（2026-09-13 严格审查补跑封存版）

作者：王沛（新疆医科大学第一附属医院）
课题：遗传预测 TSH 水平与分化型甲状腺癌（两样本 MR + MVMR + 预后模型 + 免疫浸润 + 药敏预测）

本脚本包收录 2026 年 9 月推倒重来后**全部跑通并封存**的分析脚本。所有数值均来自真实运行输出，
与论文（`D:\gwas\副课题1论文_完整版_20260913.docx`）逐字对应。
2026-09-13 晚已按严格审查意见完成三项补跑：① MR 敏感性全量重跑（含 MR-PRESSO，脚本 A）；
② MVMR 条件检验（条件 F / 条件异质性 Q / Q 最小化稳健估计，脚本 11）；③ Lasso/弹性网重跑（脚本 10）。

---

## 运行顺序与依赖

| 顺序 | 脚本 | 用途 | 关键依赖 |
|---|---|---|---|
| 1 | `01_脚本A_MR主分析（修正版）.R` | 两样本 MR 全流程（99 SNP → clump → 协调 → 五法 MR → PRESSO → 复制结局） | TwoSampleMR, MRPRESSO, readxl；**需 OpenGWAS Token（写入 ~/.Renviron）** |
| 2 | `02_脚本B_MVMR（v3方法学重写）.R` | 多变量 MR 骨架脚本（旧版，仅供追溯） | TwoSampleMR, MendelianRandomization |
| 3 | `03_脚本C_预后模型（定稿第6版）.R` | TCGA-THY 预后模型（GNG7/NFIA） | survival, survivalROC, rms, glmnet |
| 4 | `04_脚本D_GEO验证（v5）.R` | 表达矩阵构建 + GSE33630 梯度验证 | 本地 RData |
| 5 | `05_脚本E_免疫浸润（v4）.R` | ssGSEA（28 免疫细胞）+ 组间差异 + 检查点 | GSVA 2.6.6, GSEABase |
| 6 | `06_脚本F_药敏分析（v3）.R` | GDSC2 药敏预测（198 药 IC50） | oncoPredict 1.3.1, sva 3.60.0 |
| 7 | `07_脚本EF_对账诊断.R` | E/F 与历史产物逐项对账 | 本地 RData |
| 8 | `08_脚本F_论文数值提取.R` | 论文 3.6/4.5 段数值提取 | drug_sensitivity_res.RData |
| 9 | `09_脚本H_论文出图.R` | 论文全部配图（真实数据重算） | survival, survivalROC, rms, pheatmap |
| 10 | `10_脚本I_Lasso弹性网重跑.R` | Lasso-Cox / 弹性网变量选择重跑（2026-09-13） | glmnet |
| 11 | `11_脚本B_MVMR条件检验_v4.R` | MVMR 条件检验：条件 F / 异质性 Q / Q 最小化（2026-09-13） | TwoSampleMR, **MVMR 0.4.8（GitHub 安装）**；需 OpenGWAS Token |

> 注意：脚本 A 与脚本 11 需要 OpenGWAS 在线 API（JWT Token 14 天有效）。若 Token 过期，可跳过
> 在线步骤、直接复用已封存结果（见下方"封存关键值"），其余脚本全部离线可跑。

## 数据文件依赖（均位于 D:\gwas\）

- `risk_model.RData`：dat（409×30，含 riskScore / DFS_time / DFS_event / riskGroup），High 204 / Low 205
- `expr_immune.RData`：expr_mat_use（409×20387）
- `tcga_thca.RData`：GDSC2_Expr（17419×805）/ GDSC2_Res（805×198）/ ic50_pred（历史，已作废）
- `immune_ssgsea.RData`：ss（28×409）/ res（28×6）/ grp / expr_aligned
- `drug_sensitivity_res.RData`：pred_ic50（409×198）/ drug_res（198×5）/ grp
- `geo\gse33630_ext.RData`：df_ext（105×4）
- `Charoentong2017.gmt`：28 免疫细胞基因集（E 脚本用）
- `41467_2020_17718_MOESM4_ESM.xlsx` / `MOESM5_ESM.xlsx`：Zhou 2020 TSH GWAS 补充数据（工具变量来源）

## 环境（本机实测）

R 4.6.1（ucrt）；survival 3.8-12、survivalROC 1.0.3.1、rms 8.1-1、ggplot2 4.0.3、
pheatmap 1.0.13、readxl 1.5.0、TwoSampleMR 0.7.9、MendelianRandomization 0.10.0、
oncoPredict 1.3.1、GSVA 2.6.6、MRPRESSO、sva 3.60.0、glmnet、**MVMR 0.4.8**（GitHub WSpiller/mvmr，
经代理 `curl.exe -x http://127.0.0.1:7897` 安装；Bioconductor 3.23 无 mvmr 包）。

## 关键坑位（勿重犯）

1. **GSVA 2.6 的 gsva() 是单参数 S4 泛型**：正确调用 `gsva(ssgseaParam(expr, gene_sets), verbose=FALSE)`。
   双参数 `gsva(expr, sets, method="ssgsea")` 报 "param='matrix'找不到继承方法"（E v2/v3 已踩）。
2. **calcPhenotype() 1.3.1 无 pValCutoff/seed 参数**（F v2 报"参数没有用"）；seed 用 `set.seed()` 外部设定；
   GDSC IC50 是 log 尺度，`trainingPtype` 需 `exp(GDSC2_Res)` 转回；microarray→RNA-seq 用 `batchCorrect="standardize"`。
3. **harmonise_data 不删除无效行**，而是标 mr_keep=FALSE。主分析前必须 `dat <- dat[dat$mr_keep == TRUE, ]`，
   否则 mr() 与 mr_presso() 的 SNP 数不一致（71 vs 69 问题）。
4. **MR-PRESSO 流程以失真检验为准**：2026-09-13 重跑 Outlier Test 检出 4 个潜在离群 SNP，
   但 Distortion Test P = 0.124 > 0.05 → **不剔除**，主分析保留全部 69 SNP（IVW 随机效应
   b = −0.7193693，P = 1.73×10⁻⁴）。早期 README 记录"剔除 4 SNP → 65 SNP IVW b = −0.547"为
   手动操作结果，与标准 MR-PRESSO 决策不同，论文一律以 2026-09-13 重跑为准。
5. **MVMR 包名大写**：`library(MVMR)`（GitHub WSpiller/mvmr 安装后 DESCRIPTION 包名是 MVMR）。
   条件 F 用 `strength_mvmr(r_input, gencov)`；异质性 Q 用 `pleiotropy_mvmr(r_input, gencov)`；
   Q 最小化用 `qhet_mvmr(F.data, pheno_cor, CI=FALSE, iterations=1000)`（第三参数是表型相关矩阵，非 gencov）。
6. **历史 `ic50_pred`（tcga_thca 内嵌）与历史 drug_sensitivity_res_backup 已作废**（9377 负值、24898 NA、
   最大 1.8e+32），论文一律引用 F v3 新预测。

## 封存关键值（论文引用基准，2026-09-13 重跑）

### MR 主分析（脚本 A，2026-09-13 全量重跑，69 有效 SNP）
- F 统计量：最小值 30.01、中位数 70.41（全部 > 10）
- IVW 随机效应：b = −0.7193693，SE = 0.1915475，P = 1.7295×10⁻⁴，OR = 0.487（95% CI 0.335–0.709）
- MR-Egger：OR = 0.467，P = 0.0617；加权中位数：OR = 0.626，P = 0.0032；简单众数：OR = 0.405，P = 0.0115；
  加权众数：OR = 0.702，P = 0.0913
- Cochran Q：IVW Q = 279.3924，df = 68，P = 1.990×10⁻²⁷；Egger Q = 279.3341，df = 67，P = 9.885×10⁻²⁸
- MR-Egger 截距：0.0028297，SE = 0.0239，P = 0.906（无显著定向多效性）
- 留一法：b ∈ [−0.8324, −0.5873]（剔除任一 SNP 方向不变）
- MR-PRESSO（10000 次置换）：Global RSS = 289.7744，P < 1×10⁻⁴；Outlier 检出 4 个潜在离群（P < 0.05）；
  Distortion 系数 = −31.5%，P = 0.1238（> 0.05 → 保留全部 69 SNP）

### 复制结局
- 复制① ebi-a-GCST90013867（SPA 校正甲状腺癌 GWAS）：65 有效 SNP，IVW OR = 0.362（95% CI 0.233–0.562），
  P = 6.10×10⁻⁶
- 复制② finn-b-C3_THYROID_GLAND（FinnGen）：68 有效 SNP，IVW OR = 0.559（95% CI 0.384–0.816），P = 0.0025

### MVMR 条件检验（脚本 11，2026-09-13，54 有效 SNP）
- 条件 F（gencov = 0 / phenocov 近似）：TSH 32.74 / 31.81；BMI 5.66 / 5.68；睾酮 3.38 / 3.38
  （TSH 合格 > 10；BMI、睾酮 < 10，弱工具）
- MVMR-IVW：TSH b = −0.6446，SE = 0.2575，OR = 0.525（95% CI 0.317–0.869），P = 0.0156（t 检验）；
  BMI b = −0.5252，SE = 2.9446，P = 0.859；睾酮 b = 2.3938，SE = 4.6987，P = 0.613
- 条件异质性 Q（pleiotropy_mvmr）：Q = 244.76（df = 50），P = 1.81×10⁻²⁷（gencov=0）；
  Q = 244.63，P = 1.90×10⁻²⁷（phenocov）
- Q 统计量最小化（qhet_mvmr，1000 次）：TSH b = −0.6277；BMI b = −0.7455；睾酮 b = 6.3992
- 表型相关近似（phenocov_mvmr 输入）：TSH–BMI r = 0.05，TSH–睾酮 r = 0，BMI–睾酮 r = −0.10（文献近似，已注明）

### Lasso/弹性网（脚本 10，2026-09-13）
- Lasso-Cox（10 折 CV，seed 20260913，缺失值中位数填补）：lambda.min 保留 43 个变量（含 GNG7、NFIA，负向）；
  lambda.1se = 0 变量；弹性网（α = 0.5）两条路径均 0 变量 → 正则化变量选择结果不稳定
- 产物：`D:\gwas\lasso_rerun_20260913.RData`（cv1 / cv05）

### 预后模型（脚本 C 定稿第 6 版）
- riskScore = −0.0053257 × GNG7 − 0.0003788 × NFIA；C-index 0.666（95% CI 0.559–0.774，
  concordance 正态近似；Bootstrap 1000 乐观校正 0.662）；Log-rank P = 0.0412；
  时间 ROC AUC（KM）1/3/5 年 = 0.657 / 0.612 / 0.593；3 年校准 MAE = 0.037（rms::calibrate，B = 1000，
  seed 20260913）；DCA 正净获益阈值 0.01、0.02、0.04–0.08
- 诚实性红线：多因素 Cox 中 GNG7 P = 0.0717、NFIA P = 0.2413，均不显著，论文表述为"探索性"

### GEO（GSE33630，105 例）
- KW χ² = 54.147，P = 1.75×10⁻¹²；两两 P：ATC vs PTC 7.0×10⁻¹¹、ATC vs Normal 1.3×10⁻¹¹、
  PTC vs Normal 3.2×10⁻⁹；Cohen's d：3.87 / 4.37 / 1.16（95% CI 见论文）

### 免疫浸润（E v4）
- 25/28 免疫细胞 FDR < 0.05；检查点基因 9/15 显著（FDR 基于 15 基因 BH 校正：
  CD274 3.3×10⁻¹⁴、PDCD1LG2 2.8×10⁻⁹、SIGLEC15 4.1×10⁻¹¹、CD80 4.1×10⁻¹¹、TNFRSF9 1.0×10⁻⁷、
  CTLA4 1.9×10⁻⁷、TIGIT 2.9×10⁻⁷、CD86 3.7×10⁻⁷、HAVCR2 1.9×10⁻⁶）

### 药敏（F v3）
- FDR < 0.05 药物 176/198（88.9%），High 组耐药 158、敏感 18；代表药物数值见论文 3.7/4.5 段
- 全部来自 F v3 输出（`drug_sensitivity_res_full.csv`，198 行 × 8 列）

## 论文配图（D:\gwas\paper_figs\）

图 1–11 编号连续：图 1/2 MR（scatter/forest）、图 3–6 预后（KM/ROC/校准/DCA）、图 7 GEO、
图 8/9 免疫（热图/箱线图）、图 10/11 药敏（火山图/箱线图）。均由脚本 A/09 用真实数据重算生成。

## 与旧稿的区别（诚实声明）

- 论文 F 段（3.6/4.5）方向**已反转**：旧稿写"高风险对细胞毒药更敏感"，F v3 真实结果为
  **高风险组更耐药（158/176）**，论文已按真实输出改写。
- 校准 MAE：旧稿 0.036（B=200）→ 封存 0.037（rms::calibrate B=1000 真实运行值）。
- C-index CI：旧稿 0.553–0.772（来源不明）→ 实测 concordance 正态 0.559–0.774。
- MVMR：旧稿"Q = 247.92、df = 51"无脚本证据 → 2026-09-13 用 pleiotropy_mvmr 重算
  Q = 244.76、df = 50（54 SNP），并补条件 F 与 Q 最小化稳健估计。
- MR 主分析：旧稿"剔除 4 SNP 后 65 SNP IVW b = −0.547"为手动操作 → 2026-09-13 按
  MR-PRESSO 失真检验（P = 0.124）保留 69 SNP，IVW b = −0.719（P = 1.73×10⁻⁴）。
