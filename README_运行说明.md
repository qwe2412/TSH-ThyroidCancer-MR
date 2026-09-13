# TSH 与甲状腺癌 · 全流程重跑脚本包（2026-09-13 封存版）

作者：王沛（新疆医科大学第一附属医院）
课题：遗传预测 TSH 水平与分化型甲状腺癌（两样本 MR + MVMR + 预后模型 + 免疫浸润 + 药敏预测）

本脚本包收录 2026 年 9 月推倒重来后**全部跑通并封存**的分析脚本。所有数值均来自真实运行输出，
与论文（`D:\gwas\副课题1论文_结果与讨论.docx`）逐字对应。

---

## 运行顺序与依赖

| 顺序 | 脚本 | 用途 | 关键依赖 |
|---|---|---|---|
| 1 | `01_脚本A_MR主分析（修正版）.R` | 两样本 MR 全流程（99 SNP → clump → 协调 → 五法 MR → PRESSO → 复制结局） | TwoSampleMR, MRPRESSO, readxl；**需 OpenGWAS Token（写入 ~/.Renviron）** |
| 2 | `02_脚本B_MVMR（v3方法学重写）.R` | 多变量 MR（TSH 70 SNP 骨架 + BMI + 睾酮） | TwoSampleMR, MendelianRandomization |
| 3 | `03_脚本C_预后模型（定稿第6版）.R` | TCGA-THY 预后模型（GNG7/NFIA） | survival, survivalROC, rms, glmnet |
| 4 | `04_脚本D_GEO验证（v5）.R` | 表达矩阵构建 + GSE33630 梯度验证 | 本地 RData |
| 5 | `05_脚本E_免疫浸润（v4）.R` | ssGSEA（28 免疫细胞）+ 组间差异 + 检查点 | GSVA 2.6.6, GSEABase |
| 6 | `06_脚本F_药敏分析（v3）.R` | GDSC2 药敏预测（198 药 IC50） | oncoPredict 1.3.1, sva 3.60.0 |
| 7 | `07_脚本EF_对账诊断.R` | E/F 与历史产物逐项对账 | 本地 RData |
| 8 | `08_脚本F_论文数值提取.R` | 论文 3.6/4.5 段数值提取 | drug_sensitivity_res.RData |
| 9 | `09_脚本H_论文出图.R` | 论文全部配图（真实数据重算） | survival, survivalROC, rms, pheatmap |

> 注意：脚本 A 需要 OpenGWAS 在线 API（JWT Token 14 天有效）。若 Token 过期，可跳过 A、
> 直接复用已封存结果（见下方"封存关键值"），其余脚本全部离线可跑。

## 数据文件依赖（均位于 D:\gwas\）

- `risk_model.RData`：dat（409×30，含 riskScore / DFS_time / DFS_event / riskGroup），High 204 / Low 205
- `expr_immune.RData`：expr_mat_use（409×20387）
- `tcga_thca.RData`：GDSC2_Expr（17419×805）/ GDSC2_Res（805×198）/ ic50_pred（历史，已作废）
- `immune_ssgsea.RData`：ss（28×409）/ res（28×6）/ grp / expr_aligned
- `drug_sensitivity_res.RData`：pred_ic50（409×198）/ drug_res（198×5）/ grp
- `geo\gse33630_ext.RData`：df_ext（105×4）
- `Charoentong2017.gmt`：28 免疫细胞基因集（E 脚本用）

## 环境（本机实测）

R 4.6.1（ucrt）；survival 3.8-12、survivalROC 1.0.3.1、rms 8.1-1、ggplot2 4.0.3、
pheatmap 1.0.13、readxl 1.5.0、TwoSampleMR 0.7.9、MendelianRandomization 0.10.0、
oncoPredict 1.3.1、GSVA 2.6.6、MRPRESSO、sva 3.60.0、glmnet。

## 关键坑位（勿重犯）

1. **GSVA 2.6 的 gsva() 是单参数 S4 泛型**：正确调用 `gsva(ssgseaParam(expr, gene_sets), verbose=FALSE)`。
   双参数 `gsva(expr, sets, method="ssgsea")` 报 "param='matrix'找不到继承方法"（E v2/v3 已踩）。
2. **calcPhenotype() 1.3.1 无 pValCutoff/seed 参数**（F v2 报"参数没有用"）；seed 用 `set.seed()` 外部设定；
   GDSC IC50 是 log 尺度，`trainingPtype` 需 `exp(GDSC2_Res)` 转回；microarray→RNA-seq 用 `batchCorrect="standardize"`。
3. **harmonise_data 不删除无效行**，而是标 mr_keep=FALSE。主分析前必须 `dat <- dat[dat$mr_keep == TRUE, ]`，
   否则 mr() 与 mr_presso() 的 SNP 数不一致（71 vs 69 问题）。
4. **A 脚本第 12 步"剔除离群后 IVW"勿引 69 SNP 输出（b=-0.7193693）**；正确值为 Distortion Test 段
   剔除 4 SNP（rs10186921/rs116909374/rs2993047/rs925488）后 65 SNP 的随机效应 IVW。
5. **MVMR 旧版"三暴露工具并集"交集=0**，已废弃；v3 以 TSH 70 工具 SNP 为骨架提取 BMI/睾酮效应。
6. **历史 `ic50_pred`（tcga_thca 内嵌）与历史 drug_sensitivity_res_backup 已作废**（9377 负值、24898 NA、
   最大 1.8e+32），论文一律引用 F v3 新预测。

## 封存关键值（论文引用基准，勿用历史替代）

- **MR 主分析**：剔除 4 SNP 后 65 SNP，IVW 随机效应 b = −0.5469521，P = 1.44×10⁻⁵；
  主分析五法 IVW OR = 0.487（95% CI 0.335–0.709）；Egger 截距 P = 0.906（无水平多效性）；
  复制① GCST90013867 OR = 0.3621605，P = 6.10×10⁻⁶；复制② finn-b OR = 0.5594148，P = 0.00254（68 SNP）。
- **MVMR**：54 SNP，IVW OR = 0.5248554（95% CI 0.3168282–0.8694719），P = 0.0123106，b = −0.6446324；
  异质性 Q = 247.9218，df = 51。
- **预后模型**：riskScore = −0.0053257 × GNG7 − 0.0003788 × NFIA；
  C-index 0.666（95% CI 0.553–0.772，Bootstrap 1000，乐观校正 0.663）；Log-rank P = 0.0412；
  时间 ROC AUC（KM）1/3/5 年 = 0.657 / 0.612 / 0.593；3 年校准 MAE = 0.036（rms::calibrate，B=200，seed 20260913）；
  诚实性红线：mfinal 多因素 GNG7 P = 0.0717、NFIA P = 0.2413，均不显著，论文表述为"探索性"。
- **GEO（GSE33630，105 例）**：KW χ² = 54.147，P = 1.75×10⁻¹²；Cohen's d：ATC vs PTC 3.87、
  ATC vs Normal 4.37、PTC vs Normal 1.16。
- **免疫浸润（E v4）**：25/28 免疫细胞 FDR < 0.05；检查点基因 9/15 显著
  （CD274 2.22×10⁻¹⁵、PDCD1LG2 7.50×10⁻¹⁰、SIGLEC15 8.25×10⁻¹²、CD80 6.52×10⁻¹²、
  TNFRSF9 3.40×10⁻⁸、CTLA4 7.56×10⁻⁸、TIGIT 1.33×10⁻⁷、CD86 2.00×10⁻⁷、HAVCR2 1.16×10⁻⁶）。
- **药敏（F v3）**：FDR < 0.05 药物 176/198（88.9%），其中 High 组耐药 158、敏感 18；
  代表药物数值见论文 3.6/4.5 段（全部来自 F v3 输出）。

## 论文配图（D:\gwas\paper_figs\）

fig2_km / fig3_roc / fig4_calibration / fig5_dca / fig6_geo / fig7_immune_heatmap /
fig8_immune_boxplot / fig9_drug_volcano / fig10_drug_boxplot（均由脚本 09 用真实 RData 重算生成）；
fig11–14 为 MR 四图（forest/scatter/funnel/leaveoneout，脚本 A 在线运行产物）。

## 与旧稿的区别（诚实声明）

- 论文 F 段（3.6/4.5/4.8）方向**已反转**：旧稿写"高风险对细胞毒药更敏感"，F v3 真实结果为
  **高风险组更耐药（158/176）**，论文已按真实输出改写。
- 校准 MAE：旧稿 0.038 → 本包封存 0.036（rms::calibrate 真实运行值）。
- MVMR：旧稿三处误用复制② 0.559 → 本包封存 0.525（论文已同步修正）。
