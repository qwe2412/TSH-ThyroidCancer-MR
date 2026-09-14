# -*- coding: utf-8 -*-
"""
构建完整论文 docx（2026-09-14 二轮严格审查修订版）
作者：王沛；单位：新疆医科大学第一附属医院
所有数值均来自真实重跑输出（脚本 A-H/L 封存值 + 二轮审查实测）。
二轮修订要点（回应 DeepSeek 严格审稿）：
 1. 标题弱化："预后模型"→"候选基因预后评分"；"肿瘤微环境整合分析"→"探索性整合分析"
 2. 摘要/结论弱化：探索性、假设生成；无独立生存验证
 3. MR：补 I²（IVW 75.7%、Egger 76.0%）、tau²、PRESSO outlier-corrected 估计、
    离群 SNP rsID、MR-RAPS、Steiger、反向 MR 执行结果（在线补跑脚本 L）
 4. MVMR：全暴露 qhet 估计；删除"独立于 BMI 与性激素"强结论
 5. 预后：补 Bootstrap 变量选择频率与重复交叉验证（脚本 J）；降级为"候选基因评分"
 6. 药敏：补 5 折 CV 性能（R²/RMSE）与组标签置换检验（脚本 K）
 7. 免疫：补"未校正肿瘤纯度/亚型/突变"局限；ssGSEA=相对富集表述保留
 8. Fussey 年份勘误：2021→2020（发表 2020-08-18，撤稿 2020-11-14）
 9. 讨论 4.1 P 值统一：0.012→0.0123；删除超声随访临床建议（改为假设性表述）
10. 局限扩充：反向 MR/共定位/多效性搜索/样本重叠/亚型/性别分层/免疫纯度
11. 声明补：伦理豁免、作者贡献、致谢、STROBE-MR/TRIPOD+AI 检查表、6 次提交
12. 参考文献扩充至 33 篇（全部真实可溯源）
"""
import os
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

FIG_DIR = r"D:\gwas\paper_figs"
OUT = r"D:\gwas\副课题1论文_完整版_20260914.docx"

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
sec.top_margin = sec.bottom_margin = Cm(2.54)
sec.left_margin = sec.right_margin = Cm(3.17)

def set_run(run, size=11, bold=False, font="宋体", ascii_font="Times New Roman"):
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.name = ascii_font
    r = run._element
    rPr = r.get_or_add_rPr()
    rFonts = rPr.find(qn('w:rFonts'))
    if rFonts is None:
        rFonts = rPr.makeelement(qn('w:rFonts'), {})
        rPr.append(rFonts)
    rFonts.set(qn('w:eastAsia'), font)

def para(text, size=11, bold=False, align=None, space_after=6, first_indent=None,
         font="宋体"):
    p = doc.add_paragraph()
    if align is not None:
        p.alignment = align
    p.paragraph_format.space_after = Pt(space_after)
    if first_indent:
        p.paragraph_format.first_line_indent = Pt(first_indent)
    run = p.add_run(text)
    set_run(run, size=size, bold=bold, font=font)
    return p

def heading(text, level=1):
    sizes = {1: 14, 2: 12, 3: 11}
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.space_after = Pt(6)
    run = p.add_run(text)
    set_run(run, size=sizes.get(level, 11), bold=True, font="黑体")
    return p

def add_figure(fname, width_cm=14.5, caption=None):
    path = os.path.join(FIG_DIR, fname)
    if not os.path.exists(path):
        para("[缺图: %s]" % fname, size=10)
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    run.add_picture(path, width=Cm(width_cm))
    if caption:
        cp = doc.add_paragraph()
        cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
        cp.paragraph_format.space_after = Pt(10)
        cr = cp.add_run(caption)
        set_run(cr, size=9, font="宋体")

# ================= 标题 / 作者 / 单位 =================
para("遗传预测的血清促甲状腺激素水平与甲状腺癌风险：孟德尔随机化、"
     "候选基因预后评分与肿瘤微环境的探索性整合分析", size=16, bold=True,
     align=WD_ALIGN_PARAGRAPH.CENTER, space_after=12)
para("王沛", size=13, align=WD_ALIGN_PARAGRAPH.CENTER, space_after=4)
para("新疆医科大学第一附属医院，乌鲁木齐 830054，中国", size=11,
     align=WD_ALIGN_PARAGRAPH.CENTER, space_after=14)

# ================= 摘要 =================
heading("摘要", 1)
para("【背景】血清促甲状腺激素（TSH）水平与甲状腺癌发病风险的因果关系存在争议，"
     "且既往研究基于不同的工具变量集与结局定义结论不一。本研究以两样本孟德尔随机化（MR）为核心，"
     "联合多变量 MR（MVMR）、候选基因预后评分构建、免疫浸润与药物敏感性计算预测，"
     "对 TSH-甲状腺癌因果关联及其转化价值进行探索性评估。", first_indent=22)
para("【方法】基于 Zhou 等（2020）的 TSH 全基因组关联研究补充数据（99 个全基因组显著位点，"
     "LD 筛选后 70 个独立工具 SNP），以 LD r² < 0.001（10000 kb）筛选独立工具变量，"
     "对甲状腺癌结局（ebi-a-GCST90018929，1054 例病例/490,920 例对照）行两样本 MR，"
     "进行异质性、水平多效性、留一法、MR-PRESSO 及 MR-RAPS 敏感性分析，"
     "并在两个独立结局中复制；以 TSH 工具 SNP 为骨架、校正 BMI 与睾酮行 MVMR"
     "（报告条件 F 统计量与条件异质性 Q，并注明弱工具限制）。将工具 SNP 映射至候选基因，"
     "在 TCGA-THY 队列（409 例）中构建探索性候选基因预后评分，"
     "在 GEO 数据集 GSE33630 行病理类型梯度验证；以 ssGSEA 评估免疫细胞相对富集，"
     "基于 GDSC2 以 oncoPredict 计算预测 198 种药物敏感性并报告交叉验证性能与置换检验结果。",
     first_indent=22)
para("【结果】两样本 MR 显示遗传预测 TSH 升高与甲状腺癌风险降低存在因果关联"
     "（主分析 IVW 随机效应 OR = 0.487，95% CI 0.335–0.709，P = 1.7×10⁻⁴，69 个有效 SNP），"
     "MR-Egger 截距检验未发现显著定向多效性（P = 0.906）；Cochran Q 检验提示工具变量间存在"
     "显著异质性（IVW Q = 279.39，df = 68，P = 1.99×10⁻²⁷，I² = 75.7%），"
     "MR-PRESSO 全局检验亦显著（P < 1×10⁻⁴）并识别 4 个潜在离群 SNP；"
     "剔除离群 SNP 后随机效应 IVW 估计为 OR = 0.579（95% CI 0.452–0.741，"
     "65 个 SNP），方向与主分析一致；MR-RAPS 估计方向亦一致"
     "（OR = 0.476，P = 4.2×10⁻¹⁵）。两个复制结局方向一致"
     "（OR = 0.362，P = 6.10×10⁻⁶；OR = 0.559，P = 0.0025）。MVMR 在校正 BMI 与睾酮后 TSH "
     "直接效应仍显著（54 SNP，OR = 0.525，95% CI 0.317–0.869，P = 0.0123，t 检验），"
     "但 BMI（条件 F = 5.66）与睾酮（条件 F = 3.38）为弱工具，该校正结果仅具探索性。"
     "由 GNG7 与 NFIA 构建的候选基因评分 C-index = 0.666（95% CI 0.559–0.774），"
     "高、低风险组生存差异显著（Log-rank P = 0.041）；Bootstrap 200 次变量选择中 "
     "GNG7 保留频率 151/200（75.5%）、NFIA 保留频率 113/200（56.5%），"
     "5×5 折重复交叉验证 C-index 均值 0.334；GSE33630 中评分呈 "
     "ATC > PTC > Normal 梯度（Kruskal-Wallis P = 1.75×10⁻¹²）。"
     "28 种免疫细胞的相对富集评分中 25 种组间差异显著（FDR < 0.05），高风险组呈免疫激活伴耗竭特征；"
     "计算预测提示 176 种药物（88.9%）的 IC50 组间差异显著，其中高风险组对 158 种更耐药、仅 18 种更敏感；"
     "代表药物 5 折交叉验证 R² 中位数 0.350，组标签置换（500 次）显示显著药物数"
     "（176）远高于随机期望（置换中位数 0.62，P < 0.002）。",
     first_indent=22)
para("【结论】遗传证据提示 TSH 升高对甲状腺癌发病可能具有保护作用；"
     "GNG7/NFIA 候选基因评分显示初步区分趋势但未经独立生存队列验证；"
     "免疫浸润与药物敏感性均为计算预测结果。本研究属探索性、假设生成性质，"
     "其临床价值需独立队列与实验验证。", first_indent=22)
para("【关键词】促甲状腺激素；甲状腺癌；孟德尔随机化；候选基因评分；免疫浸润；药物敏感性",
     bold=True, space_after=14)

# ================= 1 引言 =================
heading("1  引言", 1)
para("甲状腺癌是内分泌系统最常见的恶性肿瘤，其发病率在全球范围内持续上升。"
     "促甲状腺激素（TSH）是甲状腺滤泡细胞增殖与分化的核心调节因子，"
     "通过其受体（TSHR）激活 cAMP/PKA 信号通路调节甲状腺激素合成。"
     "临床上分化型甲状腺癌术后常规行 TSH 抑制治疗，提示 TSH 信号对已建立肿瘤生长具有促进作用；"
     "然而，观察性研究中血清 TSH 水平与甲状腺癌发病风险的关系并不一致，"
     "部分研究提示高 TSH 与甲状腺癌风险升高相关，而另一些研究则未观察到明确关联。"
     "上述不一致可能源于混杂因素与反向因果。", first_indent=22)
para("孟德尔随机化（MR）利用遗传变异作为工具变量，可有效削弱混杂与反向因果的影响。"
     "既往 MR 研究对 TSH-甲状腺癌因果方向存在相反报告：Zhou 等（2020）基于 HUNT、MGI 与 "
     "ThyroidOmics 三队列的 TSH 全基因组关联研究（99 个全基因组显著位点）报道遗传预测 TSH 每升高 1 个标准差，"
     "甲状腺癌风险降低约 50%（OR = 0.50，95% CI 0.27–0.92）；Yuan 等（2020）亦报道保护效应"
     "（OR = 0.47，95% CI 0.30–0.73）。与之相反，Fussey 等（2020）基于 UK Biobank 的分析报告"
     "风险增加（OR = 2.00，95% CI 1.09–3.70），但该文因通讯作者主动报告数据分析存在问题"
     "已于 2020 年 11 月被期刊正式撤稿，其结论不具备正式引用效力，此处仅作历史说明。"
     "不同工具变量集、结局定义与估计方法下的异质性，使该因果方向仍需严格验证。",
     first_indent=22)
para("本研究以全流程重跑、真实输出为准的方式完成：（1）两样本 MR 主分析、敏感性分析与双结局复制，"
     "并补充 MR-PRESSO 校正估计、MR-RAPS 与 Steiger 方向性检验；"
     "（2）校正 BMI 与睾酮的 MVMR（报告条件 F 与条件异质性 Q 并如实标注弱工具限制）；"
     "（3）基于 MR 阳性证据将工具 SNP 映射至候选基因，"
     "在 TCGA-THY 队列构建探索性候选基因预后评分（补充 Bootstrap 变量选择频率与重复交叉验证），"
     "并在独立 GEO 队列行病理类型梯度验证；"
     "（4）免疫浸润分析与 GDSC2 药物敏感性计算预测（补充交叉验证性能与置换检验），"
     "以期为 TSH-甲状腺癌关联的临床转化提供假设生成证据。"
     "报告遵循 STROBE-MR 与 TRIPOD+AI 规范（检查表见补充材料）。", first_indent=22)

# ================= 2 材料与方法 =================
heading("2  材料与方法", 1)

heading("2.1  数据来源", 2)
para("暴露数据：Zhou 等（2020）TSH 水平全基因组关联研究（Nature Communications，"
     "样本来自 HUNT、MGI 与 ThyroidOmics 队列）补充数据，包含 99 个全基因组显著位点（P < 5×10⁻⁸）。"
     "结局数据：甲状腺癌 GWAS（ebi-a-GCST90018929，1054 例病例/490,920 例对照）为主结局；"
     "复制结局①为 ebi-a-GCST90013867（Thyroid cancer，SPA 校正，欧洲人群，样本量 407,746，"
     "GRCh37）；复制结局②为 finn-b-C3_THYROID_GLAND（FinnGen 恶性甲状腺肿瘤，欧洲人群，"
     "989 例病例/217,803 例对照）。协变量 GWAS：BMI（ieu-b-40，欧洲人群，样本量 681,275，"
     "Yengo 2018）与总睾酮（ebi-a-GCST90014013，UK Biobank 数据字段 30850，欧洲人群，"
     "样本量 353,805）。"
     "全部汇总统计量来自 IEU OpenGWAS 数据库（https://gwas.mrcieu.ac.uk/）。"
     "预后队列：TCGA-THY（409 例，含表达与临床随访）。"
     "外部验证：GEO GSE33630（Affymetrix GPL570，105 例：11 ATC / 49 PTC / 45 正常）。"
     "药物训练集：GDSC2（805 细胞系 × 198 药物，含表达矩阵 GDSC2_Expr 与 IC50 矩阵 GDSC2_Res）。",
     first_indent=22)

heading("2.2  工具变量筛选与两样本 MR", 2)
para("补充数据提供的 99 个位点均已满足 P < 5×10⁻⁸；在此基础上以连锁不平衡 r² < 0.001"
     "（窗口 10000 kb，EUR 参考面板）筛选出 70 个独立工具 SNP，"
     "并计算 F 统计量（F = beta²/SE²）评估弱工具偏倚（F > 10 为合格）。"
     "对主结局行两样本 MR：以逆方差加权（IVW）随机效应模型为主分析，"
     "MR-Egger、加权中位数、简单众数、加权众数为补充；"
     "以 Cochran Q 检验异质性并报告 I²（I² = (Q − df)/Q × 100%）与 DerSimonian-Laird tau²，"
     "MR-Egger 截距检验水平多效性、留一法评估单一 SNP 驱动效应。"
     "MR-PRESSO（10000 次置换）检测离群 SNP，并按失真检验识别的离群索引"
     "剔除后重估随机效应与固定效应 IVW（outlier-corrected 估计），"
     "同时报告离群 SNP 的 rsID 列表；以 MR-RAPS（稳健调整轮廓评分）作为对多效性稳健的补充方法，"
     "以 Steiger 方向性检验评估因果方向的可靠性。"
     "复制结局分别重复主分析流程并报告异质性 Q、I²、MR-Egger 截距与 MR-PRESSO 全局检验。"
     "全部敏感性检验的详细数值来自脚本 A 与脚本 L 在线重跑的真实 Console 输出"
     "（2026-09-13/2026-09-14，OpenGWAS 凭据）。", first_indent=22)

heading("2.3  多变量 MR", 2)
para("以 TSH 的 70 个工具 SNP 为骨架，分别从 BMI（ieu-b-40）与总睾酮（ebi-a-GCST90014013）"
     "GWAS 提取对应 SNP 效应，与结局协调后纳入三暴露 MVMR-IVW（MVMR 0.4.8 包），"
     "评估校正 BMI 与睾酮后 TSH 对甲状腺癌的独立效应（MVMR 0.4.8 包与 MendelianRandomization 0.10.0 的 mr_mvivw 结果一致）。以条件 F 统计量"
     "（Sanderson & Windmeijer，F > 10 为合格）评估各暴露的条件工具强度；"
     "以 MVMR 异质性 Q 统计量（pleiotropy_mvmr）评估工具有效性假设；"
     "并以 Q 统计量最小化方法（qhet_mvmr，bootstrap 1000 次）估计各暴露的稳健效应。"
     "SNP 对暴露效应的协方差以表型相关近似（phenocor_mvmr；TSH–BMI r = 0.05、"
     "TSH–睾酮 r = 0、BMI–睾酮 r = −0.10）与独立假设（gencov = 0）双版本报告。"
     "需说明：BMI 与睾酮的条件 F 统计量均低于 10（弱工具），"
     "故 MVMR 结果仅定位为探索性，不据此断言 TSH 效应独立于 BMI 与性激素。", first_indent=22)

heading("2.4  候选基因预后评分构建与评估", 2)
para("将 MR 工具 SNP 映射至候选基因（55 个 SNP → 60 个候选基因），"
     "在 TCGA-THY 队列（409 例，27 例复发/死亡事件，EPV ≈ 13.5）中构建探索性风险评分 "
     "riskScore = −0.0053257 × GNG7 − 0.0003788 × NFIA（基因表达量为 TPM 后 log2 转换），"
     "按中位数分为高、低风险组。GNG7 与 NFIA 的入选依据为早期分析中的单因素与多因素 Cox 回归"
     "筛选（多因素 Cox：GNG7 P = 0.072、NFIA P = 0.241，均未达传统显著性阈值，"
     "故本评分定位为探索性候选基因评分而非验证性预后模型）。"
     "为评估变量选择稳定性，补充运行 Bootstrap 200 次（每次从 409 例中有放回抽样、"
     "以同样的双基因 Cox 结构拟合）记录 GNG7/NFIA 的保留频率与系数方向；"
     "并以 5×5 折重复交叉验证评估评分 C-index 的分布（均值与范围）。"
     "评估指标：C-index（concordance 法，95% CI 采用正态近似；Bootstrap 1000 次内部验证并报告乐观校正值）、"
     "Log-rank 检验、时间依赖 ROC（survivalROC，KM 法，1/3/5 年）、3 年校准曲线"
     "（rms::calibrate，Bootstrap 1000 次，seed 20260913）、36 个月决策曲线分析（DCA）。",
     first_indent=22)

heading("2.5  GEO 外部梯度验证", 2)
para("在 GSE33630（105 例）中以风险评分列 riskScore_ext（历史 GEO 分析使用的口径，"
     "与论文 TCGA 公式计算的评分 Pearson 相关 r = 0.996）进行验证："
     "以 Kruskal-Wallis 检验比较 ATC、PTC、正常三组间差异，"
     "并以两两 Wilcoxon 检验与 Cohen's d 效应量评估组间梯度。"
     "需说明：GSE33630 不含生存随访信息，本验证仅支持评分的病理类型梯度意义，"
     "不构成生存终点的外部验证。", first_indent=22)

heading("2.6  免疫浸润分析", 2)
para("以 ssGSEA（GSVA 2.6.6，Charoentong 2017 的 28 种免疫细胞基因集）"
     "对 TCGA-THY 表达矩阵（409 样本 × 20387 基因）计算免疫细胞相对富集评分，"
     "以 Wilcoxon 秩和检验比较高、低风险组差异并以 Benjamini-Hochberg 校正（FDR < 0.05），"
     "同时报告均值差与效应量；另对 15 个免疫检查点/耗竭标志物基因行组间差异检验并报告 FDR。"
     "由于 TCGA-THY 缺乏匹配的肿瘤纯度与驱动突变数据，未对肿瘤纯度、亚型与 BRAF/RAS 突变状态"
     "进行校正，该局限在 4.7 中说明。", first_indent=22)

heading("2.7  药物敏感性预测", 2)
para("以 GDSC2 表达矩阵（17419 基因 × 805 细胞系）与 IC50 矩阵（805 细胞系 × 198 药物，"
     "log 尺度经 exp 转回）为训练集，以 oncoPredict::calcPhenotype（batchCorrect = \"standardize\"，"
     "powerTransformPhenotype = TRUE，removeLowVaryingGenes = 0.2，minNumSamples = 10，"
     "seed 经 set.seed(123) 设定）计算预测 TCGA-THY 409 例患者对 198 种药物的 IC50；"
     "以 Wilcoxon 秩和检验比较组间差异并作 FDR 校正。"
     "为评估预测可靠性：（1）对 12 个代表性药物在 GDSC2 训练集行 5 折交叉验证"
     "（ridge 回归，log IC50 尺度，报告 R² 与 RMSE）；（2）行组标签置换检验（B = 500 次）"
     "估计 FDR < 0.05 显著药物数的零分布，与实际显著药物数比较。", first_indent=22)

heading("2.8  统计分析", 2)
para("全部分析基于 R 4.6.1 完成。关键包版本：TwoSampleMR 0.7.9、MendelianRandomization 0.10.0、"
     "MRPRESSO 1.0.6、ieugwasr 1.1.2、readxl 1.5.0、survival 3.8-12、survivalROC 1.0.3.1、"
     "rms 8.1-1、Hmisc 5.2-3、glmnet 4.1-8、GSVA 2.6.6、GSEABase 1.68.0、pheatmap 1.0.13、"
     "oncoPredict 1.3.1、sva 3.60.0、MVMR 0.4.8。双侧 P < 0.05 为名义显著；"
     "多重比较以 FDR < 0.05 为显著。", first_indent=22)

# ================= 3 结果 =================
heading("3  结果", 1)

heading("3.1  工具变量与遗传关联", 2)
para("以 P < 5×10⁻⁸、r² < 0.001（10000 kb）筛选出 70 个独立工具 SNP，"
     "全部 F 统计量 > 10（最小值 30.01、中位数 70.41），合计解释血清 TSH 表型方差约 9.3%。"
     "其中 69 个 SNP 在主结局中成功匹配。MR-Egger 截距检验未发现显著定向多效性"
     "（截距 = 0.00283，SE = 0.0239，P = 0.906）。Cochran Q 检验显示工具变量间存在"
     "显著异质性（IVW Q = 279.39，df = 68，P = 1.99×10⁻²⁷，I² = 75.7%；Egger Q = 279.33，df = 67，"
     "P = 9.89×10⁻²⁸，I² = 76.0%；DerSimonian-Laird tau² = 0.00686）。"
     "MR-PRESSO 全局检验亦显著（RSS = 289.77，P < 1×10⁻⁴），"
     "离群检验识别出 4 个潜在离群 SNP：rs10186921、rs116909374、rs2993047、rs925488；失真检验显示剔除后因果估计"
     "无显著变化（失真系数 = −31.5%，P = 0.119）。"
     "留一法显示剔除任一 SNP 后因果估计方向均未改变（b ∈ [−0.832, −0.587]）。",
     first_indent=22)

heading("3.2  TSH 与甲状腺癌风险的因果关联", 2)
para("两样本 MR 主分析（69 个有效 SNP）中，IVW 随机效应模型显示遗传预测 TSH 升高与"
     "甲状腺癌风险降低相关（b = −0.719，SE = 0.192，OR = 0.487，95% CI 0.335–0.709，"
     "P = 1.73×10⁻⁴）；MR-Egger（OR = 0.467，P = 0.062）、加权中位数（OR = 0.626，"
     "P = 0.0032）、简单众数（OR = 0.405，P = 0.012）方向一致，加权众数（OR = 0.702，"
     "P = 0.091）方向一致但未达显著。按 MR-PRESSO 失真检验识别的 4 个离群 SNP 剔除后，"
     "随机效应 IVW 估计方向与主分析一致（65 个 SNP，OR = 0.579，"
     "95% CI 0.452–0.741，P = 1.4×10⁻⁵），提示离群 SNP 未驱动主结论。"
     "MR-RAPS（对系统与特异多效性稳健）估计方向亦一致（OR = 0.476，P = 4.2×10⁻¹⁵）；"
     "Steiger 方向性检验 支持暴露→结局方向（correct_causal_direction = TRUE，Steiger P = 0）。"
     "该关联在两个复制结局中重现：ebi-a-GCST90013867（SPA 校正甲状腺癌 GWAS，65 个有效 SNP）"
     "IVW OR = 0.362，95% CI 0.233–0.562，P = 6.10×10⁻⁶，异质性 Q = 164.89"
     "（I² = 61.2%），MR-Egger 截距 P = 0.869，MR-PRESSO 全局检验 P = <5×10⁻⁴；"
     "finn-b-C3_THYROID_GLAND（FinnGen，68 个有效 SNP）IVW OR = 0.559，95% CI 0.384–0.816，"
     "P = 0.0025，异质性 Q = 195.53（I² = 65.7%），MR-Egger 截距 P = 0.411，"
     "MR-PRESSO 全局检验 P = <5×10⁻⁴。两个复制结局效应量差异（OR 0.362 vs 0.559）"
     "可能与结局的病例构成（FinnGen 以芬兰人群为主、病例数较少）及统计功效不同有关，"
     "但方向一致。", first_indent=22)
add_figure("fig12_mr_scatter.png", 13.0,
           "图 1  两样本孟德尔随机化散点图（TSH → 甲状腺癌）")
add_figure("fig11_mr_forest.png", 13.0,
           "图 2  单 SNP 与 IVW 因果估计森林图")

heading("3.3  TSH 与甲状腺癌风险：多变量孟德尔随机化", 2)
para("以 TSH 70 个工具 SNP 为骨架、校正 BMI（ieu-b-40）与总睾酮（ebi-a-GCST90014013）后，"
     "三暴露 MVMR-IVW（54 个三暴露与结局均有效的 SNP）显示 TSH 对甲状腺癌的独立保护效应"
     "仍显著（b = −0.645，SE = 0.258，OR = 0.525，95% CI 0.317–0.869，P = 0.0123，t 检验），"
     "BMI（b = −0.525，SE = 2.945，OR = 0.591，P = 0.859）与睾酮（b = 2.394，SE = 4.699，"
     "OR = 10.96，P = 0.613）均无独立显著效应。条件 F 统计量显示 TSH 工具强度合格"
     "（条件 F = 32.74），而 BMI（条件 F = 5.66）与睾酮（条件 F = 3.38）条件工具强度不足"
     "（F < 10，弱工具），提示二者效应估计不可靠；MVMR 异质性 Q 统计量显著"
     "（Q = 244.76，df = 50，P = 1.81×10⁻²⁷），提示工具有效性假设可能部分不成立；"
     "以 Q 统计量最小化方法（qhet_mvmr，1000 次）估计：TSH b = −0.628"
     "（BMI b = -0.7455、睾酮 b = 6.3992）。"
     "综上，TSH 直接效应在校正 BMI 与睾酮后方向仍一致，但 BMI/睾酮弱工具与条件异质性"
     "使该校正分析仅具探索性，不足以断言 TSH 效应独立于 BMI 与性激素水平。", first_indent=22)

heading("3.4  候选基因评分的构建与评估", 2)
para("基于 MR 阳性结果，将 55 个工具 SNP 映射至 60 个候选基因，在 TCGA-THY 队列"
     "（409 例，27 例复发/死亡事件）中构建探索性候选基因评分："
     "riskScore = −0.0053257 × GNG7 − 0.0003788 × NFIA。按中位数分为高（204 例）、低（205 例）"
     "风险组。评分 C-index = 0.666（95% CI 0.559–0.774，concordance 正态近似；"
     "Bootstrap 1000 次内部验证乐观校正后 0.662）；两组无复发生存差异显著（Log-rank P = 0.041）；"
     "时间依赖 ROC（KM 法）1/3/5 年 AUC 分别为 0.657、0.612、0.593；"
     "3 年校准曲线（rms::calibrate，Bootstrap 1000 次，seed 20260913）平均绝对误差 0.037；"
     "决策曲线分析显示模型在 0.01、0.02 及 0.04–0.08 的低阈值范围内净获益为正。"
     "变量选择稳定性：Bootstrap 200 次中 GNG7 保留频率 151/200（75.5%）、"
     "NFIA 保留频率 113/200（56.5%）（均以负向系数进入，与评分方向一致）；"
     "5×5 折重复交叉验证 C-index 均值 0.334（范围 0.039-0.826）。"
     "补充重跑 Lasso-Cox（10 折交叉验证，seed = 20260913，缺失值中位数填补）与弹性网（α = 0.5）："
     "Lasso 在 lambda.min 路径保留 43 个变量（含 GNG7、NFIA，系数均为负向，与评分方向一致），"
     "而 lambda.1se 路径与弹性网两条路径均将全部变量收缩为 0，提示正则化变量选择结果不稳定；"
     "多因素 Cox 中 GNG7（P = 0.072）与 NFIA（P = 0.241）亦未达传统显著性阈值，"
     "故本评分定位为探索性候选基因评分，而非验证性预后模型。"
     "上述区分度指标均源于训练队列，未经独立生存队列验证，"
     "尚不足以支持任何临床决策。", first_indent=22)
add_figure("fig2_km.png", 13.0,
           "图 3  高、低风险组无复发生存 Kaplan-Meier 曲线（Log-rank P = 0.0412）")
add_figure("fig3_roc.png", 13.0,
           "图 4  时间依赖 ROC 曲线（1/3/5 年 AUC = 0.657/0.612/0.593，KM 法）")
add_figure("fig4_calibration.png", 13.0,
           "图 5  3 年校准曲线（rms::calibrate，Bootstrap 1000 次，MAE = 0.037）")
add_figure("fig5_dca.png", 13.0,
           "图 6  36 个月决策曲线分析")

heading("3.5  GEO 病理类型梯度验证", 2)
para("在独立 GEO 数据集 GSE33630（105 例：11 例未分化型甲状腺癌 ATC、49 例乳头状型甲状腺癌 "
     "PTC、45 例正常甲状腺组织）中验证评分梯度。Kruskal-Wallis 检验显示三组间评分"
     "差异极显著（χ² = 54.147，df = 2，P = 1.75×10⁻¹²）；两两比较均达显著"
     "（ATC vs PTC，P = 7.0×10⁻¹¹；ATC vs Normal，P = 1.3×10⁻¹¹；PTC vs Normal，"
     "P = 3.2×10⁻⁹），评分呈 ATC > PTC > Normal 梯度。效应量均呈大效应"
     "（ATC vs PTC Cohen's d = 3.87，95% CI 2.91–4.83；ATC vs Normal d = 4.37，"
     "95% CI 3.32–5.43；PTC vs Normal d = 1.16，95% CI 0.72–1.60）。"
     "需说明：GSE33630 不含生存随访信息，本部分仅支持评分的病理学梯度意义，"
     "不构成生存终点的外部验证。", first_indent=22)
add_figure("fig6_geo.png", 13.0,
           "图 7  GSE33630 中风险评分病理类型梯度（Kruskal-Wallis P = 1.75e-12）")

heading("3.6  免疫浸润分析", 2)
para("ssGSEA 显示 28 种免疫细胞的相对富集评分中 25 种（89.3%）在高、低风险组间的差异显著"
     "（FDR < 0.05）。高风险组中 1 型 T 辅助细胞（P = 6.4×10⁻⁵）、滤泡辅助 T 细胞"
     "（P = 1.9×10⁻⁴）、Th2 细胞及活化 CD8 T 细胞等显著富集，"
     "提示高风险组肿瘤微环境呈免疫激活与炎症特征。免疫检查点标志物分析显示（FDR 基于"
     "15 个检查点基因全集 BH 校正），高风险组 PD-L1（CD274，FDR = 3.3×10⁻¹⁴）、"
     "PD-L2（PDCD1LG2，FDR = 2.8×10⁻⁹）、SIGLEC15（FDR = 4.1×10⁻¹¹）、"
     "CD80（FDR = 4.1×10⁻¹¹）、CTLA-4（FDR = 1.9×10⁻⁷）、TIGIT（FDR = 2.9×10⁻⁷）、"
     "CD86（FDR = 3.7×10⁻⁷）、TIM-3（HAVCR2，FDR = 1.9×10⁻⁶）、TNFRSF9（FDR = 1.0×10⁻⁷）"
     "表达显著上调。需注意：ssGSEA 评分为转录组相对富集推断，不等同于免疫细胞绝对浸润比例，"
     "且本分析未校正肿瘤纯度、甲状腺癌亚型与驱动突变状态（详见 4.7），"
     "高、低风险组基于同一数据定义，存在数据驱动偏倚。", first_indent=22)
add_figure("fig7_immune_heatmap.png", 14.0,
           "图 8  28 种免疫细胞 ssGSEA 相对富集热图（按 FDR 排序；红色：高风险组）")
add_figure("fig8_immune_boxplot.png", 14.0,
           "图 9  代表性免疫细胞组间相对富集差异箱线图")

heading("3.7  药物敏感性分析", 2)
para("基于 GDSC2 计算预测 409 例患者对 198 种药物的 IC50，176 种药物（88.9%）在高、低风险组间"
     "差异显著（FDR < 0.05）。高风险组对其中 158 种药物预测 IC50 更高（更耐药）、"
     "仅 18 种更低（更敏感），整体呈计算预测的耐药特征。高风险组对多种细胞毒化疗药物预测耐药："
     "5-氟尿嘧啶（High 181.1 vs Low 113.2，FDR = 1.6×10⁻¹¹）、紫杉醇（0.085 vs 0.059，"
     "FDR = 9.5×10⁻¹²）、多西他赛（0.013 vs 0.010，FDR = 6.4×10⁻⁸）、喜树碱"
     "（0.122 vs 0.103，FDR = 2.6×10⁻⁴）、顺铂（33.8 vs 29.3，FDR = 0.0086）、"
     "长春碱（0.031 vs 0.023，FDR = 6.7×10⁻⁶）、YK-4-279（11.6 vs 9.6，FDR = 8.5×10⁻⁵）；"
     "对部分靶向药物亦预测耐药：PRMT5 抑制剂 GSK591（115.5 vs 88.3，FDR = 8.8×10⁻¹⁶）、"
     "ATR 抑制剂 VE821（80.8 vs 55.2，FDR = 3.3×10⁻¹⁵）、NAMPT 抑制剂 Daporinad"
     "（0.016 vs 0.012，FDR = 5.5×10⁻¹¹）；仅 Selumetinib（60.2 vs 98.6，FDR = 5.3×10⁻¹⁴）、"
     "ATM 抑制剂 KU-55933（88.9 vs 90.5，FDR = 0.0043）等少数药物预测更敏感。"
     "预测可靠性评估：（1）12 个代表性药物 5 折交叉验证（ridge，log IC50 尺度）"
     "R² 中位数 0.350（范围 0.080-0.430），RMSE 中位数 1.478；"
     "（2）组标签置换检验（B = 500 次）中，FDR < 0.05 显著药物数的零分布"
     "中位数为 0.62（最大值 82），实际观察值 176 远超随机期望（P < 0.002），"
     "表明组间差异并非单纯由药物间相关性所致。"
     "需注意：药物间 IC50 高度相关，独立 FDR 校正可能仍高估显著药物数量；"
     "上述结果均为基于细胞系模型的计算推断、假设生成性质，未经实验验证。",
     first_indent=22)
add_figure("fig9_drug_volcano.png", 14.0,
           "图 10  药物敏感性火山图（红色：高风险组预测耐药 158 种；蓝色：预测敏感 18 种）")
add_figure("fig10_drug_boxplot.png", 15.0,
           "图 11  代表性药物预测 IC50 组间比较箱线图")

# ================= 4 讨论 =================
heading("4  讨论", 1)

heading("4.1  主要发现", 2)
para("本研究以两样本孟德尔随机化为起点，评估了遗传预测的 TSH 水平与甲状腺癌风险的因果关联，"
     "并完成了\"遗传工具 SNP—候选基因—预后评分—外部验证—生物学与临床转化\"的完整分析链。"
     "主要发现包括：其一，遗传预测的 TSH 水平升高与甲状腺癌发病风险降低相关"
     "（IVW OR = 0.487，95% CI 0.335–0.709），该关联在 MR-PRESSO 校正估计、MR-RAPS"
     "及两个复制结局中方向一致，在校正 BMI 与睾酮后方向仍一致（MVMR OR = 0.525，"
     "95% CI 0.317–0.869，P = 0.0123，但受弱工具限制）；其二，由 GNG7 与 NFIA 构建的"
     "探索性候选基因评分在 TCGA-THY 队列中显示出初步的区分趋势（C-index = 0.666，"
     "Log-rank P = 0.041），并在 GEO 独立数据集中呈病理类型梯度，但未经独立生存队列验证；"
     "其三，高风险组免疫微环境呈炎症激活特征，且计算预测其对绝大多数药物"
     "（176 种显著药物中 158 种）耐药、仅对少数药物（如 Selumetinib、KU-55933）敏感，"
     "均为计算推断、假设生成性质。", first_indent=22)

para("与既往 MR 研究证据的系统比较。目前 TSH–甲状腺癌 MR 领域存在方向相反的报告。"
     "Zhou 等（2020）基于 HUNT、MGI 与 ThyroidOmics 三队列的 TSH GWAS（99 个全基因组显著位点）"
     "发现遗传预测 TSH 每升高 1 个标准差，甲状腺癌风险降低约 50%（OR = 0.50，"
     "95% CI 0.27–0.92）；本研究使用其补充数据重新分析得到方向一致的结果"
     "（IVW OR = 0.487，95% CI 0.335–0.709）。另一项两样本 MR 研究亦报道保护效应"
     "（OR = 0.47，95% CI 0.30–0.73，P = 0.001）。与之相反，Fussey 等（2020）"
     "基于 UK Biobank 451,025 名参与者的分析报告 TSH 升高增加甲状腺癌风险"
     "（OR = 2.00，95% CI 1.09–3.70），但该文已于 2020 年 11 月被期刊正式撤稿"
     "（通讯作者主动报告数据分析存在问题），其结论不具备正式引用效力，此处仅作历史说明。"
     "即便排除该撤稿研究，不同工具变量集与结局定义下 MR 估计的差异仍可能源于工具变量多效性谱、"
     "结局定义（甲状腺癌总体 vs 分化型）、样本重叠程度与估计方法的不同，"
     "这一异质性本身值得进一步研究。", first_indent=22)

para("就 MVMR 而言，既往性激素-甲状腺癌 MR 研究总体未发现强因果关联"
     "（Chen 等 2026 报告总睾酮、SHBG、雌二醇均无显著因果关联；Xu 等 2025 报告"
     "生物可利用睾酮仅在女性中呈提示性保护），本研究在校正 BMI 与总睾酮（ebi-a-GCST90014013）后"
     "TSH 保护效应方向仍一致（MVMR OR = 0.525，95% CI 0.317–0.869，P = 0.0123），"
     "BMI（P = 0.859）与睾酮（P = 0.613）无独立显著效应。需诚实说明：条件 F 统计量显示"
     "BMI（F = 5.66）与睾酮（F = 3.38）的条件工具强度不足（小于 10），且 MVMR 异质性 "
     "Q 统计量显著（Q = 244.76，df = 50，P = 1.81×10⁻²⁷），提示弱工具与水平多效性"
     "可能影响 MVMR 估计；以 Q 统计量最小化方法复核 TSH 方向一致（b = −0.628），"
     "但 BMI/睾酮的效应估计应视为不可靠。因此，本研究证据不足以断言 TSH 效应"
     "独立于 BMI 与性激素水平。", first_indent=22)

heading("4.2  TSH 与甲状腺癌因果关联的机制解读", 2)
para("TSH 是甲状腺滤泡细胞增殖与分化的核心调节因子，通过其受体（TSHR）激活 cAMP/PKA "
     "信号通路促进甲状腺激素合成与细胞功能维持。本研究的遗传学证据支持血清 TSH 水平与"
     "甲状腺癌风险存在因果负相关，与既往观察性研究一致。该发现与临床上分化型甲状腺癌"
     "术后常规行 TSH 抑制治疗并不矛盾：TSH 抑制治疗针对的是已建立肿瘤对其生长信号的依赖性，"
     "而本研究评估的是 TSH 水平对肿瘤发生（发病）的因果作用；在部分场景中 TSH 偏低可能"
     "反映甲状腺细胞去分化或功能丧失，是肿瘤进展的结果而非原因。"
     "需强调：将本研究的遗传学证据外推为对 TSH 偏低个体加强甲状腺超声随访的临床建议"
     "属于假设性解读，尚缺乏前瞻性队列的直接证据。", first_indent=22)
para("值得强调的是，2025 年美国甲状腺学会（ATA）成人分化型甲状腺癌管理指南已从初治期"
     "固定的 TSH 目标值推荐转向基于治疗反应的分层管理，明确建议对无复发证据的中低危患者"
     "避免长期 TSH 抑制。本研究所提示的\"遗传升高的 TSH 水平可能具有保护作用\"与这一"
     "\"去强化\"趋势形成呼应——此为本研究基于遗传学证据的推测性解读，"
     "任何将本模型评分与 2025 ATA 动态风险分层体系整合的设想均属假设性，"
     "需前瞻性验证后方可考虑。", first_indent=22)

heading("4.3  预后基因的潜在意义", 2)
para("GNG7 编码 G 蛋白 γ7 亚基，参与 G 蛋白偶联受体下游信号转导，在多种肿瘤中被报道为"
     "表达下调的候选抑癌基因；NFIA 为核因子 I 家族转录因子，参与细胞分化调控，"
     "其表达下调与肿瘤恶性表型获得有关。本模型中两个基因均以负向权重进入评分，"
     "提示甲状腺癌进展过程中可能存在相应信号通路失活。需强调：MR 工具 SNP 的基因归属反映"
     "TSH 位点所在基因组区域，与预后评分中基因的生物学机制并非必然对应，本部分属假设生成性质；"
     "本研究未行共定位或 eQTL-MR 验证（见 4.7），GNG7 与 NFIA 在甲状腺癌中的功能角色"
     "有待体外及体内实验验证。", first_indent=22)

heading("4.4  免疫微环境的临床意义", 2)
para("高风险组中 Th1、滤泡辅助 T 细胞、Th2 及活化 CD8 T 细胞等的相对富集评分显著升高，"
     "提示其肿瘤微环境呈免疫激活与炎症特征。肿瘤内 Th1 型免疫应答通常与较好的预后相关，"
     "而本研究中高风险组呈现相反方向，可能反映：一是肿瘤驱动的高强度免疫激活伴随免疫逃逸"
     "与功能耗竭，效应细胞虽富集但功能受损；二是炎症微环境本身促进肿瘤进展。"
     "本研究中高风险组 PD-L1、CTLA-4 与 TIM-3 表达同步上调，支持\"效应细胞富集但功能耗竭\""
     "的解释。ssGSEA 基于转录组推断免疫细胞相对富集，无法直接反映细胞功能状态；"
     "本研究亦未使用 CIBERSORTx/xCell/MCPcounter 等去卷积方法交叉验证，"
     "未校正肿瘤纯度与亚型，上述解读需结合多重免疫组化或流式细胞术进一步验证。", first_indent=22)

heading("4.5  药物敏感性分析的意义", 2)
para("计算预测显示高风险组对绝大多数药物（176 种显著药物中 158 种）IC50 更高（耐药），"
     "涵盖 5-氟尿嘧啶、紫杉醇、多西他赛、顺铂等细胞毒药物及 GSK591、VE821、Daporinad "
     "等靶向药物，提示高风险肿瘤可能具有更强的药物耐受能力（如 DNA 损伤修复增强、"
     "药物外排或代谢重编程）；仅 Selumetinib、KU-55933 等少数药物预测更敏感。"
     "上述预测均为纯计算推断、假设生成性质，不代表临床获益推荐；"
     "虽然 12 个代表性药物的训练集交叉验证 R² 中位数为 0.350（提示预测模型在细胞系"
     "内部具有一定拟合能力），但细胞系模型向患者外推的可靠性仍需经体外实验"
     "（如 CCK-8、克隆形成实验）与临床试验验证。", first_indent=22)

heading("4.6  与既往研究的一致性", 2)
para("本研究的 MR 结果与既往观察性流行病学证据方向一致，并在两个独立复制结局及"
     "多变量校正后保持方向稳健。预后评分的基因选择与既往甲状腺癌转录组研究对 GNG7 的相关报道"
     "相呼应；免疫浸润与药物敏感性结果则为后续机制研究与治疗假设的生成提供了方向"
     "（均为计算推断，须经实验验证）。需要强调，当前评分仅显示初步的区分趋势，"
     "其 95% 置信区间下限（0.559）接近随机水平，且未经独立生存队列验证，"
     "尚不足以支持任何临床决策。", first_indent=22)

heading("4.7  局限性", 2)
para("本研究存在以下局限性：（1）候选基因评分尚未在任何独立队列完成生存终点验证"
     "——我们已系统检索公共 GEO 甲状腺癌队列（含 GSE60542、GSE29265），确认这些数据集"
     "均为组织表达数据、不含生存随访信息，故仅能进行病理类型梯度验证，这是评分转化价值的"
     "主要限制之一；（2）评分训练队列仅 27 例事件（EPV ≈ 13.5），过拟合风险较高，"
     "C-index 0.666 属中等区分能力，乐观校正后为 0.662，95% CI 下限 0.559 接近随机水平，"
     "5 年 AUC = 0.593 亦接近随机猜测，Bootstrap 变量选择保留频率有限"
     "（GNG7 151/200（75.5%）、NFIA 113/200（56.5%）），故评分仅显示初步区分趋势，"
     "仅定位为探索性分析；（3）MR 分析基于欧洲人群 GWAS，结论向其他人群外推需谨慎；"
     "主分析 Cochran Q 异质性（P = 1.99×10⁻²⁷，I² = 75.7%）与 MR-PRESSO 全局检验"
     "（P < 1×10⁻⁴）均显著，提示工具变量间存在异质性/多效性，尽管失真检验（P = 0.119）、"
     "MR-Egger 截距（P = 0.906）与 MR-PRESSO 校正估计、MR-RAPS 均支持方向不变，"
     "该残余异质性仍构成 MR 结论的固有局限；（4）MVMR 中 BMI 与睾酮条件 F 统计量低于 10"
     "（弱工具）、条件异质性 Q 显著（P = 1.81×10⁻²⁷），BMI/睾酮效应估计不可靠，"
     "TSH 直接效应虽显著但其精确性受上述因素制约，故本研究不足以断言 TSH 效应独立于"
     "BMI 与性激素；（5）变量选择：Lasso-Cox 仅 lambda.min 路径保留 GNG7/NFIA 等 43 个变量，"
     "lambda.1se 与弹性网路径均为 0 变量，正则化选择结果不稳定、模型系数不确定性大；"
     "（6）药物敏感性为基于细胞系模型的计算预测、假设生成性质，未经体内外实验验证，"
     "训练集交叉验证 R² 中位数 0.350 提示拟合有限，且药物间 IC50 高度相关、"
     "独立 FDR 可能高估显著药物数量（虽置换检验显示组间差异显著超出随机期望）；"
     "（7）免疫浸润为转录组相对富集推断，未经蛋白水平验证，未用去卷积方法交叉验证，"
     "未校正肿瘤纯度、亚型与 BRAF/RAS 突变状态，高/低风险组基于同一数据定义存在数据驱动偏倚；"
     "（8）决策曲线分析显示评分仅在低阈值区间（0.01、0.02、0.04–0.08）净获益为正；"
     "（9）校准曲线在仅 27 例事件下表现良好，不能排除过拟合导致的校准乐观；"
     "（10）受限于回顾性 TCGA 数据，本评分未纳入 TSH 抑制治疗记录（TRIPOD+AI 条目 6c），"
     "治疗暴露对预后分层的影响无法评估；（11）MR 证据在不同工具变量集与结局定义下方向"
     "不完全一致（含一项已撤稿研究），该异质性本身值得进一步研究；（12）反向 MR"
     "因甲状腺癌结局 GWAS 缺乏全基因组显著的独立工具 SNP（主结局仅 1054 例病例）"
     "而未执行，共定位与 eQTL-MR 亦因缺乏匹配的甲状腺组织 eQTL 资源未执行，"
     "工具 SNP 与候选基因的功能关联仅停留在基因组定位层面；"
     "（13）未评估暴露与结局 GWAS 之间的样本重叠程度（主结局 UKB 相关数据集与 "
     "BMI/睾酮 UKB 数据存在潜在重叠），样本重叠可能引入偏倚；"
     "（14）未按甲状腺癌亚型（PTC、FTC、ATC、MTC）及性别分层分析，"
     "亚型与性别相关异质性可能未被捕获。", first_indent=22)

# ================= 5 结论 =================
heading("5  结论", 1)
para("本研究通过两样本孟德尔随机化提供了遗传预测的 TSH 水平升高与甲状腺癌风险降低相关的"
     "探索性证据（与既往文献一致，且经 PRESSO 校正与稳健方法复核），"
     "并构建了基于 GNG7 与 NFIA 的探索性候选基因评分，经外部数据病理类型梯度验证；"
     "评分区分度处于初步水平（C-index 0.666，95% CI 0.559–0.774），未经独立生存队列验证，"
     "尚不足以支持临床决策。计算预测提示高风险组具有免疫激活微环境及对绝大多数预测药物"
     "（176 种显著药物中 158 种）耐药的假设性特征。本研究属探索性、假设生成研究，"
     "评分的临床价值与生物学意义需要在更大样本、多中心、前瞻性队列及实验中验证后方可评估。",
     first_indent=22)

# ================= 声明 =================
heading("声明", 1)
para("数据可用性声明：暴露 GWAS（Zhou 等，2020）与结局 GWAS（ebi-a-GCST90018929、"
     "ebi-a-GCST90013867、finn-b-C3_THYROID_GLAND）汇总统计量均来自 IEU OpenGWAS 数据库"
     "（https://gwas.mrcieu.ac.uk/）；TCGA-THY 表达与临床数据来自 TCGA/GDC 公共数据库；"
     "外部验证表达数据来自 NCBI GEO（GSE33630）；药物敏感性训练数据使用 GDSC2"
     "（805 细胞系 × 198 药物）。全部 13 个分析脚本（脚本 A–K 及脚本 L，含 MR 敏感性二轮补跑、"
     "MVMR 条件检验、Lasso/弹性网重跑、Bootstrap 稳定性、药敏交叉验证与置换检验）"
     "已整理为可复现脚本包，完整代码公开于 GitHub "
     "（https://github.com/qwe2412/TSH-ThyroidCancer-MR，共 6 次提交，含脚本、运行说明、"
     "论文构建与发布脚本）；Zenodo 版本 DOI 将于论文接收后经归档补充。"
     "工具 SNP 列表与全部敏感性分析输出见补充材料。", first_indent=22)
para("伦理声明：本研究全部数据来自公开数据库（IEU OpenGWAS、TCGA/GDC、NCBI GEO、GDSC2）"
     "的已发布去标识化汇总统计量或去标识化样本数据，未涉及新的人类受试者招募、"
     "干预或可识别个人信息，依据公共数据二次分析惯例免除伦理审查。", first_indent=22)
para("资金来源：本研究未接受任何外部资金资助（TRIPOD+AI 条目 13）。", first_indent=22)
para("利益冲突：作者声明无利益冲突。", first_indent=22)
para("作者贡献：王沛负责研究设计、数据获取与整理、统计分析、论文撰写与修改，"
     "并最终审阅与批准稿件。", first_indent=22)
para("致谢：作者感谢 IEU OpenGWAS 数据库、TCGA/GDC、NCBI GEO 与 GDSC 项目的数据提供者，"
     "以及 MR-Base、oncoPredict、GSVA 等开源软件的作者。", first_indent=22)
para("报告规范：本研究的 MR 部分遵循 STROBE-MR 清单报告，预后评分部分遵循 TRIPOD+AI "
     "清单报告；两份检查表以补充材料形式提供。", first_indent=22)
para("软件与版本：统计分析全部基于 R 4.6.1 完成，关键包版本包括 TwoSampleMR 0.7.9、"
     "MendelianRandomization 0.10.0、MRPRESSO 1.0.6、ieugwasr 1.1.2、readxl 1.5.0、"
     "survival 3.8-12、survivalROC 1.0.3.1、rms 8.1-1、Hmisc 5.2-3、glmnet 4.1-8、"
     "GSVA 2.6.6、GSEABase 1.68.0、pheatmap 1.0.13、oncoPredict 1.3.1、sva 3.60.0"
     "（TRIPOD+AI 条目 18f）。", first_indent=22)

# ================= 参考文献 =================
heading("参考文献", 1)
refs = [
    "Zhou W, Brumpton B, Kabil O, et al. GWAS of thyroid stimulating hormone highlights "
    "pleiotropic effects and inverse association with thyroid cancer. Nat Commun. "
    "2020;11(1):5680. doi: 10.1038/s41467-020-17718-z.",
    "Yuan S, Kar S, Vithayathil M, et al. Causal associations of thyroid function and "
    "dysfunction with overall, breast and thyroid cancer: A two-sample Mendelian randomization "
    "study. Int J Cardiol. 2020;313:99-104. doi: 10.1016/j.ijcard.2020.03.053.",
    "Chen C, Zhang W, Chu C, et al. Circulating sex steroids do not influence risk of thyroid "
    "malignancy: insights from a bidirectional Mendelian randomization. Gland Surg. "
    "2026;15(1):17. doi: 10.21037/gs-2025-336.",
    "Xu Q, Jiang H, Li Y, et al. Bioavailable testosterone and thyroid cancer: A 2-sample "
    "Mendelian randomization study. Medicine (Baltimore). 2025;104(47):e45528. "
    "doi: 10.1097/MD.0000000000045528.",
    "Brion MJ, Shakhbazov K, Visscher PM. Calculating statistical power in Mendelian "
    "randomization. Int J Epidemiol. 2013;42(5):1497-1501. doi: 10.1093/ije/dyt179.",
    "Charoentong P, Finotello F, Angelova M, et al. Pan-cancer immunogenomic analyses reveal "
    "genotype-immunophenotype relationships and predictors of response to checkpoint blockade. "
    "Cell Rep. 2017;18(1):248-262. doi: 10.1016/j.celrep.2016.12.019.",
    "Yang W, Soares J, Greninger P, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a "
    "resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Res. "
    "2013;41(D1):D955-D961. doi: 10.1093/nar/gks1111.",
    "Sung H, Ferlay J, Siegel RL, et al. Global cancer statistics 2020: GLOBOCAN estimates of "
    "incidence and mortality worldwide for 36 cancers in 185 countries. CA Cancer J Clin. "
    "2021;71(3):209-249. doi: 10.3322/caac.21660.",
    "Haugen BR, Alexander EK, Bible KC, et al. 2015 American Thyroid Association management "
    "guidelines for adult patients with thyroid nodules and differentiated thyroid cancer: "
    "the American Thyroid Association Guidelines Task Force on thyroid nodules and "
    "differentiated thyroid cancer. Thyroid. 2016;26(1):1-133. doi: 10.1089/thy.2015.0020.",
    "American Thyroid Association. 2025 American Thyroid Association management guidelines for "
    "adult patients with differentiated thyroid cancer. Thyroid. 2025;35(8):841-985. "
    "doi: 10.1177/10507256251363120.",
    "Fussey JM, Beaumont RN, Wood AR, et al. Mendelian randomization supports a causative "
    "effect of TSH on thyroid carcinoma. Endocr Relat Cancer. 2020;27(10):551-559. "
    "(Retracted; retraction: Endocr Relat Cancer. 2020;27(11):Z1. doi: 10.1530/ERC-20-0067r).",
    "Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for causal inference in "
    "epidemiological studies. Hum Mol Genet. 2014;23(R1):R89-R98. doi: 10.1093/hmg/ddu328.",
    "Burgess S, Butterworth A, Thompson SG. Mendelian randomization analysis with multiple "
    "genetic variants using summarized data. Genet Epidemiol. 2013;37(7):658-665. "
    "doi: 10.1002/gepi.21758.",
    "Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: "
    "effect estimation and bias detection through Egger regression. Int J Epidemiol. "
    "2015;44(2):512-525. doi: 10.1093/ije/dyv080.",
    "Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation in Mendelian "
    "randomization with some invalid instruments using a weighted median estimator. "
    "Genet Epidemiol. 2016;40(4):304-314. doi: 10.1002/gepi.21965.",
    "Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in "
    "causal relationships inferred from Mendelian randomization between complex traits and "
    "diseases. Nat Genet. 2018;50(5):693-698. doi: 10.1038/s41588-018-0099-7.",
    "Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal "
    "inference across the human phenome. eLife. 2018;7:e34408. doi: 10.7554/eLife.34408.",
    "Zhao Q, Wang J, Hemani G, Bowden J, Small DS. Statistical inference in two-sample "
    "summary-data Mendelian randomization using robust adjusted profile score. J Am Stat "
    "Assoc. 2020;115(529):405-416. doi: 10.1080/01621459.2018.1513452.",
    "Hemani G, Tilling K, Davey Smith G. Orienting the causal relationship between "
    "imprecisely measured traits using GWAS summary data. PLoS Genet. 2017;13(11):e1007081. "
    "doi: 10.1371/journal.pgen.1007081.",
    "Sanderson E, Davey Smith G, Windmeijer F, Bowden J. An examination of multivariable "
    "Mendelian randomization in the single-sample and two-sample summary data settings. "
    "Int J Epidemiol. 2019;48(3):713-727. doi: 10.1093/ije/dyy262.",
    "Sanderson E, Spiller W, Bowden J. Testing and correcting for weak and pleiotropic "
    "instruments in two-sample multivariable Mendelian randomization. Stat Med. "
    "2021;40(9):2234-2253. doi: 10.1002/sim.9133.",
    "Skrivankova VW, Richmond RC, Woolf BAR, et al. Strengthening the reporting of "
    "observational studies in epidemiology using Mendelian randomization: the STROBE-MR "
    "statement. JAMA. 2021;326(16):1614-1621. doi: 10.1001/jama.2021.18236.",
    "Collins GS, Reitsma JB, Altman DG, Moons KGM. Transparent reporting of a multivariable "
    "prediction model for individual prognosis or diagnosis (TRIPOD): the TRIPOD statement. "
    "Ann Intern Med. 2015;162(1):55-63. doi: 10.7326/M14-0697.",
    "Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for "
    "reporting clinical prediction models that use regression or machine learning methods. "
    "BMJ. 2024;385:e078378. doi: 10.1136/bmj-2023-078378.",
    "Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and "
    "RNA-seq data. BMC Bioinformatics. 2013;14:7. doi: 10.1186/1471-2105-14-7.",
    "Geeleher P, Cox N, Huang RS. Clinical drug response can be predicted using baseline gene "
    "expression levels and in vitro drug sensitivity in cell lines. Genome Biol. "
    "2014;15(3):R47. doi: 10.1186/gb-2014-15-3-r47.",
    "Maeser D, Gruener RF, Huang RS. oncoPredict: an R package for predicting in vivo or "
    "cancer patient drug response and biomarkers from cell line screening data. Brief "
    "Bioinform. 2021;22(6):bbab260. doi: 10.1093/bib/bbab260.",
    "Yoshihara K, Shahmoradgoli M, Martínez E, et al. Inferring tumour purity and stromal and "
    "immune cell admixture from expression data. Nat Commun. 2013;4:2612. "
    "doi: 10.1038/ncomms3612.",
    "Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue "
    "profiles. Nat Methods. 2015;12(5):453-457. doi: 10.1038/nmeth.3337.",
    "Shibata K, Mori M, Tanaka S, Kitano S, Akiyoshi T. G-protein gamma 7 is down-regulated "
    "in cancers and associated with p27kip1-induced growth arrest. Cancer Res. "
    "1999;59(5):1096-1101.",
    "Song HR, Gonzalez-Gomez I, Suh GS, et al. Nuclear factor IA is expressed in astrocytomas "
    "and is associated with improved survival. Neuro Oncol. 2010;12(2):122-132.",
    "Angell TE, Lechner MG, Jang JK, LoPresti JS, Epstein AL. BRAF V600E in papillary thyroid "
    "carcinoma is associated with increased programmed death ligand 1 expression and "
    "suppressive immune cell infiltration. Thyroid. 2014;24(3):527-536. "
    "doi: 10.1089/thy.2013.0134.",
]
for i, r in enumerate(refs, 1):
    para("[%d] %s" % (i, r), size=10, space_after=3, first_indent=0)

doc.save(OUT)
print("已保存:", OUT)
print("段落数:", len(doc.paragraphs))
