# 完整示例：Crohn disease ↔ Psoriasis（2026-10-12 实跑）

本文件是**一次真实运行的留档**，用于示范第 1–6 步每一步长什么样。
命中数会随时间变化，重跑请以自己的 `*_query_log.md` 为准；PMID 与结论均来自本次真实检索，不得照抄到别的疾病对上。

## 第 0 步 · 术语

```bash
python scripts/pubmed_search.py --mesh "Crohn disease"
# → Crohn Disease（UID 68003424，1980 年引入）
#   同义入口词：Crohn's Disease / Crohns Disease / Crohn's Enteritis / Inflammatory Bowel Disease 1 / Regional Enteritis
```

## 第 1 步 · 检索

```bash
python scripts/pubmed_search.py \
  --a "Crohn disease" --a-terms "Crohn disease|Crohn's disease|regional enteritis" \
  --b "Psoriasis" \
  --direction both --years 2000:2026 --retmax 60 \
  --out crohn_psoriasis \
  --mr-stats "C:/Users/PC/Desktop/MR/mr_stats.js"
```

输出：

```
[general(两病共现)]   命中 463，取回 60
[A→B 方向短语]        命中   4，取回  4
[B→A 方向短语]        命中   4，取回  4
去重后候选 = 66 篇
MR A→B: Weighted median, OR=1.108, P=0.000455
MR B→A: Weighted median, OR=1.058, P=0.0831
MR 方向提示: 仅 A→B 显著 → 支持 A→B 单向
```

> MR 主方法由脚本按 CoGAP 判定树自动选取（n=1 → Wald；Egger intercept P<.05 → MR-Egger；
> Cochran's Q P<.05 → Weighted median；否则 IVW）。本例异质性 P 极小 → 走 Weighted median，
> 与网站结论句（OR = 1.108, P = 4.55×10⁻⁴, Weighted median）完全一致。

## 第 2 步 · L0 硬门槛

66 篇里，仅 8 篇命中方向短语分支；其余 58 篇多为"两病各提一句"的治疗/机制综述或无关共现，**剔 57 篇**。

## 第 3–4 步 · 打分与判定（示范 3 篇）

| PMID | 设计 | S1 | S2 | S3 | S4 | S5 | S6 | 总分 | 方向 | 级别 |
|---|---|---|---|---|---|---|---|---|---|---|
| 33556947 | 韩国全国匹配队列，IBD（2005–2008）→ 新发银屑病（2009–2018），ICD-10+处方定义 | 2 | 2 | 1 | 0 | 1 | 1 | **7** | A→B | B |
| 22941766 | NHS/NHS II 前瞻队列，银屑病 → 新发克罗恩病（女性，病历审阅确认） | 2 | 1 | 1 | 0 | 1 | 1 | **6** | B→A | B |
| 36443384 | 双向两样本 MR，银屑病/PsA ↔ CD/UC | 2 | 2 | 0 | 0 | 2 | 1 | **7** | 双向 | B |

扣分理由（写进 `note` 列，不要只写分数）：

- 33556947：CD 与 UC 合并报告，未单列克罗恩病 → 作为 A→B 证据时须注明属于 IBD 层面；
  未校正就诊频率（IBD 患者皮肤科就诊机会更多）→ S5 偏倚部分只给 1。
- 22941766：银屑病靠补充问卷自报、克罗恩病自报+病历审阅，非编码登记 → S2=1；
  护士队列本身健康素养高、筛查可及性好，检测偏倚未校正 → S5=1；仅限女性 → 外推受限。
- 36443384：MR 天然满足时序、免疫于 Berkson 与检测偏倚 → S1=2、S5=2；
  但 MR 无法做"单纯 A / 单纯 B / A+B 共病组"分层，也无交互检验 → S3=S4=0。

## 第 5 步 · 交叉验证（本例的关键发现）

| 来源 | A→B（CD→银屑病） | B→A（银屑病→CD） |
|---|---|---|
| CoGAP MR（`mr_stats.js`） | OR=1.108，P=4.55×10⁻⁴ ✅ | OR=1.058，P=0.083 ❌ |
| 36443384（双向 MR） | OR=1.425，P<0.001 ✅ | OR=1.350，P=0.013 ✅ |
| 22941766（前瞻队列） | — | 风险升高 ✅ |

**结论：冲突，必须标注待查，不许强行调和。**

CoGAP MR 只支持 A→B 单向；36443384 报告双向均显著。可能原因（写在 summary 里，不要下断言）：

1. 工具变量来源与 SNP 数量不同（本例 A→B 用 114 个 SNP，B→A 仅 31 个，反向效力不足）；
2. 结局 GWAS 人群/版本不同（FinnGen  vs  UKB/IIBDGC）；
3. 多效性处理不同（是否走 MR-Egger、是否剔除 HLA 区 SNP——银屑病与 CD 共享 HLA 位点，是最典型的多效性来源）。

处理：把 36443384 的 `consistency` 标为 `冲突`，`note` 写清冲突点；
最终方向结论写 **A→B 证据更稳健**，B→A 标为"存在正向提示但受效力与多效性限制，待查"。

## 第 6 步 · summary 骨架（照此写，不要自由发挥）

```markdown
# Crohn disease ↔ Psoriasis 方向性证据筛查

## PRISMA 计数
检索命中 471（463+4+4）→ 去重 66 → L0 剔 57 → 进入评分 9 → A 级 0 / B 级 3 / C 级 x / D 级 x

## 方向性结论（一句话）
遗传与队列证据一致支持 **克罗恩病 → 银屑病**（A→B）；
反向（银屑病 → 克罗恩病）存在双向 MR 的正向提示，但在本项目的 MR 结果中不显著，
且两方向共享 HLA 位点带来的多效性未完全排除，暂判 **证据不足以认定 B→A**。

## Top 文献
- A→B｜33556947｜全国匹配队列｜OR/HR 见原文｜B 级｜理由：时序明确、编码+处方定义、排除基线病例
- B→A｜22941766｜NHS 前瞻队列｜B 级｜理由：前瞻设计但自报表型、仅限女性
- 双向｜36443384｜双向 MR｜B 级｜理由：唯一同时报告两方向的设计；与本项目 MR 冲突，待查

## 证据空白
1. 缺"单纯 CD / 单纯银屑病 / CD+银屑病"三组比较与协同交互检验（S3、S4 普遍为 0）。
2. 缺检测偏倚校正：IBD 患者就诊频率更高，银屑病更易被发现，现有研究均未校正。
3. B→A 方向的 MR 效力不足（31 个 SNP），需换更大样本结局 GWAS 重做。
```
