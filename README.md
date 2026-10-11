# comorbidity-causal-evidence

检索并筛选「疾病 A → 疾病 B」**方向性因果共病**文献的 WorkBuddy skill：仅人类、按六项方法学门槛打分、显式判定方向、输出结构化 CSV 与 PRISMA 留痕。

## 为什么自己写

GitHub 上 `comorbidity + causal`、`disease pair causality` 这类关键词**搜不到现成项目**；
现有的 PubMed 相关工具（pubmed-mcp-server、mcp-simple-pubmed、PubMed-MCP-Server、
medsci-skills、paper-search-pro）与技能市场候选**全部只做检索层**——给你一堆文献就结束，
没有方向性判定、没有人类过滤、没有质量打分。

## 目录

| 文件 | 说明 |
|---|---|
| `SKILL.md` | 主文件：5 条硬约束 + 6 步工作流 + A/B/C/D 分级 + 方向五选一 |
| `references/rubric.md` | 六项 0/1/2 打分细则（S1 时序 / S2 表型 / S3 共病组区分 / S4 协同交互 / S5 偏倚与因果鉴别 / S6 结论尺度） |
| `references/example.md` | Crohn disease ↔ Psoriasis 实跑留档：打分表、与 MR 冲突的处理、summary 骨架 |
| `scripts/pubmed_search.py` | 单疾病对检索脚本，仅标准库，E-utilities 只读 GET |

## 用法

对任意疾病对说一句话即可，例如：
「查一下类风湿关节炎会不会导致结肠癌，要人类研究、有方向性的高质量证据。」

可选环境变量：`NCBI_EMAIL`、`NCBI_API_KEY`（10 req/s，否则 3 req/s）。

## 硬约束

1. 绝不凭记忆生成文献——PMID 与结论必须来自本次真实检索。
2. 检索式与命中数必须留痕（写进 `*_query_log.md`）。
3. `humans[Filter]` 必须出现在每条检索式里。
4. 方向必须显式判定，判不出就写「证据不足」。
5. 每条结论挂 PMID 与链接。

## 批量版（2237 对疾病）

若要对大量疾病对批量跑，用配套的服务器端工作流（`pubmed_work/`）：
共现检索 → 抓摘要 → 本地正则判方向 → DeepSeek 打分 → 与 MR 交叉验证 → 出报告。
纯标准库，可断点续跑，支持并发与限速。

> ⚠️ 踩过的坑：不要在 PubMed 查询层做方向性判断。
> `"risk of colon cancer"[tiab]` 返回 0 条（短语索引里没有），
> 而 `"risk of breast cancer"[tiab]` 又有 11846 条——同一模板时灵时不灵且静默失败。
> 方向性应在抓回来的原文上判定。
