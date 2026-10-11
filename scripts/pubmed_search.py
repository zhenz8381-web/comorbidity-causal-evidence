#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PubMed 共病（方向性）文献检索 —— 只做「召回 + 去重 + 留痕」，不做任何主观判定。

用法示例：
  python pubmed_search.py --a "Crohn disease" --b "Psoriasis" --direction both \
      --years 2000:2026 --retmax 300 --email you@example.com \
      --out out/crohn_psoriasis --log out/crohn_psoriasis_query_log.md \
      --mr-stats "C:/Users/PC/Desktop/MR/mr_stats.js"

  python pubmed_search.py --mesh "Crohn disease"      # 只查 MeSH 主题词

只依赖标准库；仅向 eutils.ncbi.nlm.nih.gov 发起只读 GET 请求。
NCBI 建议：带 --email；若要提高速率可设环境变量 NCBI_API_KEY（10 req/s，否则 3 req/s）。
"""

import argparse
import csv
import json
import math
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "comorbidity-causal-evidence"
_last = [0.0]


def _get(path, params, key=None):
    """带速率限制的只读 GET。"""
    gap = 0.11 if key else 0.34
    dt = time.time() - _last[0]
    if dt < gap:
        time.sleep(gap - dt)
    params = dict(params)
    params.setdefault("tool", TOOL)
    url = EUTILS + path + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": TOOL})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read()
    _last[0] = time.time()
    return data


def esearch(term, retmax, email, key, mindate=None, maxdate=None, db="pubmed"):
    p = {"db": db, "term": term, "retmax": retmax, "retmode": "json"}
    if email:
        p["email"] = email
    if key:
        p["api_key"] = key
    if mindate and maxdate:
        p.update({"datetype": "pdat", "mindate": mindate, "maxdate": maxdate})
    d = json.loads(_get("esearch.fcgi", p, key).decode("utf-8"))
    res = d.get("esearchresult", {})
    return int(res.get("count", 0)), res.get("idlist", [])


def efetch(pmids, email, key):
    if not pmids:
        return []
    out = []
    for i in range(0, len(pmids), 100):
        p = {"db": "pubmed", "id": ",".join(pmids[i:i + 100]), "retmode": "xml"}
        if email:
            p["email"] = email
        if key:
            p["api_key"] = key
        out.append(_get("efetch.fcgi", p, key).decode("utf-8", "replace"))
    return out


def _text(el):
    if el is None:
        return ""
    return "".join(el.itertext()).strip()


def parse(xmls):
    rows = []
    for xml in xmls:
        root = ET.fromstring(xml)
        for art in root.findall(".//PubmedArticle"):
            pmid = _text(art.find(".//MedlineCitation/PMID"))
            title = _text(art.find(".//Article/ArticleTitle"))
            abs_parts = [_text(a) for a in art.findall(".//Abstract/AbstractText")]
            abstract = " ".join(x for x in abs_parts if x)
            journal = _text(art.find(".//Journal/ISOAbbreviation")) or _text(art.find(".//Journal/Title"))
            year = _text(art.find(".//Journal/JournalIssue/PubDate/Year")) \
                or _text(art.find(".//Journal/JournalIssue/PubDate/MedlineDate"))[:4]
            pubtypes = ";".join(_text(p) for p in art.findall(".//PublicationTypeList/PublicationType"))
            doi = ""
            for aid in art.findall(".//ArticleIdList/ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = _text(aid)
            if not doi:
                for e in art.findall(".//ELocationID"):
                    if e.get("EIdType") == "doi":
                        doi = _text(e)
            rows.append({
                "pmid": pmid,
                "doi": doi,
                "year": year,
                "journal": journal,
                "pubtypes": pubtypes,
                "title": title,
                "abstract": abstract,
                "url": "https://pubmed.ncbi.nlm.nih.gov/%s/" % pmid if pmid else "",
            })
    return rows


def block(terms):
    """把一组同义词拼成 ( "x"[tiab] OR "x"[MeSH Terms] OR ... )"""
    parts = []
    for t in terms:
        t = t.strip().strip('"')
        if not t:
            continue
        parts.append('"%s"[tiab] OR "%s"[MeSH Terms]' % (t, t))
    return "(" + " OR ".join(parts) + ")" if parts else ""


DESIGN = ('("cohort"[tiab] OR "prospective"[tiab] OR "longitudinal"[tiab] OR "follow-up"[tiab] '
          'OR "followed"[tiab] OR "incidence"[tiab] OR "incident"[tiab] OR "risk of"[tiab] '
          'OR "case-control"[tiab] OR "registry"[tiab] OR "Mendelian randomization"[tiab] '
          'OR "self-controlled"[tiab] OR "causal"[tiab] OR "causality"[tiab] OR "precedes"[tiab] '
          'OR "subsequent"[tiab])')

EXCL_DEFAULT = 'NOT ("case reports"[pt] OR "comment"[pt] OR "editorial"[pt] OR "letter"[pt])'
EXCL_REVIEWS = 'NOT ("review"[pt] OR "systematic review"[pt] OR "meta-analysis"[pt])'


def directional(a, b):
    """方向性短语块：强调 A 在前、B 在后。"""
    return ('("%s"[tiab] AND ("risk of %s"[tiab] OR "incident %s"[tiab] OR "%s following %s"[tiab] '
            'OR "subsequent %s"[tiab] OR "developing %s"[tiab] OR "%s after %s"[tiab] '
            'OR "%s precedes %s"[tiab] OR "later %s"[tiab] OR "%s in patients with %s"[tiab]))'
            % (a, b, b, b, a, b, b, b, a, a, b, b, b, a))


def mesh_lookup(term, email, key):
    hits, ids = esearch(term, 10, email, key, db="mesh")
    if not ids:
        print("MeSH: 未找到 '%s'" % term)
        return
    p = {"db": "mesh", "id": ",".join(ids), "retmode": "json"}
    if email:
        p["email"] = email
    if key:
        p["api_key"] = key
    d = json.loads(_get("esummary.fcgi", p, key).decode("utf-8"))
    print("MeSH 命中 %d 条，前 %d 个主题词：" % (hits, len(ids)))
    res = d.get("result", {})
    for uid in res.get("uids", []):
        rec = res.get(uid) or {}
        names = rec.get("ds_meshterms") or []
        head = names[0] if names else uid
        scope = (rec.get("ds_scopenote") or "").replace("\n", " ")
        print("  - %s（UID %s，%s 年引入）" % (head, uid, rec.get("ds_yearintroduced", "-")))
        if scope:
            print("      %s" % scope[:160])
        if len(names) > 1:
            print("      同义入口词：%s" % " / ".join(names[1:6]))


def load_mr(path):
    """从 mr_stats.js 里取出 window.MR_STATS 字典（CoGAP 项目的静态预计算结果）。"""
    with open(path, "r", encoding="utf-8") as f:
        txt = f.read()
    i = txt.find("window.MR_STATS")
    if i < 0:
        return None
    start = txt.find("{", i)
    if start < 0:
        return None
    depth, end = 0, -1
    for j in range(start, len(txt)):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
            if depth == 0:
                end = j
                break
    if end < 0:
        return None
    try:
        return json.loads(txt[start:end + 1])
    except Exception:
        return None


def mr_lookup(stats, a, b):
    """按 CoGAP 判定树选主方法，返回 (method, OR, P, note)；找不到就返回空。

    判定树（多效性优先于异质性）：
      n==1                  → Wald ratio
      Egger intercept P<.05 → MR-Egger
      Cochran's Q P<.05     → Weighted median
      否则                  → IVW
    """
    if not stats:
        return ("", "", "", "")
    key = None
    for k in stats:
        if k == "%s||%s" % (a, b):
            key = k
            break
    if key is None:
        la, lb = a.lower(), b.lower()
        for k in stats:
            ka, kb = k.split("||")
            if ka.lower() == la and kb.lower() == lb:
                key = k
                break
    if key is None:
        return ("", "", "", "MR 无此疾病对")
    s = stats[key]
    m = s.get("m", {}) or {}
    n = s.get("n")
    pl, he = s.get("pl"), s.get("he")
    if n == 1 or "wald" in m:
        method, v = "Wald ratio", m.get("wald")
    elif pl is not None and pl < 0.05 and m.get("egger"):
        method, v = "MR-Egger", m.get("egger")
    elif he is not None and he < 0.05 and m.get("wmed"):
        method, v = "Weighted median", m.get("wmed")
    else:
        method, v = "IVW", m.get("ivw")
    if not v:
        return (method, "", "", "无可用效应量")
    beta, p = v[0], v[2]
    note = "SNPs=%s; 多效性P=%s; 异质性P=%s" % (
        n,
        "%.3g" % pl if pl is not None else "-",
        "%.3g" % he if he is not None else "-")
    return (method, "%.3f" % math.exp(beta), ("%.3g" % p), note)


def mr_hint(a2b, b2a):
    """仅按 P<0.05 给出粗略方向提示；正式判定必须由 SKILL.md 第 5 步人工/模型复核。"""

    def sig(x):
        try:
            return float(x[2]) < 0.05
        except (TypeError, ValueError, IndexError):
            return False

    sa, sb = sig(a2b), sig(b2a)
    if sa and sb:
        return "双向均显著 → 可能为双向因果或存在共同危险因素/遗传相关，需重点排查"
    if sa:
        return "仅 A→B 显著 → 支持 A→B 单向"
    if sb:
        return "仅 B→A 显著 → 支持 B→A 单向"
    return "双向均不显著 → MR 不支持因果，文献若声称因果应降级"


SCORE_COLS = ["S1_temporality", "S2_phenotype", "S3_groups", "S4_interaction", "S5_bias",
              "S6_conclusion", "total", "direction", "tier", "note"]


def main():
    ap = argparse.ArgumentParser(description="PubMed 共病方向性文献检索（仅人类）")
    ap.add_argument("--a", help="疾病 A（暴露）名称")
    ap.add_argument("--b", help="疾病 B（结局）名称")
    ap.add_argument("--a-terms", help="A 的同义词，用 | 分隔")
    ap.add_argument("--b-terms", help="B 的同义词，用 | 分隔")
    ap.add_argument("--direction", choices=["both", "a2b", "b2a"], default="both")
    ap.add_argument("--years", default="", help="如 2000:2026")
    ap.add_argument("--retmax", type=int, default=200)
    ap.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    ap.add_argument("--out", default="comorbidity_candidates", help="输出前缀（不含扩展名）")
    ap.add_argument("--log", default="", help="检索式留痕 md 路径，默认 <out>_query_log.md")
    ap.add_argument("--no-design-filter", action="store_true", help="不加队列/时序设计词过滤")
    ap.add_argument("--allow-case-reports", action="store_true")
    ap.add_argument("--no-reviews", action="store_true")
    ap.add_argument("--mr-stats", default="", help="可选：mr_stats.js 路径，用于 MR 方向交叉验证")
    ap.add_argument("--mesh", default="", help="只查 MeSH 主题词后退出")
    args = ap.parse_args()

    key = os.environ.get("NCBI_API_KEY")

    if args.mesh:
        mesh_lookup(args.mesh, args.email, key)
        return

    if not args.a or not args.b:
        ap.error("需要 --a 与 --b（或用 --mesh 只查主题词）")

    a_terms = (args.a_terms or args.a).split("|")
    b_terms = (args.b_terms or args.b).split("|")
    A, B = block(a_terms), block(b_terms)
    years = args.years.split(":") if args.years else []
    mindate = years[0] if len(years) == 2 else None
    maxdate = years[1] if len(years) == 2 else None

    tail = " AND humans[Filter]"
    if not args.no_design_filter:
        tail += " AND " + DESIGN
    if not args.allow_case_reports:
        tail += " " + EXCL_DEFAULT
    if args.no_reviews:
        tail += " " + EXCL_REVIEWS

    queries = [("general(两病共现)", "%s AND %s%s" % (A, B, tail))]
    if args.direction in ("both", "a2b"):
        queries.append(("A→B 方向短语", "%s AND %s AND %s%s" % (A, B, directional(a_terms[0], b_terms[0]), tail)))
    if args.direction in ("both", "b2a"):
        queries.append(("B→A 方向短语", "%s AND %s AND %s%s" % (A, B, directional(b_terms[0], a_terms[0]), tail)))

    print("疾病 A = %s ｜ 疾病 B = %s ｜ 物种过滤 = humans" % (args.a, args.b))
    result = {}
    branch = {}
    uniq = []
    seen = set()
    for label, q in queries:
        cnt, ids = esearch(q, args.retmax, args.email, key, mindate, maxdate)
        print("  [%s] 命中 %d，取回 %d" % (label, cnt, len(ids)))
        result[label] = (q, cnt, len(ids))
        for i in ids:
            branch.setdefault(i, []).append(label)
            if i not in seen:
                seen.add(i)
                uniq.append(i)
    print("去重后候选 = %d 篇，正在抓取题录…" % len(uniq))

    rows = parse(efetch(uniq, args.email, key))
    by_pmid = {r["pmid"]: r for r in rows}
    ordered = [by_pmid[p] for p in uniq if p in by_pmid]
    for r in ordered:
        r["hit_in"] = ";".join(sorted(set(branch.get(r["pmid"], []))))

    # MR 交叉验证（可选）
    mr_a = mr_b = ("", "", "", "")
    if args.mr_stats and os.path.exists(args.mr_stats):
        stats = load_mr(args.mr_stats)
        mr_a = mr_lookup(stats, args.a, args.b)
        mr_b = mr_lookup(stats, args.b, args.a)
        print("MR A→B: %s ｜ MR B→A: %s" % (mr_a[:3], mr_b[:3]))
        print("MR 方向提示: %s" % mr_hint(mr_a, mr_b))
    for r in ordered:
        r["mr_a2b_method"], r["mr_a2b_or"], r["mr_a2b_p"], r["mr_a2b_note"] = mr_a
        r["mr_b2a_method"], r["mr_b2a_or"], r["mr_b2a_p"], r["mr_b2a_note"] = mr_b

    cols = ["pmid", "doi", "year", "journal", "pubtypes", "title", "abstract", "url", "hit_in",
            "mr_a2b_method", "mr_a2b_or", "mr_a2b_p", "mr_a2b_note",
            "mr_b2a_method", "mr_b2a_or", "mr_b2a_p", "mr_b2a_note"] + SCORE_COLS
    out_csv = args.out + "_candidates.csv"
    os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in ordered:
            for c in SCORE_COLS:
                r.setdefault(c, "")
            w.writerow(r)
    print("候选表 →", out_csv)

    log = args.log or (args.out + "_query_log.md")
    with open(log, "w", encoding="utf-8") as f:
        f.write("# 检索留痕（PRISMA identification）\n\n")
        f.write("- 时间：%s\n" % datetime.now().strftime("%Y-%m-%d %H:%M"))
        f.write("- 疾病 A：%s（同义词：%s）\n" % (args.a, " / ".join(a_terms)))
        f.write("- 疾病 B：%s（同义词：%s）\n" % (args.b, " / ".join(b_terms)))
        f.write("- 年份：%s ｜ 物种过滤：humans ｜ retmax=%d\n\n" % (args.years or "不限", args.retmax))
        for label, (q, cnt, got) in result.items():
            f.write("## %s\n\n命中 **%d** 条，取回 %d 条。\n\n```\n%s\n```\n\n" % (label, cnt, got, q))
        f.write("## 去重合计\n\n去重后候选 **%d** 篇 → 进入 L0 硬门槛筛选。\n" % len(ordered))
        if args.mr_stats:
            f.write("\n## MR 交叉验证\n\n- A→B：%s（OR=%s，P=%s，%s）\n- B→A：%s（OR=%s，P=%s，%s）\n"
                    % (mr_a[0], mr_a[1], mr_a[2], mr_a[3], mr_b[0], mr_b[1], mr_b[2], mr_b[3]))
            f.write("- 粗略方向提示：%s（仅供交叉验证，不能替代文献自身证据）\n" % mr_hint(mr_a, mr_b))
    print("检索留痕 →", log)


if __name__ == "__main__":
    main()
