"""BED12 转录本迁移结果导出：BED12 文本与预览 JSON。

输出严格是 12 列 BED12（无 BOM、tab 分隔，与 UCSC 工具链一致）：

- 只导出 **mapped / ambiguous** 行；歧义行的每个完整候选各占一行
  （全部导出，不自动选取）；unmapped 行没有任何完整路径，**绝不**
  生成伪 BED12，仅出现在预览中并带原因。
- 每个候选按**目标坐标**重排块（blockStarts 升序、与方向无关），
  chromStart/chromEnd 为目标块的最小左端点 / 最大右端点；
- thickStart/thickEnd 依据映射后的外显子内编码部分重算，非编码记录
  保持零长 thick（thickStart = thickEnd = chromStart）；
- strand 采用迁移后的方向（必要时已翻转）；name / score / itemRgb
  原样保留。
"""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

from .bed12 import Bed12Record
from .engine12 import Bed12Candidate, Bed12Row, MappedBlock


def candidate_to_bed12(rec: Bed12Record, cand: Bed12Candidate) -> str:
    """把一个完整路径候选渲染为一行 12 列 BED12。"""
    ordered = sorted(cand.blocks, key=lambda b: b.dest_start)
    chrom_start = ordered[0].dest_start
    chrom_end = ordered[-1].dest_end

    sizes = [b.dest_end - b.dest_start for b in ordered]
    starts = [b.dest_start - chrom_start for b in ordered]

    if rec.coding:
        thick_start = cand.thick_dest_start
        thick_end = cand.thick_dest_end
    else:
        # 非编码记录：零长 thick（UCSC 惯例 thickStart=thickEnd=chromStart）
        thick_start = thick_end = chrom_start

    fields = [
        cand.dest_chrom,
        chrom_start,
        chrom_end,
        rec.name,
        rec.score,
        cand.output_strand,
        thick_start,
        thick_end,
        rec.item_rgb,
        len(ordered),
        ",".join(str(n) for n in sizes) + ",",
        ",".join(str(n) for n in starts) + ",",
    ]
    return "\t".join(str(x) for x in fields)


def build_bed12(records: list[Bed12Record], rows: list[Bed12Row]) -> str:
    """生成 BED12 下载文本：只含完整映射与歧义候选，每行 12 列。"""
    by_line = {r.line_no: r for r in rows}
    lines: list[str] = []
    for rec in records:
        row = by_line[rec.line_no]
        if not row.candidates:
            continue
        for cand in row.candidates:
            lines.append(candidate_to_bed12(rec, cand))
    return "\n".join(lines) + ("\n" if lines else "")


def _block_dict(b: MappedBlock) -> dict[str, Any]:
    d = asdict(b)
    # 预览里 mid_* 对单段为 null（dataclass 默认 None，asdict 已处理）
    return d


def _candidate_dict(cand: Bed12Candidate) -> dict[str, Any]:
    return {
        "chain1_id": cand.chain1_id,
        "chain2_id": cand.chain2_id,
        "dest_chrom": cand.dest_chrom,
        "dest_strand": cand.dest_strand,
        "output_strand": cand.output_strand,
        "mid_chrom": cand.mid_chrom,
        "mid_strand": cand.mid_strand,
        "forward": cand.forward,
        "thick_dest_start": cand.thick_dest_start,
        "thick_dest_end": cand.thick_dest_end,
        "blocks": [_block_dict(b) for b in cand.blocks],
    }


def rows_to_preview12(rows: list[Bed12Row]) -> dict[str, Any]:
    """BED12 预览 JSON。

    计数仍返回 mapped/partial/ambiguous/unmapped 四键（BED12 整体迁移
    没有 partial，恒为 0），便于前端沿用同一套统计展示。
    """
    counts = {"mapped": 0, "partial": 0, "ambiguous": 0, "unmapped": 0}
    out_rows = []
    for r in rows:
        counts[r.status] += 1
        out_rows.append({
            "line_no": r.line_no,
            "src_chrom": r.src_chrom,
            "src_start": r.src_start,
            "src_end": r.src_end,
            "name": r.name,
            "strand": r.strand,
            "block_count": r.block_count,
            "status": r.status,
            "reason": r.reason,
            "candidates": [_candidate_dict(c) for c in r.candidates],
        })
    return {"total": len(rows), "counts": counts, "rows": out_rows}
