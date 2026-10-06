"""BED12 转录本整体迁移引擎。

与 BED3 逐碱基/逐区间迁移不同，转录本（BED12）必须作为**整体**迁移：
只有一条路径把**全部外显子完整**映射到**同一目标染色体**、各外显子在
目标坐标系上**方向一致**、且按源顺序保持同向或**整体反转**、目标块
**互不重叠**时，该路径才成为候选；任何外显子落入缺口/未对齐区域、被
拆成多段、目标染色体或方向不一致、目标块顺序冲突或相互重叠，都只淘汰
该条路径，不产出“半张”BED12。

- 单段：对每个外显子求 ``LiftOverEngine.candidates_for``，只保留“恰好
  一个片段且完整覆盖”的 chain 命中；所有外显子命中同一条 chain 并通过
  整体校验才构成候选 chain。
- 双段：每个外显子先走第一段、再走第二段，同样要求每段完整（恰好一片
  段且全长覆盖）；按 (chain1, chain2) 组合归并为整条转录本路径。

候选构造时按**目标坐标**重排块、重算 chromStart/End 与 block 列；
thickStart/thickEnd 依据“映射后的外显子内编码部分”重新计算（所有编码
子区间映射目标的最小左端点 / 最大右端点），非编码记录保持零长 thick；
反向迁移翻转 strand（源链 '.' 视同 '+'）。

多个完整路径全部保留并标为歧义，**不自动选取**；没有任何完整路径时
标记 unmapped 并给出原因，绝不导出伪 BED12（由接口层负责不提供下载）。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .bed12 import Bed12Record
from .chain import Chain
from .engine import LiftOverEngine


@dataclass
class MappedBlock:
    """一个源外显子经一条路径映射后的目标块（坐标零基左闭右开）。"""

    src_start: int
    src_end: int
    dest_start: int
    dest_end: int
    mid_start: int | None = None    # 双段：中间坐标（第一段目标）
    mid_end: int | None = None


@dataclass
class Bed12Candidate:
    """一条完整的转录本迁移路径。"""

    chain1_id: str
    chain2_id: str = ""             # 空串表示单段结果
    dest_chrom: str = ""
    dest_strand: str = ""           # 路径目标链（chain 目标侧方向）
    output_strand: str = "+"        # 迁移后转录本的 strand（必要时已翻转）
    mid_chrom: str = ""
    mid_strand: str = ""
    forward: bool = True            # 目标块是否按源外显子顺序同向排列
    blocks: list[MappedBlock] = field(default_factory=list)
    # 映射后编码区（thick）范围；非编码记录保持零长（== chromStart）。
    # 依据“映射后的外显子内编码部分”的最小左端点 / 最大右端点计算，
    # 与是否整体反转无关，输出时 thickStart < thickEnd。
    thick_dest_start: int = 0
    thick_dest_end: int = 0


@dataclass
class Bed12Row:
    line_no: int
    src_chrom: str
    src_start: int
    src_end: int
    name: str
    strand: str
    block_count: int
    status: str                     # mapped / ambiguous / unmapped（无 partial）
    reason: str = ""
    candidates: list[Bed12Candidate] = field(default_factory=list)


def _flip(strand: str) -> str:
    """整体反转时翻转方向；'.'（未知）在按 '+' 处理后翻转成 '-'。"""
    return "-" if strand != "-" else "+"


def _validate_blocks(
    blocks: list[MappedBlock],
) -> tuple[bool, str]:
    """整体校验目标块：按源顺序整体同向/反转、目标块互不重叠。

    入参 blocks 按源外显子顺序排列（外显子本身在源坐标上严格有序，且同
    一 chain/路径的目标染色体与目标链恒定）。返回 (是否通过, 失败原因)。
    """
    if not blocks:
        return False, "没有任何映射块"

    # 用源顺序相邻块的目标坐标比较，确定每一对的方向符号
    signs: list[int] = []
    for i in range(1, len(blocks)):
        prev, cur = blocks[i - 1], blocks[i]
        if cur.dest_start > prev.dest_start:
            signs.append(1)
        elif cur.dest_start < prev.dest_start:
            signs.append(-1)
        else:
            return False, (
                f"源外显子 {i} 与 {i + 1} 的目标块起点相同"
                f"（{cur.dest_start}），顺序无法保持，该路径淘汰"
            )

    # 所有相邻对方向符号必须一致（整体同向或整体反转）
    if signs and any(s != signs[0] for s in signs):
        return False, (
            "目标块相对源外显子顺序既有正向又有反向排列（顺序冲突），"
            "不满足按源顺序保持或整体反转，该路径淘汰"
        )
    forward = (not signs) or signs[0] == 1

    # 按目标坐标排序后，目标块必须互不重叠（相接允许）
    ordered = sorted(blocks, key=lambda b: (b.dest_start, b.dest_end))
    for i in range(1, len(ordered)):
        if ordered[i].dest_start < ordered[i - 1].dest_end:
            return False, (
                f"目标块 {ordered[i - 1].dest_start}-{ordered[i - 1].dest_end} "
                f"与 {ordered[i].dest_start}-{ordered[i].dest_end} 在目标染色体上"
                f"重叠，该路径淘汰"
            )
    return True, ""


class TranscriptLiftoverEngine:
    """BED12 转录本整体迁移（单段或双段 chain）。"""

    def __init__(
        self,
        chains1: list[Chain],
        side1: str,
        chains2: list[Chain] | None = None,
        side2: str | None = None,
    ) -> None:
        self.engine1 = LiftOverEngine(chains1, side1)
        self.double = chains2 is not None
        self.engine2 = LiftOverEngine(chains2, side2) if self.double else None
        self.side1 = side1
        self.side2 = side2

    # ------------------------------------------------------------------
    # 单段
    # ------------------------------------------------------------------

    def _exon_candidates_single(
        self, rec: Bed12Record, lo: int, hi: int
    ) -> tuple[dict, list[str]]:
        """单个外显子的“完整命中” chain 映射。

        返回 (uid → 映射信息, 失败说明)。只保留恰好一个片段且完整覆盖
        整个外显子的命中；多片段（跨过缺口）或未覆盖全长的命中不作为
        整体迁移候选。
        """
        length = hi - lo
        cands = self.engine1.candidates_for(rec.chrom, lo, hi)
        if not cands:
            if not self.engine1.has_chrom(rec.chrom):
                return {}, [
                    f"外显子 [{lo},{hi})：染色体 {rec.chrom} 在 chain 源侧不存在"
                ]
            return {}, [
                f"外显子 [{lo},{hi})：完全落在所有 chain 的对齐范围之外"
                f"（链内缺口或未覆盖区域），外显子未完整映射"
            ]
        out = {}
        notes: list[str] = []
        for cand in cands:
            if cand.covered != length or len(cand.fragments) != 1:
                notes.append(
                    f"外显子 [{lo},{hi})：chain {cand.chain_id} 仅覆盖 "
                    f"{cand.covered}/{length} bp、命中 {len(cand.fragments)} 段，"
                    f"外显子未完整映射，该路径淘汰"
                )
                continue
            out[cand.chain_uid] = {
                "chain_id": cand.chain_id,
                "dest_chrom": cand.dest_chrom,
                "dest_strand": cand.dest_strand,
                "frag": cand.fragments[0],
            }
        if not out:
            return {}, notes
        return out, notes

    def convert_single(self, rec: Bed12Record) -> Bed12Row:
        row = Bed12Row(
            line_no=rec.line_no,
            src_chrom=rec.chrom,
            src_start=rec.start,
            src_end=rec.end,
            name=rec.name,
            strand=rec.strand,
            block_count=rec.block_count,
            status="unmapped",
        )
        exons = rec.exons()

        # 每个外显子的完整命中（uid → 片段信息）
        per_exon: list[dict] = []
        all_notes: list[str] = []
        for lo, hi in exons:
            hit, notes = self._exon_candidates_single(rec, lo, hi)
            all_notes.extend(notes)
            per_exon.append(hit)

        # 出现在每个外显子完整命中集合里的 chain 才可能是完整路径。
        # 以第一外显子的命中发现顺序（dict 插入序）为准遍历，保证多候选时
        # 路径序号稳定（set 迭代顺序不确定，不能直接遍历交集）。
        common = set(per_exon[0])
        for hit in per_exon[1:]:
            common &= set(hit)
        ordered_uids = [uid for uid in per_exon[0] if uid in common]

        candidates: list[Bed12Candidate] = []
        for uid in ordered_uids:
            info0 = per_exon[0][uid]
            dest_chrom = info0["dest_chrom"]
            dest_strand = info0["dest_strand"]
            chain_id = info0["chain_id"]
            blocks: list[MappedBlock] = []
            for i, (lo, hi) in enumerate(exons):
                info = per_exon[i][uid]
                frag = info["frag"]
                blocks.append(MappedBlock(
                    src_start=lo, src_end=hi,
                    dest_start=frag.dest_start, dest_end=frag.dest_end,
                ))
            ok, why = _validate_blocks(blocks)
            if not ok:
                all_notes.append(f"chain {chain_id}：{why}")
                continue
            forward = self.engine1.relative_forward(uid) if len(blocks) == 1 \
                else all(
                    blocks[i].dest_start < blocks[i + 1].dest_start
                    for i in range(len(blocks) - 1)
                )
            t0, t1 = self._thick_single(rec, uid, blocks)
            candidates.append(Bed12Candidate(
                chain1_id=chain_id,
                dest_chrom=dest_chrom,
                dest_strand=dest_strand,
                output_strand=rec.strand if forward else _flip(rec.strand),
                forward=forward,
                blocks=blocks,
                thick_dest_start=t0,
                thick_dest_end=t1,
            ))

        return self._finalize(row, candidates, all_notes)

    # ------------------------------------------------------------------
    # 双段
    # ------------------------------------------------------------------

    def convert_double(self, rec: Bed12Record) -> Bed12Row:
        row = Bed12Row(
            line_no=rec.line_no,
            src_chrom=rec.chrom,
            src_start=rec.start,
            src_end=rec.end,
            name=rec.name,
            strand=rec.strand,
            block_count=rec.block_count,
            status="unmapped",
        )
        exons = rec.exons()

        # 每外显子：key=(uid1, uid2) → 映射信息
        per_exon: list[dict[tuple[int, int], dict]] = []
        all_notes: list[str] = []

        for ei, (lo, hi) in enumerate(exons, start=1):
            length = hi - lo
            cands1 = self.engine1.candidates_for(rec.chrom, lo, hi)
            if not cands1:
                if not self.engine1.has_chrom(rec.chrom):
                    all_notes.append(
                        f"外显子 {ei} [{lo},{hi})：第一段无映射，染色体 "
                        f"{rec.chrom} 在 chain 文件 1 源侧不存在"
                    )
                else:
                    all_notes.append(
                        f"外显子 {ei} [{lo},{hi})：第一段未完整映射"
                        f"（落在 chain 文件 1 的缺口/未覆盖区域），路径淘汰"
                    )
                per_exon.append({})
                continue

            exon_hits: dict[tuple[int, int], dict] = {}
            exon_alive1 = 0
            for c1 in cands1:
                if c1.covered != length or len(c1.fragments) != 1:
                    all_notes.append(
                        f"外显子 {ei} [{lo},{hi})：chain 文件 1 的 chain "
                        f"{c1.chain_id} 仅覆盖 {c1.covered}/{length} bp、"
                        f"{len(c1.fragments)} 段，外显子未完整映射，路径淘汰"
                    )
                    continue
                f1 = c1.fragments[0]
                exon_alive1 += 1
                cands2 = self.engine2.candidates_for(
                    c1.dest_chrom, f1.dest_start, f1.dest_end
                )
                if not cands2:
                    if not self.engine2.has_chrom(c1.dest_chrom):
                        all_notes.append(
                            f"外显子 {ei} [{lo},{hi})：经第一段 chain "
                            f"{c1.chain_id} 到中间染色体 {c1.dest_chrom} 后，"
                            f"该染色体在 chain 文件 2 源侧不存在，路径淘汰"
                        )
                    else:
                        all_notes.append(
                            f"外显子 {ei} [{lo},{hi})：经第一段 chain "
                            f"{c1.chain_id} 到中间片段 "
                            f"{c1.dest_chrom}:[{f1.dest_start},{f1.dest_end}) "
                            f"后，第二段未完整映射（缺口/未覆盖区域），路径淘汰"
                        )
                    continue
                for c2 in cands2:
                    if c2.covered != length or len(c2.fragments) != 1:
                        all_notes.append(
                            f"外显子 {ei} [{lo},{hi})：第二段 chain {c2.chain_id} "
                            f"仅覆盖 {c2.covered}/{length} bp、{len(c2.fragments)} 段，"
                            f"外显子未完整映射，路径淘汰"
                        )
                        continue
                    f2 = c2.fragments[0]
                    exon_hits[(c1.chain_uid, c2.chain_uid)] = {
                        "chain1_id": c1.chain_id,
                        "chain2_id": c2.chain_id,
                        "mid_chrom": c1.dest_chrom,
                        "mid_strand": c1.dest_strand,
                        "mid_start": f1.dest_start,
                        "mid_end": f1.dest_end,
                        "dest_chrom": c2.dest_chrom,
                        "dest_strand": c2.dest_strand,
                        "dest_start": f2.dest_start,
                        "dest_end": f2.dest_end,
                        "uid1": c1.chain_uid,
                        "uid2": c2.chain_uid,
                    }
            if not exon_hits and exon_alive1:
                all_notes.append(
                    f"外显子 {ei} [{lo},{hi})：第一段可完整映射，但没有任何"
                    f"第二段 chain 能完整承接，未形成两段都完整的路径"
                )
            per_exon.append(exon_hits)

        common = set(per_exon[0])
        for hit in per_exon[1:]:
            common &= set(hit)
        ordered_keys = [k for k in per_exon[0] if k in common]

        candidates: list[Bed12Candidate] = []
        for key in ordered_keys:
            first = per_exon[0][key]
            blocks: list[MappedBlock] = []
            bad = False
            for i, (lo, hi) in enumerate(exons):
                info = per_exon[i][key]
                if (info["dest_chrom"] != first["dest_chrom"]
                        or info["dest_strand"] != first["dest_strand"]):
                    bad = True
                    break
                blocks.append(MappedBlock(
                    src_start=lo, src_end=hi,
                    dest_start=info["dest_start"], dest_end=info["dest_end"],
                    mid_start=info["mid_start"], mid_end=info["mid_end"],
                ))
            if bad:
                all_notes.append(
                    f"路径（chain {first['chain1_id']} × chain "
                    f"{first['chain2_id']}）各外显子目标染色体/方向不一致，淘汰"
                )
                continue
            ok, why = _validate_blocks(blocks)
            if not ok:
                all_notes.append(
                    f"路径（chain {first['chain1_id']} × chain "
                    f"{first['chain2_id']}）：{why}"
                )
                continue
            if len(blocks) == 1:
                forward = (
                    self.engine1.relative_forward(key[0])
                    == self.engine2.relative_forward(key[1])
                )
            else:
                forward = all(
                    blocks[i].dest_start < blocks[i + 1].dest_start
                    for i in range(len(blocks) - 1)
                )
            t0, t1 = self._thick_double(rec, key[0], key[1], blocks)
            candidates.append(Bed12Candidate(
                chain1_id=first["chain1_id"],
                chain2_id=first["chain2_id"],
                dest_chrom=first["dest_chrom"],
                dest_strand=first["dest_strand"],
                mid_chrom=first["mid_chrom"],
                mid_strand=first["mid_strand"],
                output_strand=rec.strand if forward else _flip(rec.strand),
                forward=forward,
                blocks=blocks,
                thick_dest_start=t0,
                thick_dest_end=t1,
            ))

        return self._finalize(row, candidates, all_notes, double=True)

    # ------------------------------------------------------------------
    # 汇总
    # ------------------------------------------------------------------

    @staticmethod
    def _summarize(notes: list[str], cap: int = 4) -> str:
        dedup = list(dict.fromkeys(notes))
        if len(dedup) > cap:
            return "；".join(dedup[:cap]) + f"；等 {len(dedup)} 处"
        return "；".join(dedup)

    def _finalize(
        self,
        row: Bed12Row,
        candidates: list[Bed12Candidate],
        notes: list[str],
        double: bool = False,
    ) -> Bed12Row:
        if not candidates:
            summary = self._summarize(notes)
            row.reason = (
                "没有任何完整路径：转录本必须全部外显子完整映射到同一目标"
                "染色体、方向一致、按源顺序保持或整体反转且目标块互不重叠；"
                + (summary or "所有 chain 均未给出满足条件的映射")
            )
            return row
        if len(candidates) > 1:
            row.status = "ambiguous"
            ids = "；".join(
                (f"chain {c.chain1_id}" + (f" × chain {c.chain2_id}" if double else ""))
                + f" → {c.dest_chrom}({c.output_strand})"
                for c in candidates
            )
            row.reason = (
                f"存在 {len(candidates)} 条完整路径（{ids}），按要求全部标为歧义"
                f"并导出，不自动选取；缺口或冲突的路径已被淘汰"
            )
            row.candidates = candidates
            return row
        row.status = "mapped"
        c = candidates[0]
        row.reason = (
            f"全部 {row.block_count} 个外显子经 chain {c.chain1_id}"
            + (f" × chain {c.chain2_id}" if double else "")
            + f" 完整映射到 {c.dest_chrom}，目标块按"
            + ("源顺序排列" if c.forward else "源顺序整体反转（strand 已翻转）")
        )
        row.candidates = candidates
        return row

    def convert(self, rec: Bed12Record) -> Bed12Row:
        return self.convert_double(rec) if self.double else self.convert_single(rec)

    def convert_all(self, records: list[Bed12Record]) -> list[Bed12Row]:
        return [self.convert(r) for r in records]

    # ------------------------------------------------------------------
    # 编码区（thick）映射
    # ------------------------------------------------------------------

    def _thick_single(
        self, rec: Bed12Record, uid: int, blocks: list[MappedBlock]
    ) -> tuple[int, int]:
        """单段：把源编码部分（各外显子与 [thickStart,thickEnd) 的交集）
        逐段映射到目标坐标，返回映射区间的并包络 (min 左端点, max 右端点)。

        路径已保证每个外显子在该 chain 上整体落在单个对齐块内，因此编码
        子区间必然完整映射，不需要再处理拆分；这里仍取所有映射片段的包络。
        """
        if not rec.coding:
            first = min(blocks, key=lambda b: b.dest_start)
            return first.dest_start, first.dest_start
        lo_dst: int | None = None
        hi_dst = 0
        for lo, hi in rec.exons():
            c_lo = max(lo, rec.thick_start)
            c_hi = min(hi, rec.thick_end)
            if c_lo >= c_hi:
                continue
            for frag in self.engine1.map_uid(uid, c_lo, c_hi):
                a, b = frag.dest_start, frag.dest_end
                lo_dst = a if lo_dst is None else min(lo_dst, a)
                hi_dst = max(hi_dst, b)
        if lo_dst is None:
            # 理论不可达：编码端点已校验在外显子内
            raise RuntimeError("编码区映射后为空")
        return lo_dst, hi_dst

    def _thick_double(
        self,
        rec: Bed12Record,
        uid1: int,
        uid2: int,
        blocks: list[MappedBlock],
    ) -> tuple[int, int]:
        """双段：源编码部分经第一段到中间坐标，再经第二段到最终坐标。"""
        if not rec.coding:
            first = min(blocks, key=lambda b: b.dest_start)
            return first.dest_start, first.dest_start
        lo_dst: int | None = None
        hi_dst = 0
        for lo, hi in rec.exons():
            c_lo = max(lo, rec.thick_start)
            c_hi = min(hi, rec.thick_end)
            if c_lo >= c_hi:
                continue
            for m1 in self.engine1.map_uid(uid1, c_lo, c_hi):
                for m2 in self.engine2.map_uid(uid2, m1.dest_start, m1.dest_end):
                    a, b = m2.dest_start, m2.dest_end
                    lo_dst = a if lo_dst is None else min(lo_dst, a)
                    hi_dst = max(hi_dst, b)
        if lo_dst is None:
            raise RuntimeError("编码区双段映射后为空")
        return lo_dst, hi_dst
