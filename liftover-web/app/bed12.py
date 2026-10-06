"""BED12（转录本）严格解析器。

BED12 的 12 列::

    chrom, start, end, name, score, strand,
    thickStart, thickEnd, itemRgb,
    blockCount, blockSizes, blockStarts

坐标一律零基、左闭右开。解析采用“收集全部错误”策略：逐行校验，坏行
不中断，调用方一次性拿到所有带物理行号的问题；错误列表非空时整次转换
被拒绝，绝不产生任何可下载结果。

校验内容：

- 必须恰好 12 列（少于 3 列的裸 BED3 行也不接受：整体迁移接口只收 BED12）；
- start < end、score ∈ [0,1000]、strand ∈ {'+', '-', '.'}；
- blockCount 与 blockSizes / blockStarts 中的整数个数一致，块长为正整数；
- 块起点严格递增、块之间互不重叠、首块贴 start、末块贴 end（外显子边界）；
- thickStart/thickEnd 落在 [start, end] 内且 thickStart <= thickEnd；
  编码记录（thickStart < thickEnd）的两个端点必须落在某个外显子内部
  （闭区间，端点允许正好落在外显子边界）；非编码记录（两者相等，
  常见写法均为 start，也允许相等的零长 thick 落在其他位置）原样保留，
  迁移后仍输出零长 thick；
- itemRgb 为 ``0`` 或 ``r,g,b``（各分量 0~255 整数）。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .bed import _parse_uint

BED12_COLUMNS = 12
VALID_STRANDS = "+-."


@dataclass(frozen=True)
class Bed12Record:
    line_no: int           # 物理文件行号，从 1 开始
    chrom: str
    start: int
    end: int
    name: str
    score: str             # 0..1000 整数（按 UCSC 语义为整数，保留原始文本）
    strand: str            # '+' / '-' / '.'
    thick_start: int
    thick_end: int
    item_rgb: str          # 原始 itemRgb 文本（"0" 或 "r,g,b"）
    block_count: int
    block_sizes: tuple[int, ...]
    block_starts: tuple[int, ...]
    raw: str

    @property
    def coding(self) -> bool:
        """编码记录：thickStart < thickEnd；相等即非编码（零长 thick）。"""
        return self.thick_start < self.thick_end

    def exons(self) -> list[tuple[int, int]]:
        """返回各外显子在转录本染色体坐标上的 [abs_start, abs_end)。"""
        return [
            (self.start + s, self.start + s + size)
            for s, size in zip(self.block_starts, self.block_sizes)
        ]


@dataclass
class Bed12ParseResult:
    records: list[Bed12Record] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)


def _parse_int(text: str, line_no: int, field_name: str) -> int:
    """解析十进制整数（允许可选负号）；BED 坐标字段必须非负。"""
    t = text.strip()
    sign = 1
    if t[:1] in "+-":
        if t[:1] == "-":
            sign = -1
        t = t[1:]
    if not t or not t.isdigit():
        raise ValueError(
            f"BED 第 {line_no} 行: {field_name} 必须是整数，得到 {text!r}"
        )
    return sign * int(t)


def _parse_block_list(
    text: str,
    line_no: int,
    field_name: str,
    expected: int,
) -> list[int]:
    """解析 blockSizes/blockStarts：UCSC 习惯在末尾多写一个逗号。

    空串、空项（连续逗号）、非整数都视为错误；返回前还须与 blockCount
    个数一致（由调用方/此处统一报错）。
    """
    items = text.split(",")
    if items and items[-1] == "":
        items = items[:-1]
    if any(x == "" for x in items):
        raise ValueError(
            f"BED 第 {line_no} 行: {field_name} 含空项，应为 {expected} 个"
            f"以逗号分隔的非负整数（可带末尾逗号），得到 {text!r}"
        )
    vals = [_parse_uint(x, line_no, field_name) for x in items]
    if len(vals) != expected:
        raise ValueError(
            f"BED 第 {line_no} 行: {field_name} 有 {len(vals)} 项，"
            f"与 blockCount={expected} 不一致"
        )
    return vals


def parse_bed12(text: str) -> Bed12ParseResult:
    result = Bed12ParseResult()

    for idx, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(("#", "browser", "track")):
            continue

        fields = line.split()
        if len(fields) != BED12_COLUMNS:
            result.errors.append(
                f"BED 第 {idx} 行: BED12 必须恰好 12 列，实际 {len(fields)} 列"
                + ("（少于 3 列；整体迁移只接受完整 BED12 转录本行）"
                   if len(fields) < 3 else "")
            )
            continue

        (chrom, start_s, end_s, name, score_s, strand,
         thick_start_s, thick_end_s, item_rgb,
         block_count_s, block_sizes_s, block_starts_s) = fields

        ok = True

        try:
            start = _parse_uint(start_s, idx, "start")
            end = _parse_uint(end_s, idx, "end")
        except ValueError as exc:
            result.errors.append(str(exc))
            ok = False
            start, end = 0, 0
        else:
            if start >= end:
                result.errors.append(
                    f"BED 第 {idx} 行: start ({start}) 必须小于 end ({end})"
                )
                ok = False

        # score：0..1000 整数
        try:
            score = _parse_int(score_s, idx, "score")
        except ValueError as exc:
            result.errors.append(str(exc))
            ok = False
        else:
            if not 0 <= score <= 1000:
                result.errors.append(
                    f"BED 第 {idx} 行: score 必须在 0..1000 之间，得到 {score_s!r}"
                )
                ok = False

        if strand not in VALID_STRANDS:
            result.errors.append(
                f"BED 第 {idx} 行: strand 只能是 '+'、'-' 或 '.'，得到 {strand!r}"
            )
            ok = False

        try:
            thick_start = _parse_uint(thick_start_s, idx, "thickStart")
            thick_end = _parse_uint(thick_end_s, idx, "thickEnd")
        except ValueError as exc:
            result.errors.append(str(exc))
            ok = False
            thick_start = thick_end = 0
        else:
            if not (start <= thick_start <= end and start <= thick_end <= end):
                result.errors.append(
                    f"BED 第 {idx} 行: thickStart/thickEnd 必须落在 "
                    f"[{start},{end}] 内，得到 thickStart={thick_start}, "
                    f"thickEnd={thick_end}"
                )
                ok = False
            if thick_start > thick_end:
                result.errors.append(
                    f"BED 第 {idx} 行: thickStart ({thick_start}) 不能大于 "
                    f"thickEnd ({thick_end})"
                )
                ok = False

        # itemRgb：0 或 r,g,b（0~255）
        rgb_parts = item_rgb.split(",")
        if item_rgb != "0":
            rgb_err = False
            if len(rgb_parts) != 3:
                rgb_err = True
            else:
                for part in rgb_parts:
                    t = part.strip()
                    if not t.isdigit() or not 0 <= int(t) <= 255:
                        rgb_err = True
                        break
            if rgb_err:
                result.errors.append(
                    f"BED 第 {idx} 行: itemRgb 必须是 0 或 r,g,b（各分量 "
                    f"0~255 整数），得到 {item_rgb!r}"
                )
                ok = False

        # blockCount 与两个块列表
        try:
            block_count = _parse_uint(block_count_s, idx, "blockCount")
            if block_count == 0:
                result.errors.append(
                    f"BED 第 {idx} 行: blockCount 必须是正整数，得到 {block_count_s!r}"
                )
                ok = False
        except ValueError as exc:
            result.errors.append(str(exc))
            ok = False
            block_count = 0

        try:
            block_sizes = _parse_block_list(
                block_sizes_s, idx, "blockSizes", block_count
            )
            block_starts = _parse_block_list(
                block_starts_s, idx, "blockStarts", block_count
            )
        except ValueError as exc:
            result.errors.append(str(exc))
            ok = False
            block_sizes, block_starts = [], []

        # 外显子边界与顺序校验（块列表解析成功时）
        if ok and block_sizes and block_starts:
            span = end - start

            for i, (bs, bst) in enumerate(zip(block_sizes, block_starts)):
                if bs <= 0:
                    result.errors.append(
                        f"BED 第 {idx} 行: 第 {i + 1} 个块长为 0 或负数，"
                        f"外显子长度必须为正"
                    )
                    ok = False
                if bst < 0 or bst + bs > span:
                    result.errors.append(
                        f"BED 第 {idx} 行: 第 {i + 1} 个块 "
                        f"[start+{bst}, start+{bst + bs}) 超出转录本范围 "
                        f"[{start},{end})"
                    )
                    ok = False

            if block_starts and block_starts[0] != 0:
                result.errors.append(
                    f"BED 第 {idx} 行: 首个外显子必须从 start 开始"
                    f"（blockStarts 首项应为 0），实际为 {block_starts[0]}"
                )
                ok = False
            if block_sizes and block_starts and \
                    block_starts[-1] + block_sizes[-1] != span:
                result.errors.append(
                    f"BED 第 {idx} 行: 末个外显子必须到 end 结束"
                    f"（最后一个 blockStart + blockSize 应等于 {span}），"
                    f"实际为 {block_starts[-1] + block_sizes[-1]}"
                )
                ok = False

            # 起点严格递增（相等即重叠）；相邻块允许相接但不允许重叠
            for i in range(1, len(block_starts)):
                prev_end = block_starts[i - 1] + block_sizes[i - 1]
                if block_starts[i] <= block_starts[i - 1]:
                    result.errors.append(
                        f"BED 第 {idx} 行: blockStarts 必须严格递增，"
                        f"第 {i} 与第 {i + 1} 个块起点分别为 "
                        f"{block_starts[i - 1]}、{block_starts[i]}"
                    )
                    ok = False
                elif block_starts[i] < prev_end:
                    result.errors.append(
                        f"BED 第 {idx} 行: 第 {i} 与第 {i + 1} 个外显子重叠"
                        f"（前一块结束于相对位置 {prev_end}，后一块从 "
                        f"{block_starts[i]} 开始）"
                    )
                    ok = False

            # 编码区端点必须在外显子内（端点允许正好落在外显子边界）
            if thick_start < thick_end:
                exons = [
                    (start + s, start + s + size)
                    for s, size in zip(block_starts, block_sizes)
                ]
                for label, pos in (("thickStart", thick_start),
                                   ("thickEnd", thick_end)):
                    if not any(lo <= pos <= hi for lo, hi in exons):
                        result.errors.append(
                            f"BED 第 {idx} 行: 编码区端点 {label}={pos} "
                            f"不在任何外显子内（外显子之间是内含子，"
                            f"编码端点不能落在缺口处）"
                        )
                        ok = False

        if ok:
            result.records.append(
                Bed12Record(
                    line_no=idx,
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    score=score_s,
                    strand=strand,
                    thick_start=thick_start,
                    thick_end=thick_end,
                    item_rgb=item_rgb,
                    block_count=block_count,
                    block_sizes=tuple(block_sizes),
                    block_starts=tuple(block_starts),
                    raw=line,
                )
            )

    if not result.records and not result.errors:
        result.errors.append("BED 文件为空（没有任何有效数据行）")

    return result
