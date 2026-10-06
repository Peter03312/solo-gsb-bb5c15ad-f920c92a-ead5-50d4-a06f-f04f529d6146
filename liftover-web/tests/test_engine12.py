"""BED12 转录本整体迁移引擎测试（单段与双段）。

直接构造 Chain 对象（与 app.chain 解析产物同构），覆盖：全部外显子完整
映射、按目标坐标重排与范围重算、thick 重算（含负链反转）、非编码零长
thick、strand 翻转、缺口/跨块淘汰、染色体/方向不一致淘汰、顺序冲突与
目标块重叠淘汰、多 chain 歧义全保留、双段路径归并等。
"""

from app.bed12 import parse_bed12
from app.chain import Block, Chain
from app.engine12 import TranscriptLiftoverEngine
from app.report12 import build_bed12


def _chain(
    cid,
    blocks,
    t_name="chr1",
    t_size=1000,
    t_strand="+",
    t_span=None,
    q_name="chrA",
    q_size=1000,
    q_strand="+",
    q_span=None,
):
    """blocks: [(t0,t1,q0,q1), ...]（存储帧；正链即正向帧）。"""
    t_span = t_span or (blocks[0][0], blocks[-1][1])
    q_span = q_span or (blocks[0][2], blocks[-1][3])
    return Chain(
        chain_id=cid,
        t_name=t_name, t_size=t_size, t_strand=t_strand,
        t_start=t_span[0], t_end=t_span[1],
        q_name=q_name, q_size=q_size, q_strand=q_strand,
        q_start=q_span[0], q_end=q_span[1],
        score=1000,
        blocks=[Block(*b) for b in blocks],
        header_line=1,
    )


def _rec(line):
    r = parse_bed12(line)
    assert r.errors == [], r.errors
    return r.records[0]


# --------------------------------------------------------------------------
# 单段
# --------------------------------------------------------------------------

def test_single_forward_multi_exon_mapped():
    # chain 覆盖 [0,100) [200,300)，正向平移 +1000
    chain = _chain("1", [(0, 100, 1000, 1100), (200, 300, 1200, 1300)],
                   t_span=(0, 300), q_span=(1000, 1300))
    eng = TranscriptLiftoverEngine([chain], "target")
    # 外显子 [0,50)、[250,300)
    rec = _rec("chr1\t0\t300\tTX\t0\t+\t20\t280\t0\t2\t50,50,\t0,250,\n")
    row = eng.convert(rec)
    assert row.status == "mapped"
    cand = row.candidates[0]
    assert cand.dest_chrom == "chrA"
    assert cand.output_strand == "+"
    assert cand.forward is True
    # 块按目标坐标排序、范围重算
    assert [(b.dest_start, b.dest_end) for b in cand.blocks] == [
        (1000, 1050), (1250, 1300)
    ]
    # thick 依据映射后的编码部分
    assert (cand.thick_dest_start, cand.thick_dest_end) == (1020, 1280)


def test_single_exon_in_gap_is_unmapped_not_partial():
    chain = _chain("1", [(0, 100, 0, 100), (200, 300, 200, 300)],
                   t_span=(0, 300), q_span=(0, 300))
    eng = TranscriptLiftoverEngine([chain], "target")
    # 外显子 2 [120,150) 落在缺口
    rec = _rec("chr1\t0\t150\tTX\t0\t+\t0\t150\t0\t2\t100,30,\t0,120,\n")
    row = eng.convert(rec)
    assert row.status == "unmapped"
    assert row.candidates == []
    assert "缺口" in row.reason or "未覆盖" in row.reason


def test_exon_straddling_gap_not_a_path():
    # 外显子 [80,220) 跨过 [100,200) 缺口：单 chain 会拆成两个片段，淘汰
    chain = _chain("1", [(0, 100, 0, 100), (200, 300, 200, 300)],
                   t_span=(0, 300), q_span=(0, 300))
    eng = TranscriptLiftoverEngine([chain], "target")
    rec = _rec("chr1\t80\t220\tTX\t0\t+\t80\t220\t0\t1\t140,\t0,\n")
    row = eng.convert(rec)
    assert row.status == "unmapped"
    assert row.candidates == []


def test_reverse_chain_flips_strand_and_block_order():
    # target [300,400) <-> query 负链，q 存储 [0,100)，qSize=500
    chain = _chain(
        "2", [(300, 400, 0, 100)],
        t_span=(300, 400), q_name="chrB", q_size=500, q_strand="-",
        q_span=(0, 100),
    )
    eng = TranscriptLiftoverEngine([chain], "target")
    rec = _rec("chr1\t300\t400\tTX\t0\t+\t310\t390\t0\t1\t100,\t0,\n")
    row = eng.convert(rec)
    cand = row.candidates[0]
    # 负链归一化后 query 正向帧 [400,500)，target 300 -> q 499
    assert (cand.blocks[0].dest_start, cand.blocks[0].dest_end) == (400, 500)
    assert cand.forward is False
    assert cand.output_strand == "-"
    # 编码部分 [310,390) -> [410,490)
    assert (cand.thick_dest_start, cand.thick_dest_end) == (410, 490)


def test_reverse_chain_source_minus_becomes_plus():
    chain = _chain(
        "2", [(300, 400, 0, 100)],
        t_span=(300, 400), q_name="chrB", q_size=500, q_strand="-",
        q_span=(0, 100),
    )
    eng = TranscriptLiftoverEngine([chain], "target")
    rec = _rec("chr1\t300\t400\tTX\t0\t-\t310\t390\t0\t1\t100,\t0,\n")
    assert eng.convert(rec).candidates[0].output_strand == "+"


def test_two_chains_ambiguous_both_exported():
    c1 = _chain("1", [(300, 400, 0, 100)], t_span=(300, 400),
                q_name="chrA", q_span=(0, 100))
    c2 = _chain("2", [(300, 400, 50, 150)], t_span=(300, 400),
                q_name="chrC", q_span=(50, 150))
    eng = TranscriptLiftoverEngine([c1, c2], "target")
    rec = _rec("chr1\t300\t400\tTX\t0\t+\t320\t380\t0\t1\t100,\t0,\n")
    row = eng.convert(rec)
    assert row.status == "ambiguous"
    assert len(row.candidates) == 2
    assert {c.dest_chrom for c in row.candidates} == {"chrA", "chrC"}


def test_chain_missing_one_exon_is_eliminated_other_survives():
    # chain1 只覆盖外显子 1；chain2 覆盖两个外显子
    c1 = _chain("1", [(0, 100, 0, 100)], t_span=(0, 100), q_name="chrA",
                q_span=(0, 100))
    c2 = _chain("2", [(0, 100, 0, 100), (200, 300, 200, 300)],
                t_span=(0, 300), q_name="chrD", q_span=(0, 300))
    eng = TranscriptLiftoverEngine([c1, c2], "target")
    rec = _rec("chr1\t0\t300\tTX\t0\t+\t0\t300\t0\t2\t100,100,\t0,200,\n")
    row = eng.convert(rec)
    # 不是歧义：chain1 缺外显子 2 被淘汰，只剩 chain2 一条完整路径
    assert row.status == "mapped"
    assert len(row.candidates) == 1
    assert row.candidates[0].chain1_id == "2"
    assert row.candidates[0].dest_chrom == "chrD"


def test_noncoding_zero_length_thick_in_output():
    chain = _chain("1", [(0, 100, 1000, 1100)], t_span=(0, 100),
                   q_span=(1000, 1100))
    eng = TranscriptLiftoverEngine([chain], "target")
    rec = _rec("chr1\t0\t100\tNCRNA\t42\t.\t0\t0\t255,0,0\t1\t100,\t0,\n")
    row = eng.convert(rec)
    cand = row.candidates[0]
    assert cand.thick_dest_start == cand.thick_dest_end == 1000
    bed = build_bed12([rec], [row])
    fields = bed.splitlines()[0].split("\t")
    assert fields[6] == fields[7] == "1000"
    # name/score/itemRgb 保留
    assert fields[3] == "NCRNA"
    assert fields[4] == "42"
    assert fields[8] == "255,0,0"


def test_output_bed12_only_contains_complete_rows():
    c_good = _chain("1", [(0, 100, 0, 100)], t_span=(0, 100),
                    q_span=(0, 100))
    # 另一条 chain 只覆盖 200 起的区域，不影响第 1 行
    c_other = _chain("2", [(200, 300, 0, 100), (400, 450, 100, 150)],
                     t_span=(200, 450), q_name="chrX", q_span=(0, 150))
    eng = TranscriptLiftoverEngine([c_good, c_other], "target")
    r = parse_bed12(
        "chr1\t0\t100\tOK\t0\t+\t0\t100\t0\t1\t100,\t0,\n"
        # 外显子 2 [350,400) 落在 c_other 的链内缺口 [300,400) => 不完整
        "chr1\t200\t400\tBAD\t0\t+\t200\t400\t0\t2\t50,50,\t0,150,\n"
    )
    rows = eng.convert_all(r.records)
    assert rows[0].status == "mapped"
    assert rows[1].status == "unmapped"
    bed = build_bed12(r.records, rows)
    assert bed.count("\n") == 1
    assert "OK" in bed and "BAD" not in bed


def test_overlapping_target_blocks_reject_path():
    # chain 把外显子 1 [0,100) 映射到 [500,600)，外显子 2 [200,300) 也
    # 映射到 [550,650)（目标块重叠）：该路径必须淘汰
    chain = _chain(
        "9", [(0, 100, 500, 600), (200, 300, 550, 650)],
        t_span=(0, 300), q_span=(500, 650),
    )
    eng = TranscriptLiftoverEngine([chain], "target")
    rec = _rec("chr1\t0\t300\tTX\t0\t+\t0\t300\t0\t2\t100,100,\t0,200,\n")
    row = eng.convert(rec)
    assert row.status == "unmapped"
    assert row.candidates == []
    assert "重叠" in row.reason


# --------------------------------------------------------------------------
# 双段
# --------------------------------------------------------------------------

def test_double_full_path_with_mid_coordinates():
    # 段1：chr1 [0,100) -> chrA [1000,1100)
    c1 = _chain("1", [(0, 100, 1000, 1100)], t_span=(0, 100),
                q_name="chrA", q_span=(1000, 1100))
    # 段2：chrA [1000,1100) -> chrZ [0,100)
    c2 = _chain("10", [(1000, 1100, 0, 100)], t_name="chrA",
                t_span=(1000, 1100), q_name="chrZ", q_span=(0, 100))
    eng = TranscriptLiftoverEngine([c1], "target", [c2], "target")
    rec = _rec("chr1\t0\t100\tTX\t0\t+\t10\t90\t0\t1\t100,\t0,\n")
    row = eng.convert(rec)
    assert row.status == "mapped"
    cand = row.candidates[0]
    assert cand.chain1_id == "1" and cand.chain2_id == "10"
    assert cand.mid_chrom == "chrA"
    assert (cand.blocks[0].mid_start, cand.blocks[0].mid_end) == (1000, 1100)
    assert (cand.blocks[0].dest_start, cand.blocks[0].dest_end) == (0, 100)
    assert (cand.thick_dest_start, cand.thick_dest_end) == (10, 90)


def test_double_exon_gap_in_segment2_unmapped():
    # 段1：[0,200) 两块；段2：chrA 只有 [0,100) 与 [200,300) 对齐
    c1 = _chain("1", [(0, 100, 0, 100), (100, 200, 100, 200)],
                t_span=(0, 200), q_name="chrA", q_span=(0, 200))
    c2 = _chain("10", [(0, 100, 0, 100), (200, 300, 200, 300)],
                t_name="chrA", t_span=(0, 300),
                q_name="chrZ", q_span=(0, 300))
    eng = TranscriptLiftoverEngine([c1], "target", [c2], "target")
    # 外显子 2 [100,200) 在第二段位于 chrA 缺口 [100,200) => 无完整路径
    rec = _rec("chr1\t0\t200\tTX\t0\t+\t0\t200\t0\t2\t100,100,\t0,100,\n")
    row = eng.convert(rec)
    assert row.status == "unmapped"
    assert row.candidates == []


def test_double_two_chain2_paths_ambiguous():
    c1 = _chain("1", [(0, 100, 0, 100)], t_span=(0, 100),
                q_name="chrA", q_span=(0, 100))
    c2a = _chain("11", [(0, 100, 0, 100)], t_name="chrA", t_span=(0, 100),
                 q_name="chrW", q_span=(0, 100))
    c2b = _chain("12", [(0, 100, 50, 150)], t_name="chrA", t_span=(0, 100),
                 q_name="chrV", q_span=(50, 150))
    eng = TranscriptLiftoverEngine([c1], "target", [c2a, c2b], "target")
    rec = _rec("chr1\t0\t100\tTX\t0\t+\t0\t100\t0\t1\t100,\t0,\n")
    row = eng.convert(rec)
    assert row.status == "ambiguous"
    assert len(row.candidates) == 2
    assert {c.dest_chrom for c in row.candidates} == {"chrW", "chrV"}


def test_double_relative_orientation_single_exon():
    # 段1：chr1 正链 [0,100) -> chrA 负链，q 存储帧 [0,100)、qSize=500，
    #       归一化后中间正向帧为 [400,500)（方向反转一次）。
    c1 = _chain("1", [(0, 100, 0, 100)], t_span=(0, 100),
                q_name="chrA", q_size=500, q_strand="-", q_span=(0, 100))
    # 段2：chrA 负链（源）存储帧 [0,100) -> chrW 正链 [700,800)。
    #       源负、目标正 => 段2 也反转；两次反转抵消 => 总同向。
    c2 = Chain(
        chain_id="11",
        t_name="chrA", t_size=500, t_strand="-", t_start=0, t_end=100,
        q_name="chrW", q_size=800, q_strand="+", q_start=700, q_end=800,
        score=100, blocks=[Block(0, 100, 700, 800)], header_line=1,
    )
    eng = TranscriptLiftoverEngine([c1], "target", [c2], "target")
    rec = _rec("chr1\t0\t100\tTX\t0\t+\t0\t100\t0\t1\t100,\t0,\n")
    row = eng.convert(rec)
    cand = row.candidates[0]
    # 段1：chr1[0,100) -> chrA 正向帧 [400,500)（反向）
    assert (cand.blocks[0].mid_start, cand.blocks[0].mid_end) == (400, 500)
    # 段2 归一化后输出 chrW 正向帧 [700,800)
    assert (cand.blocks[0].dest_start, cand.blocks[0].dest_end) == (700, 800)
    # 两次反向相互抵消 => 总同向，源 '+' 输出 '+'
    assert cand.forward is True
    assert cand.output_strand == "+"
