"""BED12 严格解析器测试：12 列、块、外显子边界、thick 范围等。"""

from app.bed12 import parse_bed12


VALID = (
    "chr1\t56\t200\tTX1\t0\t+\t100\t190\t0\t3\t44,76,21,\t0,44,123,\n"
)


def test_valid_bed12_parses():
    r = parse_bed12(VALID)
    assert r.errors == []
    rec = r.records[0]
    assert rec.block_count == 3
    assert rec.block_sizes == (44, 76, 21)
    assert rec.block_starts == (0, 44, 123)
    assert rec.exons() == [(56, 100), (100, 176), (179, 200)]
    assert rec.coding is True


def test_must_be_exactly_12_columns():
    r = parse_bed12("chr1\t0\t100\n")  # 3 列的 BED3 行不接受
    assert len(r.errors) == 1
    assert "恰好 12 列" in r.errors[0]
    assert "3 列" in r.errors[0]

    r2 = parse_bed12("chr1\t0\t100\tn\t0\t+\t0\t100\t0\t1\t100,\t0,\textra\n")
    assert any("恰好 12 列" in e and "13 列" in e for e in r2.errors)


def test_block_count_mismatch_collected_with_line():
    text = (
        "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t50,50,\t0,50,\n"    # 行1 合法
        "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t50,50,\t0,,\n"     # 行2 空项
    )
    r = parse_bed12(text)
    assert len(r.records) == 1
    assert any("第 2 行" in e and "blockStarts" in e for e in r.errors)

    # 个数不匹配（去掉一个 start 项）
    text2 = (
        "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t50,50,\t0,\n"
    )
    r2 = parse_bed12(text2)
    assert any("blockStarts 有 1 项" in e and "blockCount=2" in e
               for e in r2.errors)


def test_block_count_must_be_positive():
    line = "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t0\t100,\t0,\n"
    r = parse_bed12(line)
    assert any("blockCount 必须是正整数" in e for e in r.errors)


def test_zero_block_size_rejected():
    line = "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t100,0,\t0,100,\n"
    r = parse_bed12(line)
    assert any("块长为 0" in e for e in r.errors)


def test_overlapping_blocks_rejected():
    line = "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t60,60,\t0,50,\n"
    r = parse_bed12(line)
    assert any("重叠" in e for e in r.errors)


def test_blocks_must_be_strictly_increasing():
    # 两个起点相同（既触发严格递增，也触发块越界/重叠）
    line = "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t50,50,\t0,0,\n"
    r = parse_bed12(line)
    assert any("严格递增" in e for e in r.errors)


def test_first_block_must_start_at_zero_last_must_reach_end():
    bad_first = "chr1\t10\t110\ta\t0\t+\t10\t110\t0\t1\t100,\t10,\n"
    r = parse_bed12(bad_first)
    assert any("首个外显子" in e for e in r.errors)

    bad_last = "chr1\t0\t40\ta\t0\t+\t0\t40\t0\t2\t20,10,\t0,20,\n"
    r2 = parse_bed12(bad_last)
    assert any("末个外显子" in e for e in r2.errors)

    # 外显子相接（intron 长度 0）是合法的
    touching = "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t50,50,\t0,50,\n"
    assert parse_bed12(touching).errors == []


def test_thick_range_must_be_inside_transcript():
    line = "chr1\t10\t110\ta\t0\t+\t5\t100\t0\t1\t100,\t0,\n"
    r = parse_bed12(line)
    assert any("thickStart/thickEnd 必须落在" in e for e in r.errors)


def test_thick_start_after_end_rejected():
    line = "chr1\t0\t100\ta\t0\t+\t90\t10\t0\t1\t100,\t0,\n"
    r = parse_bed12(line)
    assert any("thickStart (90) 不能大于 thickEnd (10)" in e for e in r.errors)


def test_coding_thick_endpoints_must_be_within_exons():
    # 外显子 [0,10) [90,100)；thickStart=50 落在内含子里
    line = "chr1\t0\t100\ta\t0\t+\t50\t95\t0\t2\t10,10,\t0,90,\n"
    r = parse_bed12(line)
    assert any("thickStart=50" in e and "不在任何外显子内" in e for e in r.errors)

    # thick 端点正好落在外显子边界（10）是允许的
    ok = "chr1\t0\t100\ta\t0\t+\t10\t90\t0\t2\t10,10,\t0,90,\n"
    assert parse_bed12(ok).errors == []


def test_noncoding_zero_length_thick_preserved():
    # thickStart=thickEnd（非编码）不要求端点在外显子内，但必须在 [start,end]
    line = "chr1\t56\t176\ta\t0\t.\t56\t56\t0\t1\t120,\t0,\n"
    r = parse_bed12(line)
    assert r.errors == []
    assert r.records[0].coding is False

    line2 = "chr1\t56\t176\ta\t0\t.\t100\t100\t0\t1\t120,\t0,\n"
    assert parse_bed12(line2).errors == []


def test_score_and_strand_validated():
    assert parse_bed12(
        "chr1\t0\t100\ta\t1001\t+\t0\t100\t0\t1\t100,\t0,\n"
    ).errors
    assert parse_bed12(
        "chr1\t0\t100\ta\t0\t?\t0\t100\t0\t1\t100,\t0,\n"
    ).errors
    # '.' strand 合法
    assert parse_bed12(
        "chr1\t0\t100\ta\t0\t.\t0\t100\t0\t1\t100,\t0,\n"
    ).errors == []


def test_item_rgb_validation():
    assert parse_bed12(
        "chr1\t0\t100\ta\t0\t+\t0\t100\t256,0,0\t1\t100,\t0,\n"
    ).errors
    assert parse_bed12(
        "chr1\t0\t100\ta\t0\t+\t0\t100\tred\t1\t100,\t0,\n"
    ).errors
    assert parse_bed12(
        "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t1\t100,\t0,\n"
    ).errors == []
    assert parse_bed12(
        "chr1\t0\t100\ta\t0\t+\t0\t100\t12,34,56\t1\t100,\t0,\n"
    ).errors == []


def test_comments_tracks_blank_skipped_and_empty_is_error():
    text = "# c\n\ntrack x=y\nbrowser foo\n"
    r = parse_bed12(text)
    assert any("为空" in e for e in r.errors)


def test_multiple_errors_all_collected():
    text = (
        "chr1\t0\t10\n"
        "chr1\t0\t100\ta\t0\t+\t0\t100\t0\t2\t60,60,\t0,50,\n"
    )
    r = parse_bed12(text)
    # 两行都坏，错误都保留
    assert any("12 列" in e for e in r.errors)
    assert any("重叠" in e for e in r.errors)
    assert r.records == []
