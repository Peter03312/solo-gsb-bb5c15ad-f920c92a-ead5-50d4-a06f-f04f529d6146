"""BED12 转录本整体迁移 HTTP 接口测试（单段 / 双段）。"""

import os

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


def _post12(bed_name, chain_name="demo.chain", side="target"):
    with open(os.path.join(EXAMPLES, bed_name), "rb") as bf, \
         open(os.path.join(EXAMPLES, chain_name), "rb") as cf:
        return client.post(
            "/api/liftover12",
            files={"bed": (bed_name, bf, "text/plain"),
                   "chain": (chain_name, cf, "text/plain")},
            data={"side": side},
        )


def _post12_double(bed_name):
    with open(os.path.join(EXAMPLES, bed_name), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo_seg1.chain"), "rb") as c1f, \
         open(os.path.join(EXAMPLES, "demo_seg2.chain"), "rb") as c2f:
        return client.post(
            "/api/liftover12/double",
            files={"bed": (bed_name, bf), "chain1": ("a.chain", c1f),
                   "chain2": ("b.chain", c2f)},
            data={"side1": "target", "side2": "target"},
        )


def test_bed12_single_preview_counts_and_kind():
    resp = _post12("example_bed12.bed")
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["kind"] == "single12"
    preview = data["preview"]
    assert preview["total"] == 4
    # 2 mapped（TX1/TX4）、1 ambiguous（TX2）、1 unmapped（TX3），无 partial
    assert preview["counts"] == {
        "mapped": 2, "partial": 0, "ambiguous": 1, "unmapped": 1,
    }
    # 没有 partial 状态出现在任何行
    assert all(r["status"] != "partial" for r in preview["rows"])


def test_bed12_download_is_12_columns_and_skips_unmapped():
    data = _post12("example_bed12.bed").json()
    token = data["download_token"]
    dl = client.get(f"/api/download12/{token}")
    assert dl.status_code == 200
    assert dl.headers["content-type"].startswith("text/x-bed")
    body = dl.text
    lines = body.splitlines()
    # TX1 一行 + TX2 两个完整候选各一行 + TX4 一行；TX3 不导出
    assert len(lines) == 4
    for ln in lines:
        assert len(ln.split("\t")) == 12
    assert "TX3" not in body
    # 歧义转录本 TX2 的两条路径都导出
    tx2 = [ln for ln in lines if "\tTX2\t" in ln]
    assert len(tx2) == 2
    assert {ln.split("\t")[0] for ln in tx2} == {"chrA", "chrB"}
    # 反向路径 strand 翻转（源 '-' -> chain 2 反向 -> '+'）
    chrB = next(ln for ln in tx2 if ln.startswith("chrB"))
    assert chrB.split("\t")[5] == "+"
    # 非编码 TX4：零长 thick 保留
    tx4 = next(ln for ln in lines if "\tTX4\t" in ln).split("\t")
    assert tx4[6] == tx4[7]
    assert tx4[8] == "255,0,0"


def test_bed12_download_url_in_response():
    data = _post12("example_bed12.bed").json()
    assert data["download_url"] == f"/api/download12/{data['download_token']}"
    assert data["download_name"].endswith(".bed")


def test_bed12_bad_file_collects_line_numbers_and_blocks_export():
    resp = _post12("example_bed12_bad.bed")
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    # 6 个坏行各有错误，且带物理行号
    for n in (8, 9, 10, 11, 12, 13):
        assert any(f"第 {n} 行" in e for e in errors), (n, errors)
    # 整次转换被拒绝：没有任何 token 被生成
    assert "download_token" not in resp.json()


def test_bed12_chain_error_still_blocked():
    # BED12 合法但 chain 损坏 => 400 且不产生下载
    bad_chain = "this is not a chain\n"
    with open(os.path.join(EXAMPLES, "example_bed12.bed"), "rb") as bf:
        resp = client.post(
            "/api/liftover12",
            files={"bed": ("b.bed", bf), "chain": ("c.chain", bad_chain)},
            data={"side": "target"},
        )
    assert resp.status_code == 400
    assert resp.json()["errors"]


def test_bed12_invalid_side():
    with open(os.path.join(EXAMPLES, "example_bed12.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo.chain"), "rb") as cf:
        resp = client.post(
            "/api/liftover12",
            files={"bed": ("b.bed", bf), "chain": ("c.chain", cf)},
            data={"side": "left"},
        )
    assert resp.status_code == 400


def test_bed12_empty_download_when_all_unmapped():
    # 所有转录本都无完整路径：仍 200 返回预览（逐行原因），下载为空 BED12
    bed = "chr9\t0\t100\tNOPE\t0\t+\t0\t100\t0\t1\t100,\t0,\n"
    with open(os.path.join(EXAMPLES, "demo.chain"), "rb") as cf:
        resp = client.post(
            "/api/liftover12",
            files={"bed": ("b.bed", bed), "chain": ("c.chain", cf)},
            data={"side": "target"},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["preview"]["counts"]["unmapped"] == 1
    dl = client.get(data["download_url"])
    assert dl.status_code == 200
    assert dl.content == b""      # 无伪 BED12


def test_bed12_download_token_404():
    assert client.get("/api/download12/not-a-token").status_code == 404


def test_bed12_double_preview_and_download():
    resp = _post12_double("example_bed12_double.bed")
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["kind"] == "double12"
    assert data["preview"]["counts"] == {
        "mapped": 2, "partial": 0, "ambiguous": 1, "unmapped": 1,
    }
    dl = client.get(data["download_url"])
    lines = dl.text.splitlines()
    assert len(lines) == 4          # DX1、DX2 各 1 行，DX3 2 行，DX4 不导出
    for ln in lines:
        assert len(ln.split("\t")) == 12
    names = [ln.split("\t")[3] for ln in lines]
    assert names == ["DX1", "DX2", "DX3", "DX3"]
    # 双段歧义两目标
    assert {ln.split("\t")[0] for ln in lines if "\tDX3\t" in ln} == {
        "chrW", "chrV"
    }
    # 文件名带 double
    assert "double" in data["download_name"]


def test_bed12_double_bad_bed_rejected():
    resp = client.post(
        "/api/liftover12/double",
        files={
            "bed": ("b.bed", "chr1\t0\t10\tx\n"),
            "chain1": ("a", open(os.path.join(EXAMPLES, "demo_seg1.chain"), "rb")),
            "chain2": ("b", open(os.path.join(EXAMPLES, "demo_seg2.chain"), "rb")),
        },
        data={"side1": "target", "side2": "target"},
    )
    assert resp.status_code == 400
    assert any("12 列" in e for e in resp.json()["errors"])


def test_bed12_double_bad_chain2_prefixed():
    with open(os.path.join(EXAMPLES, "example_bed12_double.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo_seg1.chain"), "rb") as c1f:
        resp = client.post(
            "/api/liftover12/double",
            files={"bed": ("b.bed", bf), "chain1": ("a.chain", c1f),
                   "chain2": ("b.chain", "bad chain content")},
            data={"side1": "target", "side2": "target"},
        )
    assert resp.status_code == 400
    assert any(e.startswith("chain 文件 2：") for e in resp.json()["errors"])


def test_bed3_endpoint_still_accepts_3_column_bed():
    # 兼容性：原 BED3 接口行为不变
    resp = client.post(
        "/api/liftover",
        files={
            "bed": ("b.bed", "chr1\t56\t176\n"),
            "chain": ("c.chain",
                      open(os.path.join(EXAMPLES, "demo.chain"), "rb")),
        },
        data={"side": "target"},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["kind"] == "single"
    assert data["preview"]["rows"][0]["status"] == "mapped"


def test_bed12_endpoint_rejects_bed3_line():
    # BED12 接口不接受 3 列 BED3 行（与 BED3 接口的宽松解析区分开）
    with open(os.path.join(EXAMPLES, "demo.chain"), "rb") as cf:
        resp = client.post(
            "/api/liftover12",
            files={"bed": ("b.bed", "chr1\t56\t176\n"),
                   "chain": ("c.chain", cf)},
            data={"side": "target"},
        )
    assert resp.status_code == 400
    assert any("12 列" in e for e in resp.json()["errors"])
