"""Download and verify all inputs. Portable (Linux/macOS), resumable, md5-checked.

  python src/download.py kather   # Zenodo 2530835, four CRC_DX zips (21.7 GB) -> data/raw/CRC_DX_*/
  python src/download.py ref      # Zenodo 22938170, Vardar & Vardar deposit (0.6 GB) -> ref/main, ref/emb_parts
  MSI_ROOT=/mnt/slides/msimut/stad python src/download.py stad   # Zenodo 2530835 STAD zips (25.4 GB) -> $MSI_ROOT/data/raw/
  python src/download.py nct      # Zenodo 1214456, NCT-CRC-HE-100K (Macenko) + CRC-VAL-HE-7K (12.5 GB) -> data/nct/

Zenodo throttles single connections hard (we saw <1 MB/s), so each file is fetched as 8 parallel HTTP range
requests (segmented_download.py) with curl, resuming from whatever is already on disk.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "segmented_download.py")

KATHER = {  # name: (bytes, md5) from the Zenodo 2530835 API
    "CRC_DX_TRAIN_MSIMUT": (5256478484, "d98f4ee32f5923a09779d9aabd773e37"),
    "CRC_DX_TRAIN_MSS": (5259435967, "0dbb048590925920c383e423fa282f77"),
    "CRC_DX_TEST_MSIMUT": (3251312941, "006962fd57ffa71db43f3d768a63d9ec"),
    "CRC_DX_TEST_MSS": (7890181037, "3407a53dddec68e67f49b863a3520e88"),
}
STAD = {
    "STAD_TRAIN_MSIMUT": (5841008485, "c785506ccdf49d4d08b65c5e7cdd28ae"),
    "STAD_TRAIN_MSS": (5848305951, "273295ed37ebccc3a574ab9106145e32"),
    "STAD_TEST_MSIMUT": (3226442356, "77941b8914a2fac469ce4843c75ee19a"),
    "STAD_TEST_MSS": (10484110033, "8f324496c6feffe5353f493ae5e8c9ff"),
}
REF_RECORD = "22938170"
NCT = {  # Kather 2018 9-class tissue tiles (same group, 224 px, 0.5 um/px, Macenko); not TCGA, so no patient overlap
    "NCT-CRC-HE-100K": (11690284003, "6fd702d11df6292bc054397ae038a464"),
    "CRC-VAL-HE-7K": (800276929, "2fd1651b4f94ebd818ebf90ad2b6ce06"),
}


def fetch(url, out, size, md5, parts=8):
    for attempt in range(5):
        r = subprocess.run([sys.executable, SEG, url, out, str(size), md5, "--parts", str(parts)])
        if r.returncode == 0:
            return
        # md5 mismatch after a full-size file: start over
        if os.path.exists(out) and os.path.getsize(out) >= size:
            os.remove(out)
        time.sleep(30)
    raise SystemExit(f"failed to download {out}")


def unzip(zf, dest):
    with zipfile.ZipFile(zf) as z:
        for m in z.infolist():
            if m.filename.startswith("__MACOSX") or os.path.basename(m.filename).startswith("._"):
                continue
            z.extract(m, dest)


def kather(files=None, root=None):
    files = files or KATHER
    raw = os.path.join(root or ROOT, "data/raw")
    os.makedirs(raw, exist_ok=True)

    def one(item):
        name, (size, md5) = item
        dest = os.path.join(raw, name)
        if os.path.exists(os.path.join(dest, ".done")):
            print("already done", name, flush=True)
            return
        zf = os.path.join(raw, name + ".zip")
        fetch(f"https://zenodo.org/records/2530835/files/{name}.zip?download=1", zf, size, md5)
        unzip(zf, dest)
        os.remove(zf)
        open(os.path.join(dest, ".done"), "w").close()
        print("unzipped", name, flush=True)

    with ThreadPoolExecutor(4) as ex:
        list(ex.map(one, files.items()))


def nct():
    d = os.path.join(ROOT, "data/nct")
    os.makedirs(d, exist_ok=True)
    for name, (size, md5) in NCT.items():
        if os.path.exists(os.path.join(d, name, ".done")):
            continue
        zf = os.path.join(d, name + ".zip")
        fetch(f"https://zenodo.org/records/1214456/files/{name}.zip?download=1", zf, size, md5)
        unzip(zf, os.path.join(d, name))
        os.remove(zf)
        open(os.path.join(d, name, ".done"), "w").close()
        print("unzipped", name, flush=True)


def ref():
    d = os.path.join(ROOT, "ref")
    os.makedirs(d, exist_ok=True)
    meta = json.loads(subprocess.run(["curl", "-sfL", "--retry", "10", f"https://zenodo.org/api/records/{REF_RECORD}"],
                                     capture_output=True, check=True).stdout)
    for f in meta["files"]:
        out = os.path.join(d, f["key"])
        md5 = f["checksum"].split(":")[1]
        fetch(f["links"]["self"], out, f["size"], md5, parts=2)
        if f["key"].endswith(".zip"):
            sub = "main" if "EMBEDDINGS" not in f["key"] else os.path.join("emb_parts", f["key"][:-4])
            if not os.path.exists(os.path.join(d, sub)):
                unzip(out, os.path.join(d, sub))
        print("ok", f["key"], flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["kather", "ref", "nct", "stad", "all"])
    a = ap.parse_args()
    if a.what in ("ref", "all"):
        ref()
    if a.what in ("kather", "all"):
        kather()
    if a.what == "stad":
        kather(STAD, os.environ["MSI_ROOT"])
    if a.what in ("nct", "all"):
        nct()
    shutil.rmtree(os.path.join(ROOT, "data/raw/__MACOSX"), ignore_errors=True)
