"""Resume a partial download with parallel HTTP range requests, then verify md5.

  python src/segmented_download.py URL OUT SIZE MD5 [--parts 8]
"""
import argparse
import hashlib
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor


def fetch(url, out, start, end):
    """Download bytes [start, end] into out (resumable by its own size)."""
    want = end - start + 1
    for attempt in range(200):
        have = os.path.getsize(out) if os.path.exists(out) else 0
        if have == want:
            return
        if have > want:
            os.remove(out)
            continue
        subprocess.run(["curl", "-sfL", "--speed-limit", "20000", "--speed-time", "60", "-r",
                        f"{start + have}-{end}", "-o", "-", url], stdout=open(out, "ab"))
        time.sleep(5)
    raise RuntimeError(f"failed {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("url")
    ap.add_argument("out")
    ap.add_argument("size", type=int)
    ap.add_argument("md5")
    ap.add_argument("--parts", type=int, default=8)
    a = ap.parse_args()
    have = os.path.getsize(a.out) if os.path.exists(a.out) else 0
    rem = a.size - have
    step = -(-rem // a.parts)
    segs = [(i, have + i * step, min(a.size, have + (i + 1) * step) - 1) for i in range(a.parts) if have + i * step < a.size]
    if segs:  # an already-complete file has no segments left
        with ThreadPoolExecutor(len(segs)) as ex:
            list(ex.map(lambda s: fetch(a.url, f"{a.out}.seg{s[0]}", s[1], s[2]), segs))
    with open(a.out, "ab") as f:
        for i, _, _ in segs:
            with open(f"{a.out}.seg{i}", "rb") as g:
                while chunk := g.read(1 << 24):
                    f.write(chunk)
    h = hashlib.md5()
    with open(a.out, "rb") as f:
        while chunk := f.read(1 << 24):
            h.update(chunk)
    ok = h.hexdigest() == a.md5
    if ok:
        for i, _, _ in segs:
            os.remove(f"{a.out}.seg{i}")
    print("md5", "ok" if ok else "MISMATCH", a.out, flush=True)
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
