#!/usr/bin/env python3
"""
System-level (5G-LENA) experiment matrix for the UDN handover study.

Usage (from ns-3.48-nr):  python3 scratch/nr-udn-experiments.py [--sets E1,E2,E3,E4] [--procs 6]
                                                                [--seeds 5] [--ogThr -5.5]
Outputs: results/sys/raw/<tag>_r<seed>-{summary,handovers}.csv and results/sys/<set>.csv
"""

import argparse
import csv
import glob
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "build", "scratch", "ns3.48-nr-udn-handover-optimized")
# v3 adds the radio-level gap, the MAC gap, the A3 + return-block policy and the RLC options
BIN3 = os.path.join(ROOT, "build", "scratch", "ns3.48-nr-udn-handover-v3-optimized")
V3_SETS = ("E7", "E8")
# per-set run length (s); E7 is the higher-power scheme comparison of the revision
SIM_TIME = {"E7": 60.0}
MIX = "3,30,60,120"


def matrix(og_thr, og_cal):
    m = {"E1": [], "E2": [], "E3": [], "E4": [], "E6": [], "E7": [], "E8": []}
    for isd in (100, 75, 50, 35, 25):
        for tm in (40, 200):
            m["E1"].append(dict(isd=isd, tm=tm, speeds="30", policy="a3", ttt=100))
    for v in (3, 30, 60, 120):
        m["E2"].append(dict(isd=50, tm=40, speeds=str(v), policy="a3", ttt=100))
    for ttt in (40, 100, 256, 512):
        m["E3"].append(dict(isd=50, tm=40, speeds="60", policy="a3", ttt=ttt))
    for isd in (25, 50, 100):
        for ttt in (40, 100, 256):
            m["E4"].append(dict(isd=isd, tm=40, speeds=MIX, policy="a3", ttt=ttt))
        for thr in sorted({og_thr, og_cal}):
            m["E4"].append(dict(isd=isd, tm=40, speeds=MIX, policy="og", guardThr=thr))
        m["E4"].append(dict(isd=isd, tm=200, speeds=MIX, policy="a3", ttt=100))
        m["E4"].append(dict(isd=isd, tm=200, speeds=MIX, policy="og", guardThr=og_cal))
    # E6: OG-A3-lite (longer guard averaging window and larger guard offset)
    for isd in (25, 50, 100):
        m["E6"].append(dict(isd=isd, tm=40, speeds=MIX, policy="og", guardThr=og_thr, guardHys=2.0, guardTavg=0.05))
    # E7: scheme comparison with 10 seeds x 60 s, plus A3 with the OG-A3 return block
    lite = dict(guardHys=2.0, guardTavg=0.05)
    for isd in (25, 50, 100):
        for ttt in (40, 100, 256):
            m["E7"].append(dict(isd=isd, tm=40, speeds=MIX, policy="a3", ttt=ttt))
        for ttt in (40, 100):
            m["E7"].append(dict(isd=isd, tm=40, speeds=MIX, policy="a3rb", ttt=ttt))
        m["E7"].append(dict(isd=isd, tm=40, speeds=MIX, policy="og", guardThr=og_thr))
        m["E7"].append(dict(isd=isd, tm=40, speeds=MIX, policy="og", guardThr=og_thr, **lite))
        m["E7"].append(dict(isd=isd, tm=200, speeds=MIX, policy="a3", ttt=100))
        m["E7"].append(dict(isd=isd, tm=200, speeds=MIX, policy="og", guardThr=og_thr))
    # E8: RLC-UM t-Reordering of 10 ms instead of the 100 ms default (source of the
    # trigger-independent component of the pre-HO gap)
    for isd in (25, 50, 100):
        for ttt in (40, 100):
            m["E8"].append(dict(isd=isd, tm=40, speeds=MIX, policy="a3", ttt=ttt, rlcReorderMs=10))
        m["E8"].append(dict(isd=isd, tm=40, speeds=MIX, policy="og", guardThr=og_thr, rlcReorderMs=10))
        m["E8"].append(dict(isd=isd, tm=40, speeds=MIX, policy="og", guardThr=og_thr, rlcReorderMs=10, **lite))
    return m


def tag_of(p):
    if p["policy"] == "og":
        pol = "og%+.1f" % p["guardThr"]
    elif p["policy"] == "a3rb":
        pol = "a3rbt%d" % p["ttt"]
    else:
        pol = "a3t%d" % p["ttt"]
    if p.get("guardHys", 1.0) != 1.0 or p.get("guardTavg", 0.02) != 0.02:
        pol += "-lite"
    if "rlcReorderMs" in p:
        pol += "-ro%d" % p["rlcReorderMs"]
    sp = "mix" if p["speeds"] == MIX else "v" + p["speeds"]
    return "isd%d_tm%d_%s_%s" % (p["isd"], p["tm"], sp, pol)


def run_one(setname, p, seed, out, sim_time):
    prefix = os.path.join(out, "raw", "%s_%s_r%d" % (setname, tag_of(p), seed))
    if os.path.exists(prefix + "-summary.csv"):
        return prefix  # resume support
    # per-run lock so that several runner instances can share the matrix safely;
    # a stale lock (owner process gone, e.g. after a reboot) is taken over
    lock = prefix + ".lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        try:
            pid = int(open(lock).read().strip() or 0)
            os.kill(pid, 0)
            return prefix  # another live runner owns this run
        except (ValueError, ProcessLookupError, PermissionError):
            os.remove(lock)
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    sim_time = SIM_TIME.get(setname, sim_time)
    args = [BIN3 if setname in V3_SETS else BIN, "--gridN=4", "--simTime=%g" % sim_time, "--warm=2", "--run=%d" % seed,
            "--prefix=" + prefix] + ["--%s=%s" % (k, v) for k, v in p.items()]
    try:
        r = subprocess.run(args, capture_output=True, text=True)
    finally:
        if os.path.exists(lock):
            os.remove(lock)
    if r.returncode != 0:
        raise RuntimeError(" ".join(args) + "\n" + r.stderr[-1500:])
    return prefix


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", default="E1,E2,E3,E4")
    ap.add_argument("--procs", type=int, default=6)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--simTime", type=float, default=20.0)
    ap.add_argument("--ogThr", type=float, default=-5.5, help="a-priori guard threshold (TL-MC design)")
    ap.add_argument("--ogCal", type=float, default=-5.5, help="guard threshold from the calibration runs")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "sys"))
    ap.add_argument("--skipFailed", nargs="*", default=[],
                    help="runner logs whose FAILED runs are skipped (aborts are deterministic; they are "
                         "handled by nr-udn-rerun.py)")
    a = ap.parse_args()
    os.makedirs(os.path.join(a.out, "raw"), exist_ok=True)
    m = matrix(a.ogThr, a.ogCal)
    import re
    failed = set()
    for lg in a.skipFailed:
        failed |= set(re.findall(r"FAILED (\w+) (\S+) r(\d+):", open(lg).read()))
    jobs = [(s, p, seed) for s in a.sets.split(",") for p in m[s] for seed in range(1, a.seeds + 1)
            if (s, tag_of(p), str(seed)) not in failed]
    print(len(jobs), "runs", flush=True)
    done = 0
    with ThreadPoolExecutor(a.procs) as ex:
        futs = {ex.submit(run_one, s, p, seed, a.out, a.simTime): (s, p, seed) for s, p, seed in jobs}
        for f in as_completed(futs):
            done += 1
            s, p, seed = futs[f]
            try:
                f.result()
                print("[%d/%d] %s %s r%d ok" % (done, len(jobs), s, tag_of(p), seed), flush=True)
            except Exception as e:
                print("[%d/%d] FAILED %s %s r%d: %s" % (done, len(jobs), s, tag_of(p), seed, e), flush=True)
    for s in a.sets.split(","):
        rows = []
        for f in sorted(glob.glob(os.path.join(a.out, "raw", s + "_*-summary.csv"))):
            r = next(csv.DictReader(open(f)))
            r["set"] = s
            r["tag"] = os.path.basename(f).rsplit("_r", 1)[0][len(s) + 1:]
            r["handoverFile"] = os.path.relpath(f.replace("-summary", "-handovers"), a.out)
            rows.append(r)
        if rows:
            with open(os.path.join(a.out, s + ".csv"), "w", newline="") as fo:
                w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
            print("wrote", s, len(rows))


if __name__ == "__main__":
    main()
