#!/usr/bin/env python3
"""
Re-run the 5G-LENA runs of an experiment set that were aborted by a simulator assertion.

5G-LENA v5.0 occasionally aborts with an RRC/PHY assertion after many handovers (e.g. "method
unexpected in state CONNECTION_RECONFIGURATION"). The simulation is deterministic for a given seed,
so an aborted run is repeated with the same seed and a duration 0.1 s shorter than the abort time
when the abort happened late enough (at least MIN_TIME s); otherwise the seed is replaced by the
next unused seed of that configuration, run for the full duration. Replacement runs that abort
again are handled in the same way.

Usage (from ns-3.48-nr): python3 scratch/nr-udn-rerun.py --log results/sys_e7.log [more logs] --set E7 [--procs 7]
"""

import argparse
import glob
import math
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MIN_TIME = 20.0


def parse_failures(log):
    txt = open(log).read()
    out = []
    for m in re.finditer(r"FAILED (\w+) (\S+) r(\d+): (\S+ [^\n]*)\n(.*?)(?=\n\[\d+/\d+\]|\Z)", txt, re.S):
        s, tag, seed, cmd, err = m.groups()
        t = re.findall(r"\+(\d+\.\d+)s", err)
        out.append((s, tag, int(seed), cmd.split(" "), float(t[-1]) if t else float("nan")))
    return out


def arg(args, key):
    for a in args:
        if a.startswith("--%s=" % key):
            return a.split("=", 1)[1]
    return None


def with_arg(args, key, val):
    return [a for a in args if not a.startswith("--%s=" % key)] + ["--%s=%s" % (key, val)]


def run(args):
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode == 0:
        return None
    t = re.findall(r"\+(\d+\.\d+)s", r.stderr)
    return float(t[-1]) if t else float("nan")


def handle(job, used, lock_free_seed):
    args, abort_t = job
    while True:
        prefix = arg(args, "prefix")
        if os.path.exists(prefix + "-summary.csv"):
            return "exists " + os.path.basename(prefix)
        full = float(arg(args, "simTime"))
        if math.isfinite(abort_t) and abort_t - 0.1 >= MIN_TIME:
            a = with_arg(args, "simTime", "%.1f" % (math.floor((abort_t - 0.1) * 10) / 10))
            if run(a) is None:
                return "truncated %s to %s s" % (os.path.basename(prefix), arg(a, "simTime"))
            abort_t = float("nan")  # should not happen (deterministic); fall back to a new seed
        seed = lock_free_seed(prefix)
        base = prefix.rsplit("_r", 1)[0]
        args = with_arg(with_arg(with_arg(args, "run", seed), "prefix", "%s_r%d" % (base, seed)), "simTime", "%g" % full)
        abort_t = run(args)
        if abort_t is None:
            return "replaced %s by seed %d" % (os.path.basename(prefix), seed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True, nargs="+")
    ap.add_argument("--set", required=True)
    ap.add_argument("--procs", type=int, default=7)
    a = ap.parse_args()
    fails = [f for lg in a.log for f in parse_failures(lg) if f[0] == a.set]
    import threading
    lk = threading.Lock()
    taken = set()

    def free_seed(prefix):
        base = prefix.rsplit("_r", 1)[0]
        with lk:
            s = 11
            while glob.glob("%s_r%d-summary.csv" % (base, s)) or (base, s) in taken:
                s += 1
            taken.add((base, s))
            return s

    jobs = [(cmd, t) for _, _, _, cmd, t in fails]
    print(len(jobs), "aborted runs", flush=True)
    with ThreadPoolExecutor(a.procs) as ex:
        for msg in ex.map(lambda j: handle(j, None, free_seed), jobs):
            print(msg, flush=True)


if __name__ == "__main__":
    main()
