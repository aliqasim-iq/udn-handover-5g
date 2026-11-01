#!/usr/bin/env python3
"""
Aggregate the 5G-LENA runs (results/sys/E*.csv + per-HO files) into sys_stats.json and draw the
system-level figures (fig7..fig11) used in Section 7 of the manuscript.

Usage: <python with matplotlib> scratch/nr-udn-analysis.py [--res results/sys] [--figs DIR]
"""

import argparse
import csv
import json
import math
import os
import statistics as st
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#7a5cc2"]
MK = ["o", "s", "^", "D", "v", "P", "X"]
INK2, GRID = "#555555", "#e3e3e3"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.edgecolor": INK2,
    "axes.grid": True, "grid.color": GRID, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.8, "lines.markersize": 5.5, "legend.frameon": False, "savefig.dpi": 300,
    "savefig.bbox": "tight"})
T95 = {1: 12.71, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57, 6: 2.45, 7: 2.36, 8: 2.31, 9: 2.26}
METRICS = ["hoRate", "texec_ms", "tpre_ms", "tpost_ms", "tint_ms", "lostPerHo", "plr_pct",
           "hoPlr_pct", "pp1_pct", "pp3_pct", "guard_pct", "hoFail", "rlf", "outage_pct", "thr_Mbps",
           "delay_ms", "numHO"]
# v3 runs (E7, E8) additionally report the radio-level gap and the MAC gap
METRICS_V3 = ["tpreRadio_ms", "pRadioGap", "macGap_ms"]


def mci(v):
    v = [x for x in v if not math.isnan(x)]
    if not v:
        return [float("nan"), float("nan")]
    if len(v) == 1:
        return [v[0], 0.0]
    return [st.mean(v), T95.get(len(v) - 1, 1.96) * st.stdev(v) / math.sqrt(len(v))]


def load(res, s):
    """Read every finished run of set s directly from raw/ (works while the sweep is running)."""
    import glob
    rows = []
    for f in sorted(glob.glob(os.path.join(res, "raw", s + "_*-summary.csv"))):
        r = next(csv.DictReader(open(f)))
        r["tag"] = os.path.basename(f).rsplit("_r", 1)[0][len(s) + 1:]
        r["handoverFile"] = os.path.relpath(f.replace("-summary", "-handovers"), res)
        rows.append(r)
    return rows


def group(rows):
    g = defaultdict(list)
    for r in rows:
        g[r["tag"]].append(r)
    return g


def summarize(res, rows):
    out = {}
    for tag, rs in group(rows).items():
        d = {m: mci([float(r[m]) for r in rs]) for m in METRICS}
        d.update({m: mci([float(r[m]) for r in rs]) for m in METRICS_V3 if m in rs[0]})
        # RLF per UE-hour from each run's own duration (numHO / (UEs * hoRate) = simTime - warm),
        # so that runs truncated before a simulator abort are weighted correctly
        d["rlfUeH"] = mci([float(r["rlf"]) * 3600 * float(r["hoRate"]) / float(r["numHO"])
                           for r in rs if float(r["numHO"]) > 0])
        d["seeds"] = len(rs)
        d.update({k: rs[0][k] for k in ("isd", "tm", "ttt", "policy", "speeds")})
        # per-HO pooled statistics: decision delay and Tpre distribution
        tdf, tdl, prep, pre = [], [], [], []
        # system-level outage margin kappa = v K of HOs with a radio-level gap (K = to - th); only
        # K < T_D is observable, so kappa is reported per HO of the fast classes (60, 120 km/h)
        kap, nfast, nradio = [], 0, 0
        for r in rs:
            for h in csv.DictReader(open(os.path.join(res, r["handoverFile"]))):
                if h["tPre_ms"] not in ("nan", ""):
                    pre.append(float(h["tPre_ms"]))
                if "k_first_ms" in h:
                    v = float(h["speed_kmh"]) / 3.6
                    tr = float(h["tPreRadio_ms"])
                    nradio += tr > 0
                    if v >= 60 / 3.6 - 1e-9:
                        nfast += 1
                        k = float(h["k_first_ms"])
                        if tr > 0 and math.isfinite(k):
                            kap.append(v * k / 1e3)
                if h["trigger"] != "a3":
                    continue
                for key, lst in (("td_first_ms", tdf), ("td_last_ms", tdl), ("tPrep_ms", prep)):
                    x = float(h[key])
                    if math.isfinite(x) and 0 <= x < 5000:
                        lst.append(x)
        d["td_first_med"] = st.median(tdf) if tdf else float("nan")
        d["td_last_med"] = st.median(tdl) if tdl else float("nan")
        d["tprep_med"] = st.median(prep) if prep else float("nan")
        d["n_td"] = len(tdf)
        d["pgap"] = sum(1 for x in pre if x > 0) / len(pre) if pre else float("nan")
        d["tpre_p90"] = sorted(pre)[int(0.9 * (len(pre) - 1))] if pre else float("nan")
        if nfast:
            d["kappa_fast"] = sorted(round(x, 3) for x in kap)
            d["n_fast"] = nfast
            d["kneg_fast"] = sum(1 for x in kap if x < 0) / nfast
        out[tag] = d
    return out


def eb(ax, xs, ys, cis, i, lab, ls="-"):
    ax.errorbar(xs, ys, yerr=cis, color=C[i], marker=MK[i], capsize=3, label=lab, ls=ls,
                markeredgecolor="white", markeredgewidth=1)


def pick(S, **kw):
    return {t: d for t, d in S.items() if all(str(d[k]) == str(v) for k, v in kw.items())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default=os.path.join(ROOT, "results", "sys"))
    ap.add_argument("--figs", default=os.path.join(ROOT, "results", "sys"))
    a = ap.parse_args()
    S = {s: summarize(a.res, load(a.res, s)) for s in ("E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8")}
    # buffer-flush floor: Tpre of A3 HOs (T_m = 40 ms) that start with serving PDCCH SINR > 0 dB
    good, delays = [], []
    for s in ("E1", "E4"):
        for r in load(a.res, s):
            if r["policy"] != "a3" or float(r["tm"]) != 40:
                continue
            delays.append(float(r["delay_ms"]))
            for h in csv.DictReader(open(os.path.join(a.res, r["handoverFile"]))):
                try:
                    tp, ss = float(h["tPre_ms"]), float(h["sinrStart_db"])
                except ValueError:
                    continue
                if math.isfinite(tp) and math.isfinite(ss) and ss > 0:
                    good.append(tp)
    S["floor"] = {"goodMed": st.median(good) if good else float("nan"), "n": len(good),
                  "delay": st.mean(delays) if delays else float("nan")}
    S["calib"] = json.load(open(os.path.join(a.res, "calib.json"))) if os.path.exists(
        os.path.join(a.res, "calib.json")) else {}
    json.dump(S, open(os.path.join(a.res, "sys_stats.json"), "w"), indent=1)

    # ---------------- Figure 7: density (E1)
    if S["E1"]:
        fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.1))
        isds = [100, 75, 50, 35, 25]
        for j, tm in enumerate((40, 200)):
            sel = [pick(S["E1"], isd=i, tm=tm) for i in isds]
            sel = [list(x.values())[0] for x in sel if x]
            xs = [float(d["isd"]) for d in sel]
            eb(axs[0], xs, [d["tpre_ms"][0] for d in sel], [d["tpre_ms"][1] for d in sel], j,
               r"$T_{pre}$, $T_m$=%d ms" % tm)
            eb(axs[1], xs, [d["plr_pct"][0] for d in sel], [d["plr_pct"][1] for d in sel], j,
               r"$T_m$=%d ms" % tm)
            eb(axs[2], xs, [d["hoRate"][0] for d in sel], [d["hoRate"][1] for d in sel], j,
               r"$T_m$=%d ms" % tm)
            if tm == 40:
                eb(axs[0], xs, [d["texec_ms"][0] for d in sel], [d["texec_ms"][1] for d in sel], 2,
                   r"$T_{exec}$, $T_m$=40 ms", ls="--")
        geo = [4 * (30 / 3.6) / (math.pi * i) for i in isds]
        axs[2].plot(isds, geo, color=INK2, ls=":", marker="x", label=r"geometric $4v\sqrt{\lambda}/\pi$")
        axs[0].set(xlabel="ISD (m)", ylabel="ms", title="(a) Pre-HO gap and execution")
        axs[1].set(xlabel="ISD (m)", ylabel="PLR (%)", title="(b) Packet loss ratio")
        axs[2].set(xlabel="ISD (m)", ylabel="HOs per UE per s", title="(c) HO rate")
        for ax in axs:
            ax.invert_xaxis()
            ax.legend()
        fig.savefig(os.path.join(a.figs, "fig7_sys_density.png"))
        plt.close(fig)

    # ---------------- Figure 8: speed (E2) and TTT (E3)
    if S["E2"] or S["E3"]:
        fig, axs = plt.subplots(1, 4, figsize=(13, 3.1))
        sel = sorted(S["E2"].values(), key=lambda d: float(d["speeds"]))
        if sel:
            xs = [float(d["speeds"]) for d in sel]
            eb(axs[0], xs, [d["tpre_ms"][0] for d in sel], [d["tpre_ms"][1] for d in sel], 0, r"$T_{pre}$")
            eb(axs[0], xs, [d["texec_ms"][0] for d in sel], [d["texec_ms"][1] for d in sel], 2, r"$T_{exec}$", ls="--")
            eb(axs[1], xs, [d["plr_pct"][0] for d in sel], [d["plr_pct"][1] for d in sel], 0, "PLR")
            eb(axs[1], xs, [d["pp1_pct"][0] for d in sel], [d["pp1_pct"][1] for d in sel], 1, "ping-pong (1 s)")
            axs[0].set(xlabel="UE speed (km/h)", ylabel="ms", title="(a) Gap vs speed (ISD 50 m)")
            axs[1].set(xlabel="UE speed (km/h)", ylabel="%", title="(b) PLR and ping-pong vs speed")
            axs[0].legend()
            axs[1].legend()
        sel = sorted(S["E3"].values(), key=lambda d: float(d["ttt"]))
        if sel:
            xs = [float(d["ttt"]) for d in sel]
            eb(axs[2], xs, [d["tpre_ms"][0] for d in sel], [d["tpre_ms"][1] for d in sel], 0, r"$T_{pre}$ (ms)")
            eb(axs[3], xs, [d["pp1_pct"][0] for d in sel], [d["pp1_pct"][1] for d in sel], 1, "ping-pong (1 s)")
            eb(axs[3], xs, [d["plr_pct"][0] for d in sel], [d["plr_pct"][1] for d in sel], 0, "PLR")
            axs[2].set(xlabel="TTT (ms)", ylabel="ms", title="(c) Gap vs TTT (ISD 50 m, 60 km/h)")
            axs[3].set(xlabel="TTT (ms)", ylabel="%", title="(d) PLR and ping-pong vs TTT")
            axs[2].legend()
            axs[3].legend()
        fig.savefig(os.path.join(a.figs, "fig8_sys_speed_ttt.png"))
        plt.close(fig)

    # ---------------- Figure 9: scheme comparison (E4) trade-off
    SC = "E7" if S["E7"] else "E4"
    if S[SC]:
        fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.2), sharey=False)
        for k, isd in enumerate((25, 50, 100)):
            ax = axs[k]
            a3 = sorted(pick(S[SC], isd=isd, tm=40, policy="a3").values(), key=lambda d: float(d["ttt"]))
            if a3:
                ax.errorbar([d["pp3_pct"][0] for d in a3], [d["tpre_ms"][0] for d in a3],
                            xerr=[d["pp3_pct"][1] for d in a3], yerr=[d["tpre_ms"][1] for d in a3],
                            color=C[0], marker="o", capsize=2, label="A3, TTT 40/100/256")
                for d in a3:
                    ax.annotate("%d" % float(d["ttt"]), (d["pp3_pct"][0], d["tpre_ms"][0]),
                                textcoords="offset points", xytext=(4, 4), fontsize=7, color=C[0])
            marks = [("isd%d_tm40_mix_og-5.5", "OG-A3", "D"), ("isd%d_tm40_mix_og-5.5-lite", "OG-A3-lite", "X"),
                     ("isd%d_tm40_mix_a3rbt40", "A3-40 + RB", "P"), ("isd%d_tm40_mix_a3rbt100", "A3-100 + RB", "v"),
                     ("isd%d_tm200_mix_a3t100", r"A3-100, $T_m$=200", "s"),
                     ("isd%d_tm200_mix_og-5.5", r"OG-A3, $T_m$=200", "^")]
            if SC == "E4":
                marks = [("isd%d_tm40_mix_og-5.5", r"OG-A3, $\gamma_g$=−5.5 dB", "D"),
                         ("isd%d_tm40_mix_og-5.5-lite", "OG-A3-lite", "X"),
                         ("isd%d_tm40_mix_og-2.5", r"OG-A3, $\gamma_g$=−2.5 dB", "P"),
                         ("isd%d_tm200_mix_a3t100", r"A3-100, $T_m$=200", "s"),
                         ("isd%d_tm200_mix_og-2.5", r"OG-A3 (−2.5), $T_m$=200", "^")]
            for j, (tag, lab, mk) in enumerate(marks):
                d = S[SC].get(tag % isd) or S["E6"].get(tag % isd)
                if d:
                    ax.errorbar([d["pp3_pct"][0]], [d["tpre_ms"][0]], xerr=[d["pp3_pct"][1]], yerr=[d["tpre_ms"][1]],
                                color=C[j + 1], marker=mk, capsize=2, ls="none", label=lab, ms=7)
            ax.set(xlabel="ping-pong within 3 s (%)", ylabel=r"mean $T_{pre}$ (ms)", title="ISD %d m" % isd)
        h, lab = axs[0].get_legend_handles_labels()
        fig.legend(h, lab, loc="upper center", ncol=4, fontsize=7.2, bbox_to_anchor=(0.5, -0.02))
        fig.savefig(os.path.join(a.figs, "fig9_sys_schemes.png"))
        plt.close(fig)

    # ---------------- Figure 10: system-level validation of the decision delay
    pts = []
    for s in ("E1", "E2", "E3", "E4", "E7"):
        for d in S[s].values():
            if d["policy"] != "a3" or not d["n_td"]:
                continue
            tm, ttt = float(d["tm"]), float(d["ttt"])
            pred = 2 * tm + ttt + (d["tprep_med"] if math.isfinite(d["tprep_med"]) else 0)
            pts.append((pred, d["td_first_med"], d["td_last_med"], s))
    if pts:
        fig, ax = plt.subplots(figsize=(4.4, 3.6))
        lo, hi = 0, max(max(p[0], p[1], p[2]) for p in pts) * 1.08
        ax.plot([lo, hi], [lo, hi], color=INK2, ls=":", lw=1, label="y = x")
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], color=C[1], marker="^", s=22, label="first A3 entry")
        ax.scatter([p[0] for p in pts], [p[2] for p in pts], color=C[0], marker="o", s=22, label="last A3 entry")
        ax.set(xlabel=r"predicted $E[T_D] = T_m/a + T_{TTT} + T_{prep}$ (ms)",
               ylabel=r"measured median $T_D$ (ms)", xlim=(lo, hi), ylim=(lo, hi))
        ax.legend()
        fig.savefig(os.path.join(a.figs, "fig10_sys_validation.png"))
        plt.close(fig)

    # ---------------- Figure 11: radio-level gap (E7) and system-level outage margin
    if S["E7"]:
        fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.4))
        pols = [("a3t40", "A3-40"), ("a3t100", "A3-100"), ("a3t256", "A3-256"), ("a3rbt40", "A3-40+RB"),
                ("a3rbt100", "A3-100+RB"), ("og-5.5", "OG-A3"), ("og-5.5-lite", "OG-A3-lite")]
        w = 0.26
        for k, isd in enumerate((25, 50, 100)):
            ds = [S["E7"].get("isd%d_tm40_mix_%s" % (isd, p)) for p, _ in pols]
            xs = [i + (k - 1) * w for i in range(len(pols))]
            axs[0].bar(xs, [d["tpreRadio_ms"][0] if d else 0 for d in ds], w, yerr=[d["tpreRadio_ms"][1] if d else 0 for d in ds],
                       color=C[k], label="ISD %d m" % isd, capsize=2)
        axs[0].set_xticks(range(len(pols)))
        axs[0].set_xticklabels([l for _, l in pols], rotation=30, ha="right")
        axs[0].set(ylabel=r"mean radio-level gap $T_{pre,radio}$ (ms)", title="(a) Radio-level pre-HO gap (E7)")
        axs[0].legend()
        kfile = os.path.expanduser("~/ns-3-dev/paper_v2/model/res_kappa.json")
        if os.path.exists(kfile):
            for j, r in enumerate([r for r in json.load(open(kfile)) if r["shadow"] == "full" and r["D"] in (25, 50, 100)]):
                kk = sorted(x for x in r["kappa"] if math.isfinite(x))
                axs[1].plot(kk, [(i + 1) / len(kk) for i in range(len(kk))], color=C[j], ls=":", lw=1.3,
                            label="TL-MC, ISD %d m" % r["D"])
        for j, isd in enumerate((25, 50, 100)):
            d = S["E7"].get("isd%d_tm40_mix_a3t100" % isd)
            if d and d.get("n_fast"):
                kk = d["kappa_fast"]
                axs[1].step(kk, [(i + 1) / d["n_fast"] for i in range(len(kk))], where="post", color=C[j],
                            label="5G-LENA, ISD %d m" % isd)
        axs[1].set(xlim=(-3, 3.5), ylim=(0, 0.7), xlabel=r"outage-margin distance $\kappa$ (m)", ylabel="CDF",
                   title=r"(b) Lower tail of $F_\kappa$ (A3-100; 60/120 km/h)")
        axs[1].legend(fontsize=6.8, ncol=2)
        fig.savefig(os.path.join(a.figs, "fig11_sys_radio_kappa.png"))
        plt.close(fig)
    print("wrote sys_stats.json and figures to", a.figs)


if __name__ == "__main__":
    main()
