import json, numpy as np, matplotlib.pyplot as plt
from collections import defaultdict
from style import *
apply()
res = json.load(open('../model/res_scheme.json')) + json.load(open('../model/res_sg.json')) + json.load(open('../model/res_lite.json'))
A = defaultdict(list)
for r in res: A[(r['D'], r['tag'])].append(r)
def m(D, tag, k): return np.mean([x[k] for x in A[(D, tag)]])
def ci(D, tag, k):
    v = [x[k] for x in A[(D, tag)]]
    return 2.262 * np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.9), sharey=True)
pts = [('mse', 'MSE (oracle)', C[3], 'D'), ('ga', 'MB-TTT', C[4], 'v'), ('ltm_fixed100', 'LTM-like A3', C[6], 'P'),
       ('og', 'OG-A3', C[1], 's'), ('og_lite', 'OG-A3-lite', '#555555', 'h'), ('og_cho', 'OG-A3 + CHO-like', C[2], '^'), ('og_ltm', 'OG-A3 + LTM-like', C[5], 'X')]
for ax, D in zip(axs, [25, 50, 100]):
    tt = ['fixed40', 'fixed100', 'fixed256', 'fixed512']
    x = [m(D, t, 'pp3')*100 for t in tt]; y = [m(D, t, 'mc_mean')*1e3 for t in tt]
    ax.errorbar(x, y, xerr=[ci(D, t, 'pp3')*100 for t in tt], yerr=[ci(D, t, 'mc_mean')*1e3 for t in tt],
                fmt='-o', color=C[0], label='A3, fixed TTT', ms=5, mec='white', mew=0.6, capsize=2)
    for xi, yi, t in zip(x, y, ['40', '100', '256', '512']):
        ax.annotate(t, (xi, yi), textcoords='offset points', xytext=(4, 3), fontsize=6.5, color=INK2)
    for tag, lab, c, mk in pts:
        ax.errorbar(m(D, tag, 'pp3')*100, m(D, tag, 'mc_mean')*1e3, xerr=ci(D, tag, 'pp3')*100,
                    yerr=ci(D, tag, 'mc_mean')*1e3, fmt=mk, color=c, label=lab, ms=7, mec='white', mew=0.6, capsize=2)
    ax.set_title('$d_{ISD}$ = %d m' % D, fontsize=9); ax.set_yscale('log'); ax.set_ylim(5, 700); ax.set_xlim(-1, 32)
axs[0].set_ylabel('Mean pre-HO gap (ms, log)')
fig.supxlabel('Ping-pong: returns to the previous cell within 3 s (%)', fontsize=9, y=0.13)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=4, fontsize=7.2, bbox_to_anchor=(0.5, -0.07))
plt.tight_layout(rect=(0, 0.1, 1, 1)); plt.savefig('fig5_tradeoff.png')
# ---- Fig 6 per speed, D=50
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7))
tags = [('fixed100', 'A3-100', C[0]), ('mse', 'MSE', C[3]), ('ga', 'MB-TTT', C[4]), ('og', 'OG-A3', C[1])]
sp = ['3', '30', '60', '120']; w = 0.2; X = np.arange(4)
for j, (tag, lab, c) in enumerate(tags):
    for ax, key, sc in [(axs[0], 'mean', 1e3), (axs[1], 'pp', 100)]:
        vals = [np.mean([x['by_speed'][s][key] for x in A[(50, tag)]])*sc for s in sp]
        ax.bar(X + (j-1.5)*w, vals, w*0.9, color=c, label=lab)
for ax, yl in [(axs[0], 'Mean pre-HO gap (ms)'), (axs[1], 'Ping-pong ratio (%)')]:
    ax.set_xticks(X); ax.set_xticklabels(sp); ax.set_xlabel('UE speed (km/h)'); ax.set_ylabel(yl); ax.grid(axis='x', visible=False)
axs[0].legend(fontsize=7, ncol=2); axs[0].set_title('(a) Pre-HO gap per speed class', fontsize=9); axs[1].set_title('(b) Ping-pong per speed class', fontsize=9)
plt.tight_layout(); plt.savefig('fig6_speed.png'); print('ok')
