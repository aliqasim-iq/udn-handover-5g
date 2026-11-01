import json, numpy as np, matplotlib.pyplot as plt
from style import *
apply()
V = json.load(open('../model/res_val.json'))
# ---- Fig 3
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw=dict(width_ratios=[1, 1.35]))
det = [r for r in V if r['tag'] == 'val_det']
for i, D in enumerate([25, 50, 100]):
    xs = [r['mc_mean_ok']*1e3 for r in det if r['D'] == D]; ys = [r['pred_mean']*1e3 for r in det if r['D'] == D]
    a.plot(xs, ys, MK[i], color=C[i], ls='none', label='$d_{ISD}$ = %d m' % D, mec='white', mew=0.6, ms=6)
m = 650; a.plot([0, m], [0, m], color=INK2, lw=0.9, ls='--', zorder=0)
a.set_xlim(0, m); a.set_ylim(0, m); a.set_aspect('equal')
a.set_xlabel('TL-MC mean pre-HO gap (ms)'); a.set_ylabel('Model, Eq. (14) (ms)')
a.set_title('(a) Deterministic propagation, $R^2$ = 0.998', fontsize=9); a.legend(loc='upper left')
full = [r for r in V if r['tag'] == 'val_full']
keys = sorted({(r['D'], r['speeds'][0], r['Tm']) for r in full})
lo = {(r['D'], r['speeds'][0], r['Tm']): r['pred_mean'] for r in full if r['first_entry']}
hi = {(r['D'], r['speeds'][0], r['Tm']): r['pred_mean'] for r in full if not r['first_entry']}
mcv = {(r['D'], r['speeds'][0], r['Tm']): r['mc_mean_ok'] for r in full if r['first_entry']}
keys = sorted(keys, key=lambda k: (k[2], k[0], k[1]))
x = np.arange(len(keys))
b.vlines(x, [lo[k]*1e3 for k in keys], [hi[k]*1e3 for k in keys], color=C[0], lw=5, alpha=0.35, label='model bounds (first / last A3 entry)')
b.plot(x, [mcv[k]*1e3 for k in keys], 'o', color=C[1], ms=5, mec='white', mew=0.6, label='TL-MC (shadowed)')
b.set_xticks(x); b.set_xticklabels(['%d\n%d' % (k[0], k[1]) for k in keys], fontsize=6)
b.axvline(8.5, color=INK2, lw=0.8, ls=':')
b.text(4, 535, '$T_m$ = 40 ms', ha='center', fontsize=8, color=INK2); b.text(13, 535, '$T_m$ = 200 ms', ha='center', fontsize=8, color=INK2)
b.set_xlabel('$d_{ISD}$ (m) / speed (km/h)', fontsize=8); b.set_ylabel('Mean pre-HO gap (ms)')
b.set_ylim(0, 570); b.legend(loc='upper left', bbox_to_anchor=(0, 0.9), fontsize=7); b.set_title('(b) Shadowed propagation, $T_{TTT}$ = 100 ms', fontsize=9)
plt.tight_layout(); plt.savefig('fig3_validation.png')
# ---- Fig 4
K = json.load(open('../model/res_kappa.json'))
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
for ax, sh, ttl in [(axs[0], 'none', '(a) Deterministic'), (axs[1], 'full', '(b) LOS/NLOS + shadowing')]:
    for i, D in enumerate([25, 35, 50, 75, 100]):
        r = [q for q in K if q['D'] == D and q['shadow'] == sh][0]
        k = np.sort(np.array(r['kappa'])); k = k[np.isfinite(k)]
        ax.plot(k, np.arange(1, k.size+1)/k.size, color=C[i], label='%d m' % D, lw=1.6)
    ax.axvline(0, color=INK2, lw=0.8, ls=':'); ax.set_xlim(-6, 20); ax.set_xlabel('Outage-margin distance κ (m)'); ax.set_title(ttl, fontsize=9)
axs[0].set_ylabel('CDF'); axs[0].legend(title='$d_{ISD}$', loc='lower right', fontsize=7, title_fontsize=7)
plt.tight_layout(); plt.savefig('fig4_kappa.png'); print('ok')
