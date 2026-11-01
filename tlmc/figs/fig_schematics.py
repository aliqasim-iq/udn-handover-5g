import numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from style import *
apply()
# ---------------- Figure 1
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw=dict(width_ratios=[1, 1.15]))
N = 4
# central N x N tile and faint wrapped images of it (torus): every link uses the closest image
for di in (-1, 0, 1):
    for dj in (-1, 0, 1):
        ghost = (di, dj) != (0, 0)
        for i in range(N):
            for j in range(N):
                x, y = i + di * N, j + dj * N
                if not (-1.2 <= x <= N + 0.2 and -1.2 <= y <= N + 0.2):
                    continue
                ax.add_patch(Rectangle((x, y), 1, 1, fc='#f6f6f4' if ghost else '#eaf2fc', ec='#d6d6d1' if ghost else '#c9c9c4', lw=0.6))
                ax.plot(x + .5, y + .5, marker='^', color='#b9b9b4' if ghost else C[0], ms=5 if ghost else 6.5, ls='none')
ax.add_patch(Rectangle((0, 0), N, N, fc='none', ec=C[0], lw=1.6))
rng = np.random.default_rng(4)
for k in range(4):
    p = rng.uniform(0.4, N - 0.4, 2); pts = [p.copy()]
    for s_ in range(6):
        a_ = rng.uniform(0, 2*np.pi); p = (p + 0.9 * np.array([np.cos(a_), np.sin(a_)])) % N; pts.append(p.copy())
    pts = np.array(pts)
    for q in range(len(pts) - 1):
        if np.hypot(*(pts[q + 1] - pts[q])) < 1.0:
            ax.plot(pts[q:q + 2, 0], pts[q:q + 2, 1], color=C[1], lw=1.1)
    ax.plot(*pts[0], 'o', color=C[1], ms=4)
for (x0, y0, x1, y1) in [(N + 0.05, 2.5, N + 0.6, 2.5), (2.5, N + 0.05, 2.5, N + 0.6)]:
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='->', color=C[2], lw=1.3))
ax.text(N + 0.65, 2.75, 'wrap', color=C[2], fontsize=7, rotation=90)
ax.annotate('', xy=(1.0, -0.35), xytext=(0.0, -0.35), arrowprops=dict(arrowstyle='<->', color=INK2, lw=1))
ax.text(0.5, -0.55, '$d_{ISD}$', ha='center', va='top', color=INK2, fontsize=8)
ax.set_xlim(-1.25, N + 1.25); ax.set_ylim(-1.25, N + 1.25); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('(a) N×N grid on a torus (wrap-around)', fontsize=9)
ax.plot([], [], '^', color=C[0], label='gNB'); ax.plot([], [], '^', color='#b9b9b4', label='wrapped image')
ax.plot([], [], '-', color=C[1], label='UE random walk')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=7, handlelength=1.2)
bx.axis('off'); bx.set_xlim(0, 10); bx.set_ylim(0, 10)
def box(x, y, w, h, t, c):
    bx.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.08,rounding_size=0.25', fc='white', ec=c, lw=1.4))
    bx.text(x + w/2, y + h/2, t, ha='center', va='center', fontsize=7.5, color=INK)
box(2.4, 8.6, 5.2, 0.9, 'Remote host (DL UDP CBR)', INK2)
box(2.6, 6.5, 4.8, 1.0, 'EPC: SGW/PGW + MME\n(5G-LENA NR over EPC)', INK2)
for k, x in enumerate([0.3, 3.7, 7.1]):
    box(x, 3.9, 2.6, 0.9, 'NR gNB %d' % (k+1), C[0])
    bx.plot([x + 1.3, 5.0], [4.85, 6.45], color=INK2, lw=1)
bx.text(1.6, 5.75, 'S1-U', fontsize=7, color=INK2)
bx.plot([5, 5], [7.55, 8.55], color=INK2, lw=1)
for x0, x1 in [(2.95, 3.65), (6.35, 7.05)]:
    bx.plot([x0, x1], [4.35, 4.35], color=C[2], lw=1.6, ls='--')
bx.text(3.3, 4.55, 'X2', fontsize=7, color=C[2], ha='center'); bx.text(6.7, 4.55, 'X2', fontsize=7, color=C[2], ha='center')
box(3.9, 1.0, 2.2, 0.9, 'UE', C[1])
bx.plot([5.0, 1.6], [1.95, 3.85], color=C[0], ls='--', lw=1.1); bx.plot([5.0, 5.0], [1.95, 3.85], color=C[1], ls='--', lw=1.1)
bx.text(6.3, 1.45, 'L1 RSRP/SINR → L3 filter →\nA3 (Hys, TTT) + OG-A3 guard', fontsize=7, color=INK2, va='center')
bx.set_title('(b) Modelled protocol stack and interfaces', fontsize=9)
plt.tight_layout(); plt.savefig('fig1_topology.png')
# ---------------- Figure 2: decision timeline
fig, ax = plt.subplots(figsize=(7.2, 3.3))
t = np.linspace(-0.4, 0.75, 500); s = 20.0
D = s * t + 3.0
Tm = 0.08; a = 0.5
F = s * (t - Tm/2 - Tm*(1-a)/a) + 3.0
ax.plot(t, D, color=C[0], label=r'$\Delta(t)$: neighbour − serving RSRP')
ax.step(np.arange(-0.4, 0.75, Tm), s * (np.arange(-0.4, 0.75, Tm) - Tm/2 - Tm) + 3.0, where='post', color=C[2], label=r'$F(t)$: L1-averaged, L3-filtered, sampled')
dout = 7.2
ax.axhline(3.0, color=INK2, lw=0.9, ls=':'); ax.text(-0.39, 3.25, r'$Hys$', color=INK2, fontsize=8)
ax.axhline(dout, color=C[1], lw=0.9, ls=':'); ax.text(-0.39, dout + 0.25, r'$\Delta_{out}$ (serving SINR $=\gamma_{min}$)', color=C[1], fontsize=8)
t_h = 0.0; t_det = Tm/2 + Tm + 0.03; TTT = 0.16; tprep = 0.03
t_trig = t_det + TTT; t_start = t_trig + tprep; t_o = (dout - 3)/s
for x, lab, yy in [(t_h, r'$t_h$', 9.6), (t_det, 'A3 detected', 9.6), (t_trig, 'TTT expires', 9.6), (t_start, 'HandoverStart', 8.7)]:
    ax.axvline(x, color='#9a9a95', lw=0.8)
    ax.text(x + 0.008, yy, lab, fontsize=7, color=INK, va='top', bbox=dict(fc='white', ec='none', pad=0.5))
ax.axvspan(t_o, t_start, color=C[1], alpha=0.13)
ax.annotate('', xy=(t_start, 1.2), xytext=(t_o, 1.2), arrowprops=dict(arrowstyle='<->', color=C[1]))
ax.text((t_o + t_start)/2, 0.55, r'$T_{pre}=[T_D-K]^+$', color=C[1], ha='center', fontsize=8)
ax.annotate('', xy=(t_o, -0.6), xytext=(t_h, -0.6), arrowprops=dict(arrowstyle='<->', color=INK2))
ax.text(t_o/2, -1.25, r'$K=\kappa/v$', ha='center', color=INK2, fontsize=8)
ax.annotate('', xy=(t_start, -2.2), xytext=(t_h, -2.2), arrowprops=dict(arrowstyle='<->', color=C[0]))
ax.text(t_start/2 + 0.05, -2.9, r'$T_D = T_m/2 + T_m(1-a)/a + U\,T_m + T_{TTT} + T_{prep}$', ha='center', color=C[0], fontsize=8)
ax.set_xlim(-0.4, 0.75); ax.set_ylim(-3.4, 9.8)
ax.set_xlabel('time relative to the instant the A3 condition is met by the unfiltered RSRP (s)'); ax.set_ylabel('dB')
ax.legend(loc='center right', bbox_to_anchor=(1.0, 0.42), fontsize=7.5)
plt.tight_layout(); plt.savefig('fig2_timeline.png')
print('ok')
