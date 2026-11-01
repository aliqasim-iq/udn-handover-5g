import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from style import *
apply()
specs = {
 'ph_fig7_density.png': ('PLACEHOLDER — Figure 7 (5G-LENA, d_ISD sweep)', [
   '(a) x: small-cell density (100–1600 cells/km², log); y: time (ms); series: T_exec, T_pre, T_int (mean ± 95% CI)',
   '(b) x: density; y: PLR (%); series: all losses, HO-attributed losses',
   '(c) x: density; y: HO rate (1/s per UE); series: simulated, Eq. (9) reference',
   '(d) x: density; y: ping-pong ratio (%) and RLF per UE-hour',
   'Config: A3 3 dB / 100 ms, 30 km/h, T_m = 40 ms and 200 ms (two line styles), UDP CBR']),
 'ph_fig8_speed.png': ('PLACEHOLDER — Figure 8 (5G-LENA, speed sweep)', [
   'x: UE speed 3, 30, 60, 90, 120 km/h at d_ISD = 50 m',
   'panels: (a) T_pre and T_int (ms); (b) PLR (%); (c) HO rate; (d) ping-pong (%)',
   'overlay in (a): model prediction Eq. (14) using measured κ distribution (dashed)']),
 'ph_fig9_ttt.png': ('PLACEHOLDER — Figure 9 (5G-LENA, TTT × Hys)', [
   'd_ISD = 50 m, 60 km/h; x: TTT 40, 100, 256, 512 ms; series: Hys 1, 3, 5 dB',
   'panels: (a) T_int (ms); (b) PLR (%); (c) ping-pong (%)',
   'row 2 (optional): same with T_m = 200 ms']),
 'ph_fig10_schemes.png': ('PLACEHOLDER — Figure 10 (5G-LENA, scheme comparison)', [
   'mixed-speed population (3/30/60/120 km/h), d_ISD = 25, 50, 100 m (groups)',
   'policies: A3-40, A3-100, A3-256, MSE, MB-TTT, OG-A3, OG-A3 + CHO-like, LTM-like',
   'panels: (a) mean T_int (ms) with 95% CI; (b) PLR (%); (c) ping-pong (%); (d) HO rate (1/s)',
   'and a scatter: outage fraction vs ping-pong (same layout as Figure 5)']),
 'ph_fig11_transport.png': ('PLACEHOLDER — Figure 11 (5G-LENA, transport/app impact)', [
   '(a) CDF of per-HO service interruption T_int for A3-100 vs OG-A3 (UDP)',
   '(b) TCP CUBIC goodput (Mb/s) per policy, RLC AM, 60 km/h, d_ISD = 50 m',
   '(c) URLLC-like flow (32 B every 1 ms): CDF of packet delay; mark 10 ms deadline']),
}
for fn, (title, lines) in specs.items():
    fig, ax = plt.subplots(figsize=(7.2, 2.4)); ax.axis('off')
    ax.add_patch(FancyBboxPatch((0.01, 0.03), 0.98, 0.94, boxstyle='round,pad=0.0,rounding_size=0.02', fc='#fffbe6', ec='#c9a400', lw=1.4, ls='--', transform=ax.transAxes))
    ax.text(0.5, 0.86, title, ha='center', va='center', fontsize=10, weight='bold', color='#7a6200', transform=ax.transAxes)
    for i, l in enumerate(lines):
        ax.text(0.04, 0.68 - i*0.14, '• ' + l, ha='left', va='center', fontsize=7.6, color=INK, transform=ax.transAxes)
    plt.savefig(fn); plt.close()
print('ok')
