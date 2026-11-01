"""
Analytical model of the pre-handover service gap (T_pre) in ultra-dense
small-cell networks, and the gradient-aware time-to-trigger (GA-TTT) rule.

All quantities in SI units (m, s) unless stated; powers in dB.
Notation follows the manuscript (Section 4).
"""
import numpy as np
from scipy.stats import norm

LN10 = np.log(10.0)
NR_TTT_MS = np.array([0, 40, 64, 80, 100, 128, 160, 256, 320, 480, 512, 640,
                      1024, 1280, 2560, 5120])  # TS 38.331 TimeToTrigger


# ---------------------------------------------------------------- geometry
def gradient(D, y=0.0, alpha=2.1, dh=8.5):
    """Spatial gradient (dB/m) of the mean RSRP difference across the
    equal-power boundary of two sites separated by D, at offset y along
    the boundary (Eq. 10)."""
    L = D / 2.0
    return (20.0 * alpha / LN10) * L / (L ** 2 + y ** 2 + dh ** 2)


def gradient_eff(D, y=0.0, dh=8.5, fc_ghz=3.5, h=0.01):
    """Gradient (dB/m) of the mean RSRP difference for the composite
    LOS/NLOS path-loss law pl_eff, by central differences."""
    L = D / 2.0
    def delta(x):
        ds = np.hypot(L + x, y); dn = np.hypot(L - x, y)
        return pl_eff(ds, dh, fc_ghz) - pl_eff(dn, dh, fc_ghz)
    return (delta(h) - delta(-h)) / (2 * h)


def gradient_harmonic(D, alpha=2.1, dh=8.5):
    """1 / E_y[1/g] for crossings uniform along a square-grid edge of
    length D (Eq. 11)."""
    L = D / 2.0
    return 20.0 * alpha * L / (LN10 * (4.0 * L ** 2 / 3.0 + dh ** 2))


def steepest_isd(dh=8.5):
    """ISD that maximises the harmonic-mean gradient (Eq. 12)."""
    return np.sqrt(3.0) * dh


# ---------------------------------------------------------------- outage
def pl_los(d3d, fc_ghz=3.5):
    return 32.4 + 21.0 * np.log10(d3d) + 20.0 * np.log10(fc_ghz)


def pr_los(d2d):
    d = np.maximum(d2d, 1e-3)
    p = 18 / d + np.exp(-d / 36) * (1 - 18 / d)
    return np.where(d <= 18, 1.0, p)


def pl_eff(d2d, dh=8.5, fc_ghz=3.5):
    """Mean (LOS-probability-weighted, linear domain) UMi-SC path loss in dB,
    used by the deterministic model (Eq. 9)."""
    d3 = np.sqrt(d2d ** 2 + dh ** 2)
    pll = pl_los(d3, fc_ghz)
    pln = np.maximum(pll, 35.3 * np.log10(d3) + 22.4 + 21.3 * np.log10(fc_ghz))
    p = pr_los(d2d)
    return -10 * np.log10(p * 10 ** (-pll / 10) + (1 - p) * 10 ** (-pln / 10))


def rho_interference(D, y=0.0, eta=1.0, R=6, dh=8.5, fc_ghz=3.5,
                     ptx_dbm=30.0, noise_dbm=-89.0):
    """Other-cell interference + noise, relative to the serving power, at a
    point of the boundary between sites (0,0) and (D,0), offset y."""
    px, py = D / 2.0, y
    ps = None
    acc = 0.0
    for i in range(-R, R + 2):
        for j in range(-R, R + 1):
            d2 = np.hypot(px - i * D, py - j * D)
            p = ptx_dbm - pl_eff(d2, dh, fc_ghz)
            if (i, j) in ((0, 0), (1, 0)):
                if (i, j) == (0, 0):
                    ps = p
                continue
            acc += eta * 10 ** (p / 10)
    acc += 10 ** (noise_dbm / 10)
    return acc / 10 ** (ps / 10)


def delta_out(rho, gamma_min_db=-6.0, eta=1.0):
    """RSRP difference (neighbour minus serving, dB) at which the serving
    SINR falls below gamma_min (Eq. 14); the neighbour (target) cell
    interferes with activity factor eta. Returns -inf for a coverage hole."""
    arg = (10 ** (-gamma_min_db / 10.0) - rho) / eta
    return np.where(arg > 0, 10 * np.log10(np.maximum(arg, 1e-30)), -np.inf)


# ---------------------------------------------------------------- delays
def decision_delay(Tm, ttt, t_prep, k=4):
    """T_D = T_m / a + T_TTT + T_prep (Eq. 17) with a = 2^(-k/4)."""
    a = 2.0 ** (-k / 4.0)
    return Tm / a + ttt + t_prep


# ---------------------------------------------------------------- T_pre
def phi(c0):
    """Normalised mean pre-HO gap phi(c0) (Eq. 20)."""
    c0 = np.asarray(c0, dtype=float)
    out = np.zeros_like(c0)
    neg = c0 <= 0
    mid = (c0 > 0) & (c0 < 1)
    out[neg] = 1.0 - 0.5 * np.pi * c0[neg]
    cm = c0[mid]
    out[mid] = np.sqrt(1 - cm ** 2) - cm * np.arccos(cm)
    return out


def mean_tpre_det(v, D, Hys, T_D, dout, alpha=None, dh=8.5, ny=41):
    """Mean pre-HO gap of the drift-only model, averaged over the crossing
    angle (closed form, Eq. 19-20) and numerically over the crossing point
    along the edge. dout may be an array over y (same length as ny)."""
    ys = np.linspace(-D / 2, D / 2, ny)
    g = (gradient(D, ys, alpha, dh) if alpha else gradient_eff(D, ys, dh))
    dout = np.broadcast_to(dout, ys.shape)
    c0 = (dout - Hys) / (g * v * T_D)
    return float(np.mean(T_D * phi(c0)))


def critical_speed(D, Hys, T_D, dout, alpha=2.1, dh=8.5):
    """Speed above which a perpendicular crossing at the edge midpoint
    suffers a pre-HO gap (Eq. 21)."""
    return (dout - Hys) / (gradient(D, 0.0, alpha, dh) * T_D)


def tail_prob(t, T_D, K):
    """P(T_pre > t) in the margin-limited regime (Hys > Delta_out),
    K = (Hys - Delta_out)/(g v) (Eq. 22)."""
    x = (t - T_D) / K
    return np.where(x <= 1, 1.0, 1 - np.sqrt(1 - 1 / np.maximum(x, 1) ** 2))


# --------------------------------------------- drift-diffusion extension
def fp_cdf(t, m, mu, beta):
    """P(first passage of level m>0 by Brownian motion with drift mu and
    variance rate beta occurs before t) (inverse-Gaussian CDF)."""
    t = np.maximum(t, 1e-12)
    s = np.sqrt(beta * t)
    a1 = norm.logcdf((mu * t - m) / s)
    a2 = 2 * mu * m / beta + norm.logcdf((-m - mu * t) / s)
    return np.exp(a1) + np.exp(a2)


def mean_tpre_dd(v, D, Hys, T_D, dout, sigma, dc, alpha=2.1, dh=8.5,
                 ny=21, nth=41, nt=200):
    """Mean pre-HO gap of the drift-diffusion model (Eq. 24):
    E[(T_D - tau)^+] = int_0^T_D P(tau <= t) dt, averaged over the
    flux-weighted crossing angle and the crossing point."""
    beta = 4.0 * sigma ** 2 * v / dc       # variance rate of shadowing diff.
    ys = np.linspace(-D / 2, D / 2, ny)
    dout = np.broadcast_to(dout, ys.shape)
    th = np.linspace(-np.pi / 2, np.pi / 2, nth)
    w_th = np.cos(th); w_th /= w_th.sum()
    ts = np.linspace(0, T_D, nt)
    res = 0.0
    for yi, y in enumerate(ys):
        m = dout[yi] - Hys
        g = gradient(D, y, alpha, dh)
        for thi, t_ in enumerate(th):
            mu = g * v * np.cos(t_)
            if m <= 0:
                val = T_D + (-m) / max(mu, 1e-9)
            else:
                val = np.trapezoid(fp_cdf(ts, m, mu, beta), ts)
            res += w_th[thi] * val
    return res / len(ys)


# ---------------------------------------------------------------- GA-TTT
def ga_ttt(v_hat, D, Hys, Tm, t_prep, dout, margin=1.0, k=4, alpha=2.1,
           dh=8.5, ttt_set_ms=NR_TTT_MS, ttt_cap_ms=640):
    """Gradient-aware TTT (Algorithm 1): the largest standard TTT such that
    a perpendicular crossing reaches the HO command before the serving
    link enters outage, keeping a safety margin (dB)."""
    a = 2.0 ** (-k / 4.0)
    s = gradient_harmonic(D, alpha, dh) * max(v_hat, 0.1)
    ttt_max = (dout - margin - Hys) / s - (Tm / a + t_prep)
    cand = ttt_set_ms[(ttt_set_ms <= ttt_max * 1e3) & (ttt_set_ms <= ttt_cap_ms)]
    return (cand.max() if cand.size else 0) / 1e3


# ------------------------------------------- conditional (per-crossing) model
def tpre_conditional(K, Tm, ttt, t_prep, k=4):
    """Mean pre-HO gap for a crossing with outage margin time K (s),
    averaged over the uniform phase of the measurement grid (Eq. 19):
    E_U[(T_D0 + U*T_m - K)^+], T_D0 = T_m(1/a - 1/2) + T_TTT + T_prep."""
    a = 2.0 ** (-k / 4.0)
    A = Tm * (1 / a - 0.5) + ttt + t_prep - K - 0.0  # = T_D0 - K
    A = np.asarray(A, float)
    return np.where(A >= 0, A + Tm / 2,
                    np.where(A <= -Tm, 0.0, (A + Tm) ** 2 / (2 * Tm)))


def outage_margin_time(D, y, vperp, Hys, eta=0.5, gamma_min=-6.0, dh=8.5):
    """K = (Delta_out(y) - Hys) / (g(y) * v_perp) for one crossing."""
    rho = rho_interference(D, y, eta, dh=dh)
    do = float(delta_out(rho, gamma_min, eta))
    return (do - Hys) / (gradient_eff(D, y, dh) * vperp), do


def ga_ttt_kappa(v_hat, kappa_q, Tm, t_prep, k=4, l3=True,
                 ttt_set_ms=NR_TTT_MS, ttt_cap_ms=640):
    """GA-TTT (Algorithm 1): largest standard TTT such that the distance
    travelled during the worst-phase decision delay does not exceed the
    q-quantile outage-margin distance kappa_q of the cell (Eq. 27)."""
    a = 2.0 ** (-k / 4.0) if l3 else 1.0
    ttt_max = kappa_q / max(v_hat, 0.3) - (Tm / a + Tm / 2 + t_prep)
    cand = ttt_set_ms[(ttt_set_ms <= ttt_max * 1e3) & (ttt_set_ms <= ttt_cap_ms)]
    return (cand.max() if cand.size else 0) / 1e3
