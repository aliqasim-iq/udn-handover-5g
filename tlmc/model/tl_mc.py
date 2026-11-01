"""
Trigger-level Monte-Carlo (TL-MC) simulator used to validate the analytical
pre-HO gap model and to compare TTT policies.

It models: NxN square grid on a torus (wrap-around, no edge effects),
3GPP TR 38.901 UMi-SC path loss, optional LOS/NLOS state and spatially
correlated shadowing (exponential correlation, Gudmundson), L1 averaging
over the measurement period, NR L3 filtering, event A3 with hysteresis and
TTT, HO preparation and execution delays, serving-cell SINR outage, and
ping-pong detection. It does NOT model scheduling, HARQ, RLC or transport
protocols; those are left to the ns-3/5G-LENA system-level evaluation.
"""
import numpy as np
from scipy.stats import norm
import udn_model as um


def pl_umi(d3d, d2d, los, fc=3.5, hut=1.5):
    pl_l = 32.4 + 21 * np.log10(d3d) + 20 * np.log10(fc)
    pl_n = 35.3 * np.log10(d3d) + 22.4 + 21.3 * np.log10(fc) - 0.3 * (hut - 1.5)
    return np.where(los, pl_l, np.maximum(pl_l, pl_n))


def pr_los(d2d):
    d = np.maximum(d2d, 1e-3)
    p = 18 / d + np.exp(-d / 36) * (1 - 18 / d)
    return np.where(d <= 18, 1.0, p)


def run(D=50.0, speeds_kmh=(30,), ttt=0.1, hys=3.0, Tm=0.04, k=4,
        policy="fixed", N=6, n_ue=300, T=60.0, dt=0.005, warm=3.0,
        shadow="none", sigma_l1=0.0, eta=1.0, gamma_min=-6.0,
        t_prep=0.02, t_exec=0.05, ptx=30.0, fc=3.5, noise_dbm=-89.0,
        hbs=10.0, hut=1.5, pp_win=1.0, pp_win_block=None, l3=True, seed=0, ga_margin=1.0,
        speed_err=0.0, init=None, turn_dist=None, trace=None, first_entry=False,
        kappa_q=None, mse_base=0.256, mse_sf=(1.0, 0.75, 0.5),
        guard=False, g_noreturn=False, a3_noreturn=False, g_thr=-3.0, g_hys=0.0, g_tavg=0.02, g_adapt=False, g_slope=None, g_sigma=None):
    rng = np.random.default_rng(seed)
    a = 2.0 ** (-k / 4.0) if l3 else 1.0
    L = N * D
    bx, by = np.meshgrid((np.arange(N) + 0.5) * D, (np.arange(N) + 0.5) * D)
    bx, by = bx.ravel(), by.ravel()
    C = bx.size
    dh = hbs - hut

    # UE speeds (mixed population allowed)
    sp = np.array(speeds_kmh, float) / 3.6
    v = sp[rng.integers(0, sp.size, n_ue)]
    px, py = rng.uniform(0, L, n_ue), rng.uniform(0, L, n_ue)
    ph = rng.uniform(0, 2 * np.pi, n_ue)
    walked = rng.uniform(0, D, n_ue)
    if init is not None:
        px, py, ph = [np.array(q, float) for q in init]
    turn_d = D if turn_dist is None else turn_dist

    # spatial processes (per UE-cell link)
    z = rng.standard_normal((n_ue, C))   # shadowing OU, unit variance
    u = rng.standard_normal((n_ue, C))   # LOS-state OU
    dc_sf_los, dc_sf_nlos, dc_los = 10.0, 13.0, 50.0
    s_los, s_nlos = 4.0, 7.82

    # per-UE TTT
    # SINR guard (SG-A3): fast L1 SINR/RSRP averages (EMA over g_tavg)
    ga_ = min(1.0, dt / g_tavg)
    sinr_f = None; rsrp_f = None
    if guard and g_adapt:
        # speed-aware guard threshold (Eq. 30): gamma_g = gamma_min + s_hat*(T_g + T_prep)
        thr_ue = gamma_min + g_slope * v * (g_tavg + t_prep) + 0.5
    else:
        thr_ue = np.full(n_ue, g_thr)
    if policy == "fixed":
        ttt_ue = np.full(n_ue, ttt)
    elif policy == "ga":
        vh = v * (1 + speed_err * rng.standard_normal(n_ue))
        ttt_ue = np.array([um.ga_ttt_kappa(x, kappa_q, Tm, t_prep, k, l3)
                           for x in vh])
    elif policy == "mse":   # LTE-style speed-state TTT scaling (oracle state)
        vk_ = v * 3.6
        sf = np.where(vk_ < 30, mse_sf[0], np.where(vk_ <= 60, mse_sf[1], mse_sf[2]))
        ttt_ue = mse_base * sf
    else:
        raise ValueError(policy)

    mstep = int(round(Tm / dt))
    nsteps = int(round(T / dt))
    serv = np.full(n_ue, -1)
    state = np.zeros(n_ue, int)           # 0 attached, 1 preparing, 2 exec
    tgt = np.zeros(n_ue, int)
    t_cmd = np.zeros(n_ue); t_end = np.zeros(n_ue)
    enter = np.full((n_ue, C), np.inf)
    F = np.full((n_ue, C), np.nan)
    acc = np.zeros((n_ue, C)); cnt = 0
    last_good = np.zeros(n_ue)
    wait_first = np.zeros(n_ue, bool); t_end_prev = np.zeros(n_ue)
    prev_cell = np.full(n_ue, -1); prev_ho_t = np.full(n_ue, -np.inf)

    rec_tpre, rec_tpost, rec_v, rec_pp = [], [], [], []
    rec_cos, rec_y = [], []
    hraw = np.full((n_ue, C), np.inf)          # unfiltered A3-entry times
    vtrk = np.zeros(n_ue, bool); vsrc = np.zeros(n_ue, int)
    vth = np.zeros(n_ue); vrec = np.zeros(n_ue, int)
    rec_K = []; rec_ppv = []; rec_pp3 = []
    n_ho = 0; out_slots = 0; tot_slots = 0; n_guard = 0
    idx = np.arange(n_ue)

    for n in range(nsteps):
        t = n * dt
        # ---- mobility (random walk, turn every D metres)
        step = v * dt
        px = (px + step * np.cos(ph)) % L
        py = (py + step * np.sin(ph)) % L
        walked += step
        turn = walked >= turn_d
        ph[turn] = rng.uniform(0, 2 * np.pi, turn.sum()); walked[turn] = 0
        # ---- geometry (torus minimum image)
        dx = (px[:, None] - bx[None, :] + L / 2) % L - L / 2
        dy = (py[:, None] - by[None, :] + L / 2) % L - L / 2
        d2 = np.hypot(dx, dy)
        d3 = np.sqrt(d2 ** 2 + dh ** 2)
        if shadow == "none":
            pl = um.pl_eff(d2, dh, fc)
        if shadow == "full":
            r = np.exp(-step / dc_los)[:, None]
            u = r * u + np.sqrt(1 - r ** 2) * rng.standard_normal((n_ue, C))
            los = norm.cdf(u) < pr_los(d2)
        else:
            los = np.ones((n_ue, C), bool)
        if shadow != "none":
            pl = pl_umi(d3, d2, los, fc, hut)
        if shadow in ("los", "full"):
            dcs = np.where(los, dc_sf_los, dc_sf_nlos)
            r = np.exp(-step[:, None] / dcs)
            z = r * z + np.sqrt(1 - r ** 2) * rng.standard_normal((n_ue, C))
            pl = pl - np.where(los, s_los, s_nlos) * z
        rsrp = ptx - pl
        plin = 10 ** (rsrp / 10)
        if n == 0:
            serv = np.argmax(rsrp, axis=1)
        ps = plin[idx, serv]
        sinr = 10 * np.log10(ps / (eta * (plin.sum(1) - ps) + 10 ** (noise_dbm / 10)))
        good = (sinr >= gamma_min) & (state != 2)
        # unfiltered A3 entry (for the conditional-model validation)
        rs = rsrp[idx, serv][:, None]
        craw = (rsrp - rs > hys); craw[idx, serv] = False
        hraw = np.where(craw, np.minimum(hraw, t), hraw if first_entry else np.inf)
        # virtual tracking of the source link after HandoverStart
        if vtrk.any():
            vi = np.where(vtrk)[0]
            pv = plin[vi, vsrc[vi]]
            sv = 10 * np.log10(pv / (eta * (plin[vi].sum(1) - pv) + 10 ** (noise_dbm / 10)))
            hit = (sv < gamma_min) | (t - vth[vi] > 3.0)
            for q in vi[hit]:
                rec_K[vrec[q]] = (t - vth[q]) if t - vth[q] <= 3.0 else np.inf
            vtrk[vi[hit]] = False
        if trace is not None:
            j = trace['ue']
            rr = rsrp[j].copy(); rr[serv[j]] = -1e9
            fr = F[j].copy(); fr[serv[j]] = -1e9
            trace.setdefault('rows', []).append((t, px[j], py[j], serv[j], state[j], sinr[j],
                                                 rr.max() - rsrp[j, serv[j]],
                                                 np.nanmax(fr) - F[j, serv[j]] if not np.isnan(F[j,0]) else np.nan,
                                                 np.argmax(rr)))
        if t >= warm:
            tot_slots += n_ue; out_slots += (~good).sum()
        upd = good & (state != 2)
        last_good[upd] = t
        nf = wait_first & good
        if t >= warm:
            rec_tpost.extend((t - t_end_prev[nf]).tolist())
        wait_first[nf] = False

        # ---- L1 averaging and L3 filtering
        acc += rsrp; cnt += 1
        if (n + 1) % mstep == 0:
            M = acc / cnt + sigma_l1 * rng.standard_normal((n_ue, C))
            acc[:] = 0; cnt = 0
            F = np.where(np.isnan(F), M, (1 - a) * F + a * M)
            fs = F[idx, serv][:, None]
            cond = (F - fs > hys)
            cond[idx, serv] = False
            cond &= (state == 0)[:, None]
            enter = np.where(cond, np.minimum(enter, t), np.inf)
        # ---- TTT expiry -> measurement report -> preparation
        cand = (state == 0) & np.any(enter <= t - ttt_ue[:, None] + 1e-9, axis=1)
        n_guard_now = 0
        if guard:
            # guard SINR estimate: own estimation-noise level (default: same as L1 RSRP error)
            sm = sinr + (sigma_l1 if g_sigma is None else g_sigma) * rng.standard_normal(n_ue)
            rm = rsrp + sigma_l1 * rng.standard_normal((n_ue, C))
            sinr_f = sm if sinr_f is None else (1 - ga_) * sinr_f + ga_ * sm
            rsrp_f = rm if rsrp_f is None else (1 - ga_) * rsrp_f + ga_ * rm
            rn = rsrp_f.copy(); rn[idx, serv] = -np.inf
            if g_noreturn:   # no guard-triggered return to the cell left < pp_win ago
                blk = (prev_cell >= 0) & (t - prev_ho_t <= (pp_win if pp_win_block is None else pp_win_block))
                rn[idx[blk], prev_cell[blk]] = -np.inf
            best = np.argmax(rn, axis=1)
            gfire = (state == 0) & ~cand & (sinr_f < thr_ue) & \
                    (rn[idx, best] > rsrp_f[idx, serv] + g_hys)
            if gfire.any():
                gi = np.where(gfire)[0]
                tgt[gi] = best[gi]; state[gi] = 1; t_cmd[gi] = t + t_prep
                enter[gi] = np.inf
                if t >= warm:
                    n_guard += gi.size
        if cand.any():
            ci = np.where(cand)[0]
            Fm = np.where(enter[ci] <= t - ttt_ue[ci, None] + 1e-9, F[ci], -np.inf)
            if a3_noreturn:
                blk = (prev_cell[ci] >= 0) & (t - prev_ho_t[ci] <= pp_win)
                Fm[np.where(blk)[0], prev_cell[ci][blk]] = -np.inf
            okc = np.isfinite(Fm.max(axis=1))
            ci = ci[okc]; Fm = Fm[okc]
            tgt[ci] = np.argmax(Fm, axis=1)
            state[ci] = 1; t_cmd[ci] = t + t_prep
            enter[ci] = np.inf
        # ---- HO command (HandoverStart)
        hs = (state == 1) & (t >= t_cmd - 1e-9)
        if hs.any():
            hi = np.where(hs)[0]
            th_ = hraw[hi, tgt[hi]]
            gap_ = np.maximum(t - last_good[hi] - dt, 0)
            if t >= warm:
                for q, gq, tq in zip(hi, gap_, th_):
                    if gq > 0:
                        rec_K.append((last_good[q] + dt) - tq)
                    else:
                        rec_K.append(np.nan)
                        vtrk[q] = True; vsrc[q] = serv[q]; vth[q] = tq
                        vrec[q] = len(rec_K) - 1
                rec_tpre.extend(gap_.tolist())
                rec_v.extend(v[hi].tolist())
                ex = (bx[tgt[hi]] - bx[serv[hi]] + L / 2) % L - L / 2
                ey = (by[tgt[hi]] - by[serv[hi]] + L / 2) % L - L / 2
                nn = np.hypot(ex, ey); ex, ey = ex / nn, ey / nn
                rec_cos.extend((np.cos(ph[hi]) * ex + np.sin(ph[hi]) * ey).tolist())
                mx = (px[hi] - bx[serv[hi]] + L / 2) % L - L / 2
                my = (py[hi] - by[serv[hi]] + L / 2) % L - L / 2
                rec_y.extend((-mx * ey + my * ex).tolist())
            state[hi] = 2; t_end[hi] = t + t_exec
        # ---- HO completion
        he = (state == 2) & (t >= t_end - 1e-9)
        if he.any():
            ei = np.where(he)[0]
            if t >= warm:
                n_ho += ei.size
                pp = (tgt[ei] == prev_cell[ei]) & (t - prev_ho_t[ei] <= pp_win)
                rec_pp.extend(pp.tolist()); rec_ppv.extend(v[ei].tolist())
                rec_pp3.extend(((tgt[ei] == prev_cell[ei]) & (t - prev_ho_t[ei] <= 3.0)).tolist())
            prev_cell[ei] = serv[ei]; prev_ho_t[ei] = t
            serv[ei] = tgt[ei]; state[ei] = 0
            wait_first[ei] = True; t_end_prev[ei] = t
            enter[ei] = np.inf; hraw[ei] = np.inf
    dur = T - warm
    return dict(tpre=np.array(rec_tpre), tpost=np.array(rec_tpost),
                v=np.array(rec_v), pp=np.array(rec_pp, bool),
                ho_rate=n_ho / (n_ue * dur), outage=out_slots / max(tot_slots, 1),
                ttt_ue=ttt_ue, speeds=v, cos=np.array(rec_cos), y=np.array(rec_y), K=np.array(rec_K), ppv=np.array(rec_ppv), n_guard=n_guard, pp3=np.array(rec_pp3, bool))
