"""Run all TL-MC experiments used in the manuscript and save results."""
import numpy as np, json, itertools, sys, time
from multiprocessing import Pool
import os
import udn_model as um, tl_mc as mc

SEEDS = tuple(range(1, int(os.environ.get('SEEDS', '10')) + 1))

ETA = 0.5
OUT = 'results.json'


def summarize(r, extra):
    K = r['K']; t = r['tpre']; ok = ~np.isnan(K)
    d = dict(extra)
    d.update(mc_mean=float(t.mean()), mc_p_gap=float((t > 0).mean()),
             mc_p90=float(np.percentile(t, 90)), n_ho=int(t.size),
             pp=float(r['pp'].mean()) if r['pp'].size else 0.0,
             ho_rate=float(r['ho_rate']), outage=float(r['outage']),
             tpost=float(r['tpost'].mean()) if r['tpost'].size else 0.0,
             guard=float(r.get('n_guard', 0)) / max(1, t.size),
             pp3=float(r['pp3'].mean()) if r['pp3'].size else 0.0)
    if 'Tm' in extra:
        Kk = np.where(np.isinf(K[ok]), 1e9, K[ok])
        d['pred_mean'] = float(um.tpre_conditional(Kk, extra['Tm'], extra['ttt'],
                                                   extra.get('t_prep', 0.02)).mean())
        d['mc_mean_ok'] = float(t[ok].mean())
        d['kappa'] = (Kk * r['v'][ok]).tolist() if extra.get('keep_kappa') else None
        d['kneg'] = float((K[ok] < 0).mean())
    # per-speed breakdown
    vs = np.round(r['v'] * 3.6).astype(int)
    d['by_speed'] = {int(s): dict(mean=float(t[vs == s].mean()),
                                  p_gap=float((t[vs == s] > 0).mean()),
                                  pp=float(r['pp'][np.round(r['ppv'] * 3.6).astype(int) == s].mean()),
                                  n=int((vs == s).sum()))
                     for s in np.unique(vs)}
    return d


def job(p):
    p = dict(p)
    tag = p.pop('tag'); keep = p.pop('keep_kappa', False)
    r = mc.run(eta=ETA, **p)
    extra = dict(tag=tag, D=p.get('D'), Tm=p.get('Tm', 0.04), ttt=p.get('ttt', 0.1),
                 t_prep=p.get('t_prep', 0.02), speeds=list(p.get('speeds_kmh', ())),
                 policy=p.get('policy', 'fixed'), shadow=p.get('shadow', 'none'),
                 first_entry=p.get('first_entry', False), keep_kappa=keep,
                 l3=p.get('l3', True), seed=p.get('seed', 0))
    return summarize(r, extra)


def main(stage):
    jobs = []
    if stage == 'val':
        for D, v, Tm, ttt in itertools.product([25, 50, 100], [30, 60, 120],
                                               [0.04, 0.2], [0.1, 0.256]):
            jobs.append(dict(tag='val_det', D=D, speeds_kmh=(v,), Tm=Tm, ttt=ttt,
                             n_ue=300, T=30, shadow='none', turn_dist=None))
        for D, v, Tm in itertools.product([25, 50, 100], [30, 60, 120], [0.04, 0.2]):
            for fe in (True, False):
                jobs.append(dict(tag='val_full', D=D, speeds_kmh=(v,), Tm=Tm, ttt=0.1,
                                 n_ue=300, T=30, shadow='full', sigma_l1=1.0,
                                 first_entry=fe))
    elif stage == 'kappa':
        for D in [25, 35, 50, 75, 100]:
            for sh in ('none', 'full'):
                jobs.append(dict(tag='kappa', D=D, speeds_kmh=(60,), Tm=0.04, ttt=0.1,
                                 n_ue=400, T=30, shadow=sh, sigma_l1=1.0 if sh == 'full' else 0.0,
                                 first_entry=True, keep_kappa=True, seed=11))
    elif stage == 'scheme':
        kq = json.load(open('kappa_q.json'))
        mix = (3, 30, 60, 120)
        for D in [25, 50, 100]:
            base = dict(D=D, speeds_kmh=mix, n_ue=400, T=40, shadow='full',
                        sigma_l1=1.0, Tm=0.04)
            for seed in SEEDS:
                for ttt in (0.04, 0.1, 0.256, 0.512):
                    jobs.append(dict(base, tag='fixed%d' % int(ttt * 1e3), ttt=ttt, seed=seed))
                jobs.append(dict(base, tag='mse', policy='mse', seed=seed))
                k50 = kq[str(D)]['50']; k25 = kq[str(D)]['25']
                jobs.append(dict(base, tag='ga', policy='ga', kappa_q=k50, seed=seed))
                jobs.append(dict(base, tag='ga_q25', policy='ga', kappa_q=k25, seed=seed))
                jobs.append(dict(base, tag='ga_noisy', policy='ga', kappa_q=k50,
                                 speed_err=0.3, seed=seed))
                jobs.append(dict(base, tag='ga_cho', policy='ga', kappa_q=k50,
                                 t_prep=0.0, seed=seed))
                jobs.append(dict(base, tag='ltm_like', policy='ga', kappa_q=k50, Tm=0.02,
                                 t_prep=0.0, l3=False, t_exec=0.025, seed=seed))
            # 3GPP-compliant measurement period
            for seed in SEEDS:
                b2 = dict(base, Tm=0.2)
                jobs.append(dict(b2, tag='fixed100_Tm200', ttt=0.1, seed=seed))
                jobs.append(dict(b2, tag='ga_Tm200', policy='ga', kappa_q=kq[str(D)]['50'], seed=seed))
    if stage == 'sg':
        mix = (3, 30, 60, 120)
        for D in [25, 50, 100]:
            base = dict(D=D, speeds_kmh=mix, n_ue=400, T=40, shadow='full', sigma_l1=1.0, Tm=0.04)
            for seed in SEEDS:
                for ttt in (0.04, 0.1, 0.256):
                    jobs.append(dict(base, tag='fixed%d_nr' % int(ttt * 1e3), ttt=ttt, a3_noreturn=True, seed=seed))
                jobs.append(dict(base, tag='og', ttt=0.256, guard=True, g_thr=-5.5, g_hys=1.0, g_noreturn=True, seed=seed))
                jobs.append(dict(base, tag='og_adapt', ttt=0.256, guard=True, g_adapt=True, g_slope=0.8, g_hys=1.0, g_noreturn=True, seed=seed))
                jobs.append(dict(base, tag='og_nonr', ttt=0.256, guard=True, g_thr=-5.5, g_hys=1.0, g_noreturn=False, seed=seed))
                jobs.append(dict(base, tag='og_Tm200', Tm=0.2, ttt=0.256, guard=True, g_thr=-5.5, g_hys=1.0, g_noreturn=True, seed=seed))
                jobs.append(dict(base, tag='og_cho', ttt=0.256, guard=True, g_thr=-5.5, g_hys=1.0, g_noreturn=True, t_prep=0.0, seed=seed))
                jobs.append(dict(base, tag='og_ltm', ttt=0.256, guard=True, g_thr=-5.5, g_hys=1.0, g_noreturn=True, t_prep=0.0, t_exec=0.025, seed=seed))
                jobs.append(dict(base, tag='ltm_fixed100', ttt=0.1, Tm=0.02, l3=False, t_prep=0.0, t_exec=0.025, seed=seed))
    if stage == 'pp3':
        kq = json.load(open('kappa_q.json')); mix = (3, 30, 60, 120)
        for D in [25, 50, 100]:
            for seed in SEEDS:
                base = dict(D=D, speeds_kmh=mix, n_ue=400, T=40, shadow='full', sigma_l1=1.0, Tm=0.04, seed=seed)
                for ttt in (0.04, 0.1, 0.256):
                    jobs.append(dict(base, tag='fixed%d' % int(ttt * 1e3), ttt=ttt))
                jobs.append(dict(base, tag='fixed100_nr', ttt=0.1, a3_noreturn=True))
                jobs.append(dict(base, tag='mse', policy='mse'))
                jobs.append(dict(base, tag='ga', policy='ga', kappa_q=kq[str(D)]['50']))
                jobs.append(dict(base, tag='og', ttt=0.256, guard=True, g_thr=-5.5, g_hys=1.0, g_noreturn=True))
                jobs.append(dict(base, tag='og_rb3', ttt=0.256, guard=True, g_thr=-5.5, g_hys=1.0, g_noreturn=True, pp_win_block=3.0))
    if stage == 'sens':
        # robustness of OG-A3 to guard-threshold bias and SINR-estimation noise/averaging
        mix = (3, 30, 60, 120)
        og = dict(ttt=0.256, guard=True, g_hys=1.0, g_noreturn=True)
        for D in [25, 50, 100]:
            for seed in SEEDS:
                base = dict(D=D, speeds_kmh=mix, n_ue=400, T=40, shadow='full', sigma_l1=1.0, Tm=0.04, seed=seed)
                for thr in (-6.5, -5.5, -4.5):
                    jobs.append(dict(base, **og, tag='og_thr%+.1f' % thr, g_thr=thr))
                jobs.append(dict(base, **og, tag='og_sig3', g_thr=-5.5, g_sigma=3.0))
                jobs.append(dict(base, **og, tag='og_sig3_avg50', g_thr=-5.5, g_sigma=3.0, g_tavg=0.05))
    if stage == 'lite':
        # OG-A3-lite: guard averaging window 50 ms and guard offset 2 dB
        mix = (3, 30, 60, 120)
        for D in [25, 50, 100]:
            for seed in SEEDS:
                base = dict(D=D, speeds_kmh=mix, n_ue=400, T=40, shadow='full', sigma_l1=1.0, Tm=0.04, seed=seed)
                jobs.append(dict(base, tag='og_lite', ttt=0.256, guard=True, g_thr=-5.5, g_hys=2.0, g_tavg=0.05, g_noreturn=True))
    t0 = time.time()
    with Pool(int(os.environ.get('PROCS', '2'))) as pool:
        res = pool.map(job, jobs)
    json.dump(res, open('res_%s.json' % stage, 'w'))
    print(stage, len(res), 'runs in', round(time.time() - t0), 's')


if __name__ == '__main__':
    main(sys.argv[1])
