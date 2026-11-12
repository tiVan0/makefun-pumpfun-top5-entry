#!/usr/bin/env python3
"""Lifecycle forecaster v2 for make.fun pumpfun-top5-forecast-01.

What changed vs v1 (which scored -38.0, rank #370, on the dev phase;
order -19.75, size -18.21; top-10 cutoff -32.8):

  * SIZES are no longer taken from the coin simulation.  The size score
    (-18.21, our main leak: the winners score ~-13) is the distribution
    of the CAP AT EACH PLACE, so it is modelled directly: for each place,
    the log-cap time series across the marks gives a drift (recent
    increments, shrunk toward a prior) and a volatility (robust, shrunk);
    each read's cap is a Student-t (df 4, fat tails for jumps) centred on
    the drift-adjusted last cap, and the seven reads are mixed.  Place 1
    gets an extra "explosion" mixture component (a new coin exploding to
    #1).  This is the level the winners model ("Student-t + place-1
    explosions"); the coin simulation cannot see it, because a place's
    cap is the max of ~100 coins and is far more stable than any one coin.
  * NAMING fills all 100 rows per place, but with the right coins.  v1's
    p=0 rows named coins by their rank at their LAST APPEARANCE, which
    is often an old mark: most of them had already fallen off the chart,
    so they did nothing for the unnamed share (and 16 of the top 25 were
    not named at all).  v2 names, at p=0, the coins most likely to be ON
    the answer chart (last mark's rank, then slope) whose simulated place
    probability is negligible; the unnamed share is then divided among
    far fewer chart coins, which is exactly what scores when a newborn
    takes a place.  The rules themselves describe this trick.
  * HAZARDS are calibrated from the input history instead of fixed
    constants: graduation hazard by log-cap bucket (shrunk toward the
    v1 priors), pool drift/vol from pool->pool increments, and the
    newborn field (rate and cap level) from the measured influx of new
    coins into the top 5 / top 100.  With a 50-minute practice history
    the measurements are thin, so shrinkage dominates; with the final's
    11-hour history the measurements take over.  The code never hardcodes
    a horizon: every read time comes from task.json.

Standard library only. No network. No model calls. Deterministic
(fixed seed): two runs are byte-identical.
Reads $MAKEFUN_INPUT_DIR (else /input), writes $MAKEFUN_OUTPUT_DIR (else
/output).
"""

import json
import math
import os
import random
from datetime import datetime, timezone

SEED = 987654321
N_PATHS = 3000

# log10 of a market cap in SOL at which a dead curve coin sits.
FLOOR_LOG = 1.45
MIN_LOG = 1.0
MAX_LOG = 6.5

# --- priors (v1's constants); measurements shrink toward these ---
POOL_DRIFT0_PRIOR = 0.45
POOL_DRIFT_TAU = 1.5
POOL_VOL0_PRIOR = 0.45
CURVE_VOL0 = 0.25
VOL_TAU = 3.0
CURVE_DEATH_HAZARD = 0.10
POOL_DEATH_BASE_PRIOR = 0.04
POOL_DEATH_AGE = 0.03
POOL_DEATH_AGE_CAP = 4.0
GRAD_HAZARD_PRIOR = ((2.2, 1.2), (2.0, 0.4), (1.8, 0.08),
                     (float("-inf"), 0.01))
GRAD_HAZARD_PRIOR_EXPOSURE_H = 2.0
GRAD_HAZARD_COMPLETE = 2.0
GRAD_REPRICE_MEAN = -0.30
GRAD_REPRICE_SD = 0.60
NEWBORN_LAMBDA0_PRIOR = 2.0
NEWBORN_LAMBDA_PER_HOUR = 0.15
NEWBORN_SD = 0.40
# Place-cap (sizes) model priors.
SIZE_DRIFT_PRIOR = 0.55      # log10/h, shrunk toward this
SIZE_DRIFT_TAU_H = 1.5       # drift saturates over this horizon
SIZE_VOL_PRIOR = 0.35        # log10 per sqrt(h)
SIZE_DF = 4.0                # Student-t degrees of freedom
SIZE_EXPLODE_W1 = 0.18       # place-1 explosion mixture weight
SIZE_EXPLODE_D1 = 0.45       # place-1 explosion centre offset (log10)
SIZE_EXPLODE_S1 = 2.0        # place-1 explosion scale multiplier
SIZE_EXPLODE_WK = 0.08       # places 2-5 explosion weight
SIZE_EXPLODE_DK = 0.30
SIZE_EXPLODE_SK = 1.5
MIN_P = 0.003                # name at simulated p at/above this
MIN_SIZE_P = 0.002           # name a size bin at/above this


def parse_time(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def logcap_sol(lamports):
    return math.log10(int(lamports)) - 9.0


def floor6(x):
    return math.floor(x * 1000000.0 + 1e-9) / 1000000.0


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return 0.0
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def iqr_sd(xs):
    """Robust SD from the interquartile range; 0 when degenerate."""
    if len(xs) < 4:
        return 0.0
    xs = sorted(xs)
    q1 = xs[len(xs) // 4]
    q3 = xs[(3 * len(xs)) // 4]
    return max((q3 - q1) / 1.349, 1e-6)


# --- Student-t CDF via the regularized incomplete beta function ---
# (Numerical Recipes betacf, Lentz's continued fraction).  Tested in
# __main__ self-test against known t(4) values when run directly.

def _betacf(a, b, x):
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, 201):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 3e-7:
            break
    return h


def _betainc(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lb = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
          + a * math.log(x) + b * math.log(1.0 - x))
    bt = math.exp(lb)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def student_t_cdf(x, df, mu=0.0, scale=1.0):
    """P(T <= x) for T ~ Student-t(df) centred at mu with scale."""
    if scale <= 0:
        return 1.0 if x >= mu else 0.0
    z = (x - mu) / scale
    xb = df / (df + z * z)
    ib = _betainc(df / 2.0, 0.5, xb)
    if z >= 0:
        return 1.0 - 0.5 * ib
    return 0.5 * ib


def grad_hazard(logcap, complete, table):
    if complete:
        return GRAD_HAZARD_COMPLETE
    for threshold, hazard in table:
        if logcap >= threshold:
            return hazard
    return 0.01


def pool_drift(age_h, horizon_h, slope, drift0):
    d0 = 0.5 * slope + 0.5 * (drift0 * math.exp(-age_h / 1.0) + 0.05)
    return d0 * POOL_DRIFT_TAU * (1.0 - math.exp(-horizon_h / POOL_DRIFT_TAU))


def vol_sd(v0, age_h, horizon_h):
    a0 = math.exp(-2.0 * age_h / VOL_TAU)
    a1 = math.exp(-2.0 * (age_h + horizon_h) / VOL_TAU)
    var = v0 * v0 * (VOL_TAU / 2.0) * (a0 - a1)
    return math.sqrt(max(0.0, var))


def pool_death_prob(age_h, horizon_h, base):
    def lam(a):
        return base + POOL_DEATH_AGE * min(max(0.0, a - 1.0),
                                           POOL_DEATH_AGE_CAP)
    avg = 0.5 * (lam(age_h) + lam(age_h + horizon_h))
    return 1.0 - math.exp(-avg * horizon_h)


def exp_time(rng, hazard):
    if hazard <= 0.0:
        return float("inf")
    u = rng.random()
    if u <= 0.0:
        u = 1e-12
    return -math.log(1.0 - u) / hazard


def poisson(rng, lam):
    limit = math.exp(-lam)
    k = 0
    p = 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


class Coin(object):
    __slots__ = ("mint", "logcap", "venue", "age_h", "slope", "complete",
                 "last_rank")

    def __init__(self, mint, logcap, venue, age_h, slope, complete,
                 last_rank):
        self.mint = mint
        self.logcap = logcap
        self.venue = venue
        self.age_h = age_h
        self.slope = slope
        self.complete = complete
        self.last_rank = last_rank


def build_coins(marks, cut):
    """One Coin per mint seen in the history, with a shrunk slope."""
    traj = {}
    for m in marks:
        t = parse_time(m["at"])
        for r in m["rows"]:
            mint = r["mint"]
            entry = traj.get(mint)
            if entry is None:
                entry = {"obs": [], "created": r["created_at"],
                         "last": None}
                traj[mint] = entry
            entry["obs"].append((t, logcap_sol(r["market_cap_lamports"]),
                                 r["venue"]))
            entry["last"] = r
    coins = []
    for mint, entry in traj.items():
        obs = entry["obs"]
        last = entry["last"]
        created = parse_time(entry["created"])
        age_h = max(0.0, (cut - created).total_seconds() / 3600.0)
        slope = 0.0
        pts = obs[-6:]
        if len(pts) >= 2:
            t0 = pts[0][0]
            xs = [(p[0] - t0).total_seconds() / 3600.0 for p in pts]
            ys = [p[1] for p in pts]
            n = len(pts)
            mx = sum(xs) / n
            my = sum(ys) / n
            den = sum((x - mx) ** 2 for x in xs)
            if den > 1e-9:
                raw = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
                raw = max(-0.6, min(0.6, raw))
                slope = raw * (0.6 if n >= 4 else (0.4 if n == 3 else 0.25))
        venue = last["venue"]
        if venue == "migrating":
            venue = "pool"
        coins.append(Coin(mint, logcap_sol(last["market_cap_lamports"]),
                          venue, age_h, slope, bool(last.get("complete")),
                          int(last["rank"])))
    coins.sort(key=lambda c: (-c.logcap, c.mint))
    return coins, traj


def calibrate(marks, coins, traj):
    """Estimate hazards / drift / vol / newborn field from the history.

    Every estimate is shrunk toward the v1 prior, so a 50-minute
    practice history (thin data) stays near the priors while the final's
    11-hour history lets the measurements dominate.
    Returns a dict of calibrated parameters.
    """
    times = [parse_time(m["at"]) for m in marks]
    # --- graduation hazard by log-cap bucket ---
    buckets = [(2.2, 1.2), (2.0, 0.4), (1.8, 0.08), (float("-inf"), 0.01)]
    exposure = [0.0] * 4
    grads = [0] * 4

    def bucket_of(lc):
        for i, (thr, _r) in enumerate(buckets):
            if lc >= thr:
                return i
        return 3

    for i in range(len(marks) - 1):
        dt = (times[i + 1] - times[i]).total_seconds() / 3600.0
        if dt <= 0:
            continue
        cur = {r["mint"]: r for r in marks[i]["rows"]}
        nxt = {r["mint"]: r for r in marks[i + 1]["rows"]}
        for mint, r in cur.items():
            if r["venue"] != "curve":
                continue
            b = bucket_of(logcap_sol(r["market_cap_lamports"]))
            exposure[b] += dt
            n = nxt.get(mint)
            if n is not None and n["venue"] in ("pool", "migrating"):
                grads[b] += 1
    grad_table = []
    for i, (thr, prior) in enumerate(buckets):
        rate = (grads[i] + GRAD_HAZARD_PRIOR_EXPOSURE_H * prior) / (
            exposure[i] + GRAD_HAZARD_PRIOR_EXPOSURE_H)
        grad_table.append((thr, min(3.0, max(0.005, rate))))
    # --- pool drift / vol from pool->pool increments ---
    pinc = []
    for i in range(len(marks) - 1):
        dt = (times[i + 1] - times[i]).total_seconds() / 3600.0
        if dt <= 0:
            continue
        cur = {r["mint"]: r for r in marks[i]["rows"]}
        for r in marks[i + 1]["rows"]:
            c = cur.get(r["mint"])
            if c is not None and c["venue"] == "pool" \
                    and r["venue"] == "pool":
                pinc.append((logcap_sol(r["market_cap_lamports"])
                             - logcap_sol(c["market_cap_lamports"])) / dt)
    if pinc:
        m = sum(pinc) / len(pinc)
        # prior weight: 20 observations of drift 0.45/h
        drift0 = (len(pinc) * m + 20.0 * POOL_DRIFT0_PRIOR) / (len(pinc) + 20.0)
        drift0 = max(0.0, min(0.9, drift0))
        # vol per hour from robust sd of the per-hour increments
        sd = iqr_sd(pinc)
        if sd <= 0:
            sd = (sum((x - m) ** 2 for x in pinc) / len(pinc)) ** 0.5
        vol0 = (len(pinc) * sd + 8.0 * POOL_VOL0_PRIOR) / (len(pinc) + 8.0)
        vol0 = max(0.15, min(0.9, vol0))
    else:
        drift0, vol0 = POOL_DRIFT0_PRIOR, POOL_VOL0_PRIOR
    # --- pool death base from pool crashes (increment < -0.3 in a step) ---
    crash = 0
    pool_steps = 0
    for i in range(len(marks) - 1):
        cur = {r["mint"]: r for r in marks[i]["rows"]}
        for r in marks[i + 1]["rows"]:
            c = cur.get(r["mint"])
            if c is not None and c["venue"] == "pool" \
                    and r["venue"] == "pool":
                pool_steps += 1
                if logcap_sol(r["market_cap_lamports"]) - logcap_sol(
                        c["market_cap_lamports"]) < -0.3:
                    crash += 1
    if pool_steps:
        crash_rate = crash / pool_steps * 6.0  # per hour
        pool_death_base = (pool_steps * crash_rate + 10.0
                           * POOL_DEATH_BASE_PRIOR) / (pool_steps + 10.0)
        pool_death_base = max(0.01, min(0.25, pool_death_base))
    else:
        pool_death_base = POOL_DEATH_BASE_PRIOR
    # --- newborn field: new coins entering the top 5 / top 100 ---
    # The first mark is the initial state, not a measurement: every coin
    # in it is "new" by construction, so it must not inflate the rate.
    seen = set()
    new_top5 = []
    new_top100 = []
    for idx, m in enumerate(marks):
        rows = sorted(m["rows"], key=lambda r: r["rank"])
        t5 = {r["mint"] for r in rows[:5]}
        t100 = {r["mint"] for r in rows[:100]}
        if idx == 0:
            new_top5.append(0)
            new_top100.append(0)
        else:
            new_top5.append(len(t5 - seen))
            new_top100.append(len(t100 - seen))
        seen |= {r["mint"] for r in rows}
    # recent rate (last two marks) shrunk toward the overall rate
    measured = new_top5[1:]
    if len(measured) >= 2:
        recent = sum(measured[-2:]) / 2.0
        overall = sum(measured) / len(measured)
        rate5 = 0.6 * recent + 0.4 * overall
    elif measured:
        rate5 = float(measured[-1])
    else:
        rate5 = 1.0
    # lambda0: Poisson mean of newborns drawn per read.  A newborn takes
    # a top-5 place only if its cap beats the 5th incumbent (~1 in 3 of
    # the draws at the top-5 cap level), so scale the measured rate up.
    lambda0 = min(6.0, max(1.0, rate5 * 2.0))
    lambda0 = 0.5 * lambda0 + 0.5 * NEWBORN_LAMBDA0_PRIOR
    # newborn cap level: median of the last mark's top-5 caps
    last_rows = sorted(marks[-1]["rows"], key=lambda r: r["rank"])
    top5_caps = sorted(logcap_sol(r["market_cap_lamports"])
                       for r in last_rows[:5])
    newborn_mu = top5_caps[len(top5_caps) // 2] if top5_caps else 3.0
    # --- place-cap drift / vol per place (for the size model) ---
    # Drift: the most recent increments, extended backwards only while
    # they stay consistent (a structural break -- the early graduation
    # wave, where place caps jump +0.5..+0.9 in one step -- must not
    # contaminate the settled drift).  Vol: pooled residuals of those
    # windows around each place's drift, shrunk toward the prior.
    place_stats = []
    windows = []
    for k in range(5):
        caps = []
        for m in marks:
            rows = sorted(m["rows"], key=lambda r: r["rank"])
            if len(rows) > k:
                caps.append(logcap_sol(rows[k]["market_cap_lamports"]))
        if len(caps) < 2:
            place_stats.append((caps[-1] if caps else 3.0,
                                SIZE_DRIFT_PRIOR, SIZE_VOL_PRIOR, 0.30))
            windows.append([])
            continue
        incs = [(b - a) * 6.0 for a, b in zip(caps, caps[1:])]
        window = [incs[-1]]
        for x in reversed(incs[:-1]):
            if len(window) >= 12:
                break
            if abs(x - median(window)) <= 0.35:
                window.append(x)
            else:
                break
        windows.append(window)
        raw_drift = median(window)
        n = len(window)
        drift = (n * raw_drift + 2.0 * SIZE_DRIFT_PRIOR) / (n + 2.0)
        drift = max(0.0, min(0.6, drift))
        # standard error of the drift estimate (parameter uncertainty,
        # not process noise): the window's dispersion / sqrt(n).
        if n >= 2:
            mw = sum(window) / n
            dse = (sum((x - mw) ** 2 for x in window) / n) ** 0.5
            dse = dse / math.sqrt(n)
        else:
            dse = 0.30
        dse = max(0.10, min(0.30, dse))
        place_stats.append((caps[-1], drift, None, dse))
    # pooled vol from the window residuals, winsorized: a regime-change
    # jump must not inflate the per-read noise (it is a level shift,
    # not recurring dispersion).
    res = []
    for k, window in enumerate(windows):
        d = place_stats[k][1]
        res.extend(max(-0.5, min(0.5, x - d)) for x in window)
    if len(res) >= 4:
        mres = sum(res) / len(res)
        pooled = (sum((x - mres) ** 2 for x in res) / len(res)) ** 0.5
        vol = 0.5 * pooled + 0.5 * 0.25
    else:
        vol = 0.25
    vol = max(0.12, min(0.35, vol))
    place_stats = [(c, d, vol, se) for (c, d, _v, se) in place_stats]
    return {
        "grad_table": grad_table,
        "grad_exposure": exposure,
        "grad_counts": grads,
        "pool_drift0": drift0,
        "pool_vol0": vol0,
        "pool_death_base": pool_death_base,
        "pool_increments": len(pinc),
        "newborn_lambda0": lambda0,
        "newborn_mu": newborn_mu,
        "new_top5_per_mark": new_top5,
        "place_stats": place_stats,
    }


def simulate_coin(rng, coin, horizon_h, par):
    if coin.venue == "pool":
        if rng.random() < pool_death_prob(coin.age_h, horizon_h,
                                          par["pool_death_base"]):
            return FLOOR_LOG + rng.gauss(0.0, 0.05)
        drift = pool_drift(coin.age_h, horizon_h, coin.slope,
                           par["pool_drift0"])
        cap = coin.logcap + drift + vol_sd(
            par["pool_vol0"], coin.age_h, horizon_h) * rng.gauss(0.0, 1.0)
        return max(MIN_LOG, min(MAX_LOG, cap))
    lam_d = CURVE_DEATH_HAZARD
    lam_g = grad_hazard(coin.logcap, coin.complete, par["grad_table"])
    t_die = exp_time(rng, lam_d)
    t_grad = exp_time(rng, lam_g)
    if t_die < t_grad and t_die < horizon_h:
        return FLOOR_LOG + rng.gauss(0.0, 0.05)
    if t_grad < horizon_h:
        pre = vol_sd(CURVE_VOL0, coin.age_h, t_grad) * rng.gauss(0.0, 1.0)
        reprice = GRAD_REPRICE_MEAN + GRAD_REPRICE_SD * rng.gauss(0.0, 1.0)
        rest = horizon_h - t_grad
        age_at_grad = coin.age_h + t_grad
        if rng.random() < pool_death_prob(age_at_grad, rest,
                                          par["pool_death_base"]):
            return FLOOR_LOG + rng.gauss(0.0, 0.05)
        drift = pool_drift(age_at_grad, rest, 0.0, par["pool_drift0"])
        cap = coin.logcap + pre + reprice + drift + vol_sd(
            par["pool_vol0"], age_at_grad, rest) * rng.gauss(0.0, 1.0)
        return max(MIN_LOG, min(MAX_LOG, cap))
    cap = coin.logcap + vol_sd(CURVE_VOL0, coin.age_h,
                               horizon_h) * rng.gauss(0.0, 1.0)
    return max(MIN_LOG, min(MAX_LOG, cap))


def bin_of(logcap, edges):
    if logcap < edges[0]:
        return 0
    if logcap >= edges[-1]:
        return len(edges)
    lo, hi = 0, len(edges)
    while lo < hi:
        mid = (lo + hi) // 2
        if edges[mid] <= logcap:
            lo = mid + 1
        else:
            hi = mid
    return lo


def size_distribution(place_stats, k, horizons, edges):
    """P(bin) at place k+1: Student-t mixture over the read horizons.

    Each read's cap is Student-t(df) centred on the drift-adjusted last
    cap; place 1 mixes in an explosion component.  Returns a list of
    (bin, p) with p >= MIN_SIZE_P.
    """
    last_cap, drift, vol, drift_se = place_stats[k]
    n_bins = len(edges) + 1
    probs = [0.0] * n_bins
    w_exp = SIZE_EXPLODE_W1 if k == 0 else SIZE_EXPLODE_WK
    d_exp = SIZE_EXPLODE_D1 if k == 0 else SIZE_EXPLODE_DK
    s_exp = SIZE_EXPLODE_S1 if k == 0 else SIZE_EXPLODE_SK
    for h in horizons:
        mu = last_cap + drift * SIZE_DRIFT_TAU_H * (
            1.0 - math.exp(-h / SIZE_DRIFT_TAU_H))
        # process noise grows like sqrt(h); drift-estimate error grows
        # like h (it shifts the whole centre).  Both widen the read.
        scale = math.sqrt((vol * math.sqrt(max(h, 0.05))) ** 2
                          + (drift_se * h) ** 2)
        scale = max(0.08, min(0.55, scale))
        for b in range(n_bins):
            lo = float("-inf") if b == 0 else edges[b - 1]
            hi = float("inf") if b == len(edges) else edges[b]
            p_main = (student_t_cdf(hi, SIZE_DF, mu, scale)
                      - student_t_cdf(lo, SIZE_DF, mu, scale))
            p_exp = (student_t_cdf(hi, SIZE_DF, mu + d_exp, scale * s_exp)
                     - student_t_cdf(lo, SIZE_DF, mu + d_exp, scale * s_exp))
            probs[b] += ((1.0 - w_exp) * p_main + w_exp * p_exp) / len(horizons)
    out = []
    for b in range(n_bins):
        if probs[b] >= MIN_SIZE_P:
            out.append((b, floor6(probs[b])))
    return out


def main():
    in_dir = os.environ.get("MAKEFUN_INPUT_DIR", "/input")
    out_dir = os.environ.get("MAKEFUN_OUTPUT_DIR", "/output")
    with open(os.path.join(in_dir, "task.json"), "r", encoding="utf-8") as f:
        task = json.load(f)
    places = int(task["places"])
    edges = [float(e) for e in task["size"]["edges"]]
    reads = [parse_time(s) for s in task["reads"]]
    marks = []
    with open(os.path.join(in_dir, "history.ndjson"), "r",
              encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = json.loads(line)
            if m.get("ok") and isinstance(m.get("rows"), list) and m["rows"]:
                marks.append(m)
    if not marks:
        raise SystemExit("the history holds no ok mark with rows")
    marks.sort(key=lambda m: m["at"])
    cut = parse_time(marks[-1]["at"])
    horizons = []
    for r in reads:
        h = (r - cut).total_seconds() / 3600.0
        horizons.append(max(h, 0.0))
    if not horizons:
        horizons = [0.0]

    coins, traj = build_coins(marks, cut)
    par = calibrate(marks, coins, traj)
    last_rows = sorted(marks[-1]["rows"], key=lambda r: r["rank"])
    last_chart = last_rows
    last_caps = [logcap_sol(r["market_cap_lamports"])
                 for r in last_chart[:3]]
    newborn_mu = sorted(last_caps)[len(last_caps) // 2] if last_caps \
        else par["newborn_mu"]
    # The newborn cap level is the calibrated top-5 median; the +0.45
    # selection offset and the drift tracking are applied per read in
    # the simulation loop below (see there).
    newborn_mu = par["newborn_mu"]

    # Candidates: top 120 by last cap + every pool/migrating coin + every
    # fast grower; capped so the simulation stays fast.
    cand = {}
    for c in coins[:120]:
        cand[c.mint] = c
    for c in coins:
        if c.venue == "pool":
            cand[c.mint] = c
        elif c.slope > 0.2:
            cand[c.mint] = c
    cand_list = sorted(cand.values(), key=lambda c: (-c.logcap, c.mint))
    cand_list = cand_list[:200]
    cand_index = {c.mint: i for i, c in enumerate(cand_list)}

    # The newborn cap level: measured from the practice history, the
    # coins that take a top-5 place after the cut have caps ~+0.2..+0.9
    # (centred +0.45) above the cut's top-5 median -- they are the
    # selected upper tail of the newborn cohort, not typical newborns.
    # The level tracks the market's own growth with the horizon (the
    # place-cap drift), so a newborn at a late read must beat the grown
    # incumbents, not the cut's caps.
    drift_top5 = median([par["place_stats"][k][1] for k in range(places)])
    rng = random.Random(SEED)
    n_reads = len(reads)
    total = N_PATHS * n_reads
    place_counts = [[0] * len(cand_list) for _ in range(places)]

    for _ in range(N_PATHS):
        for h in horizons:
            sims = []
            for i, c in enumerate(cand_list):
                sims.append((simulate_coin(rng, c, h, par), i))
            lam_new = par["newborn_lambda0"] + NEWBORN_LAMBDA_PER_HOUR * h
            # Selection offset: the coins that take a top-5 place after
            # the cut are the selected upper tail of the newborn cohort
            # (measured +0.2..+0.9 above the cut's top-5 median).  Kept
            # constant across horizons.  An age-ramped variant (each
            # newborn's cap ramping with its own age) was tried and
            # REJECTED: it fixed the two short-horizon origins but
            # cost -1.51 overall (cut 20:10 -3.48, cut 20:20 -2.06) —
            # the strong newborn field is what drives the gains at the
            # long horizons the dev (10-70 min) and final (7-13 h)
            # reads live in.  A horizon-ramped offset was likewise
            # rejected (overall +1.96 vs +2.93).
            mu_new = (newborn_mu + 0.45 + drift_top5 * SIZE_DRIFT_TAU_H
                      * (1.0 - math.exp(-h / SIZE_DRIFT_TAU_H)))
            k_new = poisson(rng, lam_new)
            for _ in range(k_new):
                sims.append((mu_new + NEWBORN_SD * rng.gauss(0.0, 1.0), -1))
            sims.sort(key=lambda t: (-t[0], t[1]))
            for k in range(places):
                _cap, idx = sims[k]
                if idx >= 0:
                    place_counts[k][idx] += 1

    # Simulated place probabilities.
    sim_p = [[place_counts[k][i] / total for i in range(len(cand_list))]
             for k in range(places)]
    # Chart-likelihood order for the p=0 fill: last mark's rank (coins
    # absent from the last mark rank last), then slope (risers first).
    last_rank = {}
    for r in last_chart:
        last_rank[r["mint"]] = int(r["rank"])
    for c in coins:
        if c.mint not in last_rank:
            last_rank[c.mint] = 999
    # All coins that could be on the answer chart, ordered by how likely:
    # low last rank first; among equals, higher slope first.
    chart_pool = sorted(coins, key=lambda c: (last_rank[c.mint],
                                              -c.slope, -c.logcap, c.mint))

    probs = []
    named = [set() for _ in range(places)]
    for k in range(places):
        rows_here = 0
        for i, c in enumerate(cand_list):
            p = sim_p[k][i]
            if p >= MIN_P:
                probs.append((c.mint, k + 1, floor6(p)))
                named[k].add(c.mint)
                rows_here += 1
        # p=0 fill: coins most likely on the answer chart whose simulated
        # place probability is negligible.  They take no probability mass
        # (the unnamed share is untouched) but they leave the unnamed
        # share to the few chart coins not named -- the newborns.
        for c in chart_pool:
            if rows_here >= 100:
                break
            if c.mint in named[k]:
                continue
            idx = cand_index.get(c.mint)
            if idx is not None and sim_p[k][idx] >= MIN_P:
                continue  # already named above; unreachable, kept for safety
            if idx is not None and sim_p[k][idx] >= 0.001:
                continue  # a real chance: do not zero it; it stays unnamed
            probs.append((c.mint, k + 1, 0.0))
            named[k].add(c.mint)
            rows_here += 1
    probs.sort(key=lambda t: (t[1], -t[2], t[0]))
    probabilities = [{"mint": m, "place": pl, "p": p}
                     for (m, pl, p) in probs]

    # top5: rank by expected place score.
    score = {}
    for k in range(places):
        for i, c in enumerate(cand_list):
            p = sim_p[k][i]
            if p > 0.0:
                score[c.mint] = score.get(c.mint, 0.0) + p * (
                    places + 1 - (k + 1))
    ranked = sorted(
        score.items(),
        key=lambda kv: (-kv[1],
                        -sim_p[0][cand_index[kv[0]]],
                        kv[0]))
    top5 = [m for m, _ in ranked[:places]]
    if len(top5) < places:
        seen = set(top5)
        for r in last_chart:
            if len(top5) >= places:
                break
            if r["mint"] not in seen:
                top5.append(r["mint"])
                seen.add(r["mint"])

    # Sizes from the place-cap model (not the coin simulation).
    sizes = []
    for k in range(places):
        for (b, p) in size_distribution(par["place_stats"], k, horizons,
                                        edges):
            sizes.append({"place": k + 1, "bin": b, "p": p})
    sizes.sort(key=lambda s: (s["place"], -s["p"], s["bin"]))

    forecast = {"probabilities": probabilities, "top5": top5, "sizes": sizes}
    out_path = os.path.join(out_dir, task["output"])
    os.makedirs(out_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(json.dumps(forecast) + "\n")


if __name__ == "__main__":
    # Self-test the Student-t CDF against known t(4) values.
    checks = [(0.0, 0.5), (1.0, 0.8129), (2.0, 0.9410), (-1.0, 0.1871)]
    ok = True
    for x, want in checks:
        got = student_t_cdf(x, 4.0)
        if abs(got - want) > 0.002:
            print("t-cdf self-test FAIL: x=%s got %.4f want %.4f"
                  % (x, got, want))
            ok = False
    if not ok:
        raise SystemExit("student-t CDF self-test failed")
    main()
