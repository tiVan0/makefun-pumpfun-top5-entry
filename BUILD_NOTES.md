# BUILD NOTES — pumpfun-top5-forecast-01 entry, version 2

## Why v2 exists

v1 scored **-38.0 (rank #370)** on the dev-phase leaderboard
(order -19.75, size -18.21; top-10 cutoff -32.8; median -37.5;
house baseline practice -40.408). Verified against
`~/workspace/makefun/.leaderboard/pumpfun-top5-forecast-01.json`
(645 entries). The winners' pattern, from the leaderboard itself:

| rank | log | order | size | title (approach) |
|---|---|---|---|---|
| #1 crave | -31.7 | -18.18 (r13) | -13.56 (r3) | Scenario simulation with history-scaled newcomer rates |
| #2 chaboss-work | -31.8 | -17.61 | -14.22 | MC of new-coin caps: day-calibrated crash hazard + newcomer field |
| #3 wahyu-friyadi-1 | -32.0 | -19.03 (r87) | -13.01 (r2) | Momentum Monte-Carlo v5 |
| #6 coral-lark-73 | -32.4 | -19.47 (r200) | -12.95 (r1) | MC v9 Student-t + place-1 explosions |

Two facts drive this revision:

1. **Size is the leak.** The winners score ~-13 on size; we scored
   -18.21 — a ~5-point gap, nearly the whole distance to the top 10
   (5.2 points). Order is only ~1–2 points behind the best.
   The winners model size with a fat-tailed (Student-t) distribution
   and a place-1 explosion component, not with a coin simulation.
2. **Order and size are scored independently** (RULES, Scoring), so
   each is modelled at its own natural level: the coin simulation
   for *which coin* takes each place, the place-cap process for
   *what cap* each place holds.

## What v2 changes (three fixes)

### 1. Sizes from a direct place-cap model (the main fix)

For each place k = 1..5, the log-cap series across the marks is
measured directly:

- **Drift** from the most recent increments, extended backwards
  only while they stay within 0.35 log10/h of the running median
  (a structural break — the early graduation wave, where place
  caps jump +0.5..+0.9 in one 10-minute step — must not
  contaminate the settled drift). Shrunk toward a +0.55/h prior
  (weight 2 increments), capped at 0.6/h. Measured on the
  practice: place 1 drift +0.79/h overall, +0.91/h from the last
  two increments; places 2–5 +0.39..+0.56/h recent (their
  all-history means of +1.8..+2.2/h are the graduation wave, not
  a sustainable drift).
- **Volatility**: pooled residuals of those windows around each
  place's drift, winsorized to ±0.5 (a regime jump is a level
  shift, not recurring noise), shrunk 50% toward 0.25, floored
  0.12, capped 0.35.
- **Per read**, the place-k cap is Student-t (df 4) centred on
  `last_cap + drift * 1.5 * (1 - exp(-h/1.5))` (drift saturates;
  a market cannot compound +0.9/h for 13 hours), scale
  `sqrt((vol*sqrt(h))^2 + (drift_se*h)^2)` (process noise plus
  drift-estimate error), floored 0.08, capped 0.55. The seven
  reads are mixed with equal weights.
- **Place 1** mixes in an explosion component (weight 0.18,
  centre +0.45, scale 2x); places 2–5 a smaller one (0.08, +0.30,
  1.5x). This is the winners' "place-1 explosions": a new coin
  exploding to #1.
- The Student-t CDF is implemented in pure stdlib via the
  regularized incomplete beta function (continued fraction);
  self-tested at import against known t(4) values
  (CDF(0)=0.5, CDF(1)=0.8129, CDF(2)=0.9410).

### 2. Naming fills all 100 rows per place with the right coins

v1's p=0 rows named coins by their rank at their **last
appearance** — often an old mark. Measured on the v1 forecast:
at place 1, 91 p=0 rows, but their median last-mark rank was
"not in the last mark at all" (rank 999), and **16 of the top 25
were not named at all**. Those rows did nothing for the unnamed
share (the denominator counts chart coins not named), so a
newborn taking a place scored from a share divided by ~80 coins
instead of ~35.

v2 names, at p=0, the coins most likely to be **on the answer
chart** (last mark's rank ascending, then slope descending)
whose simulated place probability is < 0.001, up to 100 rows
per place. A p=0 row takes no probability mass (the unnamed
share is untouched) but removes the coin from the denominator —
exactly the trick the rules describe ("naming at p 0 the coins
you are sure will not hold a place ... leaves that share to the
rest of the chart's coins, the coins not yet born among them").
Coins with simulated p in [0.001, 0.003) stay unnamed rather
than being zeroed (a real, if small, chance must not score the
floor). MIN_P (name at simulated p) is 0.003, as in v1.

Measured retention behind the choice (practice marks): top-25
coins stay in the top 100 for 23/25 over 10 min and 19/25 over
20 min; flat curve coins ranked 25–90 keep only 35/64 over
10 min and 14/66 over 20 min — v1's p=0 set was the wrong
population.

### 3. Hazards calibrated from the input history

Every estimate is shrunk toward the v1 prior, so a 50-minute
practice history stays near the priors while the final's
11-hour history lets the measurements dominate:

- **Graduation hazard by log-cap bucket** (≥2.2, 2.0–2.2,
  1.8–2.0, <1.8): graduations per hour of curve-coin exposure,
  shrunk with 2.0 h of prior exposure toward v1's
  (1.2, 0.4, 0.08, 0.01)/h. Practice measured 6 graduations.
- **Pool drift/vol** from pool→pool 10-minute increments
  (practice n=19: mean +0.081, sd 0.135 per step), shrunk
  toward drift 0.45/h (weight 20 obs) and vol 0.45 (weight 8).
- **Pool death base** from the pool crash rate (increments
  < −0.3), shrunk toward 0.04/h.
- **Newborn field**: the rate of new coins entering the top 5
  (not in any earlier mark's rows) — practice 4, 4, 1, 1 per
  mark, recent rate shrunk toward the overall rate; λ0 =
  0.5·min(6, 2·rate) + 0.5·2.0. The first mark is excluded
  (it is the initial state, not a measurement). **Cap level**:
  the coins that take a top-5 place after the cut have caps
  ~+0.2..+0.9 (centred +0.45) above the cut's top-5 median —
  they are the selected upper tail, not typical newborns — so
  μ(h) = top-5 median + 0.45 + drift·1.5·(1−exp(−h/1.5)),
  σ = 0.40. (v1 used the top-3 median with no offset, which
  made the simulation overconfident on the incumbent #1: at
  cut 20:10 it gave the leader p≈0.94 while newborns took
  place 1 at 2 of the next 4 reads.)

Other parameters kept from v1: seed 987654321, lifecycle
(curve death 0.10/h, graduation repricing −0.30 ± 0.60,
pool death 0.04 + 0.03·age, vol decaying with age tau 3 h),
candidates top 120 + all pool + fast growers (cap 200, was
100/150), N_PATHS 3000 (was 2500).

## Backtest (rolling-origin, exact competition scorer)

Harness `~/workspace/makefun/.verify-pumpfun-v2/backtest.py`
implements the rules' log score (floor ln(0.0001), unnamed
shares for order and sizes, history-mint dropping for
probabilities, top-5 agreement). Origins: cut after each
non-empty practice mark except the last; subsequent marks are
the answer reads (10 reads total). The scorer is validated: it
reproduces v1's own reported backtest exactly at cut 20:30
(v1 −25.53, baseline −33.48) and the rules' baseline practice
score structure.

| origin | reads | v1 order | v1 size | v1 log | v2 order | v2 size | v2 log |
|---|---|---|---|---|---|---|---|
| cut 20:10 | 4 | −27.37 | −35.38 | −62.76 | −24.80 | −30.46 | **−55.26** |
| cut 20:20 | 3 | −24.80 | −19.77 | −44.57 | −20.96 | −22.07 | **−43.03** |
| cut 20:30 | 2 | −14.95 | −10.59 | −25.53 | −15.40 | −11.03 | −26.44 |
| cut 20:40 | 1 | −12.48 | −9.99 | −22.48 | −15.76 | −10.26 | −26.01 |
| **mean** | **10** | **−22.63** | **−23.20** | **−45.83** | **−20.86** | **−22.04** | **−42.90** |

**v2 − v1: log +2.93 (order +1.76, size +1.16); v2 − baseline
+5.88.** The gain concentrates at the longest horizons
(cut 20:10, up to 40 min: +7.50), which is the regime the dev
(10–70 min) and final (7–13 h) reads live in. The two
short-horizon origins (cut 20:30/20:40, 10–20 min) favour v1
by −0.91/−3.53 — the second is a single read, and at 10-minute
horizons the newborn field matters least.

Rejected variants (all measured on the same 10-read backtest,
same harness):

- **Newborn offset +0.25** (instead of +0.45): overall +1.59
  (order +1.11, size +0.49); better on the short origins but
  worse where the horizon is longest.
- **Horizon-ramped newborn offset** (0.45·(1−exp(−h/2))):
  overall +1.96 (order +0.79, size +1.16); fixed the short
  origins (−0.19 vs −1.79 on cut 20:30+20:40) but cost −0.97
  overall against the constant offset.
- **Age-ramped newborn caps** (each newborn's cap ramping
  with its own age, uniform over the horizon): overall +1.42
  (order +0.26, size +1.16); fixed the short origins
  (cut 20:40 −23.09 vs −26.01) but cost −1.51 overall
  (cut 20:10 −3.48, cut 20:20 −2.06).
- **Size drift prior 0.40** (instead of 0.55): overall +2.18
  (order +1.70, size +0.49); the 0.55 prior is worth +0.75,
  almost all of it in size.

The constant +0.45 offset + 0.55 drift prior is the final
configuration; it is the best overall and the best at the
longest horizons.

## Format + determinism (full practice input)

- Exactly 3 keys; place sums ≤ 1 (0.70–0.85); ≤100 rows/place
  (100 each); top5 = 5 distinct history mints; size bins 0–81,
  place size sums ≤ 1 (0.97–0.98); 53,778 bytes ≤ 1 MB.
- Two runs byte-identical (fixed seed). Runtime ≈ 17 s on the
  practice input (limit 300 s); the final's 66-mark history
  only adds cheap parsing — the simulation cost is unchanged
  (paths × reads × candidates).

## Known weaknesses (read before trusting it)

- **The backtest cannot see the dev reads.** The dev-phase
  score (−38.0 for v1) is computed on answer reads 21:00–22:00
  on 2026-11-11, which are not in the kit; the longest
  backtest horizon is 40 minutes. The +2.93 backtest gain is
  real at 10–40 min, and the structural fixes (drift-tracking
  sizes, chart-likely p=0 naming, calibrated newborn level)
  are aimed at the 10–70 min dev reads and the 7–13 h final
  reads, but the 7–13 h behaviour is untested by data — it
  rests on the drift saturation, the Student-t tails and the
  explosion component, not on measurements.
- **Long-horizon drift is a guess.** The place-cap drift is
  estimated from at most a few increments in the practice
  (the dev history has 6 marks); if the real drift after the
  cut is much higher or lower than measured, the narrow size
  distributions miss. The drift-se term widens the reads
  with horizon, but it cannot fix a wrong centre.
- **Newborn rate is regime-dependent.** It was 4 per mark
  early in the practice and 1 per mark late; the estimator
  (0.6·recent + 0.4·overall) lags a regime change. If the
  final's newborn rate is much lower than measured, the
  simulation over-weights newborns and dilutes the incumbents.
- **p=0 rows are a calculated risk.** A coin named at p=0
  that takes a place scores the floor (−9.21) there. The
  fill excludes coins with simulated p ≥ 0.001 at that place,
  but a true 0.1%-chance coin can still land. The expected
  cost is small against the denominator gain; the backtest
  showed no p=0 coin taking a place.
- **Only coins in the history's rows can be named.** A coin
  below the top 200 that explodes is covered only through the
  unnamed share (in the practice, one such coin took a top-5
  place).
- **Sizes are pooled over reads.** One distribution per place
  must cover all seven reads; if the cap drifts far across
  the reads, the mixture is wide by necessity and the
  per-read peak probability falls.
