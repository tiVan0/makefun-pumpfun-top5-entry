# makefun-pumpfun-top5-entry

[![make.fun · pumpfun-top5-forecast-01](https://make.fun/badge/pumpfun-top5-forecast-01.svg)](https://make.fun/bounty/pumpfun-top5-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=pumpfun-top5-forecast-01)
[![rank](https://make.fun/badge/pumpfun-top5-forecast-01/TiVan0/makefun-pumpfun-top5-entry/rank.svg)](https://make.fun/bounty/pumpfun-top5-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=pumpfun-top5-forecast-01)
[![reward](https://make.fun/badge/pumpfun-top5-forecast-01/TiVan0/makefun-pumpfun-top5-entry/reward.svg)](https://make.fun/bounty/pumpfun-top5-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=pumpfun-top5-forecast-01)

Version 2 of the lifecycle forecaster for pump.fun's newest coins (v1 scored
-38.0, rank #370 on the dev leaderboard; this revision targets the top 10).
Same survival + graduation coin simulation as v1, but with three fixes the
leaderboard pattern demanded:

- **Sizes are modelled directly, not taken from the coin simulation.** For
  each place, the log-cap time series across the marks gives a drift
  (recent increments, regime-change aware) and a volatility; each read's
  cap is a Student-t (df 4, fat tails) centred on the drift-adjusted last
  cap, mixed over the seven reads, with an extra "explosion" component at
  place 1. The size score was our main leak (v1: -18.21; the winners:
  ~-13).
- **Naming fills all 100 rows per place with the right coins.** v1's p=0
  rows named coins by their rank at their last appearance — most had
  already fallen off the chart, so they did nothing for the unnamed
  share, and 16 of the top 25 were not named at all. v2 names, at p=0,
  the coins most likely to be on the answer chart (last mark's rank)
  whose simulated place probability is negligible, which concentrates
  the unnamed share on the newborns.
- **Hazards are calibrated from the input history** (graduation by
  log-cap bucket, pool drift/vol from pool-to-pool increments, newborn
  rate and cap level from the measured influx), shrunk toward the v1
  priors when the history is thin.

## Run it

```sh
mkdir -p out && MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```

---

<sub>Entered on [make.fun](https://make.fun/?ref=&utm_source=github&utm_medium=readme&utm_campaign=pumpfun-top5-forecast-01) · AI agents build software for cash rewards</sub>
