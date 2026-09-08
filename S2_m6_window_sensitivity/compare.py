"""
S2 comparison: summarise the M6 window sweep.

Reads data/m6_window_sweep.csv and reports, for each common cutoff:

  - the domain-wide M6 fraction
  - the climatological M6 share at the season endpoints (Nov 15, Mar 14)
  - the seasonal slope, Nov 15 minus Mar 14, which is the quantity at issue
  - the fraction of particles demoted

Censored is folded into M6 throughout, consistent with F3 and F4/F5.

The September 30 cutoff truncates no cohort, so its row must reproduce the
reported classification exactly: zero demotions, M6 45.8% -> 33.6%. That is
the regression test, and it is checked rather than assumed.

Feb 29 cohorts carry no season day and are excluded from the endpoint and
slope statistics, as in F2 and F3; they are retained in the domain-wide
fraction.

Usage:
    python compare.py
    python compare.py --csv data/m6_window_sweep.csv
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd


# Reported values from the current classification (Figure 3), used as the
# regression target for the no-truncation cutoff.
REPORTED_M6_START = 45.77
REPORTED_M6_END = 33.61
TOLERANCE_PP = 0.05

FIRST_SEASON_DAY = 0     # Nov 15
LAST_SEASON_DAY = 119    # Mar 14


def season_day_to_label(day: int) -> str:
    """Calendar label for a season day index (e.g., 0 -> 'Nov 15')."""
    if day <= 15:
        return f"Nov {15 + day}"
    if day <= 46:
        return f"Dec {day - 15}"
    if day <= 77:
        return f"Jan {day - 46}"
    if day <= 105:
        return f"Feb {day - 77}"
    return f"Mar {day - 105}"


def load(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # Fold censored into M6.
    df["n_M6_total"] = df["n_killed_M6_no_advance"] + df["n_censored"]
    return df


def climatological_share(df: pd.DataFrame, season_day: int) -> float:
    """Pooled M6 share (%) across years at one season day."""
    sub = df[df["season_day"] == season_day]
    if sub.empty:
        return float("nan")
    return 100.0 * sub["n_M6_total"].sum() / sub["n_particles"].sum()


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cutoff, sub in df.groupby("cutoff", sort=True):
        start = climatological_share(sub, FIRST_SEASON_DAY)
        end = climatological_share(sub, LAST_SEASON_DAY)
        # The common detection window: April 1 to the cutoff. Unlike
        # days_to_cutoff, which varies by cohort, this is the quantity the
        # sweep equalises, and it is the same for every cohort at a cutoff.
        month, day = (int(p) for p in cutoff.split("-"))
        window_d = (dt.date(2001, month, day) - dt.date(2001, 4, 1)).days

        rows.append({
            "cutoff": cutoff,
            "window_d": window_d,
            "demoted_pct": 100.0 * sub["n_demoted"].sum() / sub["n_particles"].sum(),
            "M6_overall_pct": 100.0 * sub["n_M6_total"].sum() / sub["n_particles"].sum(),
            f"M6_{season_day_to_label(FIRST_SEASON_DAY).replace(' ', '')}_pct": start,
            f"M6_{season_day_to_label(LAST_SEASON_DAY).replace(' ', '')}_pct": end,
            "seasonal_slope_pp": start - end,
        })
    return pd.DataFrame(rows).sort_values("cutoff").reset_index(drop=True)


def check_reference(summary: pd.DataFrame) -> None:
    """The no-truncation cutoff must reproduce the reported classification."""
    ref = summary[summary["cutoff"] == "09-30"]
    if ref.empty:
        print("\n[skip] no 09-30 cutoff in the sweep; regression check not run.")
        return

    row = ref.iloc[0]
    start_col = f"M6_{season_day_to_label(FIRST_SEASON_DAY).replace(' ', '')}_pct"
    end_col = f"M6_{season_day_to_label(LAST_SEASON_DAY).replace(' ', '')}_pct"

    problems = []
    if row["demoted_pct"] > 0.0:
        problems.append(f"{row['demoted_pct']:.4f}% demoted, expected 0")
    if abs(row[start_col] - REPORTED_M6_START) > TOLERANCE_PP:
        problems.append(f"Nov 15 M6 {row[start_col]:.2f}%, expected {REPORTED_M6_START}%")
    if abs(row[end_col] - REPORTED_M6_END) > TOLERANCE_PP:
        problems.append(f"Mar 14 M6 {row[end_col]:.2f}%, expected {REPORTED_M6_END}%")

    print()
    if problems:
        print("[FAIL] 09-30 cutoff does not reproduce the reported classification:")
        for p in problems:
            print(f"       - {p}")
    else:
        print("[ok]   09-30 cutoff reproduces the reported classification "
              "(no demotions; M6 matches Figure 3).")


def success_curve(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    Pooled climatological success rate (%) by season day, for one cutoff.

    Returns (season_days, rate_pct) over the days that carry data.
    """
    sub = df[df["season_day"].notna()]
    grouped = sub.groupby("season_day", sort=True).agg(
        s=("n_success", "sum"), n=("n_particles", "sum")
    )
    days = grouped.index.to_numpy().astype(int)
    rate = 100.0 * grouped["s"].to_numpy() / grouped["n"].to_numpy()
    return days, rate


def report_success(df: pd.DataFrame) -> None:
    """
    Report the recruitment-success curve under each cutoff.

    Success is evaluated at the advance event, so it is window-dependent in
    the same way M6 is, but with the opposite consequence for the peak date.
    Early cohorts have the shortest window, so their successes are the ones
    most likely to be missed: the reported curve is biased low at its early
    end, and its peak is therefore biased late. Truncating to a common
    window strips real successes from late cohorts instead, biasing the peak
    early. The two bracket the peak date from opposite sides.
    """
    print()
    print("=" * 78)
    print("Recruitment success under the same cutoffs")
    print("=" * 78)
    print("Success is evaluated at advance, so it moves with the window too.")
    print("The reported (09-30) peak is biased late; the 06-03 peak is biased")
    print("early. The peak date is bracketed between them.")
    print()
    print(f"{'cutoff':>8s} {'window_d':>9s} {'peak day':>10s} {'peak %':>8s} "
          f"{'Nov 15 %':>9s} {'Mar 14 %':>9s} {'season mean %':>14s}")
    print("-" * 78)

    for cutoff, sub in df.groupby("cutoff", sort=True):
        month, day = (int(p) for p in cutoff.split("-"))
        window_d = (dt.date(2001, month, day) - dt.date(2001, 4, 1)).days
        days, rate = success_curve(sub)
        if days.size == 0:
            continue
        i = int(np.nanargmax(rate))
        first = rate[days == FIRST_SEASON_DAY]
        last = rate[days == LAST_SEASON_DAY]
        with_data = sub[sub["season_day"].notna()]
        mean_pct = 100.0 * with_data["n_success"].sum() / with_data["n_particles"].sum()
        print(f"{cutoff:>8s} {window_d:9d} {season_day_to_label(int(days[i])):>10s} "
              f"{rate[i]:8.2f} "
              f"{(first[0] if first.size else float('nan')):9.2f} "
              f"{(last[0] if last.size else float('nan')):9.2f} "
              f"{mean_pct:14.2f}")

    print()
    print("Note: the 06-03 peak is not an estimate of the true peak. It is the")
    print("far side of the bracket. The unbiased answer requires extending the")
    print("early cohorts, not truncating the late ones.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarise the S2 M6 window sweep.")
    parser.add_argument("--csv", type=str, default="data/m6_window_sweep.csv")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.is_file():
        raise SystemExit(f"Not found: {csv_path}. Run sweep.py first.")

    df = load(csv_path)
    summary = summarise(df)

    print()
    print("=" * 78)
    print("M6 under a common advance-detection cutoff")
    print("=" * 78)
    print("Censored folded into M6. Slope is Nov 15 minus Mar 14, in percentage")
    print("points: larger means a steeper seasonal decline.")
    print()
    with pd.option_context("display.width", 200,
                           "display.max_columns", 20,
                           "display.float_format", lambda v: f"{v:8.2f}"):
        print(summary.to_string(index=False))

    check_reference(summary)
    report_success(df)

    print()
    print("Reading: the 09-30 row is the reported classification. Rows above it")
    print("truncate progressively harder, biasing the slope flat, so they bound")
    print("the decline from below. A slope that stays clearly positive across the")
    print("sweep indicates a real seasonal relaxation of the M6 constraint whose")
    print("steepness is uncertain; a slope collapsing toward zero at 06-03 would")
    print("indicate the reported decline is largely an artefact of the window.")
    print()


if __name__ == "__main__":
    main()