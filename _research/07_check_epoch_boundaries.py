
from pathlib import Path
import logging

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]

SYNCED_DIR = Path("/scratch/data/ROAMM/derivatives/synced")
EVENTS_PATH = PROJECT_ROOT / "outputs" / "05_roamm_b2b_events.csv"
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "07_epoch_boundary_check.csv"

LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = LOG_DIR / f"{Path(__file__).stem}.log"

print("SCRIPT:", Path(__file__).resolve())
print("PROJECT ROOT:", PROJECT_ROOT)
print("LOG PATH:", LOG_PATH.resolve())
print("OUTPUT PATH:", OUTPUT_PATH.resolve())

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(LOG_PATH, mode="w"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)


SFREQ = 256
TMIN = -0.1
TMAX = 1.0

pre_samples = int(np.ceil(abs(TMIN) * SFREQ))
post_samples = int(np.ceil(TMAX * SFREQ))

logger.info("Starting epoch boundary sanity check")
logger.info("Sampling rate: %d Hz", SFREQ)
logger.info("Epoch window: [%.1f, %.1f] s", TMIN, TMAX)
logger.info("Pre/post samples: %d / %d", pre_samples, post_samples)

# Load B2B events
events = pd.read_csv(EVENTS_PATH)
logger.info("Loaded %d B2B fixation events", len(events))

results = []

files = sorted(SYNCED_DIR.glob("sub-*/*_mldata.pkl"))
logger.info("Found %d synced files", len(files))

for file in files:
    df = pd.read_pickle(file)

    subject_id = file.parent.name
    run_num = int(df["run_num"].dropna().iloc[0])

    n_samples = len(df)

    run_events = events.loc[
        (events["subject_id"] == subject_id)
        & (events["run_num"] == run_num)
    ].copy()

    if run_events.empty:
        logger.warning(
            "%s run %d: no B2B fixation events",
            subject_id, run_num
        )
        continue

    # Julia uses 1-based sample indexing
    latency = run_events["latency"].to_numpy()

    too_early = latency - pre_samples < 1
    too_late = latency + post_samples > n_samples
    invalid = too_early | too_late

    results.append({
        "subject_id": subject_id,
        "run_num": run_num,
        "n_samples": n_samples,
        "n_fixations": len(run_events),
        "too_early": int(too_early.sum()),
        "too_late": int(too_late.sum()),
        "invalid": int(invalid.sum()),
        "invalid_pct": 100 * invalid.mean(),
    })

    logger.info(
        "%s run %d: %d fixations, %d cross-boundary (%.3f%%)",
        subject_id,
        run_num,
        len(run_events),
        invalid.sum(),
        100 * invalid.mean()
    )


# Save run-level results
results = pd.DataFrame(results)

if results.empty:
    raise ValueError("No matching B2B fixation events were found.")

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
results.to_csv(OUTPUT_PATH, index=False)

logger.info("Saved boundary check results to %s", OUTPUT_PATH)


# Overall summary
total = results["n_fixations"].sum()
n_early = results["too_early"].sum()
n_late = results["too_late"].sum()
n_invalid = results["invalid"].sum()

logger.info("========== Overall summary ==========")
logger.info("Total B2B fixations: %d", total)
logger.info("Too early: %d", n_early)
logger.info("Too late: %d", n_late)
logger.info("Cross-boundary epochs: %d", n_invalid)
logger.info("Percentage: %.3f%%", 100 * n_invalid / total)
logger.info("Within-run complete epochs: %d", total - n_invalid)

# Sanity checks
if total != len(events):
    logger.warning(
        "Checked %d events, but CSV contains %d events",
        total, len(events)
    )

if n_invalid == 0:
    logger.info(
        "All checked B2B epochs fit within their original runs."
    )
else:
    logger.warning(
        "%d epochs extend beyond their original run boundaries.",
        n_invalid
    )
