from pathlib import Path
import logging
import numpy as np 
import pandas as pd 


SYNCED_DIR = Path("/scratch/data/ROAMM/derivatives/synced")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = PROJECT_ROOT / "outputs" / "roamm_fixation_events.csv"

LOG_DIR = PROJECT_ROOT / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_PATH = LOG_DIR / f"{Path(__file__).stem}.log"

logging.basicConfig(
    level = logging.INFO,
    format = "%(asctime)s - %(levelname)s - %(message)s",
    handlers = [
        logging.FileHandler(LOG_PATH, mode="w"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

def nearest_sample_indices(time, onsets):
    """Return nearest EEG sample index for each fixation onset."""
    idx = np.searchsorted(time, onsets)
    idx = np.clip(idx, 1, len(time) - 1)

    left = idx - 1
    right = idx

    use_right = (np.abs(time[right] - onsets) < np.abs(time[left] - onsets))
    return np.where(use_right, right, left)

files = sorted(file for file in SYNCED_DIR.glob("sub-*/*_mldata.pkl") if file.exists())
all_events = []

for file in files:
    df = pd.read_pickle(file)
    subject_id = file.parent.name
    run_num = int(df["run_num"].dropna().iloc[0])

    # First-pass right-eye fixations
    events = df.loc[
        df["first_pass_reading"].eq(True) 
        & df["fix_R_tStart"].notna(),
        [
            "page_num",
            "story_name",
            "fix_R_tStart",
            "fix_R_tEnd",
            "fix_R_duration",
            "fix_R_fixed_word",
            "fix_R_fixed_word_key"
        ]
    ].copy()

    # One row per fixation
    events = (
        events
        .drop_duplicates(subset=["fix_R_tStart"])
        .sort_values("fix_R_tStart")
        .reset_index(drop=True)
    )
    events = events.rename(
        columns={
            "page_num": "page",
            "fix_R_tStart": "onset_s",
            "fix_R_tEnd": "offset_s",
            "fix_R_duration": "duration_ms",
            "fix_R_fixed_word": "word",
            "fix_R_fixed_word_key": "word_key"
        }
    )

    # Map fixation onsets to EEF sample
    time = df["time"].to_numpy()
    sample_idx = nearest_sample_indices(time, events["onset_s"].to_numpy())

    events["subject_id"] = subject_id
    events["run_num"] = run_num
    events["eye"] = "R"

    # Julia uses 1-based indexing
    events["latency"] = sample_idx + 1

    sync_error_ms = (time[sample_idx] - events["onset_s"].to_numpy()) * 1000

    logger.info(
        "%s run %d eye R: %d fixation, max sync error %.3f ms",
        subject_id,
        run_num,
        len(events),
        np.abs(sync_error_ms).max()
    )

    all_events.append(events)





# Combine + checks
events = pd.concat(all_events, ignore_index=True)
events = events[
    [
        "subject_id",
        "run_num",
        "story_name",
        "page",
        "eye",
        "onset_s",
        "offset_s",
        "duration_ms",
        "word",
        "word_key",
        "latency",
    ]
]

assert not events.duplicated(subset=["subject_id", "run_num", "eye", "onset_s"]).any()

# Save
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
events.to_csv(OUTPUT_PATH, index=False)
logger.info("Saved %d fixation events to %s", len(events), OUTPUT_PATH)