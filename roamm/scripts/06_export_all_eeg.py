from pathlib import Path
import glob
import re

import numpy as np
import pandas as pd


SYNCED_DIR = Path("/scratch/data/ROAMM/derivatives/synced")
OUT_DIR = Path("/scratch/data/ROAMM/outputs/eeg")
OUT_DIR.mkdir(parents=True, exist_ok=True)



EEG_CHANNELS = [
    "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7",
    "FC5", "FC3", "FC1", "C1", "C3", "C5", "T7", "TP7",
    "CP5", "CP3", "CP1", "P1", "P3", "P5", "P7", "P9",
    "PO7", "PO3", "O1", "Iz", "Oz", "POz", "Pz", "CPz",
    "Fpz", "Fp2", "AF8", "AF4", "Afz", "Fz", "F2", "F4",
    "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz",
    "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2",
    "P2", "P4", "P6", "P8", "P10", "PO8", "PO4", "O2",
]

pattern = re.compile(
    r"(sub-\d+)_task-ReMind_run-(\d+)_mldata\.pkl$"
)

# glob works through the subject symlink directories
files = sorted(
    glob.glob(
        str(SYNCED_DIR / "sub-*" / "*_task-ReMind_run-*_mldata.pkl")
    )
)

print(f"Found {len(files)} synced run files.")

for i, filepath in enumerate(files, start=1):
    filepath = Path(filepath)

    m = pattern.match(filepath.name)
    if m is None:
        print("Skipping unrecognized filename:", filepath.name)
        continue

    subject = m.group(1)
    run = int(m.group(2))

    outfile = OUT_DIR / f"{subject}_run{run}_eeg.npy"

    print(
        f"[{i}/{len(files)}] "
        f"{subject} run {run}"
    )

    # allows restarting the script without doing everything again
    if outfile.exists():
        print("    already exists -> skip")
        continue

    df = pd.read_pickle(filepath)

    missing_channels = [
        ch for ch in EEG_CHANNELS
        if ch not in df.columns
    ]

    if missing_channels:
        print(
            "    ERROR: missing EEG channels:",
            missing_channels
        )
        continue

    # dataframe = samples × channels
    # Julia B2B = channels × samples
    eeg = (
        df[EEG_CHANNELS]
        .to_numpy(dtype=np.float64)
        .T
    )

    print("    EEG shape:", eeg.shape)

    np.save(outfile, eeg)

    print("    saved:", outfile.name)

print("Finished EEG export.")