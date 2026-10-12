from pathlib import Path
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PREPROCESSING_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = PREPROCESSING_ROOT.parent
SFREQ = 256
TMIN, TMAX = -0.1, 1.0
OFFSETS = np.arange(round(TMIN * SFREQ), round(TMAX * SFREQ) + 1)
TIMES_MS = OFFSETS / SFREQ * 1000


def make_logger(name, path):
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    for handler in (logging.FileHandler(path, mode='w', encoding='utf-8'), logging.StreamHandler()):
        handler.setFormatter(fmt)
        logger.addHandler(handler)
    return logger

def qc_run(path, run_events):
    eeg = np.load(path, mmap_mode='r')
    if eeg.ndim != 2 or 0 in eeg.shape:
        raise ValueError(f'Expected nonempty (channels, samples), got {eeg.shape}')
    n_channels, n_samples = eeg.shape
    eeg_uv = np.asarray(eeg, dtype=np.float64) * 1e6
    finite = np.isfinite(eeg_uv)
    finite_abs = np.abs(eeg_uv[finite])
    channel_std = np.nanstd(np.where(finite, eeg_uv, np.nan), axis=1) # if one channel is all NaN, std will be NaN
    flat_channels = np.where(channel_std < 0.1)[0] + 1 # 1-based indexing

    latencies = pd.to_numeric(run_events['latency'], errors='coerce').to_numpy(dtype=float)
    valid_latency = np.isfinite(latencies) & (latencies == np.floor(np.nan_to_num(latencies)))
    valid_latency &= (latencies >= 1) & (latencies <= n_samples)
    starts = np.zeros(len(latencies), dtype=np.int64)
    starts[valid_latency] = latencies[valid_latency].astype(np.int64) - 1
    valid_epoch = valid_latency & (starts + OFFSETS[0] >= 0) & (starts + OFFSETS[-1] < n_samples)

    epoch_sum = np.zeros((n_channels, len(OFFSETS)), dtype=np.float64)
    used = 0
    for start in starts[valid_epoch]:
        epoch = eeg_uv[:, start + OFFSETS]
        if np.isfinite(epoch).all():
            epoch_sum += epoch
            used += 1
    metrics = dict(
        n_channels=n_channels, n_samples=n_samples, duration_min=n_samples / SFREQ / 60, 
        n_events=len(run_events), n_invalid_latency=int((~valid_latency).sum()),
        n_cross_boundary=int((valid_latency & ~valid_epoch).sum()),
        n_epochs_nonfinite=int(valid_epoch.sum()) - used, n_epochs_used=used,
        nan_count = int(np.isnan(eeg_uv).sum()), inf_count=int(np.isinf(eeg_uv).sum()), 
        channels_with_nonfinite=int((~finite).any(axis=1).sum()),
        p99_abs_uv=float(np.percentile(finite_abs, 99)) if finite_abs.size else np.nan,
        max_abs_uv=float(finite_abs.max()) if finite_abs.size else np.nan,
        n_flat_channels=len(flat_channels),
        flat_channels=','.join(map(str, flat_channels)),
    )
    return metrics, (epoch_sum / used if used else None)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--events', type=Path, default=PREPROCESSING_ROOT / 'outputs/05_roamm_b2b_events.csv')
    parser.add_argument('--eeg-dir', type=Path, default=Path('/scratch/data/ROAMM/outputs/eeg'))
    parser.add_argument('--out-dir', type=Path, default=PREPROCESSING_ROOT / 'outputs')
    parser.add_argument('--fig-dir', type=Path, default=PROJECT_ROOT / 'plots/2026-10-12-quality_check')
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.fig_dir.mkdir(parents=True, exist_ok=True)
    log_dir = PREPROCESSING_ROOT / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    eeg_log = make_logger('roamm.eeg_qc', log_dir / '07_quality_check.log')
    fixation_log = make_logger('roamm.fixation_qc', log_dir / '07_quality_check_fixation.log')

    events = pd.read_csv(args.events, usecols=['subject_id', 'run_num', 'latency'])
    events['subject_id'] = events['subject_id'].astype(str)
    events['run_num'] = pd.to_numeric(events['run_num'], errors='raise').astype(int)
    eeg_files = sorted(args.eeg_dir.glob('sub-*_run*_eeg.npy'))
    eeg_log.info('Found %d EEG files in %s', len(eeg_files), args.eeg_dir)
    fixation_log.info('Loadded %d fixation events from %s', len(events), args.events)

    rows = []
    for subject, sub_events in events.groupby('subject_id', sort=True):
        eeg_log.info('Processing %s', subject)
        frp_sum = None
        total_epochs = 0
        for run, run_events in sub_events.groupby('run_num', sort=True):
            path = args.eeg_dir / f'{subject}_run{run}_eeg.npy'
            if not path.is_file():
                eeg_log.warning('Missing EEG %s', path)
                fixation_log.warning('%s run %s: EEG missing; %d events unchecked', subject, run, len(run_events))
                rows.append(dict(subject_id=subject, run_num=run, status='missing_eeg', n_events=len(run_events)))
                continue
            try:
                metrics, frp = qc_run(path, run_events) 
            except Exception as exc:
                eeg_log.exception('Failed %s run %s', subject, run)
                fixation_log.exception('%s run %s: EEG QC failed; %s', subject, run, exc)
                rows.append(dict(subject_id=subject, run_num=run, status=f'error:{exc}', n_events=len(run_events)))
                continue
            rows.append(dict(subject_id=subject, run_num=run, status='ok', **metrics))
            eeg_log.info('%s run %s: %d channels, %d samples, NaN=%d, Inf=%d, p99=%.2f, µV, flat=%d',
                        subject, run, metrics['n_channels'], metrics['n_samples'], metrics['nan_count'],
                        metrics['inf_count'], metrics['p99_abs_uv'], metrics['n_flat_channels'])
            fixation_log.info('%s run %s: events=%d, invalid_latency=%d, cross-boundary=%d, nonfinite_epochs=%d, used_epochs=%d',
                        subject, run, metrics['n_events'], metrics['n_invalid_latency'], metrics['n_cross_boundary'],
                        metrics['n_epochs_nonfinite'], metrics['n_epochs_used'])
            if frp is not None:
                if frp_sum is None:
                    frp_sum = np.zeros_like(frp)
                if frp_sum.shape != frp.shape:
                    eeg_log.warning('%s run %s: channel count mismatch; omitting run from plot', subject, run)
                    continue
                frp_sum += frp * metrics['n_epochs_used']
                total_epochs += metrics['n_epochs_used']

        if total_epochs:
            avg = frp_sum / total_epochs
            fig, ax = plt.subplots(figsize=(10, 5))
            colors = plt.get_cmap("viridis")(np.linspace(0, 1, avg.shape[0]))
            for ch, color in enumerate(colors):
                ax.plot(TIMES_MS, avg[ch], color=color, alpha=0.7, linewidth=0.8)
            ax.axvline(0, color='k', linestyle='--', linewidth=0.9)
            ax.axhline(0, color='k', linewidth=0.5)
            ax.set(xlabel='Time from fixation onset (ms)', ylabel='Amplitude (µV)', title=f'{subject} | raw FRP butterfly | {total_epochs} fixation')
            fig.tight_layout()
            output_fig = args.fig_dir / f'{subject}_raw_frp.png'
            fig.savefig(output_fig, dpi=150)
            plt.close(fig)
        else:
            eeg_log.warning('%s: no valid epochs for FRP plot', subject)

    summary = pd.DataFrame(rows)
    csv_path = args.out_dir / '07_eeg_quality_check_summary.csv'
    summary.to_csv(csv_path, index=False)
    eeg_log.info('Saved summary: %s (%d run rows)', csv_path, len(summary))
    if not summary.empty:
        eeg_log.info('Status counts: %s', summary['status'].value_counts().to_dict())
        fixation_log.info('Events in input: %d; events in processed run rows: %d',
                        len(events), int(summary['n_events'].fillna (0).sum()))
        good = summary[summary['status'] == 'ok']
        if not good.empty:
            fixation_log.info('Total invalid latencies=%d; cross-boundary=%d; nonfinite epochs=%d',
                                int(good['n_invalid_latency'].sum()), int(good['n_cross_boundary'].sum()), 
                                int(good['n_epochs_nonfinite'].sum()))
    eeg_log.info('Figures: %s', args.fig_dir)
    eeg_log.info('Logs: %s', log_dir)

if __name__ == '__main__':
    main()