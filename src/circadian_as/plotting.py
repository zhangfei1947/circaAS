"""PSI and depth together; zero-depth PSI is missing, never zero."""
import json
from pathlib import Path


def plot_event(event_path, output):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    obj = json.loads(Path(event_path).read_text())
    record = obj['record']
    time = np.asarray(record['times'])
    inc, skip = np.asarray(record['I']) / record['LI'], np.asarray(record['S']) / record['LS']
    psi = np.divide(inc, inc + skip, out=np.full(len(time), np.nan), where=(inc + skip) > 0)
    reps = np.asarray(record['reps']).astype(str)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    predictions = (obj.get('analysis') or {}).get('predictions')
    for block in np.unique(reps):
        sel = reps == block
        points = axes[0].scatter(time[sel], psi[sel], label=block, s=24)
        color = points.get_facecolor()[0]
        axes[1].plot(time[sel][np.argsort(time[sel])], np.asarray(record['N'])[sel][np.argsort(time[sel])],
                     '.-', color=color, label=block)
        if predictions:
            curve = next((p for p in predictions['blocks'] if p['block'] == block), None)
            if curve:
                grid = np.linspace(time.min(), time.max(), 500)
                axes[0].plot(grid, np.interp(grid % 24, predictions['times'], curve['psi']), color=color)
    axes[0].set(ylabel='Effective-length-adjusted PSI', ylim=(-.04, 1.04))
    axes[1].set(ylabel='Informative reads (I + S)', xlabel='Time (hours)')
    axes[0].legend(title='Experiment')
    for cell in obj['coverage']['time_depth']:
        if cell['total_count'] == 0:
            axes[0].axvline(cell['time'], color='.7', linestyle=':', alpha=.6)
    p = (obj.get('analysis') or {}).get('p_value')
    fig.suptitle(f"{record['dataset']} | {record['condition']} | {record['event_id']}\n"
                 f"{obj['status']} | raw p={p if p is not None else 'NA'}")
    fig.text(.5, .01, 'Exploratory fit. Missing depth is not zero PSI. Lines are fitted means, not confidence bands.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .03, 1, .94))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
