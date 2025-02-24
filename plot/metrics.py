import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, auc


def calculate_metrics(labels, scores):
    fpr, tpr, _ = roc_curve(labels, scores)
    roc_auc = roc_auc_score(labels, scores)

    precision, recall, _ = precision_recall_curve(labels, scores)
    pr_auc = auc(recall, precision)

    return {'fpr': fpr, 'tpr': tpr, 'roc_auc': roc_auc,
            'precision': precision, 'recall': recall, 'pr_auc': pr_auc}


metrics = {
    'gene': calculate_metrics(labels_gene, scores['gene']),
    'prot': calculate_metrics(labels_prot, scores['prot']),
    'combined': calculate_metrics(labels_combined, scores['combined'])
}


def configure_plot_style():
    plt.rcParams.update({
        'font.sans-serif': 'Arial',
        'font.size': 8,
        'axes.labelsize': 8,
        'legend.fontsize': 7,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'lines.linewidth': 1.5,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'figure.dpi': 300,
        'figure.figsize': (8.5 / 2.54, 8.5 / 2.54)
    })


# --------------------------
# 绘图模块
# --------------------------
def plot_roc_curves(metrics_dict):
    configure_plot_style()

    fig, ax = plt.subplots()
    color_map = {'gene': 'green', 'prot': 'blue', 'combined': 'red'}
    line_style = (0, (3, 2))

    for model in ['gene', 'prot', 'combined']:
        ax.plot(metrics_dict[model]['fpr'], metrics_dict[model]['tpr'],
                color=color_map[model],
                alpha=0.5,
                linestyle=line_style,
                label=f'{model.capitalize()} ({metrics_dict[model]["roc_auc"]:.3f})')

    ax.plot([0, 1], [0, 1],
            color='black',
            linestyle=line_style,
            linewidth=1,
            zorder=1)

    ax.set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.02),
           xlabel='False Positive Rate (FPR)',
           ylabel='True Positive Rate (TPR)',
           aspect='equal')

    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.linspace(0, 1, 6))

    ax.legend(frameon=False,
              loc='lower right',
              handlelength=1.5,
              handletextpad=0.4,
              borderaxespad=0.5)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15)
    plt.show()


plot_roc_curves(metrics)
