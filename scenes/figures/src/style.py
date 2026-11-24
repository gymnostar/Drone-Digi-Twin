"""Shared palette, fonts and paths for the figures in scenes/figures."""
import os
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INK = '#0b0b0b'; INK2 = '#52514e'; MUTED = '#8a8984'; GRID = '#e4e3df'; SURF = '#ffffff'
BLUE = '#2a78d6'; ORANGE = '#eb6834'; AQUA = '#1baf7a'; YELLOW = '#eda100'
MAGENTA = '#e87ba4'; GREEN = '#008300'; VIOLET = '#6250d6'; RED = '#e34948'
SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
# Inputs are read from the repository's scenes/ folder; figures are written to
# scenes/figures/ unless FIG_OUT names another folder (useful for comparisons).
SCENES = Path(__file__).resolve().parents[2]
FIGURES = Path(os.environ.get('FIG_OUT', SCENES / 'figures'))
OUT = str(FIGURES) + '/'

plt.rcParams.update({
    'figure.facecolor': SURF, 'axes.facecolor': SURF, 'savefig.facecolor': SURF,
    'font.family': ['Helvetica', 'DejaVu Sans'], 'font.size': 17, 'axes.titlesize': 21,
    'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.titlepad': 12,
    'axes.labelsize': 18, 'axes.labelcolor': INK2, 'xtick.labelsize': 16,
    'ytick.labelsize': 16, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.edgecolor': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
    'legend.fontsize': 16, 'legend.frameon': False, 'text.color': INK,
    'axes.grid': False, 'grid.color': GRID, 'grid.linewidth': 1,
})


def save(fig, name):
    fig.savefig(OUT + name, dpi=150)
    plt.close(fig)
    print('wrote', OUT + name)
