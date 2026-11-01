import matplotlib as mpl
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#7a5cc2']
MK = ['o', 's', '^', 'D', 'v', 'P', 'X']
INK, INK2, GRID = '#222222', '#555555', '#e3e3e3'
def apply():
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 10,
        'axes.labelsize': 9, 'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8,
        'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.8,
        'lines.markersize': 5.5, 'legend.frameon': False, 'savefig.dpi': 300,
        'savefig.bbox': 'tight', 'figure.dpi': 150})
