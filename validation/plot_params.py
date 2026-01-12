import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load data
# data_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/validation/synaptic_params.csv"
data_path = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/validation/synaptic_params_L23PC_L5TTPC.csv"
df = pd.read_csv(data_path)

# Rename columns for plotting
rename_dict = {
    'u_syn': r'$U_{SE}$',
    'depression_time': r'$\tau_{dep}$ (ms)',
    'facilitation_time': r'$\tau_{fac}$ (ms)',
    'n_rrp_vesicles': r'$N_{RRP}$',
    'conductance': r'$\hat{G}_{AMPAR}$ (nS)',
    'g_nmdar': r'$\hat{G}_{NMDAR}$ (nS)',
    'volume_CR': r'$X$ ($\mu m^3$)'
}

# Select and order columns
cols = [
    'u_syn',
    'depression_time',
    'facilitation_time',
    'n_rrp_vesicles',
    'conductance',
    'g_nmdar',
    'volume_CR'
]

plot_df = df[cols].rename(columns=rename_dict)

# Set style
sns.set_theme(style="ticks")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

# Create PairGrid
# Decrease height to decrease overall size (default is 2.5)
g = sns.PairGrid(plot_df, diag_sharey=False, corner=True, height=1.5)

# Map lower triangle to scatter plot
g.map_lower(sns.scatterplot, s=10, alpha=0.6, color="#9C002C", edgecolor="none")

# Map diagonal to histogram
def plot_diag(x, **kwargs):
    ax = plt.gca()
    # Plot histogram
    sns.histplot(x, ax=ax, color="#9C002C", edgecolor="white", linewidth=0.5)
    
    # Calculate stats
    mu = x.mean()
    sigma = x.std()
    
    # Add annotation
    # Position: top right of the axes
    # We need to be careful with y-axis limits as histplot sets them.
    # Let's get current ylim
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    
    # Draw error bar for mean +/- std
    # We'll place it at some height, maybe 50% of max count?
    # Or just use the style from the image: orange dot with horizontal line
    
    # Find max count to position the annotation
    # We can get the patches from the histogram
    max_h = 0
    for rect in ax.patches:
        if rect.get_height() > max_h:
            max_h = rect.get_height()
            
    y_pos = max_h * 0.6
    
    # Plot mean point
    ax.plot(mu, y_pos, 'o', color='#fdae61', markersize=5, zorder=10)
    # Plot std line
    ax.plot([mu - sigma, mu + sigma], [y_pos, y_pos], '-', color='#fdae61', linewidth=1.5, zorder=10)
    
    # Add text
    # Connect text to point with a line if needed, but the image shows text floating with a line pointing to the dot?
    # Actually the image shows a bracket-like line or just a line pointing to the dot.
    # Let's just put text near it.
    
    text_x = xlim[1] * 0.6 # A bit arbitrary, need to adjust based on data range
    if text_x < mu: text_x = mu + sigma * 1.5
    
    # Better: use relative coordinates
    # ax.text(0.6, 0.8, f'$\mu = {mu:.1f}$\n$\sigma = {sigma:.1f}$', transform=ax.transAxes)
    
    # But we want to match the look. The text is outside the plot area sometimes or top right.
    # Let's try to place it at (mu + sigma, y_pos) and see.
    
    # Using annotation with arrow
    ax.annotate(f'$\mu = {mu:.1f}$\n$\sigma = {sigma:.1f}$', 
                xy=(mu, y_pos), xycoords='data',
                xytext=(20, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="-|>", connectionstyle="angle,angleA=0,angleB=90,rad=10", color='black'),
                fontsize=12)

g.map_diag(plot_diag)

# Adjust layout
plt.subplots_adjust(hspace=0.1, wspace=0.1)

# Save
# output_plot = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/validation/synaptic_params_plot.png"
output_plot = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/validation/synaptic_params_plot_L23PC_L5TTPC.png"
g.savefig(output_plot, dpi=300)
print(f"Saved plot to {output_plot}")
