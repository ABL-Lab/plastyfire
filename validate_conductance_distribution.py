"""
Validate conductance values in SONATA edges.h5 against recipe parameters
for 4 (pre_mtype, post_mtype) pathways.
"""
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gamma as gamma_dist
from bluepysnap import Circuit

# --- Config ---
CIRCUIT_CONFIG = '/project/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/circuit_config.json'
EDGES_FILE = ('/project/ctb-emuller/datasets/SSCx_plasticity/O1_2023a_Ecker/'
              'S1nonbarrel_neurons__S1nonbarrel_neurons__chemical/edges.h5')
RECIPE_CSV = '/project/ctb-emuller/dhuruva/plastyfire/biodata/recipe_mod.csv'
EDGE_POP = 'S1nonbarrel_neurons__S1nonbarrel_neurons__chemical'
NODE_POP = 'S1nonbarrel_neurons'
OUTPUT_FIG = '/lustre06/project/6077694/dhuruva/plastyfire/conductance_validation.png'

MAX_POST_GIDS = 200
PATHWAYS = [
    ('L5_TPC:A', 'L5_TPC:A'),
    ('L5_TPC:A', 'L5_TPC:B'),
    ('L5_TPC:B', 'L5_TPC:A'),
    ('L5_TPC:B', 'L5_TPC:B'),
]

# Recipe reference values (same for all 4 pathways)
RECIPE_GSYN = 1.94
RECIPE_SD = 1.0
RECIPE_SHAPE = RECIPE_GSYN**2 / RECIPE_SD**2   # 3.7636
RECIPE_SCALE = RECIPE_SD**2 / RECIPE_GSYN       # 0.5154

def load_node_mtypes():
    print("Loading node mtypes via bluepysnap...")
    c = Circuit(CIRCUIT_CONFIG)
    df = c.nodes[NODE_POP].get(properties=['mtype'])
    print(f"  Loaded {len(df)} nodes")
    return df  # index=node_id, column='mtype'

def get_conductances_for_pathway(f, node_df, pre_mtype, post_mtype, rng):
    pop_grp = f[f'edges/{EDGE_POP}']
    src_ids = pop_grp['source_node_id'][:]
    conductances_all = pop_grp['0/conductance'][:]
    try:
        syn_type_ids = pop_grp['0/syn_type_id'][:]
        has_syn_type = True
    except KeyError:
        has_syn_type = False
        print("  Warning: syn_type_id not found, skipping excitatory filter")

    tgt_to_src_grp = pop_grp['indices/target_to_source']
    node_id_to_ranges = tgt_to_src_grp['node_id_to_ranges'][:]
    range_to_edge_id = tgt_to_src_grp['range_to_edge_id'][:]

    pre_gids = set(node_df[node_df['mtype'] == pre_mtype].index.tolist())
    post_gids = node_df[node_df['mtype'] == post_mtype].index.tolist()
    print(f"  pre_mtype={pre_mtype}: {len(pre_gids)} cells | post_mtype={post_mtype}: {len(post_gids)} cells")

    if len(post_gids) > MAX_POST_GIDS:
        post_gids = rng.choice(post_gids, MAX_POST_GIDS, replace=False).tolist()
        print(f"  Sampled {MAX_POST_GIDS} post_gids")

    conductances = []
    for pgid in post_gids:
        if pgid >= len(node_id_to_ranges):
            continue
        r_start, r_end = node_id_to_ranges[pgid]
        if r_start == r_end:
            continue
        edge_ranges = range_to_edge_id[r_start:r_end]
        for e_start, e_end in edge_ranges:
            chunk_src = src_ids[e_start:e_end]
            mask = np.isin(chunk_src, list(pre_gids))
            if has_syn_type:
                mask &= (syn_type_ids[e_start:e_end] >= 100)
            cond_vals = conductances_all[e_start:e_end][mask]
            conductances.extend(cond_vals.tolist())

    return np.array(conductances, dtype=float)

def main():
    import pandas as pd
    rng = np.random.default_rng(42)

    node_df = load_node_mtypes()

    print("\nOpening edges file...")
    f = h5py.File(EDGES_FILE, 'r')

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()

    summary_rows = []

    for idx, (pre_mtype, post_mtype) in enumerate(PATHWAYS):
        pathway_label = f'{pre_mtype} -> {post_mtype}'
        print(f"\n{'='*60}")
        print(f"Pathway: {pathway_label}")

        conds = get_conductances_for_pathway(f, node_df, pre_mtype, post_mtype, rng)
        print(f"  Collected {len(conds)} conductance values")

        if len(conds) < 10:
            print("  WARNING: Too few values to fit!")
            axes[idx].set_title(f'{pathway_label}\n(insufficient data)', fontsize=9)
            continue

        # Filter out non-positive values for gamma fit
        conds_pos = conds[conds > 0]
        print(f"  Positive conductance values: {len(conds_pos)}")

        # Fit gamma distribution
        fit_shape, fit_loc, fit_scale = gamma_dist.fit(conds_pos, floc=0)
        fit_mean = fit_shape * fit_scale
        fit_sd = np.sqrt(fit_shape) * fit_scale

        print(f"\n  --- Fitted vs Recipe ---")
        print(f"  {'Param':<12} {'Fitted':>10} {'Recipe':>10}")
        print(f"  {'mean':<12} {fit_mean:>10.4f} {RECIPE_GSYN:>10.4f}")
        print(f"  {'SD':<12} {fit_sd:>10.4f} {RECIPE_SD:>10.4f}")
        print(f"  {'shape (a)':<12} {fit_shape:>10.4f} {RECIPE_SHAPE:>10.4f}")
        print(f"  {'scale':<12} {fit_scale:>10.4f} {RECIPE_SCALE:>10.4f}")
        print(f"  {'loc':<12} {fit_loc:>10.4f} {'0':>10}")

        summary_rows.append({
            'pathway': pathway_label,
            'n_synapses': len(conds_pos),
            'fitted_mean': fit_mean,
            'recipe_mean': RECIPE_GSYN,
            'fitted_sd': fit_sd,
            'recipe_sd': RECIPE_SD,
            'fitted_shape': fit_shape,
            'recipe_shape': RECIPE_SHAPE,
            'fitted_scale': fit_scale,
            'recipe_scale': RECIPE_SCALE,
        })

        # Plot
        ax = axes[idx]
        x_max = max(conds_pos.max(), RECIPE_GSYN + 4 * RECIPE_SD)
        x = np.linspace(0, x_max, 500)

        ax.hist(conds_pos, bins=60, density=True, alpha=0.5, color='steelblue',
                label=f'Data (n={len(conds_pos)})')
        ax.plot(x, gamma_dist.pdf(x, fit_shape, loc=fit_loc, scale=fit_scale),
                'b-', lw=2, label=f'Fitted γ (a={fit_shape:.2f}, s={fit_scale:.3f})')
        ax.plot(x, gamma_dist.pdf(x, RECIPE_SHAPE, loc=0, scale=RECIPE_SCALE),
                'r--', lw=2, label=f'Recipe γ (a={RECIPE_SHAPE:.2f}, s={RECIPE_SCALE:.3f})')

        ax.set_title(f'{pathway_label}\n'
                     f'Fitted: μ={fit_mean:.3f}, σ={fit_sd:.3f} | '
                     f'Recipe: μ={RECIPE_GSYN}, σ={RECIPE_SD}', fontsize=8)
        ax.set_xlabel('Conductance (nS)')
        ax.set_ylabel('Density')
        ax.legend(fontsize=7)
        ax.set_xlim(0, min(x_max, 15))

    f.close()

    fig.suptitle('Conductance Distribution Validation: Fitted vs Recipe (gamma)\n'
                 'Recipe: gsyn=1.94 nS, gsynSD=1.0 nS, dist=gamma', fontsize=11)
    fig.tight_layout()
    fig.savefig(OUTPUT_FIG, dpi=150, bbox_inches='tight')
    print(f"\nFigure saved to {OUTPUT_FIG}")

    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70)
    print(f"{'Pathway':<30} {'N':>7} {'μ_fit':>7} {'μ_rec':>7} {'σ_fit':>7} {'σ_rec':>7} {'a_fit':>7} {'a_rec':>7}")
    print("-"*70)
    for r in summary_rows:
        print(f"{r['pathway']:<30} {r['n_synapses']:>7} "
              f"{r['fitted_mean']:>7.3f} {r['recipe_mean']:>7.3f} "
              f"{r['fitted_sd']:>7.3f} {r['recipe_sd']:>7.3f} "
              f"{r['fitted_shape']:>7.3f} {r['recipe_shape']:>7.3f}")

if __name__ == '__main__':
    main()
