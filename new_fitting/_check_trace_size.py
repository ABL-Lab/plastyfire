"""Quick diagnostic: how big are the traces?"""
import pickle, numpy as np, os

trace_dir = '/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/trace_results/Chindemi_params'
dirs = sorted(os.listdir(trace_dir))[:3]
for d in dirs:
    for proto in ['2Hz_5ms', '10Hz_10ms', '50Hz_10ms']:
        pkl = os.path.join(trace_dir, d, proto, 'simulation_traces.pkl')
        if not os.path.exists(pkl):
            continue
        with open(pkl, 'rb') as f:
            data = pickle.load(f)
        t = np.asarray(data['t'])
        cai = np.asarray(data['cai_CR'])
        dt_arr = np.diff(t)
        print(f'{d}/{proto}: n_time={len(t)}, n_syn={cai.shape[0]}, '
              f't_range=[{t[0]:.1f}, {t[-1]:.1f}]ms, '
              f'dt=[{dt_arr.min():.4f}, {np.mean(dt_arr):.4f}, {dt_arr.max():.4f}]ms')
