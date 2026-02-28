import pickle
import os
import glob

workdir = "/home/dhuruva/projects/ctb-emuller/dhuruva/plastyfire/refitting_results/fitting/n100/seed19091997/L5TTPC_L5TTPC/simulations/205559-199162/10Hz_5ms"
pkl_files = glob.glob(os.path.join(workdir, "*.pkl"))

if pkl_files:
    pkl_file = pkl_files[0]
    print(f"Inspecting {pkl_file}")
    try:
        with open(pkl_file, "rb") as f:
            data = pickle.load(f)
            print("Keys:", data.keys())
            if "synprop" in data:
                print("Synprop keys:", data["synprop"].keys())
                if "Cpre" in data["synprop"]:
                    print("Cpre sample:", data["synprop"]["Cpre"][:5])
                if "Cpost" in data["synprop"]:
                    print("Cpost sample:", data["synprop"]["Cpost"][:5])
    except Exception as e:
        print(f"Error loading pickle: {e}")
else:
    print("No pickle files found")
