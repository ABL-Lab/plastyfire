"""Closed-loop check of the open-loop spine VDCC replay: rerun local_t/sims.py with ljp_VDCC = 25 (variant s25).
    python cell_audit/vdcc_decode/closed_loop.py --pair P --variants s25 --protos ap1,burst,epsp,sj20@-10 --out X.npz"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "local_t"))
import sims  # noqa: E402
sims.SYNG["s25"] = {"ljp_VDCC": 25.0}
if __name__ == "__main__":
    sims.main()
