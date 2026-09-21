"""
Restore the paper's n = 100 by dropping the pairs whose post cell has no valid
c_post and letting `find_pairs` supply replacements.

Why this leaves the other 98 untouched
--------------------------------------
`OptSimWriter.find_pairs` is fully deterministic:

    np.random.seed(self.seed); np.random.shuffle(post_gids)     # fixed order
    args_list = [(post_gid, i, ...) for i, post_gid in enumerate(post_gids)]
    ...  np.random.seed(seed + i); np.random.choice(pre_gids, 1)[0]

The pre partner is drawn from `seed + i`, where `i` is the position in the
shuffled array -- NOT the count of pairs accepted so far. Rejecting a cell
therefore shifts nothing: every surviving pair keeps its exact pre_gid, and
find_pairs just walks further down the same list until it has `npairs` again.

So the swap is only: make the two bad cells fail `check_electrical_constraint`
(the single-AP gate now re-validates stale stimulus pkls, see simwriter.py),
re-run find_pairs, and write sim files for whatever is new. Cells already
calibrated short-circuit on their pkl, so the 98 survivors cost nothing.

    python swap_bad_pairs.py                  # diff only, writes NOTHING
    python swap_bad_pairs.py --apply          # + write sim files for new pairs
    python swap_bad_pairs.py --apply --retire # + move dead pair dirs aside

Needs a compute node (find_pairs runs NEURON threshold searches).
"""

import os
import sys
import glob
import shutil
import argparse

from plastyfire.simwriter import OptSimWriter

CONFIG = "configs/L5TTPC_L5TTPC_STDP.yaml"
BASIS_DIR = "basis_results_edges_ion_channels"


def pairs_on_disk(out_dir):
    """The pairs simwriter already materialised, read back from the dir names."""
    out = []
    for d in sorted(glob.glob(os.path.join(out_dir, "*-*"))):
        pre, post = os.path.basename(d).split("-")
        out.append((int(pre), int(post)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=CONFIG)
    ap.add_argument("--apply", action="store_true",
                    help="write sim files for the pairs find_pairs added")
    ap.add_argument("--retire", action="store_true",
                    help="move the dropped pairs' workdirs and basis csv out of the way")
    args = ap.parse_args()

    w = OptSimWriter(args.config)
    out_dir = w.out_dir
    old = pairs_on_disk(out_dir)
    print("on disk : %d pairs in %s" % (len(old), out_dir))

    print("re-running find_pairs (npairs=%d, seed=%d) ..." % (w.npairs, w.seed))
    new = [(int(a), int(b)) for a, b in w.find_pairs()]
    print("find_pairs returned %d pairs" % len(new))

    old_s, new_s = set(old), set(new)
    kept    = sorted(old_s & new_s)
    dropped = sorted(old_s - new_s)
    added   = sorted(new_s - old_s)

    print("\nkept    : %d" % len(kept))
    print("dropped : %d  %s" % (len(dropped), dropped))
    print("added   : %d  %s" % (len(added), added))

    # The whole point of the exercise: n stays at npairs, and the swap is
    # confined to the pairs we meant to swap.
    problems = []
    if len(new) != w.npairs:
        problems.append("find_pairs returned %d, not npairs=%d" % (len(new), w.npairs))
    if len(dropped) != len(added):
        problems.append("dropped %d but added %d" % (len(dropped), len(added)))
    if len(kept) != len(old) - len(dropped):
        problems.append("kept count inconsistent")
    if problems:
        sys.exit("\nREFUSING TO APPLY:\n  " + "\n  ".join(problems))

    # Same post cell re-paired to a different pre cell would mean the seeding
    # argument above is wrong -- catch it rather than silently half-swapping.
    reused = {p for _, p in dropped} & {p for _, p in added}
    if reused:
        sys.exit("\nREFUSING TO APPLY: post gid(s) %s appear on both sides" % sorted(reused))

    if not args.apply:
        print("\n(dry run -- nothing written; re-run with --apply)")
        return

    if added:
        print("\nwriting sim files for %d new pair(s) ..." % len(added))
        w.write_sim_files(added)          # only these -- the 98 are not touched
        for pre, post in added:
            n = len(glob.glob(os.path.join(out_dir, "%d-%d" % (pre, post), "*Hz_*ms")))
            print("  %d-%d -> %d workdirs" % (pre, post, n))

    if args.retire and dropped:
        graveyard = os.path.join(os.path.dirname(out_dir.rstrip("/")), "retired_pairs")
        os.makedirs(graveyard, exist_ok=True)
        for pre, post in dropped:
            src = os.path.join(out_dir, "%d-%d" % (pre, post))
            if os.path.isdir(src):
                shutil.move(src, os.path.join(graveyard, "%d-%d" % (pre, post)))
                print("  retired %s" % src)
            b = os.path.join(BASIS_DIR, "basis_%d_%d.csv" % (pre, post))
            if os.path.isfile(b):
                shutil.move(b, os.path.join(graveyard, os.path.basename(b)))
                print("  retired %s" % b)

    print("\non disk now: %d pairs" % len(pairs_on_disk(out_dir)))


if __name__ == "__main__":
    main()
