"""Physical constants, read from GluSynapse.mod and the run's fit_params.

Anything here that disagrees with the mod file or the launcher makes every
downstream number wrong, so each value carries its provenance.
"""

# GluSynapse.mod
MIN_CA_CR   = 70e-6   # mM      line 161
RHO_STAR_GB = 0.5     # 1       line 164
TAU_IND_GB  = 70.0    # s       line 165
RHO_BINARY_THRESHOLD = 0.5      # evaluator_edges.compute_epsp_from_basis

# Set per-run by pairrunner_edges_fit.py, NOT the mod defaults
# (mod defaults are gamma_d=100, gamma_p=450 — both overridden).
TAU_EFFCA_GB = 278.3177658387   # ms
GAMMA_D_GB   = 101.5
GAMMA_P_GB   = 199.773931

# Markram et al. 1997, 10 Hz — the fit targets
INVITRO = {-10: (0.7922, 0.0259), 5: (1.2038, 0.0644), 10: (1.2013, 0.0626)}

PROTOCOLS = ["10Hz_-50ms", "10Hz_-30ms", "10Hz_-10ms",
             "10Hz_5ms", "10Hz_10ms", "10Hz_30ms", "10Hz_50ms"]
PROTO_DT  = {p: int(p.split("_")[1].replace("ms", "")) for p in PROTOCOLS}

A_KEYS = ["a00", "a01", "a10", "a11", "a20", "a21", "a30", "a31"]

# The run being validated against
TIED_IC4 = dict(a00=1.002, a01=2.254638, a10=1.209857, a11=2.396290,
                a20=1.002, a21=2.254638, a30=1.209857, a31=2.396290)
TIED_IC4_HASH = "cdf3a1e1db98"
