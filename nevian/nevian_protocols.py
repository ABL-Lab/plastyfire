"""
Nevian & Sakmann 2006 — Protocol Definitions
=============================================
Defines all STDP induction protocols from:
  Nevian T, Sakmann B (2006) Spine Ca2+ signaling in spike-timing-dependent plasticity.
  J Neurosci 26(43):11001-11013.

Preparation: L2/3 pyramidal neurons, somatosensory cortex (barrel field), Wistar P13-P15.
Each protocol: 60 pairings at 0.1 Hz. EPSP amplitude measured 20-40 min post-induction.

Timing convention (matching the existing plastyfire simwriter/simulator):
  dt > 0 → APs follow EPSP   (LTP direction)
  dt < 0 → APs precede EPSP  (LTD direction)
  dt here always refers to time from the EPSP onset to the *closest* AP (Δt in paper).

Burst timing convention:
  For a burst of 3 APs at 50 Hz preceding the EPSP:
    - Paper reports Δt' = -50 ms (time from first AP to EPSP onset)
    - Equivalent to Δt = -10 ms (last AP to EPSP onset, since ISI=20ms)
    - We store dt = -10 ms (the closest-AP convention used by the simulator)
  For APs following the EPSP:
    - Δt = Δt' (first AP = closest AP since they all follow)
"""

from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class NevianProtocol:
    """
    Represents one Nevian & Sakmann 2006 STDP pairing protocol.

    Fields
    ------
    protocol_id : str
        Unique identifier (used as simulation folder name suffix)
    description : str
        Human-readable description
    n_aps : int
        Number of postsynaptic APs in the burst
    freq_hz : float
        Intra-burst AP frequency (Hz). Ignored if n_aps == 1.
    dt_ms : float
        Time from EPSP onset to the AP *closest* in time to the EPSP (ms).
        Positive = APs follow EPSP (LTP side), Negative = APs precede EPSP (LTD side).
    n_pairings : int
        Number of pre-post pairing repetitions (default: 60, same as paper)
    pairing_period_ms : float
        Period between pairings = 1 / 0.1 Hz = 10,000 ms
    expected_epsp_ratio : float
        Mean EPSP ratio observed in paper (post/pre amplitude, 20-40 min after induction)
    expected_sem : float
        SEM of EPSP ratio from paper
    expected_n : int
        Sample size from paper
    expected_plasticity : str
        Qualitative outcome: 'LTP', 'LTD', or 'none'
    condition : str
        Pharmacological condition: 'control' or name of blocker
    blocker : Optional[str]
        Blocker drug name if any (e.g. 'D-APV', 'MK-801', 'MCPG', 'AM251', etc.)
    blocker_route : Optional[str]
        Route of blocker delivery: 'bath' or 'intracellular'
    figure_ref : str
        Figure reference in paper
    notes : str
        Any additional notes
    """
    protocol_id: str
    description: str
    n_aps: int
    freq_hz: float
    dt_ms: float
    expected_epsp_ratio: float
    expected_sem: float
    expected_n: int
    expected_plasticity: str
    n_pairings: int = 60
    pairing_period_ms: float = 10000.0   # 0.1 Hz = 10 s between pairings
    condition: str = "control"
    blocker: Optional[str] = None
    blocker_route: Optional[str] = None
    figure_ref: str = ""
    notes: str = ""


# =============================================================================
# GROUP 1: Core Timing Sweep (3 APs @ 50 Hz)  — Fig. 1, 2
# =============================================================================
TIMING_SWEEP = [
    NevianProtocol(
        protocol_id="pre_only_far",
        description="3 APs 50Hz, burst precedes EPSP by Δt'=-90ms (Δt=-50ms last AP)",
        n_aps=3, freq_hz=50, dt_ms=-50,
        expected_epsp_ratio=1.00, expected_sem=0.09, expected_n=6,
        expected_plasticity="none",
        figure_ref="Fig 2B",
        notes="Far pre: no change. Δt'=-90ms → last AP at Δt=-50ms from EPSP.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms",
        description="3 APs 50Hz, last AP precedes EPSP by Δt=-10ms (Δt'=-50ms) → LTD",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=0.68, expected_sem=0.05, expected_n=10,
        expected_plasticity="LTD",
        figure_ref="Fig 1C, 2B",
        notes="CANONICAL LTD protocol. Δt'=-50ms burst onset, Δt=-10ms last AP to EPSP.",
    ),
    NevianProtocol(
        protocol_id="neutral_3ap_50hz_dt-30ms",
        description="3 APs 50Hz, Δt'=-30ms (2 APs before + 1 after EPSP) → no change",
        n_aps=3, freq_hz=50, dt_ms=-30,
        expected_epsp_ratio=0.98, expected_sem=0.12, expected_n=8,
        expected_plasticity="none",
        figure_ref="Fig 2B",
        notes="EPSP is evoked within the burst (2 APs before, 1 after). No plasticity.",
    ),
    NevianProtocol(
        protocol_id="mild_LTP_3ap_50hz_dt-10ms_burst",
        description="3 APs 50Hz, Δt'=-10ms (1 AP before + 2 after EPSP) → mild LTP",
        n_aps=3, freq_hz=50, dt_ms=10,   # The AP closest to EPSP is +10ms (2 APs after)
        expected_epsp_ratio=1.42, expected_sem=0.19, expected_n=12,
        expected_plasticity="LTP",
        figure_ref="Fig 2B",
        notes="EPSP within burst: 1 AP before at -10ms, 2 APs after at +10ms, +30ms. "
              "dt_ms=+10 used because the closest AP (2nd in burst) follows EPSP by 10ms.",
    ),
    NevianProtocol(
        protocol_id="LTP_3ap_50hz_dt+10ms",
        description="3 APs 50Hz, all APs follow EPSP, first AP at Δt=+10ms → strong LTP",
        n_aps=3, freq_hz=50, dt_ms=10,
        expected_epsp_ratio=2.01, expected_sem=0.22, expected_n=11,
        expected_plasticity="LTP",
        figure_ref="Fig 1B, 2B",
        notes="CANONICAL LTP protocol. All 3 APs follow EPSP by +10, +30, +50ms.",
    ),
    NevianProtocol(
        protocol_id="post_only_far",
        description="3 APs 50Hz, all APs far after EPSP, Δt=+50ms → no change",
        n_aps=3, freq_hz=50, dt_ms=50,
        expected_epsp_ratio=0.92, expected_sem=0.11, expected_n=4,
        expected_plasticity="none",
        figure_ref="Fig 2B",
        notes="Far post: no change.",
    ),
]

# =============================================================================
# GROUP 2: Number-of-APs Dependence (50 Hz, Δt = ±10 ms)  — Fig. 3A,B
# =============================================================================
N_APS_SWEEP = [
    NevianProtocol(
        protocol_id="epsp_only",
        description="EPSP only (no APs) — control for pairing",
        n_aps=0, freq_hz=0, dt_ms=0,
        expected_epsp_ratio=0.97, expected_sem=0.08, expected_n=5,
        expected_plasticity="none",
        figure_ref="Fig 3A,B",
        notes="Synaptic stimulation without postsynaptic APs. No plasticity.",
    ),
    NevianProtocol(
        protocol_id="LTD_1ap_50hz_dt-10ms",
        description="1 AP precedes EPSP by Δt=-10ms → LTD",
        n_aps=1, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=0.80, expected_sem=0.07, expected_n=5,
        expected_plasticity="LTD",
        figure_ref="Fig 3A,B",
        notes="Single AP LTD. Requires T-VDCC Ca2+ influx.",
    ),
    NevianProtocol(
        protocol_id="nochange_1ap_50hz_dt+10ms",
        description="1 AP follows EPSP by Δt=+10ms → no change (LTP threshold not met)",
        n_aps=1, freq_hz=50, dt_ms=10,
        expected_epsp_ratio=1.04, expected_sem=0.08, expected_n=10,
        expected_plasticity="none",
        figure_ref="Fig 3A,B",
        notes="Single AP post is insufficient for LTP. Need ≥2 APs.",
    ),
    NevianProtocol(
        protocol_id="LTD_2ap_50hz_dt-10ms",
        description="2 APs 50Hz, last AP precedes EPSP by Δt=-10ms → LTD",
        n_aps=2, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=0.72, expected_sem=0.12, expected_n=5,
        expected_plasticity="LTD",
        figure_ref="Fig 3A,B",
        notes="Minimal burst LTD.",
    ),
    NevianProtocol(
        protocol_id="LTP_2ap_50hz_dt+10ms",
        description="2 APs 50Hz, first AP follows EPSP by Δt=+10ms → LTP",
        n_aps=2, freq_hz=50, dt_ms=10,
        expected_epsp_ratio=1.95, expected_sem=0.31, expected_n=9,
        expected_plasticity="LTP",
        figure_ref="Fig 3A,B",
        notes="Minimal burst required for LTP. 2 APs at ≥50Hz.",
    ),
    # 3 APs protocols already in TIMING_SWEEP — LTD_3ap_50hz_dt-10ms, LTP_3ap_50hz_dt+10ms
]

# =============================================================================
# GROUP 3: AP Frequency Dependence (3 APs, Δt = ±10 ms)  — Fig. 3C,D
# =============================================================================
FREQ_SWEEP = [
    NevianProtocol(
        protocol_id="nochange_3ap_20hz_dt+10ms",
        description="3 APs 20Hz, first AP follows EPSP by Δt=+10ms → no change",
        n_aps=3, freq_hz=20, dt_ms=10,
        expected_epsp_ratio=1.09, expected_sem=0.27, expected_n=5,
        expected_plasticity="none",
        figure_ref="Fig 3C,D",
        notes="20Hz burst insufficient for LTP. Minimum frequency >20Hz needed.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_20hz_dt-10ms",
        description="3 APs 20Hz, last AP precedes EPSP by Δt=-10ms → LTD",
        n_aps=3, freq_hz=20, dt_ms=-10,
        expected_epsp_ratio=0.72, expected_sem=0.14, expected_n=7,
        expected_plasticity="LTD",
        figure_ref="Fig 3C,D",
        notes="LTD is frequency-independent for pre-EPSP bursts.",
    ),
    NevianProtocol(
        protocol_id="LTP_3ap_100hz_dt+10ms",
        description="3 APs 100Hz, first AP follows EPSP by Δt=+10ms → strong LTP",
        n_aps=3, freq_hz=100, dt_ms=10,
        expected_epsp_ratio=2.29, expected_sem=0.48, expected_n=7,
        expected_plasticity="LTP",
        figure_ref="Fig 3C,D",
        notes="Highest LTP in the paper. 100Hz burst post.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_100hz_dt-10ms",
        description="3 APs 100Hz, last AP precedes EPSP by Δt=-10ms → strong LTD",
        n_aps=3, freq_hz=100, dt_ms=-10,
        expected_epsp_ratio=0.52, expected_sem=0.12, expected_n=4,
        expected_plasticity="LTD",
        figure_ref="Fig 3C,D",
        notes="Strong LTD at 100Hz burst pre. MCPG flips this to LTP (see BLOCKER group).",
    ),
]

# =============================================================================
# GROUP 4: NMDA Receptor Blocker Experiments  — Fig. 6
# =============================================================================
NMDAR_BLOCKERS = [
    # --- D-APV (bath, blocks all NMDARs) ---
    NevianProtocol(
        protocol_id="LTP_3ap_50hz_dt+10ms_dapv",
        description="3 APs 50Hz dt=+10ms + bath D-APV (50µM) → LTP abolished",
        n_aps=3, freq_hz=50, dt_ms=10,
        expected_epsp_ratio=1.08, expected_sem=0.11, expected_n=3,
        expected_plasticity="none",
        condition="D-APV_bath",
        blocker="D-APV", blocker_route="bath",
        figure_ref="Fig 6D,H",
        notes="Bath D-APV abolishes LTP. NMDAR activation required for LTP.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_dapv",
        description="3 APs 50Hz dt=-10ms + bath D-APV (50µM) → LTD abolished",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=0.98, expected_sem=0.11, expected_n=6,
        expected_plasticity="none",
        condition="D-APV_bath",
        blocker="D-APV", blocker_route="bath",
        figure_ref="Fig 6E,I",
        notes="Bath D-APV also abolishes LTD (removes Ca2+ entirely via NMDAR+VDCC both).",
    ),
    # --- MK-801 (intracellular = postsynaptic NMDAR block only) ---
    NevianProtocol(
        protocol_id="LTP_3ap_50hz_dt+10ms_mk801",
        description="3 APs 50Hz dt=+10ms + intracell MK-801 (1mM) → LTP abolished",
        n_aps=3, freq_hz=50, dt_ms=10,
        expected_epsp_ratio=0.95, expected_sem=0.19, expected_n=6,
        expected_plasticity="none",
        condition="MK-801_intracellular",
        blocker="MK-801", blocker_route="intracellular",
        figure_ref="Fig 6F,H",
        notes="Postsynaptic NMDAR block abolishes LTP. LTP requires postsynaptic NMDARs.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_mk801",
        description="3 APs 50Hz dt=-10ms + intracell MK-801 (1mM) → LTD NOT affected",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=0.60, expected_sem=0.11, expected_n=7,
        expected_plasticity="LTD",
        condition="MK-801_intracellular",
        blocker="MK-801", blocker_route="intracellular",
        figure_ref="Fig 6G,I",
        notes="KEY RESULT: postsynaptic NMDAR block does NOT affect LTD. "
              "LTD is independent of postsynaptic NMDARs; VDCC Ca2+ is sufficient.",
    ),
]

# =============================================================================
# GROUP 5: VDCC Blocker Experiments  — Fig. 7
# =============================================================================
VDCC_BLOCKERS = [
    NevianProtocol(
        protocol_id="LTP_3ap_50hz_dt+10ms_nimodipine",
        description="3 APs 50Hz dt=+10ms + nimodipine (10µM, L-VDCC) → LTP unaffected",
        n_aps=3, freq_hz=50, dt_ms=10,
        expected_epsp_ratio=1.92, expected_sem=0.31, expected_n=7,
        expected_plasticity="LTP",
        condition="nimodipine_bath",
        blocker="nimodipine", blocker_route="bath",
        figure_ref="Fig 7D",
        notes="L-VDCC block has no effect on LTP.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_nimodipine",
        description="3 APs 50Hz dt=-10ms + nimodipine (10µM, L-VDCC) → LTD unaffected",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=0.70, expected_sem=0.10, expected_n=8,
        expected_plasticity="LTD",
        condition="nimodipine_bath",
        blocker="nimodipine", blocker_route="bath",
        figure_ref="Fig 7E",
        notes="L-VDCC alone is not essential for burst LTD.",
    ),
    NevianProtocol(
        protocol_id="LTD_1ap_50hz_dt-10ms_ni2",
        description="1 AP dt=-10ms + NiCl2 (50µM, T-VDCC) → LTD abolished",
        n_aps=1, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=1.05, expected_sem=0.16, expected_n=3,
        expected_plasticity="none",
        condition="NiCl2_bath",
        blocker="NiCl2", blocker_route="bath",
        figure_ref="Fig 7F,I",
        notes="Single-AP LTD requires T-VDCC Ca2+ influx.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_ni2",
        description="3 APs 50Hz dt=-10ms + NiCl2 (50µM, T-VDCC) → LTD NOT abolished",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=0.70, expected_sem=0.03, expected_n=4,
        expected_plasticity="LTD",
        condition="NiCl2_bath",
        blocker="NiCl2", blocker_route="bath",
        figure_ref="Fig 7G,I",
        notes="Burst LTD survives T-VDCC block by recruiting other VDCC subtypes.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_nimo_plus_ni2",
        description="3 APs 50Hz dt=-10ms + nimodipine + NiCl2 (combined L+T VDCC block) → LTD abolished",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=1.00, expected_sem=0.05, expected_n=3,
        expected_plasticity="none",
        condition="nimodipine+NiCl2_bath",
        blocker="nimodipine+NiCl2", blocker_route="bath",
        figure_ref="Fig 7I",
        notes="Combined L+T VDCC block abolishes burst LTD completely. "
              "VDCCs are essential for LTD.",
    ),
]

# =============================================================================
# GROUP 6: mGluR Blocker Experiments  — Fig. 9
# =============================================================================
MGLUR_BLOCKERS = [
    NevianProtocol(
        protocol_id="LTP_3ap_50hz_dt+10ms_mcpg",
        description="3 APs 50Hz dt=+10ms + MCPG (500µM, mGluR antagonist) → LTP unaffected",
        n_aps=3, freq_hz=50, dt_ms=10,
        expected_epsp_ratio=1.79, expected_sem=0.33, expected_n=9,
        expected_plasticity="LTP",
        condition="MCPG_bath",
        blocker="MCPG", blocker_route="bath",
        figure_ref="Fig 9D",
        notes="mGluR block does not affect LTP pathway.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_mcpg",
        description="3 APs 50Hz dt=-10ms + MCPG (500µM) → LTD abolished",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=1.06, expected_sem=0.16, expected_n=7,
        expected_plasticity="none",
        condition="MCPG_bath",
        blocker="MCPG", blocker_route="bath",
        figure_ref="Fig 9E",
        notes="mGluR activation is required for LTD induction.",
    ),
    NevianProtocol(
        protocol_id="LTP_from_LTD_3ap_100hz_dt-10ms_mcpg",
        description="3 APs 100Hz dt=-10ms + MCPG (500µM) → LTD FLIPPED TO LTP",
        n_aps=3, freq_hz=100, dt_ms=-10,
        expected_epsp_ratio=1.53, expected_sem=0.24, expected_n=11,
        expected_plasticity="LTP",
        condition="MCPG_bath",
        blocker="MCPG", blocker_route="bath",
        figure_ref="Fig 9F",
        notes="KEY RESULT: Blocking mGluRs converts high-Ca2+ LTD protocol to LTP. "
              "Large Ca2+ (from 100Hz burst) is sufficient for LTP when mGluR switch is blocked. "
              "Demonstrates mGluR as the LTD/LTP direction switch.",
    ),
]

# =============================================================================
# GROUP 7: LTD Downstream Cascade Blockers  — Fig. 10
# =============================================================================
LTD_CASCADE_BLOCKERS = [
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_u73122",
        description="3 APs 50Hz dt=-10ms + U73122 (5µM, PLC inhibitor) → LTD abolished",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=1.11, expected_sem=0.07, expected_n=4,
        expected_plasticity="none",
        condition="U73122_bath",
        blocker="U73122", blocker_route="bath",
        figure_ref="Fig 10A,F",
        notes="PLC inhibitor blocks LTD. mGluR→PLC→DAG→2-AG cascade required.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_am251",
        description="3 APs 50Hz dt=-10ms + AM251 (2µM, CB1 antagonist) → LTD abolished",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=1.09, expected_sem=0.14, expected_n=5,
        expected_plasticity="none",
        condition="AM251_bath",
        blocker="AM251", blocker_route="bath",
        figure_ref="Fig 10E,F",
        notes="CB1 receptor block abolishes LTD. Endocannabinoid retrograde signaling required.",
    ),
    NevianProtocol(
        protocol_id="LTD_3ap_50hz_dt-10ms_heparin",
        description="3 APs 50Hz dt=-10ms + heparin intracell (400 U/ml, IP3R blocker) → LTD NOT affected",
        n_aps=3, freq_hz=50, dt_ms=-10,
        expected_epsp_ratio=0.68, expected_sem=0.12, expected_n=9,
        expected_plasticity="LTD",
        condition="heparin_intracellular",
        blocker="heparin", blocker_route="intracellular",
        figure_ref="Fig 10C,F",
        notes="IP3-mediated Ca2+ release from internal stores is NOT required for LTD. "
              "Ca2+ from VDCCs alone is sufficient.",
    ),
]

# =============================================================================
# MASTER LIST — all protocols in a flat dict keyed by protocol_id
# =============================================================================
ALL_PROTOCOLS: List[NevianProtocol] = (
    TIMING_SWEEP
    + N_APS_SWEEP
    + FREQ_SWEEP
    + NMDAR_BLOCKERS
    + VDCC_BLOCKERS
    + MGLUR_BLOCKERS
    + LTD_CASCADE_BLOCKERS
)

PROTOCOLS_BY_ID = {p.protocol_id: p for p in ALL_PROTOCOLS}

# Convenience: just the control protocols (no pharmacology)
CONTROL_PROTOCOLS: List[NevianProtocol] = [
    p for p in ALL_PROTOCOLS if p.condition == "control"
]

# The 2 canonical protocols used for fitting / validation
CANONICAL_LTP = PROTOCOLS_BY_ID["LTP_3ap_50hz_dt+10ms"]
CANONICAL_LTD = PROTOCOLS_BY_ID["LTD_3ap_50hz_dt-10ms"]


if __name__ == "__main__":
    print(f"Total protocols: {len(ALL_PROTOCOLS)}")
    print(f"Control protocols: {len(CONTROL_PROTOCOLS)}")
    print(f"\nAll protocol IDs:")
    for p in ALL_PROTOCOLS:
        tag = "★" if p.expected_plasticity != "none" else "·"
        print(f"  {tag} [{p.condition:30s}] {p.protocol_id:55s} → ratio={p.expected_epsp_ratio:.2f} ({p.expected_plasticity})")
