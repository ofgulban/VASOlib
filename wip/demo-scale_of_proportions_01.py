"""Demonstrate scale of proportions problem of compositional data.

1. Simulate a single voxel, noiseless VASO signal (BOLD and Nulled).
2. Impose percent signal change activity activity on BOCO (Nulled/BOLD).
3. Demonstrate 'scale of proportions' problem of evolution of ratios which
    results in incosistent statistical inferences.
4. Demonstrate symmetrical and coherent inferences arise from compositional
    treatment of the sampling space.
"""

import numpy as np
import compoda.core as coda


# =============================================================================
def Mz(time, M0_equi, M0_init, FA_rad, T1):
    """Longitudinal magnetization."""
    return M0_equi - (M0_equi - M0_init * np.cos(FA_rad)) * np.exp(-time/T1)


def compute_SI_SS_VASO_Mz_signal(time, T1, Tr, Ti1, Ti2, mode_nonblood=False):
    """Compute VASO Mz signal."""
    signal = np.zeros(time.shape)
    M0_equi = 1.  # This never changes
    M0_init = 1.

    # Prepare condition array
    cond = np.full(time.shape, 3)
    idx2 = time % (Tr*2) < (Ti2+Tr)  # Stages after 180 deg pulse
    idx1 = time % (Tr*2) < Ti1  # Stages after the first 90 deg pulse
    cond[idx2] = 2
    cond[idx1] = 1
    for i, t in enumerate(time):
        t %= 2*Tr
        # ---------------------------------------------------------------------
        # Handle first signal separately
        # ---------------------------------------------------------------------
        if i == 0:
            signal[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(180), T1=T1)
        # ---------------------------------------------------------------------
        # After 180 degree pulse
        # ---------------------------------------------------------------------
        elif cond[i] == 1:
            if mode_nonblood:
                if cond[i] != cond[i-1]:  # Update M0 upon condition switch
                    M0_init = Mz(time=Tr-Ti2, M0_equi=M0_equi, M0_init=M0_init,
                                 FA_rad=np.deg2rad(90), T1=T1)
            signal[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(180), T1=T1)
        # ---------------------------------------------------------------------
        # After the first 90 degree pulse
        # ---------------------------------------------------------------------
        elif cond[i] == 2:
            signal[i] = Mz(time=t-Ti1, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(90), T1=T1)
        # ---------------------------------------------------------------------
        # After the second 90 degree pulse
        # ---------------------------------------------------------------------
        else:
            signal[i] = Mz(time=t-Tr-Ti2, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(90), T1=T1)
    return signal


# =============================================================================
# Initial parameters
# =============================================================================
Ti1 = 1.45561 * 1000
Ti2 = 1.7 * 1000
Tr = 2. * 1000
max_time = 4 * Tr

# Compute readout times
time_readout1 = np.arange(Ti1, max_time, Tr*2)
time_readout2 = np.arange(Tr+Ti2, max_time, Tr*2)

nr_timepoints = time_readout1.size + time_readout2.size
time_readouts = np.zeros(nr_timepoints)
time_readouts[0::2] = time_readout1
time_readouts[1::2] = time_readout2

# =============================================================================
# Compute VASO signal
# =============================================================================
T1gm = 2.1 * 1000
T1blood = 1.9 * 1000

data = compute_SI_SS_VASO_Mz_signal(
    time_readouts-1, T1=T1gm, Tr=Tr, Ti1=Ti1, Ti2=Ti2, mode_nonblood=True)

bold = data[1]
nonblood = data[2]  # aka nulled
boco_rest = nonblood / bold

# Compute activity signal
boco_act = boco_rest - (boco_rest * 0.01)

# Signal change
delta1 = (boco_act - boco_rest) / boco_rest
print("BOCO ~ (1-CBV) signal change: {:.2f}%".format(delta1*100))

# -----------------------------------------------------------------------------
# Let's have a look at CBV (1-BOCO)
# -----------------------------------------------------------------------------
counter_boco_rest = 1 - boco_rest
counter_boco_act = 1 - boco_act

# Signal change
delta2 = (counter_boco_act - counter_boco_rest) / counter_boco_rest
print("(1-BOCO) ~ CBV signal change: {:.2f}%".format(delta2*100))

# -----------------------------------------------------------------------------
# Let's consider the sampling space as simplex
# -----------------------------------------------------------------------------
comp_rest = coda.closure(np.vstack([boco_rest, 1-boco_rest]).T)
clr_rest = coda.clr_transformation(comp_rest)

comp_act = coda.closure(np.vstack([boco_act, 1-boco_act]).T)
clr_act = coda.clr_transformation(comp_act)

# Signal change
delta3 = (clr_act[0, 0] - clr_rest[0, 0]) / clr_rest[0, 0]
delta4 = (clr_act[0, 1] - clr_rest[0, 1]) / clr_rest[0, 1]
print("COCO ~ (1-CBV) signal change: {:.2f}%".format(delta3*100))
print("1-COCO ~ (CBV) signal change: {:.2f}%".format(delta4*100))
