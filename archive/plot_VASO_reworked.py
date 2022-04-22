"""Reproduce Renzo's thesis Fig. 3.2 Panel A Page 47."""

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Gray matter (direct translation of gnuplot conditions)
# =============================================================================
def stage1_gm(time, T1):
    """Compute Mz for condition: time < Ti."""
    t = time
    return 1 - 2 * np.exp(-t / T1)


def stage2_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr."""
    t = time - Ti
    return 1 - np.exp(-t / T1)


def stage3_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Ti."""
    t = time - Tr - Tr
    return 1 - np.exp(-t / T1) * (1 + (1 - np.exp((-(Tr+Tr)+Ti) / T1)))


def stage4_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr."""
    t = time - 2*Tr - Ti
    return 1 - np.exp(-t / T1)


def stage5_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Ti+Tr+Tr."""
    t = time - 4*Tr
    return 1 - np.exp(-t / T1) * (1 + (1 - np.exp((-(4*Tr)+Ti+Tr+Tr) / T1gm)))


def stage6_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr+Tr."""
    t = time - 4*Tr - Ti
    return 1 - np.exp(-t / T1)


# =============================================================================
# Blood (direct translation of gnuplot conditions)
# =============================================================================
# NOTE[Faruk]: Renzo seems to skip two stages here because 90 deg pulse has no
# effect on the longitudinal recovery as it happens on zero crossing for blood
def stage1_blood(time, T1):
    """Compute Mz for condition: time < Ti."""
    t = time
    return 1 - 2 * np.exp(-t / T1)


def stage2_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr."""
    t = time - Ti
    return 1 - np.exp(-t / T1)


def stage3_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr."""
    t = time - 2*Tr
    return 1 - np.exp(-t / T1) * (1 + (1 - np.exp((-(2*Tr)+Ti) / T1)))


def stage4_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr+Tr."""
    t = time - 4*Tr
    return 1 - np.exp(-t / T1) * (1 + 1 - np.exp((-(2*Tr)) / T1) * (1 + (1 - 2 * np.exp(-(2*Tr) / T1))))


# =============================================================================
# Initial parameters (direct translation of gnuplot conditions)
# =============================================================================
T1csf = 5.
T1gm = 1.9
T1b = 2.1  # steady state blood
Tr = 2.
T1 = T1b  # für Tr= 3
Ti = (np.log(2) - np.log(1 + np.exp(-2 * Tr / T1))) * T1

x = 5 * Tr + 0.06
time = np.linspace(0, x, 1001)

# -----------------------------------------------------------------------------
# Compute plot data
# -----------------------------------------------------------------------------
signal1 = np.zeros(time.shape)
idx = time < Tr+Tr+Tr+Tr+Tr+Tr
signal1[idx] = stage6_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Ti+Tr+Tr
signal1[idx] = stage5_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Tr+Tr
signal1[idx] = stage4_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Ti
signal1[idx] = stage3_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr
signal1[idx] = stage2_gm(time[idx], T1gm, Ti)
idx = time < Ti
signal1[idx] = stage1_gm(time[idx], T1gm)

signal2 = np.zeros(time.shape)
idx = time < Tr+Tr+Tr+Tr+Tr+Tr
signal2[idx] = stage4_blood(time[idx], T1b, Ti)
idx = time < Tr+Tr+Tr+Tr
signal2[idx] = stage3_blood(time[idx], T1b, Ti)
idx = time < Tr+Tr
signal2[idx] = stage2_blood(time[idx], T1b, Ti)
idx = time < Ti
signal2[idx] = stage1_blood(time[idx], T1b)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Faruk's rework start from here
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def Mz(time, M0_equi, M0_init, FA_rad, T1):
    """Longitudinal magnetization."""
    return M0_equi - (M0_equi - M0_init * np.cos(FA_rad)) * np.exp(-time/T1)


def compute_VASO_Mz_signal(time, T1, Tr, Ti):
    """Compute VASO Mz signal."""
    signal = np.zeros(time.shape)
    M0_equi = 1.  # This never changes
    M0_init = 1.

    # Prepare condition array
    cond = np.full(time.shape, 2)
    idx1 = time % (Tr*2) <= Ti  # Stages after 180 deg pulse
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
        # 180 degree pulse
        # ---------------------------------------------------------------------
        elif cond[i] == 1:
            if cond[i] != cond[i-1]:  # Update M0 upon condition switch
                M0_init = Mz(time=Tr+Tr-Ti, M0_equi=M0_equi, M0_init=M0_init,
                             FA_rad=np.deg2rad(90), T1=T1)
            signal[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(180), T1=T1)
        # ----------------------------------------------------------------------
        # 90 degree pulse
        # ----------------------------------------------------------------------
        else:
            signal[i] = Mz(time=t-Ti, M0_equi=M0_equi, M0_init=M0_init,
                           FA_rad=np.deg2rad(90), T1=T1)
    return signal


signal3 = compute_VASO_Mz_signal(time, T1gm, Tr, Ti)
signal4 = compute_VASO_Mz_signal(time, T1b, Tr, Ti)

# =============================================================================
# Plotting
# =============================================================================
plt.plot(time, signal1, linewidth=3, color="blue")
plt.plot(time, signal2, linewidth=3, color="red")

plt.plot(time, signal3, linewidth=2, color="cyan")
plt.plot(time, signal4, linewidth=2, color="orange")

plt.title("Original VASO")
plt.xlabel("Time [s]")
plt.ylabel(r"$M_z$")
plt.ylim([-1, 1])
plt.grid(True)
plt.legend(["R: Gray Matter", "R: Steady State Blood",
            "F: Gray Matter", "F: Steady State Blood"])

plt.show()
