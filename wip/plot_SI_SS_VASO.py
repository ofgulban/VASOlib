"""Replicate Renzo's thesis Fig. 3.2 Panel C Page 47 from gnuplot source."""

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
    """Compute Mz for condition: time < Tr+Ti2."""
    t = time - Ti
    return 1 - np.exp(-t / T1)


def stage3_gm(time, T1, Ti2):
    """Compute Mz for condition: time < Tr+Tr."""
    t = time - Tr - Ti2
    return 1 - np.exp(-t / T1)


def stage4_gm(time, T1, Ti2):
    """Compute Mz for condition: time < Tr+Tr+Ti."""
    t = time - 2*Tr
    return 1 - np.exp(-t / T1) * (1 - np.exp((-Tr+Ti2) / T1)+1)


def stage5_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Ti2."""
    t = time - 2*Tr - Ti
    return 1 - np.exp(-t / T1)


def stage6_gm(time, T1, Ti2):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr."""
    t = time - 3*Tr - Ti2
    return 1 - np.exp(-t / T1)


def stage7_gm(time, T1, Ti2):
    """Compute Mz for condition: time < Tr+Tr+Ti+Tr+Tr."""
    t = time - 4*Tr
    return 1 - np.exp(-t / T1) * (1 - np.exp((-Tr+Ti2) / T1) + 1)


def stage8_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr+Ti."""
    t = time - 4*Tr - Ti
    return 1 - np.exp(-t / T1)


def stage9_gm(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr+Tr."""
    t = time - 5*Tr - Ti
    return 1 - np.exp(-t / T1)


# =============================================================================
# Blood (direct translation of gnuplot conditions)
# =============================================================================
def stage1_blood(time, T1):
    """Compute Mz for condition: time < Ti2+Tr."""
    t = time
    return 1 - 2 * np.exp(-t / T1)


def stage2_blood(time, T1, Ti2):
    """Compute Mz for condition: time < Tr+Tr."""
    t = time - Tr - Ti2
    return 1 - 1 * np.exp(-t / T1)


def stage3_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Ti2."""
    t = time - 2*Tr
    return 1 - 2 * np.exp(-t / T1)


def stage4_blood(time, T1, Ti2):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr."""
    t = time - 3*Tr - Ti2
    return 1 - 1 * np.exp(-t / T1)


def stage5_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr."""
    t = time - 4*Tr
    return 1 - 2 * np.exp(-t / T1)


def stage6_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr+Ti."""
    t = time - 4*Tr
    return 1 - 2 * np.exp(-t / T1)


def stage7_blood(time, T1, Ti):
    """Compute Mz for condition: time < Tr+Tr+Tr+Tr+Tr+Tr."""
    t = time - 5*Tr
    return 1 - 2 * np.exp(-t / T1)


# =============================================================================
# Initial parameters (direct translation of gnuplot conditions)
# =============================================================================
T1csf = 5
T1gm = 1.9
T1b = 2.1  # steady state blood

Ti = 1.45561
Ti2 = 1.7
Tr = 2.

x = 5 * Tr + 0.06
time = np.linspace(0, x, 1001)

# -----------------------------------------------------------------------------
# Compute plot data
# -----------------------------------------------------------------------------
signal1 = np.zeros(time.shape)
idx = time < Tr+Tr+Tr+Tr+Tr+Tr
signal1[idx] = stage9_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Tr+Tr+Tr+Ti
signal1[idx] = stage8_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Tr+Tr+Ti
signal1[idx] = stage7_gm(time[idx], T1gm, Ti2)
idx = time < Tr+Tr+Tr+Tr
signal1[idx] = stage6_gm(time[idx], T1gm, Ti2)
idx = time < Tr+Tr+Tr+Ti2
signal1[idx] = stage5_gm(time[idx], T1gm, Ti)
idx = time < Tr+Tr+Ti
signal1[idx] = stage4_gm(time[idx], T1gm, Ti2)
idx = time < Tr+Tr
signal1[idx] = stage3_gm(time[idx], T1gm, Ti2)
idx = time < Tr+Ti2
signal1[idx] = stage2_gm(time[idx], T1gm, Ti)
idx = time < Ti
signal1[idx] = stage1_gm(time[idx], T1gm)

signal2 = np.zeros(time.shape)
idx = time < Tr+Tr+Tr+Tr+Tr+Tr
signal2[idx] = stage7_blood(time[idx], T1b, Ti)
idx = time < Tr+Tr+Tr+Tr+Tr+Ti
signal2[idx] = stage6_blood(time[idx], T1b, Ti)
idx = time < Tr+Tr+Tr+Tr+Tr
signal2[idx] = stage5_blood(time[idx], T1b, Ti)
idx = time < Tr+Tr+Tr+Tr
signal2[idx] = stage4_blood(time[idx], T1b, Ti2)
idx = time < Tr+Tr+Tr+Ti2
signal2[idx] = stage3_blood(time[idx], T1b, Ti)
idx = time < Tr+Tr
signal2[idx] = stage2_blood(time[idx], T1b, Ti2)
idx = time < Ti2+Tr
signal2[idx] = stage1_blood(time[idx], T1b)


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Faruk's rework start from here
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def Mz(time, M0_equi, M0_init, FA_rad, T1):
    """Longitudinal magnetization."""
    return M0_equi - (M0_equi - M0_init * np.cos(FA_rad)) * np.exp(-time/T1)


# Prepare condition array
cond = np.full(time.shape, 3)
idx2 = time % (Tr*2) < (Ti2+Tr)  # Stages after 180 deg pulse
idx1 = time % (Tr*2) < Ti  # Stages after the first 90 deg pulse
cond[idx2] = 2
cond[idx1] = 1

# =============================================================================
# Compute signal array
# =============================================================================
signal3 = np.zeros(time.shape)

# Prepare initial parameters
M0_equi = 1.  # This never changes
M0_init = 1.
for i, t in enumerate(time):
    t %= 2*Tr

    # -------------------------------------------------------------------------
    # Handle first signal separately
    # -------------------------------------------------------------------------
    if i == 0:
        signal3[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(180), T1=T1gm)

    # -------------------------------------------------------------------------
    # After 180 degree pulse
    # -------------------------------------------------------------------------
    elif cond[i] == 1:
        if cond[i] != cond[i-1]:  # Update M0 upon condition switch
            M0_init = Mz(time=Tr-Ti2, M0_equi=M0_equi, M0_init=M0_init,
                         FA_rad=np.deg2rad(90), T1=T1gm)
        signal3[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(180), T1=T1gm)

    # -------------------------------------------------------------------------
    # After the first 90 degree pulse
    # -------------------------------------------------------------------------
    elif cond[i] == 2:
        signal3[i] = Mz(time=t-Ti, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(90), T1=T1gm)

    # -------------------------------------------------------------------------
    # After the second 90 degree pulse
    # -------------------------------------------------------------------------
    else:
        signal3[i] = Mz(time=t-Tr-Ti2, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(90), T1=T1gm)

# =============================================================================
# Compute signal array
# =============================================================================
signal4 = np.zeros(time.shape)

# Prepare initial parameters
M0_equi = 1.  # This never changes
M0_init = 1.
for i, t in enumerate(time):
    t %= 2*Tr

    # -------------------------------------------------------------------------
    # Handle first signal separately
    # -------------------------------------------------------------------------
    if i == 0:
        signal4[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(180), T1=T1b)

    # -------------------------------------------------------------------------
    # After 180 degree pulse
    # -------------------------------------------------------------------------
    elif cond[i] == 1:
        signal4[i] = Mz(time=t, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(180), T1=T1b)

    # -------------------------------------------------------------------------
    # After the first 90 degree pulse
    # -------------------------------------------------------------------------
    elif cond[i] == 2:
        signal4[i] = Mz(time=t-Ti, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(90), T1=T1b)

    # -------------------------------------------------------------------------
    # After the second 90 degree pulse
    # -------------------------------------------------------------------------
    else:
        signal4[i] = Mz(time=t-Tr-Ti2, M0_equi=M0_equi, M0_init=M0_init,
                        FA_rad=np.deg2rad(90), T1=T1b)


# =============================================================================
# Plotting
# =============================================================================
plt.plot(time, signal1, linewidth=2, color="blue")
plt.plot(time, signal2, linewidth=2, color="red")

plt.plot(time, signal3, linewidth=2, color="cyan")
plt.plot(time, signal4, linewidth=2, color="orange")

plt.title("SI-SS-VASO")
plt.xlabel("Time [s]")
plt.ylabel(r"$M_z$")
plt.ylim([-1, 1])
plt.grid(True)
plt.legend(["R: Gray Matter", "R: Steady State Blood",
            "F: Gray Matter", "F: Steady State Blood"])


plt.show()
